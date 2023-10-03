#ifndef PROC_BIGWIGS_H
#define PROC_BIGWIGS_H

#include <vector>
#include <map>
#include <algorithm>
#include <execution>
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <torch/torch.h>
#include <bigWig.h>
#include <bigWigs2tensors/util.h>

namespace constants {
    static const torch::TensorOptions tensor_opts = torch::TensorOptions()
                                                    .dtype(torch::kFloat64)
                                                    .requires_grad(false);
};

/*!
Opens a set of bigWig files whose paths are given by `bw_paths` and
returns a vector of pointers to them.
*/
std::vector<bigWigFile_t*> open_bigWigs(const std::vector<std::string>& bw_paths);

/*!
Reads a whitespace-delimited chrom_sizes file into
a map of chromosome names to their sizes.
*/
std::map<std::string, int> parse_chrom_sizes(const std::string& chrom_sizes_path);

/*!
Given a bin size, bins a vector of doubles into
a vector of their mean averages.
If any NaNs in a bin, sets result for that bin to NaN.
*/
std::vector<double> bin_vec_NaNmeans(const std::vector<double>& in_vec, size_t bin_size);

class BWBinner
/*!
a "manager" for binning a set of "alignable" bigWig files together
*/
{
public:
    /*!
    Constructs a BWBinner object that will bin the specified intervals of the bigWig files
    Args:
        bw_paths: a NULL-terminated array of paths to bigWig files
        chrom_sizes_path: path to a whitespace-delimited file of chromosome sizes
        coords_bed_path: optional path to a bed file of coordinates to bin over
    */
    BWBinner(const std::vector<std::string>& bigWig_paths, const std::string& chrom_sizes_path, const std::string& coords_bed_path)
        : bw_files(open_bigWigs(bigWig_paths)),
        num_bws(bw_files.size()),
        chrom_sizes(parse_chrom_sizes(chrom_sizes_path)),
        tens_opts(constants::tensor_opts)
    {
        // Parse the coordinates bed file
        parse_coords_bed(spec_coords, coords_bed_path, chrom_sizes);
        for (const auto& [chr, interv] : spec_coords) {
            std::cout << chr <<": "<< interv->l <<" intervals"<< std::endl;
            std::cout << '\t';
            for (int i = 0; i < interv->l; i++) {
                std::cout << ' ' << interv->start[i] << '-' << interv->end[i];
            };
            std::cout << std::endl;
        }

        std::transform(bigWig_paths.cbegin(), bigWig_paths.cend(),
                       std::back_inserter(bw_paths),
                       [](const std::string& path) {
                           return std::filesystem::path(path);
                       });
    }

    /*!
    Constructs a BWBinner object that will bin the bigWig files
    Args:
        bw_paths: a NULL-terminated array of paths to bigWig files
        chrom_sizes_path: path to a whitespace-delimited file of chromosome sizes
    */
    BWBinner(const std::vector<std::string>& bigWig_paths, const std::string& chrom_sizes_path)
        : bw_files(open_bigWigs(bigWig_paths)),
        num_bws(bw_files.size()),
        chrom_sizes(parse_chrom_sizes(chrom_sizes_path)),
        tens_opts(constants::tensor_opts)
    {
        // Will use entire chromosomes <=> "full" coordinates <=> 0 to chrom_size
        std::transform(chrom_sizes.begin(), chrom_sizes.end(),
                            std::inserter(spec_coords, spec_coords.end()),
                            [this](const std::pair<std::string, int>& chrom_size) {
                                std::string chrom = chrom_size.first;
                                unsigned size = chrom_size.second;

                                bbOverlappingEntries_t* entries = new bbOverlappingEntries_t;
                        
                                uint32_t* starts = new uint32_t[1];
                                starts[0] = 0;
                                entries->start = starts;

                                uint32_t* ends = new uint32_t[1];
                                starts[0] = size;
                                entries->end = ends;

                                entries->l = 1;
                                entries->m = 1;

                                return std::make_pair(chrom, entries);
                            });
        for (const auto& [chr, interv] : spec_coords) {
            std::cout << chr <<": "<< interv->l <<" intervals"<< std::endl;
            std::cout << '\t';
            for (int i = 0; i < interv->l; i++) {
                std::cout << ' ' << interv->start[i] << '-' << interv->end[i];
            };
            std::cout << std::endl;
        }

        std::transform(bigWig_paths.cbegin(), bigWig_paths.cend(),
                       std::back_inserter(bw_paths),
                       [](const std::string& path) {
                           return std::filesystem::path(path);
                       });
    }

    ~BWBinner() {
        std::for_each(spec_coords.begin(), spec_coords.end(),
                      [](std::pair<std::string, bbOverlappingEntries_t*> entry) {
                          delete[] entry.second->start;
                          delete[] entry.second->end;
                          delete entry.second;
                      });
        for (auto& bw : bw_files) {
            bwClose(bw);
        }
        bwCleanup();
        chrom_binneds.clear();
    }

    /*!
    Loads the binned data for all chromosomes into a map of torch Tensors,
    by calling `load_bin_chrom_tensor` on each chromosome.
    \args bin_size The size of the bins to use.
    */
    std::map<std::string, torch::Tensor> load_bin_all_chroms(unsigned bin_size);

    /*!
    Data getter for the binned data f0.5000or all chromosomes.
    \note Before binning, this will be empty.
    */
    std::map<std::string, torch::Tensor> binned_chroms() const;

    /*!
    Saves the binned data for all chromosomes (each a Tensor) and
    a text file with the bigWig filename stems in their order in the Tensors,
    one per line.
    Each Tensor is saved to a separate file, named by the chromosome name.
    \arg out_dir the directory to save to, without a trailing '/'.
    */
    void save_binneds(const std::string& out_dir) const; 

private:
    std::vector<std::filesystem::path> bw_paths;
    std::vector<bigWigFile_t*> bw_files;
    unsigned int num_bws;
    std::map<std::string, int> chrom_sizes;
    // entries are {chrom: bbOverlappingEntries_t*}
    //  '-> key struct members are array of starts, array of ends, and number of entries
    std::map<std::string, bbOverlappingEntries_t*> spec_coords;
    std::map<std::string, torch::Tensor> chrom_binneds;
    torch::TensorOptions tens_opts;

    /*!
    Loads all the data (binned series of values) for chromosome `chrom`
    into a torch Tensor, one row per bigWig file.
    */
    torch::Tensor load_bin_chrom_tensor(const std::string& chrom, unsigned bin_size);
};

#endif