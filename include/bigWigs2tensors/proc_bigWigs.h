#ifndef PROC_BIGWIGS_H
#define PROC_BIGWIGS_H

#include <vector>
#include <map>
#include <algorithm>
#include <ranges>
#include <execution>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <torch/torch.h>
#include <bigWig.h>

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

class BWBinner
/*!
a "manager" for binning a set of "alignable" bigWig files together
*/
{
public:
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
        std::transform(bigWig_paths.cbegin(), bigWig_paths.cend(),
                       std::back_inserter(bw_paths),
                       [](const std::string& path) {
                           return std::filesystem::path(path);
                       });
    }

    ~BWBinner() {
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
    unsigned num_bws;
    std::map<std::string, int> chrom_sizes;
    std::map<std::string, torch::Tensor> chrom_binneds;
    torch::TensorOptions tens_opts;

    /*!
    Loads all the data (binned series of values) for chromosome `chrom`
    into a torch Tensor, one row per bigWig file.
    */
    torch::Tensor load_bin_chrom_tensor(const std::string& chrom, unsigned bin_size);
};

#endif