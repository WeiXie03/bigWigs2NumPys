#include <vector>
#include <map>
#include <algorithm>
#include <ranges>
#include <execution>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <torch/torch.h>
#include "../libs/libBigWig/include/bigWig.h"
#include "../include/bigWigs2tensors/proc_bigWigs.h"

std::vector<bigWigFile_t*> open_bigWigs(const std::vector<std::string>& bw_paths) {
    std::vector<bigWigFile_t*> bw_files;
    for (const auto& path : bw_paths) {
        bw_files.push_back(bwOpen(path.c_str(), NULL, "r"));
    }
    return bw_files;
}

std::map<std::string, int> parse_chrom_sizes(const std::string& chrom_sizes_path) {
    std::map<std::string, int> chrom_sizes;
    std::ifstream chrom_sizes_file(chrom_sizes_path);
    std::string chrom;
    int size;
    while (chrom_sizes_file >> chrom >> size) {
        chrom_sizes[chrom] = size;
    }
    return chrom_sizes;
}

torch::Tensor BWBinner::load_bin_chrom_tensor(const std::string& chrom, unsigned bin_size) {
    // ceiling division: ceil(chrom_size / bin_size)
    // credit: https://stackoverflow.com/a/2745086
    unsigned num_bins = 1 + ((chrom_sizes[chrom] - 1) / bin_size);
    std::cout << "Binning " << chrom << " into " << num_bins << " bins" << std::endl;

    torch::Tensor chrom_tensor = torch::empty({num_bws, num_bins}, tens_opts);
    std::cout << "Created " << chrom_tensor.sizes() << " Tensor `chrom_tensor`. To fill..." << std::endl;
    
    // parallelize across bigWigs' indices within the bw_files vector
    // credit: https://stackoverflow.com/a/62829166
    std::ranges::iota_view bw_idxs((size_t)0, bw_files.size());
    std::for_each(std::execution::par_unseq,
                    bw_idxs.begin(), bw_idxs.end(),
                    [this, chrom, num_bins, &chrom_tensor](size_t bw_idx) {
                        // the whole chrom for this bigWig
                        double* binned_vals = bwStats(bw_files[bw_idx], chrom.c_str(),
                                                    0, chrom_sizes[chrom], num_bins,
                                                    bwStatsType::mean);

                        std::cout << "Inserting Tensor\n";
                        std::cout << "  heap array from libBigWig: [";
                        for (int i = 0; i < num_bins; i++)
                            std::cout << ' ' << binned_vals[i];
                        std::cout << " ]" << std::endl;

                        std::cout << "  into Tensor:\n";
                        std::cout << torch::from_blob(binned_vals, {num_bins}, torch::dtype(torch::kFloat64)) << std::endl;

                        using namespace torch::indexing;
                        // set bw_idx'th row of chrom_tensor to binned_vals
                        chrom_tensor.index_put_({(int)bw_idx, "..."},
                                                torch::from_blob(binned_vals, {num_bins},
                                                                torch::dtype(torch::kFloat64)));
                    });

    return chrom_tensor;
}

void BWBinner::load_bin_all_chroms(unsigned bin_size) {
    // parallelize over chromosomes
    std::for_each(std::execution::par_unseq,
        chrom_sizes.begin(), chrom_sizes.end(),
        //std::pair<std::string, int>
        [this, bin_size](auto& chr_entry) {
            chrom_binneds[chr_entry.first] = BWBinner::load_bin_chrom_tensor(chr_entry.first, bin_size);
        }
    );
}

std::map<std::string, torch::Tensor> BWBinner::binned_chroms() const {
    return chrom_binneds;
}

void BWBinner::save_binneds(const std::string& out_dir) const {
    std::filesystem::path out_dir_p{out_dir};
    if (!std::filesystem::exists(out_dir_p)) {
        std::cout << "Creating directory " << out_dir_p << "\n";
        std::filesystem::create_directory(out_dir_p);
    }
    
    std::filesystem::path chrom_names_path{out_dir_p / "chrom_names.txt"};
    std::ofstream chrom_names_file(chrom_names_path);
    for (const auto& [chrom, size] : chrom_sizes) {
        // add this chrom name to file
        chrom_names_file << chrom << "\n";
        // save this chrom's binned tensor
        auto bytes = torch::pickle_save(chrom_binneds.at(chrom));
        std::ofstream chr_stream{out_dir_p / (chrom + ".pt")};
        chr_stream.write(bytes.data(), bytes.size());
        chr_stream.close();
    }
    chrom_names_file.close();
    // torch::save(chrom_binneds_vec, out_dir_p/"binned_tracks_all_chroms.pt");
}