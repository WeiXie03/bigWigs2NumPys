#include <vector>
#include <map>
#include <algorithm>
#include <ranges>
#include <execution>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <torch/torch.h>
#include <bigWigs2tensors/proc_bigWigs.h>
#include <bigWig.h>

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

std::map<std::string, torch::Tensor> BWBinner::load_bin_all_chroms(unsigned bin_size) {
    // parallelize over chromosomes
    std::for_each(std::execution::par_unseq,
        chrom_sizes.begin(), chrom_sizes.end(),
        //std::pair<std::string, int>
        [this, bin_size](auto& chr_entry) {
            chrom_binneds[chr_entry.first] = BWBinner::load_bin_chrom_tensor(chr_entry.first, bin_size);
        }
    );
    return chrom_binneds;
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

    for (const auto& [chrom, size] : chrom_sizes) {
        // save this chrom's binned tensor
        auto bytes = torch::pickle_save(chrom_binneds.at(chrom));
        std::ofstream chr_stream{out_dir_p / (chrom + ".pt")};
        chr_stream.write(bytes.data(), bytes.size());
        chr_stream.close();
    }
    
    // save the indices of the bigWigs in the tensor
    std::filesystem::path bw_idx_p{out_dir_p / "tensor_bigWigs_inds.csv"};
    std::ofstream bw_idx_F(bw_idx_p);
    // headers
    bw_idx_F << "column" << ',' << "name" << '\n';
    for (size_t i = 0; i < bw_paths.size(); i++) {
        bw_idx_F << i << ',' << bw_paths[i].stem().string() << '\n';
    }
    bw_idx_F.close();
}