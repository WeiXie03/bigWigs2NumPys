#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <execution>
#include <cmath>
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

// Note: half-open interval-based, i.e. [start, end)
double NaNmean_range(std::vector<double>::iterator start, std::vector<double>::iterator end) {
    // If any NaNs in bin, set result to NaN
    if (std::any_of(start, end, [](double val) { return std::isnan(val); }))
        return std::nan("");

    return std::accumulate(start, end, 0.0) / std::distance(start, end);
}

std::vector<double> bin_vec_NaNmeans(std::vector<double>& in_vec, size_t bin_size) {
    // Need to handle last bin separately if it's not full,
    // so just go up to last full bin first
    size_t n_bins = ceil(double(in_vec.size()) / double(bin_size));
    std::vector<double> means(n_bins - 1);
    std::vector<size_t> bindx (n_bins - 1);
    std::iota(bindx.begin(), bindx.end(), 0);

    // Calculate mean means in parallel using C++17 parallel algorithms
    std::for_each(std::execution::par_unseq,
                    bindx.begin(), bindx.end(),
                    [bin_size, &in_vec, &means](size_t i_bin) {
                        std::vector<double>::iterator bin_start = in_vec.begin() + i_bin*bin_size;
                        std::vector<double>::iterator bin_end = bin_start + bin_size;
                        means[i_bin] = NaNmean_range(bin_start, bin_end);
                    });
    // std::transform(std::execution::par_unseq, bin_starts.begin(), bin_starts.begin() + n_bins-1 - 1,
    //                std::back_inserter(means), [bin_size, &in_vec](size_t indx) -> double {
    //                     std::vector<double>::iterator bin_start = in_vec.begin() + indx;
    //                     std::vector<double>::iterator bin_end = in_vec.begin() + indx + bin_size;
    //                     return NaNmean_range(bin_start, bin_end);
    //                });

    // Last bin
    auto bin_start = in_vec.begin() + (n_bins-1)*bin_size;
    auto bin_end = in_vec.end();
    means.push_back(NaNmean_range(bin_start, bin_end));

    return means;
}

torch::Tensor BWBinner::load_bin_chrom_tensor(const std::string& chrom, unsigned bin_size) {
    // ceiling division: ceil(chrom_size / bin_size)
    // credit: https://stackoverflow.com/a/2745086
    unsigned num_bins = 1 + ((chrom_sizes[chrom] - 1) / bin_size);
    std::cout << "Binning " << chrom << " into " << num_bins << " bins" << std::endl;

    torch::Tensor chrom_tensor = torch::empty({num_bins, num_bws}, tens_opts);
    std::cout << "Created " << chrom_tensor.sizes() << " tensor for "<< num_bws <<" tracks." << std::endl;
    
    // parallelize across bigWigs' indices within the bw_files vector
    // credit: https://stackoverflow.com/a/62829166
    std::vector<size_t> bw_idxs(bw_files.size());
    std::iota(bw_idxs.begin(), bw_idxs.end(), 0);

    std::cout << "bw_idxs: [";
    for (auto& idx : bw_idxs) {
        std::cout << idx << ", ";
    }
    std::cout << "]" << std::endl;

    std::for_each(std::execution::par_unseq,
                    bw_idxs.begin(), bw_idxs.end(),
                    [this, chrom, bin_size, num_bins, &chrom_tensor](size_t bw_idx) {
                        // the whole chrom for this bigWig, including missing vals
                        double* vals_arr = bwStats(bw_files[bw_idx], chrom.c_str(),
                                                        0, chrom_sizes[chrom], chrom_sizes[chrom],
                                                        bwStatsType::doesNotExist);
                        std::vector<double> chrom_vals(vals_arr, vals_arr + chrom_sizes[chrom]);

                        // std::cout << "Loaded " << chrom_vals.size() << " values for " << bw_paths[bw_idx].stem().string() << ": [";
                        // for (double val : chrom_vals) {
                        //     std::cout << " " << val;
                        // }
                        // std::cout << " ]" << std::endl;

                        std::vector<double> binned_vals = bin_vec_NaNmeans(chrom_vals, bin_size);
                        // std::cout << "Binned result: [";
                        // for (double val : binned_vals) {
                        //     std::cout << " " << val;
                        // }
                        // std::cout << "... ]"  << std::endl;

                        using namespace torch::indexing;
                        // each bigWig is a column in the tensor
                        // set bw_idx'th column of chrom_tensor to binned_vals
                        chrom_tensor.index_put_({"...", (int)bw_idx},
                                                torch::from_blob(binned_vals.data(), {num_bins},
                                                                torch::dtype(torch::kFloat64)));
                    });

    return chrom_tensor;
}

std::map<std::string, torch::Tensor> BWBinner::load_bin_all_chroms(unsigned bin_size) {
    // parallelize over chromosomes
    std::for_each(std::execution::par_unseq,
        chrom_sizes.begin(), chrom_sizes.end(),
        [this, bin_size](auto& chr_entry) {
            // std::cout << "Binning " << chr_entry.first << std::endl;
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