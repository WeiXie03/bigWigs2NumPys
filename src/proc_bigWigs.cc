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

// Note: half-open 0-based, i.e. [start, end)
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

BWBinner::BWBinner(const std::vector<std::string>& bigWig_paths, const std::string& chrom_sizes_path, const std::map<std::string, bbOverlappingEntries_t*>& coords_map)
    : bw_files(open_bigWigs(bigWig_paths)),
    num_bws(bw_files.size()),
    chrom_sizes(parse_chrom_sizes(chrom_sizes_path)),
    tens_opts(constants::tensor_opts),
    spec_coords(coords_map)
{
    // std::cout << "\n=========================\n";

    // Parse the coordinates bed file
    //spec_coords = parse_coords_bigBed(coords_bed_path, chrom_sizes);
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

BWBinner::BWBinner(const std::vector<std::string>& bigWig_paths, const std::string& chrom_sizes_path)
    : BWBinner(bigWig_paths, chrom_sizes_path,
                make_full_chroms_coords_map(parse_chrom_sizes(chrom_sizes_path))) {}

BWBinner::~BWBinner() {
    // std::cout << "BWBinner shutting down" << std::endl;
    for (auto& bw : bw_files) {
        bwClose(bw);
    }
    // coordinates specification map
    for (auto& [chr, interv] : spec_coords) {
        bbDestroyOverlappingEntries(interv);
    }

    bwCleanup();
    // final Torch tensors
    chrom_binneds.clear();
}

void BWBinner::load_bin_chrom_bigWig_tensor(const std::string& chrom, size_t bw_idx, const std::vector<unsigned>& start_bindxs, unsigned bin_size, unsigned num_bins) {
    // check that the tensor for the chrom was created
    if (!chrom_binneds.contains(chrom)) {
        throw std::invalid_argument("BWBinner::load_bin_chrom_bigWig_tensor: no tensor for chrom " + chrom);
    }

    bbOverlappingEntries_t* chrom_coords = spec_coords[chrom];
    // parallelize across intervals' indices within the spec_coords map
    // credit: https://stackoverflow.com/a/62829166
    std::vector<size_t> interv_idxs (chrom_coords->l);
    std::iota(interv_idxs.begin(), interv_idxs.end(), 0);

    // TODO: check all values in spec_coords with debugger
    std::for_each(std::execution::par_unseq,
                    interv_idxs.begin(), interv_idxs.end(),
                    [this, chrom, bw_idx, &chrom_coords, &start_bindxs, bin_size, num_bins](size_t interv_idx) {
                        // each bigWig is a column in the tensor
                        // set interv_idx'th column of chrom_tensor to binned_vals
                        unsigned start_bin = start_bindxs[interv_idx];
                        unsigned end_bin;
                        // 0-based half-open
                        if (interv_idx != chrom_coords->l - 1)
                            end_bin = start_bindxs[interv_idx+1];
                        else // last interval
                            end_bin = num_bins;

                        if (end_bin > start_bin) {
                            unsigned start = chrom_coords->start[interv_idx];
                            unsigned end = chrom_coords->end[interv_idx];
                            // libBigWig, including chrom_coords, uses 0-based half-open intervals
                            unsigned interv_len = end - start;
                            double* vals_arr = bwStats(bw_files[bw_idx], chrom.c_str(),
                                                            start, end, interv_len,
                                                            bwStatsType::doesNotExist);
                            std::vector<double> chrom_vals(vals_arr, vals_arr + interv_len);

                            // std::cout << "Loaded " << chrom_vals.size() << " values for " << bw_paths[interv_idx].stem().string()// << ": [";
                            // for (double val : chrom_vals) {
                            //     std::cout << " " << val;
                            // }
                            // std::cout << " ]" << std::endl;

                            std::vector<double> binned_vals = bin_vec_NaNmeans(chrom_vals, bin_size);

                            // std::cout << "Binned result: [";
                            // for (double val : binned_vals) {
                            //     std::cout << " " << val;
                            // }
                            // std::cout << "]" << std::endl;

                            using namespace torch::indexing;

                            std::vector<double> overlapping_binned_vals = std::vector<double>(binned_vals.begin(),
                                                                                            binned_vals.begin() + (end_bin - start_bin));
                            // std::cout << "interval "<< interv_idx <<": ["<< start_bin <<", "<< end_bin <<"), "<< end_bin - start_bin << " bins." << std::endl;

                            // 0-based half-open
                            chrom_binneds[chrom].index_put_({Slice(start_bin, end_bin, None), (int)bw_idx},
                                                    torch::from_blob(overlapping_binned_vals.data(), {num_bins},
                                                                    torch::dtype(torch::kFloat64)));
                        }
                        else {
                            // std::cout << "Interval "<< interv_idx << " is empty, skipping." << std::endl;
                            return;
                        }
                    });
    // return chrom_binneds[chrom];
}

void BWBinner::load_bin_chrom_tensor(const std::string& chrom, unsigned bin_size) {
    // calculate number of bins first
    // // ceiling division: ceil(chrom_size / bin_size)
    // // credit: https://stackoverflow.com/a/2745086
    // unsigned num_bins = 1 + ((chrom_sizes[chrom] - 1) / bin_size);

    unsigned num_intervs = spec_coords[chrom]->l;
    // cache starting indices of all intervals within tensor
    // std::cout << "Starting indices of intervals for " << chrom << ": [";
    std::vector<unsigned> start_bindxs(spec_coords[chrom]->l);
    start_bindxs[0] = 0;
    // std::cout << start_bindxs[0];
    for (int i = 1; i < num_intervs; i++) {
        // start of this interval is 1 + end of previous interval
        start_bindxs[i] = start_bindxs[i-1] + num_bins_intersect_interval(spec_coords[chrom]->start[i-1], spec_coords[chrom]->end[i-1], bin_size);
        // std::cout << ", " << start_bindxs[i];
    }
    // std::cout << "]" << std::endl;
    // total num bins just the last interval's end
    // remember, always 0-based [start, end)
    unsigned num_bins = start_bindxs.back() + num_bins_intersect_interval(spec_coords[chrom]->start[num_intervs-1], spec_coords[chrom]->end[num_intervs-1], bin_size);

    // use emplace to not copy into a temporary
    chrom_binneds.emplace(chrom, torch::empty({num_bins, num_bws}, tens_opts));
    std::cout << "Created " << chrom_binneds[chrom].sizes() << " tensor for "<< num_bws <<" tracks." << std::endl;
    
    // parallelize across bigWigs' indices within the bw_files vector
    // credit: https://stackoverflow.com/a/62829166
    std::vector<size_t> bw_idxs(bw_files.size());
    std::iota(bw_idxs.begin(), bw_idxs.end(), 0);

    // std::cout << "bw_idxs: [";
    // for (auto& idx : bw_idxs) {
    //     std::cout << idx << ", ";
    // }
    // std::cout << "]" << std::endl;

    std::for_each(std::execution::par_unseq,
                    bw_idxs.begin(), bw_idxs.end(),
                    [this, chrom, bin_size, num_bins, &start_bindxs](size_t bw_idx) {
                        load_bin_chrom_bigWig_tensor(chrom, bw_idx, start_bindxs, bin_size, num_bins);
                    });
    // return chrom_binneds[chrom];
}

const std::map<std::string, torch::Tensor>& BWBinner::load_bin_all_chroms(unsigned bin_size) {
    // parallelize over chromosomes
    std::for_each(std::execution::par_unseq,
        chrom_sizes.begin(), chrom_sizes.end(),
        [this, bin_size](auto& chr_entry) {
            // std::cout << "Binning " << chr_entry.first << std::endl;
            //chrom_binneds[chr_entry.first] = BWBinner::load_bin_chrom_tensor(chr_entry.first, bin_size);
            BWBinner::load_bin_chrom_tensor(chr_entry.first, bin_size);
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

    //std::cout << "in save_binneds(): chrom_binneds.size() = " << chrom_binneds.size() << std::endl;
    for (const auto& [chrom, size] : chrom_sizes) {
        // save this chrom's binned tensor
        //std::cout << "in save_binneds(): saving " << chrom << std::endl;
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