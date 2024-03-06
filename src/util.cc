#include <vector>
#include <map>
#include <execution>
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <bigWig.h>
#include <bigWigs2tensors/util.h>

std::vector<std::string> find_paths_filetype(const std::string& search_dir, const std::string& file_ext) {
    // collect all the bigWig files in the data directory
    std::vector<std::string> match_paths;
    std::filesystem::path search_dir_p{search_dir};

    for (auto const& entry : std::filesystem::directory_iterator(search_dir_p)) {
        if (entry.path().extension() == file_ext)
            match_paths.push_back(entry.path().string());
    }
    return match_paths;
}

std::map<std::string, int> parse_chrom_sizes(const std::string& chrom_sizes_path) {
    std::map<std::string, int> chrom_sizes;
    std::ifstream chrom_sizes_file(chrom_sizes_path);
    if (!chrom_sizes_file.is_open()) {
        std::cerr << "ERROR: could not open " << chrom_sizes_path << std::endl;
        exit(1);
    }
    std::string chrom;
    int size;
    while (chrom_sizes_file >> chrom >> size) {
        chrom_sizes[chrom] = size;
    }
    chrom_sizes_file.close();
    return chrom_sizes;
}

chroms_coords_map_t parse_coords_bigBed(const std::string& coords_bed_path, const std::map<std::string, int>& chrom_sizes) {
    bigWigFile_t* coords_bed = bbOpen(const_cast<char*>(coords_bed_path.c_str()), NULL);
    if (bbIsBigBed(const_cast<char*>(coords_bed_path.c_str()), NULL) != 1) {
        std::cerr << "ERROR: " << coords_bed_path << " is not a bigBed file" << std::endl;
        exit(1);
    }

    std::map<std::string, bbOverlappingEntries_t*> chroms_coords;
    // get a bbOverlappingEntry from the coordinates specification bigBed for each chromosome
    std::for_each(std::execution::par_unseq,
                    chrom_sizes.begin(), chrom_sizes.end(),
                    [&coords_bed, &chroms_coords](const auto& chr_entry) {
                        std::string chrom = chr_entry.first;
                        int chrom_size = chr_entry.second;
                        // REMEMBER to bwDestroyOverlappingIntervals() to free
                        bbOverlappingEntries_t* interv = bbGetOverlappingEntries(coords_bed, const_cast<char*>(chrom.c_str()), 0, chrom_size, 0);

                        // if (!interv) {
                        //     // std::cout << "Warning: no intervals found for " << chrom << " in specified coordinates BED file, skipping..." << std::endl;
                        // }
                        // else {
                        if (interv) {
                            std::cout << interv->l <<" intervals for "<< chrom <<": {";
                            std::cout << ' ';

                            for (int i = 0; i < interv->l; i++) {
                                std::cout << '[' << interv->start[i] << ',' << interv->end[i] << ')';
                                std::cout << ' ';
                            };
                            std::cout << "}" << std::endl;

                            chroms_coords.emplace(chrom, interv);
                        }
                    });
    /*
    for (const auto& chr_entry : chrom_sizes) {
        std::string chrom = chr_entry.first;
        int chrom_size = chr_entry.second;
        // REMEMBER to bwDestroyOverlappingIntervals() to free
        bbOverlappingEntries_t* interv = bbGetOverlappingEntries(coords_bed, const_cast<char*>(chrom.c_str()), 0, chrom_size, 0);
        chroms_coords.emplace(chrom, interv);

        std::cout << interv->l <<" intervals for "<< chrom <<": {";
        std::cout << ' ';

        for (int i = 0; i < interv->l; i++) {
            std::cout << '[' << interv->start[i] << ',' << interv->end[i] << ')';
            std::cout << ' ';
        };
        std::cout << "}" << std::endl;
    }
    */

    bwClose(coords_bed);

    return chroms_coords;
}

std::pair< std::map<std::string,int>, chroms_coords_map_t > parse_chrom_sizes_coords(const std::string& chrom_sizes_path, const std::string& coords_bed_path) {
    std::map<std::string, int> chrom_sizes;
    std::map<std::string, bbOverlappingEntries_t*> chroms_coords;

    std::ifstream chrom_sizes_file(chrom_sizes_path);
    if (!chrom_sizes_file.is_open()) {
        std::cerr << "ERROR: could not open " << chrom_sizes_path << std::endl;
        exit(1);
    }
    bigWigFile_t* coords_bed = bbOpen(const_cast<char*>(coords_bed_path.c_str()), NULL);
    if (bbIsBigBed(const_cast<char*>(coords_bed_path.c_str()), NULL) != 1) {
        std::cerr << "ERROR: " << coords_bed_path << " is not a bigBed file" << std::endl;
        exit(1);
    }

    std::string chrom;
    int size;
    bbOverlappingEntries_t* interv;
    // for each chromsome in the chrom_sizes file, get the intervals from the bigBed file
    while (chrom_sizes_file >> chrom >> size) {
        bbOverlappingEntries_t* interv = bbGetOverlappingEntries(coords_bed, const_cast<char*>(chrom.c_str()), 0, size, 0);
        if (!interv) {
            std::cout << "Warning: no intervals found for " << chrom << " in specified coordinates BED file, skipping..." << std::endl;
        }
        else {
            // only add to the map if chrom specified in coords_bed
            chrom_sizes[chrom] = size;
            chroms_coords.emplace(chrom, interv);
        }
    }

    chrom_sizes_file.close();
    bwClose(coords_bed);

    return std::make_pair(chrom_sizes, chroms_coords);
}

unsigned num_bins_intersect_interval(unsigned start, unsigned end, unsigned bin_size) {
    return std::floor(end / static_cast<float>(bin_size)) - std::ceil(start / static_cast<float>(bin_size));
}

bwOverlappingIntervals_t* bin_coords(const bbOverlappingEntries_t* coords, unsigned bin_size) {
    bwOverlappingIntervals_t* binned_coords = new bwOverlappingIntervals_t;
    binned_coords->l = coords->l;
    binned_coords->start = new uint32_t[coords->l];
    binned_coords->end = new uint32_t[coords->l];
    binned_coords->value = new float[coords->l];

    std::transform(coords->start, coords->start + coords->l, binned_coords->start,
                   [bin_size](uint32_t start) { return std::ceil(start / static_cast<float>(bin_size)); });

    std::transform(coords->end, coords->end + coords->l, binned_coords->end,
                   [bin_size](uint32_t end) { return std::floor(end / static_cast<float>(bin_size)); });

    std::transform(binned_coords->start, binned_coords->start + coords->l, binned_coords->end, binned_coords->value,
                   [](uint32_t start, uint32_t end) { return end - start; });

    return binned_coords;
}

chroms_coords_map_t make_full_chroms_coords_map(const std::map<std::string, int>& chrom_sizes) {
    std::map<std::string, bbOverlappingEntries_t*> chroms_coords;

    std::transform(chrom_sizes.begin(), chrom_sizes.end(),
                    std::inserter(chroms_coords, chroms_coords.end()),
                    [](const std::pair<std::string, int>& chrom_size) {
                        std::string chrom = chrom_size.first;
                        unsigned size = chrom_size.second;

                        bbOverlappingEntries_t* entries = new bbOverlappingEntries_t;

                        // Will use entire chromosomes <=> "full" coordinates <=> 0 to chrom_size
                        uint32_t* starts = new uint32_t[1];
                        starts[0] = 0;
                        entries->start = starts;

                        uint32_t* ends = new uint32_t[1];
                        // 0-based half open
                        ends[0] = size;
                        entries->end = ends;

                        entries->l = 1;
                        entries->m = 1;
                        entries->str = NULL;

                        return std::make_pair(chrom, entries);
                    });
    
    return chroms_coords;
}