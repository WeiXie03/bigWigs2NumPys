#include <vector>
#include <map>
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

void parse_coords_bed(std::map<std::string, bbOverlappingEntries_t*>& coords_map, const std::string& coords_bed_path, const std::map<std::string, int>& chrom_sizes) {
    bigWigFile_t* coords_bed = bbOpen(coords_bed_path.c_str(), NULL);
    if (coords_bed == NULL) {
        std::cerr << "Error opening coordinates bed file " << coords_bed_path << std::endl;
        exit(1);
    }
    // Read the coordinates bed file into a map of {chrom: bbOverlappingEntries_t*}
    std::transform(chrom_sizes.begin(), chrom_sizes.end(),
                        std::inserter(coords_map, coords_map.end()),
                        [&coords_bed](const std::pair<std::string, int>& chrom_size) {
                            std::string chrom = chrom_size.first;
                            int size = chrom_size.second;
                            bbOverlappingEntries_t* entries = new bbOverlappingEntries_t;
                            entries = bbGetOverlappingEntries(coords_bed, chrom.c_str(), 0, size, 0);
                            return std::make_pair(chrom, entries);
                        });
}