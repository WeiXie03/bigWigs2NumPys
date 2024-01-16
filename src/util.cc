#include <vector>
#include <map>
#include <execution>
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
    std::string chrom;
    int size;
    while (chrom_sizes_file >> chrom >> size) {
        chrom_sizes[chrom] = size;
    }
    return chrom_sizes;
}

chroms_coords_map_t parse_coords_bigBed(const std::string& coords_bed_path, const std::map<std::string, int>& chrom_sizes) {
    bigWigFile_t* coords_bed = bbOpen(coords_bed_path.c_str(), NULL);

    std::map<std::string, bbOverlappingEntries_t*> chroms_coords;
    // get a bbOverlappingEntry from the coordinates specification bigBed for each chromosome
    std::transform(chrom_sizes.begin(), chrom_sizes.end(),
                    std::inserter(chroms_coords, chroms_coords.end()),
                    [&coords_bed](const auto& chr_entry) {
                        std::string chrom = chr_entry.first;
                        int chrom_size = chr_entry.second;
                        // REMEMBER to bwDestroyOverlappingIntervals() to free
                        return std::make_pair(chrom, bbGetOverlappingEntries(coords_bed, chrom.c_str(), 0, chrom_size, 0));
                    });

    return chroms_coords;
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
                        ends[0] = size;
                        entries->end = ends;

                        entries->l = 1;
                        entries->m = 1;
                        entries->str = NULL;

                        return std::make_pair(chrom, entries);
                    });
    
    return chroms_coords;
}