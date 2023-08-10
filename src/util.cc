#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>
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