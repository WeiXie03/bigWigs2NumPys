#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <tclap/CmdLine.h>
#include <tclap/Arg.h>
#include <bigWigs2tensors/util.h>
#include <bigWigs2tensors/proc_bigWigs.h>

/*!
Merge all given paths and matching paths within given directories into a single vector and return it.
tracks_list is a list of paths to bigWig files and/or directories containing bigWig files to bin over,
directly from the TCLAP interface.
*/
std::vector<std::string> parse_tracks_list(const std::vector<std::string>& tracks_list, const std::string& file_ext) {
    std::vector<std::string> paths;
    for (auto& path : tracks_list) {
        if (std::filesystem::is_directory(path)) {
            std::vector<std::string> dir_paths = find_paths_filetype(path, file_ext);
            paths.insert(paths.end(), dir_paths.begin(), dir_paths.end());
        }
        else if (std::filesystem::is_regular_file(path)) {
            paths.push_back(path);
        }
        else std::cerr << "Invalid track path specified: " << path << std::endl;
    }
    return paths;
}

/*
 * Usage *
 bigWigs2tensors -r <resolution: unsigned int> [-c <coords BED: str>] -s <chrom sizes> -o <out-dir: str> -t <track: str> [-t <track: str> ...]
    * Options *
    track: a path to one of the bigWig files and/or directories containing bigWig files to bin over
    chrom sizes: one -s <path: str>
    _OPTIONAL_
    coords BED: one -c <path: str> to a bed file specifying the genomic coordinates _within __every__ bigWig_ to bin over
*/

int main(int argc, char** argv) {
    try {
        TCLAP::CmdLine cmd("bigWigs binner that converts bigWigs to (pickled) PyTorch tensor files",
                        ' ', "0.1");
        TCLAP::ValueArg<unsigned> res("r", "resolution", "resolution (bin size) in base pairs", false, 100, "unsigned (>= 0) int", cmd);
        // TCLAP::UnlabeledMultiArg<std::string> tracks("tracks", "bigWig files to bin", true, "name: str path: str", cmd);
        TCLAP::MultiArg<std::string> tracks_list("t", "tracks-list", "a list of paths of bigWig files and/or directories containing bigWig files to bin over", true, "path (string)", cmd);
        TCLAP::ValueArg<std::string> chrom_sizes("s", "chrom-sizes", "chromosome sizes file", true, "", "path (string)", cmd);
        TCLAP::ValueArg<std::string> coords_bed("c", "coords-bigBed", "bigBed file of genomic coordinates within all bigWigs to bin over", false, "", "path (string)", cmd);
        TCLAP::UnlabeledValueArg<std::string> out_dir("out-dir", "directory to write binned tensors to", true, "", "path (string)", cmd);
        TCLAP::SwitchArg verbose("v" , "verbose", "print verbose output", cmd, false);
        cmd.parse(argc, argv);

        if (verbose.getValue()) {
            std::cout << "resolution: " << res.getValue() << std::endl;
            std::cout << "chrom sizes: " << chrom_sizes.getValue() << std::endl;
            std::cout << "coords bed: " << coords_bed.getValue() << std::endl;
            std::cout << "out dir: " << out_dir.getValue() << std::endl;

            std::cout << tracks_list.getValue().size() << " track arguments given." << std::endl;
            std::cout << "tracks list: [\n";
            for (auto& track : tracks_list.getValue()) {
                std::cout << '\t' << track << '\n';
            }
            std::cout << "]" << std::endl;

            std::cout << "Note: assuming all bigWig files' extension is '.bigWig'" << std::endl;
        }

        // TODO: all given paths in tracks_list
        std::vector<std::string> bw_paths = parse_tracks_list(tracks_list.getValue(), ".bigWig");
        if (bw_paths.size() == 0) {
            std::cerr << "No bigWig files found in given tracks list." << std::endl;
            return 1;
        }
        if (verbose.getValue()) {
            std::cout << "Found " << bw_paths.size() << " bigWig files:\n";
            for (auto& path : bw_paths) {
                std::cout << '\t' << path << std::endl;
            }
        }

        // std::vector<std::string> bw_names;
        // for (auto& track : tracks.getValue()) {
        //     auto [name, path] = split(track, ':');
        //     bw_paths.push_back(path);
        //     bw_names.push_back(name);
        // }

        //Initialize enough space to hold 128KiB (1<<17) of data at a time
        if(bwInit(1<<17) != 0) {
            fprintf(stderr, "Received an error in bwInit\n");
            return 1;
        }

        chroms_coords_map_t coords_map;
        std::string chrom_sizes_path = chrom_sizes.getValue();
        std::map<std::string,int> chr_sizes_map = parse_chrom_sizes(chrom_sizes_path);

        BWBinner* bwb = nullptr;
        if (coords_bed.isSet()) {
            std::cout << "Using specified coordinates bigBed..." << std::endl;
            bwb = new BWBinner(bw_paths, chrom_sizes_path, coords_bed.getValue());
            // std::cout << "Parsing coordinates bigBed..." << std::endl;
            // coords_map = parse_coords_bigBed(coords_bed.getValue(), chr_sizes_map);
        }
        else {
            // coords_map = make_full_chroms_coords_map(chr_sizes_map);
            bwb = new BWBinner(bw_paths, chrom_sizes_path);
        }

        std::cout << "Binning bigWigs..." << std::endl;
        auto binned = bwb->load_bin_all_chroms(res.getValue());
        std::cout << "Done binning bigWigs: " << binned.size() << " chromosomes" << std::endl;
        for (auto& [chrom, tensor] : binned) {
            std::cout << chrom << ": " << tensor.sizes() << std::endl;
        }

        std::string save_path = out_dir.getValue();
        std::cout << "Writing tensors to "<< save_path << "..." << std::endl;
        bwb->save_binneds(save_path);
        std::cout << "Done writing tensors to disk." << std::endl;

        delete bwb;
    }
    catch (TCLAP::ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}