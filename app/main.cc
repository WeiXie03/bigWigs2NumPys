#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <tclap/CmdLine.h>
#include <bigWigs2tensors/util.h>
#include <bigWigs2tensors/proc_bigWigs.h>

/*
 * Usage *
tracks: EITHER...
    - unlimited -t <name: str> <path: str>'s, or
    - one -T <path: str> to a directory with all the bigWigs
chrom sizes: one -s <path: str>
*/

int main(int argc, char** argv) {

    try {
        TCLAP::CmdLine cmd("bigWigs binner that converts bigWigs to (pickled) PyTorch tensor files",
                        ' ', "0.1");
        TCLAP::ValueArg<unsigned> res("r", "resolution", "resolution (bin size) in base pairs", false, 100, "unsigned int");
        // TCLAP::UnlabeledMultiArg<std::string> tracks("tracks", "bigWig files to bin", true, "name: str path: str");
        TCLAP::ValueArg<std::string> tracks_dir("T", "tracks-dir", "directory containing all bigWig files to bin", true, "", "string");
        TCLAP::ValueArg<std::string> chrom_sizes("s", "chrom-sizes", "chromosome sizes file", true, "", "path: str");
        TCLAP::UnlabeledValueArg<std::string> out_dir("out-dir", "directory to write binned tensors to", true, "", "string");
        cmd.add(tracks_dir);
        cmd.add(chrom_sizes);
        cmd.parse(argc, argv);

        // For now, we'll just assume that the user has given us a directory
        std::vector<std::string> bw_paths = find_paths_filetype(tracks_dir.getValue(), ".bw");

        // std::vector<std::string> bw_names;
        // for (auto& track : tracks.getValue()) {
        //     auto [name, path] = split(track, ':');
        //     bw_paths.push_back(path);
        //     bw_names.push_back(name);
        // }

        std::string chrom_sizes_path = chrom_sizes.getValue();

        BWBinner bwb(bw_paths, chrom_sizes_path);

        std::cout << "Binning bigWigs..." << std::endl;
        auto binned = bwb.load_bin_all_chroms(res.getValue());
        std::cout << "Done binning bigWigs: " << binned.size() << " chromosomes" << std::endl;
        for (auto& [chrom, tensor] : binned) {
            std::cout << chrom << ": " << tensor.sizes() << std::endl;
        }

        std::string save_path = out_dir.getValue();
        std::cout << "Writing tensors to "<< save_path << "..." << std::endl;
        bwb.save_binneds(save_path);
        std::cout << "Done writing tensors to disk." << std::endl;
    }
    catch (TCLAP::ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}