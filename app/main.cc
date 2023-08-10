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
        TCLAP::ValueArg<unsigned> res("r", "resolution", "resolution (bin size) in base pairs", false, 100, "unsigned (>= 0) int");
        // TCLAP::UnlabeledMultiArg<std::string> tracks("tracks", "bigWig files to bin", true, "name: str path: str");
        TCLAP::ValueArg<std::string> tracks_dir("T", "tracks-dir", "directory containing all bigWig files to bin", true, "", "path (string)");
        TCLAP::ValueArg<std::string> chrom_sizes("s", "chrom-sizes", "chromosome sizes file", true, "", "path (string)");
        TCLAP::UnlabeledValueArg<std::string> out_dir("out-dir", "directory to write binned tensors to", true, "", "path (string)");
        TCLAP::SwitchArg verbose("v", "verbose", "print verbose output", false);
        cmd.add(tracks_dir);
        cmd.add(chrom_sizes);
        cmd.add(res);
        cmd.add(out_dir);
        cmd.add(verbose);
        cmd.parse(argc, argv);

        if (verbose.getValue()) {
            std::cout << "resolution: " << res.getValue() << std::endl;
            std::cout << "tracks dir: " << tracks_dir.getValue() << std::endl;
            std::cout << "chrom sizes: " << chrom_sizes.getValue() << std::endl;
            std::cout << "out dir: " << out_dir.getValue() << std::endl;

            std::cout << "Note: assuming all bigWig files' extension is '.bigWig'" << std::endl;
        }

        // For now, we'll just assume that the user has given us a directory
        std::vector<std::string> bw_paths = find_paths_filetype(tracks_dir.getValue(), ".bigWig");
        if (bw_paths.size() == 0) {
            std::cerr << "No bigWig files found in " << tracks_dir.getValue() << std::endl;
            return 1;
        }
        if (verbose.getValue()) {
            std::cout << "Found " << bw_paths.size() << " bigWig files in " << tracks_dir.getValue() <<":\n";
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