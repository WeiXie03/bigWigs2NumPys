#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "../libs/doctest/doctest.h"
#include "../include/bigWigs2tensors/proc_bigWigs.h"

const std::filesystem::path DATA_DIR = std::filesystem::current_path() / "data";

/*!
Utility function to find all files of a given type in a directory.
_Note_: __Not__ recursive.
\param search_dir The path of the directory to search.
\param file_ext The file extension to search for, including '.'.
\return A vector of the paths found.
*/
std::vector<std::string> find_paths_filetype(const std::string& search_dir, const std::string& file_ext) {
    // collect all the bigWig files in the data directory
    std::vector<std::string> match_paths;
    std::filesystem::path test_dir{search_dir};
    for (auto const& entry : std::filesystem::directory_iterator(test_dir)) {
        if (entry.path().extension() == file_ext)
            match_paths.push_back(entry.path().string());
    }
    return match_paths;
}

TEST_CASE("open_bigWigs") {
    std::filesystem::path test_dir = std::filesystem::current_path() / "data";
    std::vector<std::string> bw_paths = find_paths_filetype(test_dir, ".bw");

    std::cout << "bigWig files found: \n";
    for (auto const& path : bw_paths)
        std::cout <<"  "<< path << std::endl;

    std::vector<bigWigFile_t*> bw_files = open_bigWigs(bw_paths);
    CHECK(bw_files.size() == 2);
    for (int i = 0; i < bw_files.size(); i++) {
        CHECK(bw_files[i] != NULL && bw_files[i] != nullptr);
        bwClose(bw_files[i]);
    }
    bwCleanup();
}

TEST_CASE("mono alternate 0/1 toy data") {
    std::string bwFstem {"test_mono_alt0_1"};
    std::vector<std::string> bw_paths ({ (DATA_DIR / (bwFstem + ".bw")).string() });
    std::filesystem::path chrom_sizes_path = DATA_DIR / (bwFstem + ".chrom.sizes");

    BWBinner binner(bw_paths, chrom_sizes_path.string());
    REQUIRE(binner.binned_chroms().size() == 0);

    /* Data in
    =========================
    chr1: 0 1 0 1 0 1 0 1 0 1
    chr2: 1 0 1 0
    chr3: 0 1 0 1 0 1 0
    */

    SUBCASE("bin_size = 2") {
        /* expected binned (out)
        chr1: 0.5 0.5 0.5 0.5 0.5
        chr2: 0.5 0.5
        note: libbigwig weird, bins from end,
        leaving remainders at __beginning__
        chr3: 0 0.5 0.5 0.5
        */

        binner.load_bin_all_chroms(2);
        std::map<std::string, torch::Tensor> binned_chroms = binner.binned_chroms();

        // 3 chromosomes
        CHECK(binned_chroms.size() == 3);

        // chr1
        std::cout << "binned chr1:\n" << binned_chroms["chr1"] << std::endl;
        // check(binned_chroms["chr1"].equal( torch::full({5}, 0.5, constants::tensor_opts) ));
        // chr2
        std::cout << "binned chr2:\n" << binned_chroms["chr2"] << std::endl;
        // check(binned_chroms["chr2"].equal( torch::full({2}, 0.5, constants::tensor_opts) ));
        // chr3
        std::cout << "binned chr3:\n" << binned_chroms["chr3"] << std::endl;
        // check(binned_chroms["chr3"].equal( torch::tensor({0, 0.5, 0.5, 0.5}, constants::tensor_opts) ));

        binner.save_binneds(bwFstem+"_out");
    }
}

TEST_CASE("sequential with missings toy data") {
    std::string bwFstem {"test_sequential_missing"};
    std::vector<std::string> bw_paths ({ (DATA_DIR / (bwFstem + ".bw")).string() });
    std::filesystem::path chrom_sizes_path = DATA_DIR / (bwFstem + ".chrom.sizes");

    BWBinner binner(bw_paths, chrom_sizes_path.string());
    REQUIRE(binner.binned_chroms().size() == 0);

    /* Data in, '-' denotes missing
    =========================
    chr1: 0 - 0 1 2 3 - - 0 1
    chr2: 0 1 2 3
    chr3: - 0 - 0 1 2 -
    */

    SUBCASE("bin_size = 2") {
        /* expected binned (out)
        chr1: - 0.5 2.5 - 0.5
        chr2: 0.5 2.5
        chr3: - - 0.5 -
        */

        binner.load_bin_all_chroms(2);
        std::map<std::string, torch::Tensor> binned_chroms = binner.binned_chroms();

        // 3 chromosomes
        CHECK(binned_chroms.size() == 3);

        // chr1
        std::cout << "binned chr1:\n" << binned_chroms["chr1"] << std::endl;
        // CHECK(binned_chroms["chr1"].equal( torch::full({5}, 0.5, constants::tensor_opts) ));
        // chr2
        std::cout << "binned chr2:\n" << binned_chroms["chr2"] << std::endl;
        // CHECK(binned_chroms["chr2"].equal( torch::full({2}, 0.5, constants::tensor_opts) ));
        // chr3
        std::cout << "binned chr3:\n" << binned_chroms["chr3"] << std::endl;
        // CHECK(binned_chroms["chr3"].equal( torch::tensor({0, 0.5, 0.5, 0.5}, constants::tensor_opts) ));

        binner.save_binneds(bwFstem+"_out");
    }
}