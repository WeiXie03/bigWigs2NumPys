#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <doctest/doctest.h>
#include <bigWigs2tensors/util.h>
#include <bigWigs2tensors/proc_bigWigs.h>

const std::filesystem::path DATA_DIR = std::filesystem::current_path() / "data";

TEST_CASE("open_bigWigs") {
    std::vector<std::string> bw_paths = find_paths_filetype(DATA_DIR, ".bw");

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
    std::filesystem::path chrom_sizes_path = DATA_DIR / "toy.chrom.sizes";

    // check chrom sizes parsing
    //for (const auto& [chr, size] : parse_chrom_sizes(chrom_sizes_path.string())) {
    //    std::cout << chr << ": " << size << std::endl;
    //}

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
        chr3: 0.5 0.5 0.5 0
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
        // check(binned_chroms["chr3"].equal( torch::tensor({0.5, 0.5, 0.5, 0}, constants::tensor_opts) ));

        binner.save_binneds(bwFstem+"_out");
    }
}

TEST_CASE("sequential with missings toy data") {
    std::string bwFstem {"test_sequential_missing"};
    std::vector<std::string> bw_paths ({ (DATA_DIR / (bwFstem + ".bw")).string() });
    std::filesystem::path chrom_sizes_path = DATA_DIR / "toy.chrom.sizes";

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
        chr3: - - 1.5 -
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

TEST_CASE("combined toy datas") {
    std::vector<std::string> bw_paths = find_paths_filetype(DATA_DIR, ".bw");
    std::filesystem::path chrom_sizes_path = DATA_DIR / ("toy.chrom.sizes");

    BWBinner binner(bw_paths, chrom_sizes_path.string());
    REQUIRE(binner.binned_chroms().size() == 0);

    binner.load_bin_all_chroms(2);
    std::map<std::string, torch::Tensor> binned_chroms = binner.binned_chroms();
    // chr1
    std::cout << "binned chr1:\n" << binned_chroms["chr1"] << std::endl;
    // chr2
    std::cout << "binned chr2:\n" << binned_chroms["chr2"] << std::endl;
    // chr3
    std::cout << "binned chr3:\n" << binned_chroms["chr3"] << std::endl;

    binner.save_binneds("combined_out");
}

TEST_CASE("parse bigBed for specifying coordinates") {
    std::cout << "\n========================\n";

    std::filesystem::path coords_bed_path = DATA_DIR / ("toy.coords.bigBed");
    std::cout << "coords_bed_path: " << coords_bed_path << std::endl;

    std::map<std::string,int> seq_miss_chrom_sizes = parse_chrom_sizes((DATA_DIR / ("toy.chrom.sizes")).string());
    std::map<std::string, bbOverlappingEntries_t*> spec_coords = parse_coords_bigBed(coords_bed_path.string(), seq_miss_chrom_sizes);
    for (auto& [chr, interv] : spec_coords) {
        std::cout << chr <<": "<< interv->l <<" intervals"<< std::endl;
        std::cout << '\t';
        for (int i = 0; i < interv->l; i++) {
            std::cout << ' ' << interv->start[i] << '-' << interv->end[i];
        };
        std::cout << std::endl;
        if (interv) bbDestroyOverlappingEntries(interv);
    }
    bwCleanup();
}

TEST_CASE("specify coords on sequential missing toy data") {
    std::vector<std::string> bw_paths = find_paths_filetype(DATA_DIR, ".bw");
    std::map<std::string,int> seq_miss_chrom_sizes = parse_chrom_sizes((DATA_DIR / ("toy.chrom.sizes")).string());
    std::filesystem::path chrom_sizes_path = DATA_DIR / ("toy.chrom.sizes");
    std::filesystem::path coords_bed_path = DATA_DIR / ("toy.coords.bigBed");

    BWBinner binner(bw_paths, chrom_sizes_path.string(), parse_coords_bigBed(coords_bed_path.string(), seq_miss_chrom_sizes));
    REQUIRE(binner.binned_chroms().size() == 0);

    binner.load_bin_all_chroms(2);
    std::map<std::string, torch::Tensor> binned_chroms = binner.binned_chroms();

    // chr1
    std::cout << "binned chr1:\n" << binned_chroms["chr1"] << std::endl;
    // CHECK(binned_chroms["chr1"].equal( torch::full({5}, 0.5, constants::tensor_opts) ));
    // chr2
    std::cout << "binned chr2:\n" << binned_chroms["chr2"] << std::endl;
    // CHECK(binned_chroms["chr2"].equal( torch::full({2}, 0.5, constants::tensor_opts) ));
    // chr3
    std::cout << "binned chr3:\n" << binned_chroms["chr3"] << std::endl;

    binner.save_binneds("subset_seq_miss_out");

    bwCleanup();
}