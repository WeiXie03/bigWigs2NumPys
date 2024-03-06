#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <map>
#include <execution>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <bigWig.h>

// Map for storing any user-specified coordinates, one bbOverlappingEntries_t per chromosome
using chroms_coords_map_t = std::map<std::string, bbOverlappingEntries_t*>;

/*!
Utility function to find all files of a given type in a directory.
_Note_: __Not__ recursive.
\param search_dir The path of the directory to search.
\param file_ext The file extension to search for, including '.'.
\return A vector of the paths found.
*/
std::vector<std::string> find_paths_filetype(const std::string& search_dir, const std::string& file_ext);

/*!
Reads a whitespace-delimited chrom_sizes file into
a map of chromosome names to their sizes.
*/
std::map<std::string, int> parse_chrom_sizes(const std::string& chrom_sizes_path);

/*!
Reads a whitespace-delimited bed file of coordinates into a map of chromosome names to
references to all the intervals specified, within a bbOverlappingEntries_t per chromosome.
*/
chroms_coords_map_t parse_coords_bigBed(const std::string& coords_bed_path, const std::map<std::string, int>& chrom_sizes);

/*!
Combined function to parse chrom_sizes and coords_bed files. Extracts sizes only for chromosomes in coords_bed.
*/
std::pair< std::map<std::string,int>, chroms_coords_map_t > parse_chrom_sizes_coords(const std::string& chrom_sizes_path, const std::string& coords_bed_path);

/*!
Returns the number of bins fully covered by the interval [start, end) with given bin size.
*/
unsigned num_bins_intersect_interval(unsigned start, unsigned end, unsigned bin_size);

/*!
Generates a "binned copy" of the given bbOverlappingEntries_t, where each interval is the
starting bin and ending bin of the original interval after binning with given bin size.
The interval values are set to the number of bins in the interval.
*/
bwOverlappingIntervals_t* bin_coords(const bbOverlappingEntries_t* coords, unsigned bin_size);

/*
Returns a map of the same type as parse_coords_bigBed() for all chromosomes given in chrom_sizes
but with all intervals.
*/
chroms_coords_map_t make_full_chroms_coords_map(const std::map<std::string, int>& chrom_sizes);

#endif