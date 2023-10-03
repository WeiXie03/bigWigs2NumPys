#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <bigWig.h>

/*!
Utility function to find all files of a given type in a directory.
_Note_: __Not__ recursive.
\param search_dir The path of the directory to search.
\param file_ext The file extension to search for, including '.'.
\return A vector of the paths found.
*/
std::vector<std::string> find_paths_filetype(const std::string& search_dir, const std::string& file_ext);

/*!
Reads a bigBed file into a map of chromosome names to their overlapping entries.
*/
void parse_coords_bed(std::map<std::string, bbOverlappingEntries_t*>& coords_map, const std::string& coords_bed_path, const std::map<std::string, int>& chrom_sizes);

#endif