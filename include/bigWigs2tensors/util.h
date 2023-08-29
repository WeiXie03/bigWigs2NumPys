#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>

/*!
Utility function to find all files of a given type in a directory.
_Note_: __Not__ recursive.
\param search_dir The path of the directory to search.
\param file_ext The file extension to search for, including '.'.
\return A vector of the paths found.
*/
std::vector<std::string> find_paths_filetype(const std::string& search_dir, const std::string& file_ext);

#endif