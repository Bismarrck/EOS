#ifndef EOS_CHECKSUMS_H
#define EOS_CHECKSUMS_H

#include <map>
#include <string>

// This map will hold "relative_filepath_from_eos_data_dir" -> "expected_md5_hex_string"
// Populate this with actual MD5 sums of your official data files.
// Example:
const std::map<int, std::string> official_file_checksums = {
    {-1000, "829c62a5f6ffbcf753e900a5a991c54f"}, // Example: MD5 of an empty file
    {-2000, "8578a4e97c78f41d4e0d835c40f07b7a"},
    {0, "bde004687fb05daf07f8e26a9a3104f5"},
    {1, "bde004687fb05daf07f8e26a9a3104f5"},
    {10000, "fd66cfe62fdb2e50f12dd5fcf2585e86"}
};

#endif