#include "checksum_utils.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream> // For ostringstream to format hex
#include <iomanip> // For std::setw, std::setfill
#include <digestpp/digestpp.hpp> // For MD5
#include <digestpp/algorithm/md5.hpp>

#include "../../../third_party/digestpp/algorithm/md5.hpp"


namespace EOSCheckSumUtils {

    std::string calculate_file_md5(const std::string& filepath) {
        std::ifstream file(filepath, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error (calculate_file_md5): Could not open file: " << filepath << std::endl;
            return "";
        }
        std::ostringstream oss;
        oss << digestpp::md5().absorb(file).hexdigest();
        file.close();
        return oss.str();
    }

} // namespace EOSUtils
