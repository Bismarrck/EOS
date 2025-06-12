#ifndef EOS_CHECKSUM_UTILS_H
#define EOS_CHECKSUM_UTILS_H

#include <string>

namespace EOSCheckSumUtils {

    std::string calculate_file_md5(const std::string& filepath);

}

#endif