#ifndef PARSER_UTILS_H
#define PARSER_UTILS_H

#include <fstream>
#include <string>
#include <vector>

namespace EOSParserUtils {

int parse_complicated_eos_params(std::ifstream& eos_file,
                                 const std::string& file_path_for_error,
                                 const std::vector<std::string>& param_order,
                                 std::vector<double>& out_flat_params);

}

#endif  // PARSER_UTILS_H
