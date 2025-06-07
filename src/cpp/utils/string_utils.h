#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <vector> // Include if you plan to add more utils that might need it

namespace EOSUtils {

    // Trims whitespace from both ends of a string.
    // Returns the trimmed string.
    std::string trim_string(const std::string& str);

    // Converts a string potentially containing Fortran 'D' or 'd' exponent
    // characters to a double value.
    // Returns true on success, false on failure (e.g., invalid format, out of range).
    // The converted double is stored in out_val.
    bool string_to_double_fortran_compat(const std::string& s, double& out_val);

    // You can add more string or general utilities here later.

} // namespace EOSUtils

#endif // STRING_UTILS_H