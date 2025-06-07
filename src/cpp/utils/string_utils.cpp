#include "string_utils.h"
#include <algorithm>     // For std::replace, std::transform (though transform not used in trim here)
#include <stdexcept>     // For std::stod exceptions
#include <iostream>      // For potential error messages in string_to_double

namespace EOSUtils {

    std::string trim_string(const std::string& str) {
        const std::string whitespace = " \t\n\r\f\v";
        size_t start = str.find_first_not_of(whitespace);
        if (start == std::string::npos) {
            return ""; // String is empty or all whitespace
        }
        size_t end = str.find_last_not_of(whitespace);
        return str.substr(start, end - start + 1);
    }

    bool string_to_double_fortran_compat(const std::string& s, double& out_val) {
        if (s.empty()) {
            return false;
        }
        std::string temp_s = s;
        // Replace 'D' or 'd' with 'E' for exponent
        std::replace(temp_s.begin(), temp_s.end(), 'D', 'E');
        std::replace(temp_s.begin(), temp_s.end(), 'd', 'E');
        try {
            out_val = std::stod(temp_s);
            return true;
        } catch (const std::invalid_argument& ia) {
            return false;
        } catch (const std::out_of_range& oor) {
            return false;
        }
        // Catch any other std::exception as a fallback, though stod usually throws the above two.
        catch (const std::exception& e) {
            return false;
        }
    }

} // namespace EOSUtils