#include "string_utils.h"
#include <algorithm>     // For std::replace, std::transform (though transform not used in trim here)
#include <stdexcept>     // For std::stod exceptions
#include <iostream>      // For potential error messages in string_to_double
#include <sstream>
#include <iomanip>


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

    std::string simple_path_join(const std::string& p1, const std::string& p2) {
        char sep = '/'; // Common separator
        std::string tmp = p1;

#ifdef _WIN32 // Or other Windows checks
        sep = '\\';
#endif

        if (!p1.empty() && p1.back() != sep && !p2.empty() && p2.front() != sep) {
            tmp += sep;
        } else if (!p1.empty() && p1.back() == sep && !p2.empty() && p2.front() == sep) {
            // Remove one separator if both have it
            tmp.pop_back();
        }
        tmp += p2;
        return tmp;
    }

    std::string get_eos_file_path(const std::string& base_dir, int eos_id_full) {
        std::ostringstream eos_id_ss;
        eos_id_ss << std::setw(5) << std::setfill('0') << eos_id_full;
        std::string eos_id_str = eos_id_ss.str();

        if (eos_id_str.length() < 3) return "";

        std::string material_id_str = eos_id_str.substr(0, 3);

        std::string path = base_dir;
        path = simple_path_join(path, "mat" + material_id_str);
        path = simple_path_join(path, "eos" + eos_id_str + ".dat");
        return path;
    }

} // namespace EOSUtils