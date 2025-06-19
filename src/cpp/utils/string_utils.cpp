#include "string_utils.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace EOStringUtils {

std::string trim_string(const std::string& str) {
  const std::string whitespace = " \t\n\r\f\v";
  size_t start = str.find_first_not_of(whitespace);
  if (start == std::string::npos) {
    return "";  // String is empty or all whitespace
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
  // Catch any other std::exception as a fallback, though stod usually throws
  // the above two.
  catch (const std::exception& e) {
    return false;
  }
}

std::string simple_path_join(const std::string& p1, const std::string& p2) {
  char sep = '/';  // Common separator
  std::string tmp = p1;

#ifdef _WIN32  // Or other Windows checks
  sep = '\\';
#endif

  if (!p1.empty() && p1.back() != sep && !p2.empty() && p2.front() != sep) {
    tmp += sep;
  } else if (!p1.empty() && p1.back() == sep && !p2.empty() &&
             p2.front() == sep) {
    // Remove one separator if both have it
    tmp.pop_back();
  }
  tmp += p2;
  return tmp;
}

int parse_eos_id_from_filename_stub(
    const std::string& filename_stub) {
  // Assuming filename_stub is like "eos10000"
  if (filename_stub.rfind("eos", 0) == 0 &&
      filename_stub.length() == 8) {  // "eos" + 5 digits
    try {
      return std::stoi(filename_stub.substr(3));
    } catch (const std::exception&) {
      return -1;  // Invalid format
    }
      }
  // Fallback for other naming or if parsing from full relative path
  // This part needs to be robust based on your actual relative_path format.
  // If relative_path is "mat100/eos10000.dat", extract "10000".
  size_t last_slash = filename_stub.find_last_of("/\\");
  std::string basename = (last_slash == std::string::npos)
                             ? filename_stub
                             : filename_stub.substr(last_slash + 1);
  if (basename.rfind("eos", 0) == 0) {
    size_t dot_pos = basename.find_first_of('.');
    std::string id_part = (dot_pos == std::string::npos)
                              ? basename.substr(3)
                              : basename.substr(3, dot_pos - 3);
    if (id_part.length() == 5) {
      try {
        return std::stoi(id_part);
      } catch (const std::exception&) {
      }
    }
  }
  return -1;  // Indicate failure to parse ID
}

std::string get_eos_relative_path_stub(
    int eos_id_full, std::string& out_material_group_subpath) {
  std::ostringstream eos_id_ss;
  eos_id_ss << std::setw(5) << std::setfill('0') << eos_id_full;
  std::string eos_id_str = eos_id_ss.str();

  if (eos_id_str.length() < 3) {
    out_material_group_subpath = "matunknown";  // Or handle error
    return "eos" + eos_id_str;  // Or return empty string on error
  }
  out_material_group_subpath = "mat" + eos_id_str.substr(0, 3);
  return out_material_group_subpath + "/eos" + eos_id_str;
}

}  // namespace EOSUtils