#include "parser_utils.h"

#include <algorithm>  // For std::transform
#include <iostream>
#include <map>
#include <sstream>

#include "eos_error_codes.h"
#include "string_utils.h"

int EOSParserUtils::parse_complicated_eos_params(
    std::ifstream& eos_file, const std::string& file_path_for_error,
    const std::vector<std::string>& param_order,
    std::vector<double>& out_flat_params) {
  out_flat_params.clear();
  std::map<std::string, std::vector<double>> parsed_keyed_params;
  std::string line;
  int line_number = 0;

  bool in_header_separator_block = false;
  std::string current_active_key =
      "";  // To track the key for value continuations

  while (std::getline(eos_file, line)) {
    line_number++;
    std::string processed_line = line;

    // Handle comments (full line '#' or '!')
    size_t comment_char_pos = processed_line.find_first_of("#!");
    if (comment_char_pos != std::string::npos) {
      processed_line = processed_line.substr(0, comment_char_pos);
    }
    processed_line = EOStringUtils::trim_string(processed_line);

    if (processed_line.empty()) {
      // An empty line usually resets the "current_active_key" for
      // continuations, unless specific rules say otherwise. For now, let's
      // assume it does. However, if a key expects many values and they are
      // sparse with blank lines, this might be too strict. Let's keep
      // current_active_key unless a new key is defined.
      continue;
    }

    // Handle "=========" separator blocks
    if (processed_line.find("====") == 0) {
      in_header_separator_block = !in_header_separator_block;
      current_active_key =
          "";  // Reset active key when exiting/entering separator
      continue;
    }
    if (in_header_separator_block) {
      current_active_key = "";  // Reset active key while inside separator
      continue;
    }

    size_t eq_pos = processed_line.find('=');
    std::string key_part = "";
    std::string value_part = "";

    if (eq_pos !=
        std::string::npos) {  // This line defines a new key or redefines one
      key_part = EOStringUtils::trim_string(processed_line.substr(0, eq_pos));
      value_part =
          EOStringUtils::trim_string(processed_line.substr(eq_pos + 1));

      if (key_part.empty()) {
        std::cerr << "Error (" << file_path_for_error << ":" << line_number
                  << "): Missing key for line with '=': \"" << line << "\""
                  << std::endl;
        return EOS_ERROR_FILE_PARSE;
      }
      std::transform(key_part.begin(), key_part.end(), key_part.begin(),
                     ::tolower);      // Normalize key
      current_active_key = key_part;  // This is now the active key

      // If this key was seen before, clear its old values if this is a
      // re-definition Or decide on append vs overwrite. For simplicity, let's
      // assume re-definition clears. However, typical continuation means first
      // "key=" line starts values, subsequent are appended. So, if key_part is
      // new or different from previous non-empty current_active_key, it's a new
      // list. If parsed_keyed_params.count(current_active_key), we will append
      // to it.
      if (!parsed_keyed_params.count(current_active_key)) {
        parsed_keyed_params[current_active_key] = {};  // Ensure vector exists
      }

    } else {  // No '=', so this line must be a continuation of values for
              // current_active_key
      if (current_active_key.empty()) {
        // This line has no '=' and there's no active key (e.g., after a
        // separator or at the start) It could be an error, or just an ignorable
        // line of values without a key.
        std::cerr << "Warning (" << file_path_for_error << ":" << line_number
                  << "): Line with values but no '=' and no active key: \""
                  << line << "\". Ignoring." << std::endl;
        continue;
      }
      value_part = processed_line;  // The whole line is values
    }

    // Now parse values from value_part
    if (!value_part.empty()) {
      std::istringstream val_ss(value_part);
      std::string token;
      bool values_found_on_this_line = false;
      while (val_ss >> token) {
        double val_converted;
        if (!EOStringUtils::string_to_double_fortran_compat(token,
                                                            val_converted)) {
          std::cerr << "Error (" << file_path_for_error << ":" << line_number
                    << "): Failed to convert value token '" << token
                    << "' for key '" << current_active_key << "'" << std::endl;
          return EOS_ERROR_FILE_PARSE;
        }
        // Append to the current_active_key's vector
        if (current_active_key
                .empty()) { /* Should not happen due to check above */
        } else {
          parsed_keyed_params[current_active_key].push_back(val_converted);
          values_found_on_this_line = true;
        }
      }
      // If value_part was not empty, but no numeric tokens were extracted
      // (e.g., "key = non_numeric_stuff")
      if (!values_found_on_this_line && !value_part.empty()) {
        std::cerr << "Error (" << file_path_for_error << ":" << line_number
                  << "): No valid numeric values found for key '"
                  << current_active_key << "' from value string: \""
                  << value_part << "\"" << std::endl;
        return EOS_ERROR_FILE_PARSE;
      }
    }
    // If value_part was empty (e.g. "key = #comment" or "key ="), do nothing
    // for values.
  }

  // Assemble flat_params based on param_order (this part remains the same)
  for (const std::string& ordered_key_orig : param_order) {
    std::string ordered_key = ordered_key_orig;
    std::transform(ordered_key.begin(), ordered_key.end(), ordered_key.begin(),
                   ::tolower);

    auto it = parsed_keyed_params.find(ordered_key);
    if (it == parsed_keyed_params.end()) {
      std::cerr << "Error (" << file_path_for_error
                << "): Required parameter key '" << ordered_key_orig
                << "' not found in file." << std::endl;
      return EOS_ERROR_FILE_PARSE;
    }
    if (it->second.empty()) {
      std::cerr << "Warning (" << file_path_for_error
                << "): Required parameter key '" << ordered_key_orig
                << "' was found but has no associated values." << std::endl;
      // Depending on requirements, this could be an error.
      // For now, we'll insert an empty vector's worth of data (i.e., nothing),
      // which might cause issues later if the Fortran code expects a certain
      // number of parameters for this key.
    }
    out_flat_params.insert(out_flat_params.end(), it->second.begin(),
                           it->second.end());
  }

  return EOS_SUCCESS;
}
