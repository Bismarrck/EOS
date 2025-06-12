// src/cpp/EquationOfStateV1.cpp
#include "EquationOfStateV1.h"
#include <iostream>     // For debug prints
#include <fstream>      // for std::ifstream
#include <sstream>      // For std::istringstream
#include <iomanip>      // For std::fixed, std::setprecision in main, can be useful here too
#include <algorithm> // For std::transform for string tolower
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string>
#include <cstring>
#include <functional>
#include "utils/string_utils.h"


// Helper for basic path joining (add to EquationOfStateV1.cpp or a common utils file)
// WARNING: This is a very basic path joiner, doesn't handle all edge cases
// like multiple separators or ensuring correct separator for OS.
// For robust cross-platform paths without C++17, a small utility or careful manual
// construction is needed.
static std::string simple_path_join(const std::string& p1, const std::string& p2) {
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

// Update get_eos_file_path in EquationOfStateV1.cpp
static std::string get_eos_file_path(const std::string& base_dir, int eos_id_full) {
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

static int parse_complicated_eos_params(std::ifstream& eos_file,
                                        const std::string& file_path_for_error,
                                        const std::vector<std::string>& param_order,
                                        std::vector<double>& out_flat_params) {
    out_flat_params.clear();
    std::map<std::string, std::vector<double>> parsed_keyed_params;
    std::string line;
    int line_number = 0;

    bool in_header_separator_block = false;
    std::string current_active_key = ""; // To track the key for value continuations

    while (std::getline(eos_file, line)) {
        line_number++;
        std::string processed_line = line;

        // Handle comments (full line '#' or '!')
        size_t comment_char_pos = processed_line.find_first_of("#!");
        if (comment_char_pos != std::string::npos) {
            processed_line = processed_line.substr(0, comment_char_pos);
        }
        processed_line = EOSUtils::trim_string(processed_line);

        if (processed_line.empty()) {
            // An empty line usually resets the "current_active_key" for continuations,
            // unless specific rules say otherwise. For now, let's assume it does.
            // However, if a key expects many values and they are sparse with blank lines,
            // this might be too strict. Let's keep current_active_key unless a new key is defined.
            continue;
        }

        // Handle "=========" separator blocks
        if (processed_line.find("====") == 0) {
            in_header_separator_block = !in_header_separator_block;
            current_active_key = ""; // Reset active key when exiting/entering separator
            continue;
        }
        if (in_header_separator_block) {
            current_active_key = ""; // Reset active key while inside separator
            continue;
        }

        size_t eq_pos = processed_line.find('=');
        std::string key_part = "";
        std::string value_part = "";

        if (eq_pos != std::string::npos) { // This line defines a new key or redefines one
            key_part = EOSUtils::trim_string(processed_line.substr(0, eq_pos));
            value_part =EOSUtils:: trim_string(processed_line.substr(eq_pos + 1));

            if (key_part.empty()) {
                std::cerr << "Error (" << file_path_for_error << ":" << line_number
                          << "): Missing key for line with '=': \"" << line << "\"" << std::endl;
                return EOS_ERROR_FILE_PARSE;
            }
            std::transform(key_part.begin(), key_part.end(), key_part.begin(), ::tolower); // Normalize key
            current_active_key = key_part; // This is now the active key

            // If this key was seen before, clear its old values if this is a re-definition
            // Or decide on append vs overwrite. For simplicity, let's assume re-definition clears.
            // However, typical continuation means first "key=" line starts values, subsequent are appended.
            // So, if key_part is new or different from previous non-empty current_active_key, it's a new list.
            // If parsed_keyed_params.count(current_active_key), we will append to it.
            if (!parsed_keyed_params.count(current_active_key)) {
                 parsed_keyed_params[current_active_key] = {}; // Ensure vector exists
            }

        } else { // No '=', so this line must be a continuation of values for current_active_key
            if (current_active_key.empty()) {
                // This line has no '=' and there's no active key (e.g., after a separator or at the start)
                // It could be an error, or just an ignorable line of values without a key.
                std::cerr << "Warning (" << file_path_for_error << ":" << line_number
                          << "): Line with values but no '=' and no active key: \"" << line << "\". Ignoring." << std::endl;
                continue;
            }
            value_part = processed_line; // The whole line is values
        }

        // Now parse values from value_part
        if (!value_part.empty()) {
            std::istringstream val_ss(value_part);
            std::string token;
            bool values_found_on_this_line = false;
            while (val_ss >> token) {
                double val_converted;
                if (!EOSUtils::string_to_double_fortran_compat(token, val_converted)) {
                    std::cerr << "Error (" << file_path_for_error << ":" << line_number
                              << "): Failed to convert value token '" << token << "' for key '" << current_active_key << "'" << std::endl;
                    return EOS_ERROR_FILE_PARSE;
                }
                // Append to the current_active_key's vector
                if (current_active_key.empty()) { /* Should not happen due to check above */ }
                else {
                    parsed_keyed_params[current_active_key].push_back(val_converted);
                    values_found_on_this_line = true;
                }
            }
            // If value_part was not empty, but no numeric tokens were extracted (e.g., "key = non_numeric_stuff")
            if (!values_found_on_this_line && !value_part.empty()) {
                 std::cerr << "Error (" << file_path_for_error << ":" << line_number
                           << "): No valid numeric values found for key '" << current_active_key
                           << "' from value string: \"" << value_part << "\"" << std::endl;
                 return EOS_ERROR_FILE_PARSE;
            }
        }
        // If value_part was empty (e.g. "key = #comment" or "key ="), do nothing for values.
    }

    // Assemble flat_params based on param_order (this part remains the same)
    for (const std::string& ordered_key_orig : param_order) {
        std::string ordered_key = ordered_key_orig;
        std::transform(ordered_key.begin(), ordered_key.end(), ordered_key.begin(), ::tolower);

        auto it = parsed_keyed_params.find(ordered_key);
        if (it == parsed_keyed_params.end()) {
            std::cerr << "Error (" << file_path_for_error
                      << "): Required parameter key '" << ordered_key_orig << "' not found in file." << std::endl;
            return EOS_ERROR_FILE_PARSE;
        }
        if (it->second.empty()) {
            std::cerr << "Warning (" << file_path_for_error
                      << "): Required parameter key '" << ordered_key_orig << "' was found but has no associated values." << std::endl;
            // Depending on requirements, this could be an error.
            // For now, we'll insert an empty vector's worth of data (i.e., nothing), which might cause issues later
            // if the Fortran code expects a certain number of parameters for this key.
        }
        out_flat_params.insert(out_flat_params.end(), it->second.begin(), it->second.end());
    }

    return EOS_SUCCESS;
}

// --- Shim Implementations ---
int EquationOfStateV1::call_air_eos_2000_shim(double rho, double T, double& P, double& E) {
    int istat;
    // Outputs for dPdT etc. are not provided by this analytic EOS, so pass dummies or handle if API changes
    c_air_eos_2000(rho, T, &P, &E, &istat);
    return istat;
}

int EquationOfStateV1::call_carbon_eos_2001_shim(double rho, double T, double& P, double& E) {
    int istat;
    c_carbon_eos_2001(rho, T, &P, &E, &istat);
    return istat;
}

EquationOfStateV1::EquationOfStateV1() : tfd_use_ver1_setting_(true) { // Default TFD version
    std::cout << "C++: EquationOfStateV1 instance created." << std::endl;

    // Populate the analytic EOS registry
    // The key is the `eos_id` that would be found in the filename (e.g., from materials_list)
    // For this example, let's say analytic EOS "air" has ID 1, "carbon" has ID 2
    // These IDs need to be consistent with how you identify them in your eos_data_dir/#TYPE lines
    // or in the `eos_id_list` passed to initialize.
    // If analytic EOS files also follow the matXXX/eosXXXXX.dat pattern, then those full XXXXX IDs are used.
    // For simplicity, let's assume special IDs for these hardcoded analytic types for now.
    // If eos10000.dat said #TYPE: ANALYTIC and was for air_eos_2000:
    // analytic_eos_registry_[10000] = {InternalEOSType::ANALYTIC_AIR, &EquationOfStateV1::call_air_eos_2000_shim};
    // Or, if you have designated IDs for them:
    analytic_eos_registry_[1] = {InternalEOSType::ANALYTIC_AIR, &EquationOfStateV1::call_air_eos_2000_shim}; // User must provide '1' in eos_id_list
    analytic_eos_registry_[2] = {InternalEOSType::ANALYTIC_CARBON, &EquationOfStateV1::call_carbon_eos_2001_shim}; // User must provide '2'

    // Initialize TFD data struct
    tfd_data.initialized = false;
}

EquationOfStateV1::~EquationOfStateV1() {
    std::cout << "C++: EquationOfStateV1 instance destroyed." << std::endl;
    free_resources();
}

// --- Control Variable Management ---
int EquationOfStateV1::setComplicatedEOSUseOldColdTerm(bool value) {
    complicated_eos_use_old_cold_term_setting_ = value;
    int istat;
    c_set_complicated_eos_use_old_cold_term(value, &istat);
    if (istat != 0) {
        std::cerr << "C++ Error: Failed to set complicated_eos_use_old_cold_term. Fortran istat: " << istat << std::endl;
    }
    return istat;
}

int EquationOfStateV1::setUseTFDDataVer1(bool value) {
    tfd_use_ver1_setting_ = value;
    int istat;
    c_set_use_tfd_data_ver1(value, &istat);
    if (istat != 0) {
        std::cerr << "C++ Error: Failed to set use_tfd_data_ver1. Fortran istat: " << istat << std::endl;
    }
    return istat;
}

// --- HDF5 Helper Implementations ---
herr_t EquationOfStateV1::read_hdf5_scalar_int(hid_t group_id_or_file_id, const char* dset_name, int& out_val) {
    hid_t dset_id = H5Dopen2(group_id_or_file_id, dset_name, H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "HDF5 Error: Could not open integer dataset '" << dset_name << "'" << std::endl;
        return -1; // Indicate error
    }
    herr_t status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &out_val);
    if (status < 0) {
        std::cerr << "HDF5 Error: Could not read integer dataset '" << dset_name << "'" << std::endl;
    }
    H5Dclose(dset_id);
    return status;
}

herr_t EquationOfStateV1::read_hdf5_dataset_double(hid_t file_id, const char* dset_name,
                                                 int expected_n1, int expected_n2,
                                                 std::vector<double>& out_data) {
    hid_t dset_id = H5Dopen2(file_id, dset_name, H5P_DEFAULT);
    if (dset_id < 0) {
        std::cerr << "HDF5 Error: Could not open double dataset '" << dset_name << "'" << std::endl;
        return -1;
    }

    hid_t dataspace_id = H5Dget_space(dset_id);
    if (dataspace_id < 0) {
        std::cerr << "HDF5 Error: Could not get dataspace for '" << dset_name << "'" << std::endl;
        H5Dclose(dset_id);
        return -1;
    }

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 2) { // Assuming 2D matrices
        std::cerr << "HDF5 Error: Dataset '" << dset_name << "' is not 2D (rank=" << rank << ")" << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dset_id);
        return -1;
    }

    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    if (static_cast<int>(dims[0]) != expected_n1 || static_cast<int>(dims[1]) != expected_n2) {
        std::cerr << "HDF5 Error: Dimension mismatch for dataset '" << dset_name << "'. "
                  << "Expected " << expected_n1 << "x" << expected_n2
                  << ", Got " << dims[0] << "x" << dims[1] << std::endl;
        H5Sclose(dataspace_id);
        H5Dclose(dset_id);
        return -1;
    }

    out_data.resize(static_cast<size_t>(expected_n1) * expected_n2);
    herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_data.data());
    if (status < 0) {
        std::cerr << "HDF5 Error: Could not read double dataset '" << dset_name << "'" << std::endl;
    }

    H5Sclose(dataspace_id);
    H5Dclose(dset_id);
    return status;
}

// Helper for reading string attributes
static bool read_hdf5_string_attribute(hid_t object_id, const char* attr_name, std::string& out_str) {
    out_str.clear();
    if (H5Aexists(object_id, attr_name) <= 0) {
        return true; // Attribute doesn't exist, not an error for this helper
    }

    hid_t attr_id = H5Aopen_name(object_id, attr_name);
    if (attr_id < 0) { /* error */ return false; }

    hid_t file_type_id = H5Aget_type(attr_id);
    if (file_type_id < 0) { /* error */ H5Aclose(attr_id); return false; }

    // Check if it's a variable-length string in the file
    htri_t is_vlen_str_in_file = H5Tis_variable_str(file_type_id);
    if (is_vlen_str_in_file < 0) { /* error */ H5Tclose(file_type_id); H5Aclose(attr_id); return false; }

    hid_t mem_type_id = H5Tcopy(H5T_C_S1); // Start with a C string type
    if (mem_type_id < 0) { /* error */ H5Tclose(file_type_id); H5Aclose(attr_id); return false; }

    if (is_vlen_str_in_file > 0) { // Variable-length string in file
        if (H5Tset_size(mem_type_id, H5T_VARIABLE) < 0) { /* error */ H5Tclose(mem_type_id); H5Tclose(file_type_id); H5Aclose(attr_id); return false; }
    } else if (H5Tget_class(file_type_id) == H5T_STRING) { // Fixed-length string in file
        size_t fixed_size = H5Tget_size(file_type_id);
        if (fixed_size == 0 || fixed_size == H5T_VARIABLE) { /* error: unexpected */ H5Tclose(mem_type_id); H5Tclose(file_type_id); H5Aclose(attr_id); return false; }
        if (H5Tset_size(mem_type_id, fixed_size) < 0) { /* error */ H5Tclose(mem_type_id); H5Tclose(file_type_id); H5Aclose(attr_id); return false; }
    } else { /* error: not a string type */ H5Tclose(mem_type_id); H5Tclose(file_type_id); H5Aclose(attr_id); return false; }

    // *** Crucial: Match character set if reading UTF-8 from file ***
    // h5dump shows CSET H5T_CSET_UTF8 for the attributes.
    if (H5Tget_cset(file_type_id) == H5T_CSET_UTF8) {
        if (H5Tset_cset(mem_type_id, H5T_CSET_UTF8) < 0) {
            std::cerr << "HDF5 Warning: Failed to set memory type character set to UTF-8 for attribute '" << attr_name << "'." << std::endl;
            // Continue anyway, but conversion might still fail or be lossy.
        }
    }
    // Ensure null termination for memory type (good practice)
    // H5Tset_strpad(mem_type_id, H5T_STR_NULLTERM); // Usually H5T_C_S1 is already null-terminated

    char* str_read_buffer_vlen = nullptr;
    std::vector<char> str_read_buffer_fixed;
    void* actual_read_ptr;

    if (is_vlen_str_in_file > 0) {
        actual_read_ptr = &str_read_buffer_vlen; // H5Aread expects char**
    } else {
        size_t mem_size = H5Tget_size(mem_type_id);
        str_read_buffer_fixed.resize(mem_size + 1, '\0'); // +1 for safety null term
        actual_read_ptr = str_read_buffer_fixed.data();
    }

    herr_t read_status = H5Aread(attr_id, mem_type_id, actual_read_ptr);

    if (read_status < 0) {
        std::cerr << "HDF5 Error: H5Aread failed for attribute '" << attr_name << "'." << std::endl;
        H5Eprint(H5E_DEFAULT, stderr);
        if (is_vlen_str_in_file > 0 && str_read_buffer_vlen) H5free_memory(str_read_buffer_vlen);
        H5Tclose(mem_type_id); H5Tclose(file_type_id); H5Aclose(attr_id); return false;
    }

    if (is_vlen_str_in_file > 0) {
        if (str_read_buffer_vlen) {
            out_str = str_read_buffer_vlen;
            H5free_memory(str_read_buffer_vlen);
        } else { out_str = ""; }
    } else {
        out_str = str_read_buffer_fixed.data();
    }

    H5Tclose(mem_type_id);
    H5Tclose(file_type_id);
    H5Aclose(attr_id);
    return true;
}


// Helper for reading scalar variable-length string dataset
static bool read_hdf5_scalar_vlen_string_dataset(hid_t file_id, const char* dset_name, std::string& out_str) {
    out_str.clear();
    hid_t dset_id = H5Dopen2(file_id, dset_name, H5P_DEFAULT);
    if (dset_id < 0) { /* ... */ return false; }

    hid_t file_type_id = H5Dget_type(dset_id);
    if (file_type_id < 0) { /* ... */ H5Dclose(dset_id); return false; }

    // Verify it's a variable-length string in the file
    htri_t is_vlen_str_in_file = H5Tis_variable_str(file_type_id);
    if (is_vlen_str_in_file <= 0) { /* ... error message ... */ H5Tclose(file_type_id); H5Dclose(dset_id); return false; }

    // Verify it's scalar
    hid_t dspace_id = H5Dget_space(dset_id);
    if (dspace_id < 0 || H5Sget_simple_extent_type(dspace_id) != H5S_SCALAR) { /* ... */ if (dspace_id >= 0) H5Sclose(dspace_id); H5Tclose(file_type_id); H5Dclose(dset_id); return false; }
    H5Sclose(dspace_id);

    // Prepare memory type
    hid_t mem_type_id = H5Tcopy(H5T_C_S1);
    if (mem_type_id < 0) { /* error */ H5Tclose(file_type_id); H5Dclose(dset_id); return false; }
    if (H5Tset_size(mem_type_id, H5T_VARIABLE) < 0) { /* error */ H5Tclose(mem_type_id); H5Tclose(file_type_id); H5Dclose(dset_id); return false; }

    // *** Match character set for dataset as well ***
    if (H5Tget_cset(file_type_id) == H5T_CSET_UTF8) { // From h5dump, dataset is also UTF8
        if (H5Tset_cset(mem_type_id, H5T_CSET_UTF8) < 0) {
            std::cerr << "HDF5 Warning: Failed to set memory type character set to UTF-8 for dataset '" << dset_name << "'." << std::endl;
        }
    }

    char* str_read_buffer = nullptr;
    herr_t read_status = H5Dread(dset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &str_read_buffer);

    if (read_status < 0) {
        std::cerr << "HDF5 Error: H5Dread failed for dataset '" << dset_name << "'." << std::endl;
        H5Eprint(H5E_DEFAULT, stderr);
        if (str_read_buffer) H5free_memory(str_read_buffer);
        H5Tclose(mem_type_id); H5Tclose(file_type_id); H5Dclose(dset_id); return false;
    }

    if (str_read_buffer) {
        out_str = str_read_buffer;
        H5free_memory(str_read_buffer);
    } else { out_str = ""; }

    H5Tclose(mem_type_id);
    H5Tclose(file_type_id);
    H5Dclose(dset_id);
    return true;
}

int EquationOfStateV1::loadTFDDataInternal(const std::string& hdf5_filepath) {
    tfd_data.initialized = false; // Reset
    std::cout << "C++ (loadTFDDataInternal): Attempting to load TFD data from HDF5 file: " << hdf5_filepath << std::endl;

    // Suppress HDF5 default error printing, we'll handle errors via return codes
    H5E_auto2_t old_func;
    void *old_client_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL); // Turn off HDF5 auto error printing

    hid_t file_id = H5Fopen(hdf5_filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "HDF5 Error: Could not open TFD HDF5 file: " << hdf5_filepath << std::endl;
        H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data); // Restore error printing
        return EOS_ERROR_FILE_NOT_FOUND; // Or a specific HDF5 error code
    }

    // Read file-level attributes
    std::string creation_date, tfd_version_attr;
    if (read_hdf5_string_attribute(file_id, "file_creation_date", creation_date)) {
        std::cout << "HDF5 Attr: File Creation Date = " << creation_date << std::endl;
    }
    if (read_hdf5_string_attribute(file_id, "tfd_model_version", tfd_version_attr)) {
        std::cout << "HDF5 Attr: TFD Model Version = " << tfd_version_attr << std::endl;
    }

    // Read model description dataset
    if (read_hdf5_scalar_vlen_string_dataset(file_id, "/model_description", tfd_data.tfd_file_description)) { // Assuming tfd_file_description is a member
        std::cout << "HDF5 Dataset: Model Description loaded (length " << tfd_data.tfd_file_description.length() << ")." << std::endl;
        // For brevity, don't print long descriptions here.
    } else {
        std::cout << "HDF5 Info: No '/model_description' dataset found or failed to read." << std::endl;
    }

    herr_t status;
    hid_t common_dims_group_id = H5Gopen2(file_id, "/common_dimensions", H5P_DEFAULT);
    if (common_dims_group_id < 0) {
        std::cerr << "HDF5 Error: Could not open group '/common_dimensions' in " << hdf5_filepath << std::endl;
        H5Fclose(file_id);
        H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
        return EOS_ERROR_FILE_PARSE; // Or HDF5 specific
    }

    status = read_hdf5_scalar_int(common_dims_group_id, "N1", tfd_data.N1);
    if (status < 0) { H5Gclose(common_dims_group_id); H5Fclose(file_id); H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data); return EOS_ERROR_FILE_PARSE; }
    status = read_hdf5_scalar_int(common_dims_group_id, "N2", tfd_data.N2);
    if (status < 0) { H5Gclose(common_dims_group_id); H5Fclose(file_id); H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data); return EOS_ERROR_FILE_PARSE; }
    H5Gclose(common_dims_group_id);

    if (tfd_data.N1 <= 0 || tfd_data.N2 <= 0) {
        std::cerr << "HDF5 Error: Invalid dimensions N1=" << tfd_data.N1 << ", N2=" << tfd_data.N2
                  << " read from " << hdf5_filepath << std::endl;
        H5Fclose(file_id);
        H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
        return EOS_ERROR_INVALID_DIMENSIONS;
    }
    std::cout << "HDF5: Read TFD dimensions: " << tfd_data.N1 << "x" << tfd_data.N2 << std::endl;

    // Determine which TFD version group to read from
    // std::string tfd_version_group = (tfd_use_ver1_setting_ ? "/tfd_version_1" : "/tfd_version_2");
    // For simplicity, let's assume datasets are at root or a fixed path for now,
    // and the file name itself (tfd_ver1.h5 vs tfd_ver2.h5) distinguishes versions.
    // Or, if the file contains multiple versions:
    // hid_t version_group_id = H5Gopen2(file_id, tfd_version_group.c_str(), H5P_DEFAULT);
    // if (version_group_id < 0) { /* error */ }
    // status = read_hdf5_dataset_double(version_group_id, "MatrixA", ...);
    // H5Gclose(version_group_id);

    // Assuming MatrixA and MatrixB are top-level datasets in the HDF5 file for now
    status = read_hdf5_dataset_double(file_id, "/matrix_A", tfd_data.N1, tfd_data.N2, tfd_data.matrix_A);
    if (status < 0) { H5Fclose(file_id); H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data); return EOS_ERROR_TFD_LOAD_A; }
    std::cout << "HDF5: Successfully read Matrix A." << std::endl;

    status = read_hdf5_dataset_double(file_id, "/matrix_B", tfd_data.N1, tfd_data.N2, tfd_data.matrix_B);
    if (status < 0) { H5Fclose(file_id); H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data); return EOS_ERROR_TFD_LOAD_B; }
    std::cout << "HDF5: Successfully read Matrix B." << std::endl;

    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data); // Restore HDF5 auto error printing

    tfd_data.initialized = true;
    // tfd_data.loaded_tfd_filepath = hdf5_filepath; // Store which file was loaded
    std::cout << "C++: TFD data loaded successfully from HDF5 file: " << hdf5_filepath << std::endl;
    return EOS_SUCCESS;
}

// --- Initialization ---
int EquationOfStateV1::initialize(const std::vector<int>& eos_id_list, const std::string& eos_data_dir) {
    // 0. Ensure Fortran side knows which TFD version we prefer.
    //    The C++ tfd_use_ver1_setting_ is the master; set it if user calls the setter.
    //    Then ensure Fortran is synced.
    int f_istat;
    c_set_use_tfd_data_ver1(tfd_use_ver1_setting_, &f_istat);
    if (f_istat != 0) {
        std::cerr << "Error setting Fortran TFD version preference." << std::endl;
        // Decide if this is a fatal error for initialize
    }

    // 1. Load TFD data
    std::string tfd_filename_to_load = (tfd_use_ver1_setting_ ? "tfd_ver1.h5" : "tfd_ver2.h5");
    std::string full_tfd_hdf5_path = eos_data_dir + "/" + tfd_filename_to_load; // Basic join
    int tfd_load_status = loadTFDDataInternal(full_tfd_hdf5_path);
    if (tfd_load_status != EOS_SUCCESS) {
        std::cerr << "Error initializing: Failed to load TFD data. Code: " << tfd_load_status << std::endl;
        return tfd_load_status;
    }

    // 2. Load data for each EOS ID
    for (int id : eos_id_list) {
        std::string file_path_str = get_eos_file_path(eos_data_dir, id);
        if (file_path_str.empty()) {
            std::cerr << "Error: Could not construct path for EOS ID " << id << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND; // Or collect errors and report all
        }

        std::ifstream eos_file(file_path_str);
        if (!eos_file.is_open()) {
            // Special check: is this an ID for a known, hardcoded analytic EOS?
            auto it_analytic_reg = analytic_eos_registry_.find(id);
            if (it_analytic_reg != analytic_eos_registry_.end()) {
                MaterialData md;
                md.eos_id_key = id;
                md.type_from_file = EOSTypeFromFile::ANALYTIC; // By virtue of being in registry
                md.internal_type = it_analytic_reg->second.internal_type;
                md.analytic_func = it_analytic_reg->second.func_ptr;
                loaded_materials_[id] = md;
                std::cout << "Initialized hardcoded Analytic EOS ID " << id << " (" << file_path_str << " not strictly needed but implies registration)." << std::endl;
                continue; // Successfully processed this analytic ID
            }
            std::cerr << "Error: Could not open EOS data file: " << file_path_str << " for ID " << id << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND;
        }

        std::string line;
        std::string type_str_from_file;
        if (std::getline(eos_file, line)) {
            std::string prefix = "#TYPE:";
            size_t prefix_pos = line.find(prefix);
            if (prefix_pos != std::string::npos) {
                type_str_from_file = line.substr(prefix_pos + prefix.length());
                // Trim whitespace
                type_str_from_file.erase(0, type_str_from_file.find_first_not_of(" \t\n\r\f\v"));
                type_str_from_file.erase(type_str_from_file.find_last_not_of(" \t\n\r\f\v") + 1);
                std::transform(type_str_from_file.begin(), type_str_from_file.end(), type_str_from_file.begin(), ::tolower);
            } else {
                std::cerr << "Error: Missing or malformed #TYPE: line in " << file_path_str << std::endl;
                return EOS_ERROR_FILE_PARSE;
            }
        } else {
            std::cerr << "Error: Could not read #TYPE: line from " << file_path_str << std::endl;
            return EOS_ERROR_FILE_PARSE;
        }

        MaterialData md;
        md.eos_id_key = id;

        if (type_str_from_file == "analytic") {
            md.type_from_file = EOSTypeFromFile::ANALYTIC;
            auto it_analytic_reg = analytic_eos_registry_.find(id);
            if (it_analytic_reg != analytic_eos_registry_.end()) {
                md.internal_type = it_analytic_reg->second.internal_type;
                md.analytic_func = it_analytic_reg->second.func_ptr;
            } else {
                std::cerr << "Error: Analytic EOS ID " << id << " from file " << file_path_str
                          << " is not registered in C++ analytic_eos_registry_." << std::endl;
                return EOS_ERROR_UNKNOWN_EOS_ID;
            }
            // For analytic, the .dat file might contain no other parameters.
        } else if (type_str_from_file == "complicated") {
            md.type_from_file = EOSTypeFromFile::COMPLICATED;
            md.internal_type = InternalEOSType::COMPLICATED_V1; // Assuming one type for now

            // Define the expected order of parameter keys for this *specific* Complicated EOS model
            // This order MUST match what the Fortran routine `Complicated_EOS` expects for its flat `params_in` array.
            // This might need to be configured per complicated EOS ID if different models have different param orders.
            // For now, assume a global order for all COMPLICATED_V1 types.
            std::vector<std::string> expected_param_order = {"key1", "key2", "key3", "key4"}; // EXAMPLE ORDER

            int parse_stat = parse_complicated_eos_params(eos_file, file_path_str, expected_param_order, md.params);
            if (parse_stat != EOS_SUCCESS) {
                std::cerr << "Error parsing parameters for complicated EOS ID " << id << " from " << file_path_str << std::endl;
                return parse_stat;
            }

            if (md.params.empty() && !expected_param_order.empty()) { // If order expects params but none were loaded
                std::cerr << "Warning: No parameters loaded for complicated EOS ID " << id << " from " << file_path_str
                          << " despite expected order." << std::endl;
                // This could be an error depending on whether all ordered params are mandatory.
            }
            std::cout << "Read " << md.params.size() << " params for complicated EOS ID " << id
                      << " based on key-value format and order." << std::endl;

        } else {
            std::cerr << "Error: Unknown #TYPE: '" << type_str_from_file << "' in " << file_path_str << std::endl;
            return EOS_ERROR_INVALID_EOS_TYPE;
        }
        loaded_materials_[id] = md;
        std::cout << "Successfully loaded EOS ID " << id << " from " << file_path_str << std::endl;
    }
    return EOS_SUCCESS;
}

int EquationOfStateV1::check_eos_data_dir(const std::string& eos_data_dir, const std::vector<int>& eos_ids_to_check) {
    // Check TFD files
    std::string tfd_files_to_check[2] = { "tfd_ver1.dat", "tfd_ver2.dat" };

    for (const std::string& tfd_suffix : tfd_files_to_check) {
        std::string tfd_path_str = simple_path_join(eos_data_dir, tfd_suffix);
        std::ifstream test_file(tfd_path_str);
        if (!test_file.is_open()) {
            std::cerr << "Check Error: TFD file not found or not accessible: " << tfd_path_str << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND;
        }
        test_file.close();
    }

    // Check specific EOS ID files
    for (int id : eos_ids_to_check) {
        if (analytic_eos_registry_.count(id)) {
            std::cout << "Check Info: Analytic EOS ID " << id << " is registered, file existence check skipped/optional." << std::endl;
            continue;
        }
        std::string file_path_str = get_eos_file_path(eos_data_dir, id);
        if (file_path_str.empty()){
            std::cerr << "Check Error: Could not construct path for EOS ID " << id << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND; // Or a different error code
        }
        std::ifstream eos_file(file_path_str);
        if (!eos_file.is_open()) {
            std::cerr << "Check Error: EOS data file not found for ID " << id << " at " << file_path_str << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND;
        }
        eos_file.close();
    }
    std::cout << "check_eos_data_dir: All checked files appear accessible." << std::endl;
    return EOS_SUCCESS;
}

// --- Computation ---
int EquationOfStateV1::compute(int eos_id, double rho, double T,
                             double& P, double& E,
                             double& dPdT, double& dEdT, double& dPdrho) {
    auto it = loaded_materials_.find(eos_id);
    if (it == loaded_materials_.end()) {
        std::cerr << "Error: EOS ID " << eos_id << " not initialized." << std::endl;
        return EOS_ERROR_UNKNOWN_EOS_ID;
    }

    const MaterialData& md = it->second;
    P = 0.0; E = 0.0; dPdT = 0.0; dEdT = 0.0; dPdrho = 0.0; // Initialize outputs

    int istat = EOS_SUCCESS;

    switch (md.internal_type) {
        case InternalEOSType::ANALYTIC_AIR:
        case InternalEOSType::ANALYTIC_CARBON:
            if (md.analytic_func) {
                // Current analytic shims only set P, E. Derivatives would be zero.
                // If analytic functions also give derivatives, the AnalyticEOSFunc signature
                // and shims need to be updated.
                istat = md.analytic_func(rho, T, P, E);
                // dPdT, dEdT, dPdrho remain 0 for these specific analytic stubs
            } else {
                std::cerr << "Error: No analytic function registered for EOS ID " << eos_id << std::endl;
                istat = EOS_ERROR_ANALYTIC_DISPATCH;
            }
            break;

        case InternalEOSType::COMPLICATED_V1:
            if (!tfd_data.initialized) {
                std::cerr << "Error: TFD data not initialized for complicated EOS ID " << eos_id << std::endl;
                return EOS_ERROR_TFD_NOT_INIT;
            }
            if (md.params.empty() && md.type_from_file == EOSTypeFromFile::COMPLICATED) {
                 // This check depends on whether complicated EOS can ever have zero params
                std::cerr << "Warning: Complicated EOS ID " << eos_id << " has no parameters loaded." << std::endl;
                // return EOS_ERROR_PARAMS_NOT_LOADED; // Or proceed if allowed
            }
            // Call the C wrapper from Phase 3 directly for now
            c_complicated_eos(md.params.data(), static_cast<int>(md.params.size()), rho, T,
                              tfd_data.matrix_A.data(), tfd_data.N1, tfd_data.N2,
                              tfd_data.matrix_B.data(), tfd_data.N1, tfd_data.N2,
                              &P, &E, &dPdT, &dEdT, &dPdrho,
                              &istat);
            break;

        default:
            std::cerr << "Error: Unknown internal EOS type for ID " << eos_id << std::endl;
            istat = EOS_ERROR_UNKNOWN_EOS_ID;
            break;
    }

    if (istat != 0) {
        // Fortran istat might be directly returned, or mapped to C++ error codes
        std::cerr << "C++ compute: Fortran calculation for EOS ID " << eos_id << " failed with istat: " << istat << std::endl;
    }
    return istat; // Return Fortran's istat or a mapped C++ error code
}


void EquationOfStateV1::free_resources() {
    loaded_materials_.clear();
    tfd_data.matrix_A.clear();
    tfd_data.matrix_B.clear();
    tfd_data.initialized = false;
    std::cout << "C++: EquationOfStateV1 resources freed." << std::endl;
}


// Helper for packing data into buffer and advancing offset
template<typename T>
void pack_item(char*& current_ptr, const T& item) {
    std::memcpy(current_ptr, &item, sizeof(T));
    current_ptr += sizeof(T);
}

template<typename T>
void pack_vector_data(char*& current_ptr, const std::vector<T>& vec_data) {
    if (!vec_data.empty()) {
        std::memcpy(current_ptr, vec_data.data(), vec_data.size() * sizeof(T));
        current_ptr += vec_data.size() * sizeof(T);
    }
}

// Helper for unpacking data from buffer and advancing offset
template<typename T>
void unpack_item(const char*& current_ptr, T& item) {
    std::memcpy(&item, current_ptr, sizeof(T));
    current_ptr += sizeof(T);
}

template<typename T>
void unpack_vector_data(const char*& current_ptr, std::vector<T>& vec_data, size_t num_elements) {
    vec_data.resize(num_elements);
    if (num_elements > 0) {
        std::memcpy(vec_data.data(), current_ptr, num_elements * sizeof(T));
        current_ptr += num_elements * sizeof(T);
    }
}


int EquationOfStateV1::pack_data(char*& buffer, int& buffer_size) const {
    // --- Calculate total size needed ---
    int required_size = 0;
    required_size += sizeof(tfd_use_ver1_setting_);
    required_size += sizeof(complicated_eos_use_old_cold_term_setting_);

    required_size += sizeof(tfd_data.N1);
    required_size += sizeof(tfd_data.N2);
    if (tfd_data.initialized) {
        required_size += tfd_data.matrix_A.size() * sizeof(double);
        required_size += tfd_data.matrix_B.size() * sizeof(double);
    } else {
        // If TFD not initialized, we still pack N1=0, N2=0. No matrix data.
    }

    required_size += sizeof(int); // For num_loaded_materials_
    for (const auto& pair : loaded_materials_) {
        const MaterialData& md = pair.second;
        required_size += sizeof(md.eos_id_key);
        required_size += sizeof(InternalEOSType); // Pack the enum value
        required_size += sizeof(int); // For num_params
        required_size += md.params.size() * sizeof(double);
    }

    // --- Mode 1: Calculate size only ---
    if (buffer == nullptr) {
        buffer_size = required_size;
        return EOS_SUCCESS;
    }

    // --- Mode 2: Pack data if buffer is sufficient ---
    if (buffer_size < required_size) {
        std::cerr << "Pack Error: Provided buffer too small. Need " << required_size
                  << ", got " << buffer_size << std::endl;
        return -1; // Or a specific error code for buffer too small
    }

    char* current_ptr = buffer;

    pack_item(current_ptr, tfd_use_ver1_setting_);
    pack_item(current_ptr, complicated_eos_use_old_cold_term_setting_);

    if (tfd_data.initialized) {
        pack_item(current_ptr, tfd_data.N1);
        pack_item(current_ptr, tfd_data.N2);
        pack_vector_data(current_ptr, tfd_data.matrix_A);
        pack_vector_data(current_ptr, tfd_data.matrix_B);
    } else {
        int zero_dim = 0;
        pack_item(current_ptr, zero_dim); // N1 = 0
        pack_item(current_ptr, zero_dim); // N2 = 0
        // No matrix data to pack
    }

    int num_loaded = static_cast<int>(loaded_materials_.size());
    pack_item(current_ptr, num_loaded);

    for (const auto& pair : loaded_materials_) {
        const MaterialData& md = pair.second;
        pack_item(current_ptr, md.eos_id_key);
        pack_item(current_ptr, md.internal_type); // Packing enum directly

        int num_params = static_cast<int>(md.params.size());
        pack_item(current_ptr, num_params);
        pack_vector_data(current_ptr, md.params);
    }

    buffer_size = static_cast<int>(current_ptr - buffer); // Actual size used
    if (buffer_size != required_size) {
        // This should ideally not happen if logic is correct
        std::cerr << "Pack Error: Mismatch in calculated vs packed size. Packed: " << buffer_size
                  << ", Required: " << required_size << std::endl;
        return -2; // Packing logic error
    }

    std::cout << "C++ (pack_data): Successfully packed " << buffer_size << " bytes." << std::endl;
    return EOS_SUCCESS;
}


int EquationOfStateV1::unpack_data(const char* buffer, int buffer_size) {
    // Clear existing state before unpacking
    free_resources(); // Important to reset state

    const char* current_ptr = buffer;
    const char* buffer_end = buffer + buffer_size; // For boundary checks

    // Helper lambda for boundary check
    auto check_boundary = [&](size_t needed) -> bool {
        if (current_ptr + needed > buffer_end) {
            std::cerr << "Unpack Error: Buffer too small or data corrupted. "
                      << "Attempting to read " << needed << " bytes with "
                      << (buffer_end - current_ptr) << " bytes remaining." << std::endl;
            return false;
        }
        return true;
    };

    if (!check_boundary(sizeof(tfd_use_ver1_setting_))) return -1;
    unpack_item(current_ptr, tfd_use_ver1_setting_);

    if (!check_boundary(sizeof(complicated_eos_use_old_cold_term_setting_))) return -1;
    unpack_item(current_ptr, complicated_eos_use_old_cold_term_setting_);

    // Update Fortran side with unpacked control variables
    int f_istat;
    c_set_use_tfd_data_ver1(tfd_use_ver1_setting_, &f_istat); // Error check f_istat?
    c_set_complicated_eos_use_old_cold_term(complicated_eos_use_old_cold_term_setting_, &f_istat); // Error check?

    std::cout << "C++ (unpack_data): Unpacked tfd_use_ver1_setting=" << tfd_use_ver1_setting_
              << ", complicated_eos_use_old_cold_term_setting=" << complicated_eos_use_old_cold_term_setting_ << std::endl;


    if (!check_boundary(sizeof(tfd_data.N1))) return -1;
    unpack_item(current_ptr, tfd_data.N1);
    if (!check_boundary(sizeof(tfd_data.N2))) return -1;
    unpack_item(current_ptr, tfd_data.N2);

    if (tfd_data.N1 > 0 && tfd_data.N2 > 0) { // Check if TFD data was actually packed
        size_t tfd_a_bytes = static_cast<size_t>(tfd_data.N1) * tfd_data.N2 * sizeof(double);
        if (!check_boundary(tfd_a_bytes)) return -1;
        unpack_vector_data(current_ptr, tfd_data.matrix_A, static_cast<size_t>(tfd_data.N1) * tfd_data.N2);

        size_t tfd_b_bytes = static_cast<size_t>(tfd_data.N1) * tfd_data.N2 * sizeof(double);
        if (!check_boundary(tfd_b_bytes)) return -1;
        unpack_vector_data(current_ptr, tfd_data.matrix_B, static_cast<size_t>(tfd_data.N1) * tfd_data.N2);
        tfd_data.initialized = true;
        std::cout << "C++ (unpack_data): Unpacked TFD matrices (" << tfd_data.N1 << "x" << tfd_data.N2 << ")." << std::endl;
    } else {
        tfd_data.initialized = false; // N1 or N2 was 0, so no TFD data
        std::cout << "C++ (unpack_data): TFD dimensions were zero, no TFD matrices unpacked." << std::endl;
    }


    int num_loaded_materials;
    if (!check_boundary(sizeof(num_loaded_materials))) return -1;
    unpack_item(current_ptr, num_loaded_materials);
    std::cout << "C++ (unpack_data): Expecting " << num_loaded_materials << " materials." << std::endl;

    for (int i = 0; i < num_loaded_materials; ++i) {
        MaterialData md;
        if (!check_boundary(sizeof(md.eos_id_key))) return -1;
        unpack_item(current_ptr, md.eos_id_key);

        if (!check_boundary(sizeof(md.internal_type))) return -1;
        unpack_item(current_ptr, md.internal_type); // Unpacking enum directly

        int num_params;
        if (!check_boundary(sizeof(num_params))) return -1;
        unpack_item(current_ptr, num_params);

        if (num_params > 0) {
            size_t params_bytes = static_cast<size_t>(num_params) * sizeof(double);
            if (!check_boundary(params_bytes)) return -1;
            unpack_vector_data(current_ptr, md.params, num_params);
        }

        // Re-establish analytic_func pointer based on internal_type and eos_id_key
        // This relies on analytic_eos_registry_ being populated (e.g. in constructor)
        // Or, find it in the registry by the unpacked md.eos_id_key.
        auto reg_it = analytic_eos_registry_.find(md.eos_id_key);
        if (reg_it != analytic_eos_registry_.end()) {
            // Cross-check if internal_type matches what's in registry for this ID
            if (reg_it->second.internal_type == md.internal_type) {
                md.analytic_func = reg_it->second.func_ptr;
            } else {
                 std::cerr << "Unpack Warning: Mismatch in unpacked internal_type for analytic EOS ID "
                           << md.eos_id_key << ". Using registry's type." << std::endl;
                 md.internal_type = reg_it->second.internal_type; // Prefer registry
                 md.analytic_func = reg_it->second.func_ptr;
            }
        } else if (md.internal_type == InternalEOSType::ANALYTIC_AIR || md.internal_type == InternalEOSType::ANALYTIC_CARBON) {
            // This case should ideally not happen if eos_id_key was used for registry lookup
            std::cerr << "Unpack Error: Unpacked analytic EOS ID " << md.eos_id_key
                      << " not found in registry, but type suggests it's analytic." << std::endl;
            return -1; // Data integrity issue
        }

        // Set type_from_file based on internal_type (approximate reverse)
        if (md.internal_type == InternalEOSType::ANALYTIC_AIR || md.internal_type == InternalEOSType::ANALYTIC_CARBON) {
            md.type_from_file = EOSTypeFromFile::ANALYTIC;
        } else if (md.internal_type == InternalEOSType::COMPLICATED_V1) {
            md.type_from_file = EOSTypeFromFile::COMPLICATED;
        } else {
            md.type_from_file = EOSTypeFromFile::NOT_SET; // Or UNKNOWN
        }

        loaded_materials_[md.eos_id_key] = md;
        std::cout << "C++ (unpack_data): Unpacked material ID " << md.eos_id_key
                  << ", type " << static_cast<int>(md.internal_type)
                  << ", num_params " << num_params << std::endl;
    }

    // Final check: did we consume the whole buffer?
    if (current_ptr != buffer_end) {
        std::cerr << "Unpack Warning: " << (buffer_end - current_ptr)
                  << " bytes remaining in buffer after unpack. Expected to consume all "
                  << buffer_size << " bytes." << std::endl;
        // This might be an error depending on strictness.
    }

    std::cout << "C++ (unpack_data): Successfully unpacked data." << std::endl;
    return EOS_SUCCESS;
}
