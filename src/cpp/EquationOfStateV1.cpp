// src/cpp/EquationOfStateV1.cpp
#include "EquationOfStateV1.h"

#include <hdf5.h>

#include <algorithm>  // For std::transform for string tolower
#include <cstring>
#include <fstream>  // for std::ifstream
#include <functional>
#include <iomanip>  // For std::fixed, std::setprecision in main, can be useful here too
#include <iostream>  // For debug prints
#include <sstream>   // For std::istringstream
#include <string>

#include "eos_checksums.h"
#include "materials/AnalyticEOS.h"
#include "materials/ComplicatedLegacyEOS.h"
#include "materials/MaterialEOS.h"  // Base class for EOS models
#include "materials/PolyEOS.h"
#include "materials/TFDMatrices.h"  // For TFDMatrices definition
#include "utils/checksum_utils.h"
#include "utils/hdf5_utils.h"
#include "utils/string_utils.h"

std::string EquationOfStateV1::get_eos_relative_path_stub(
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

// Parses "XXXXX" from "eosXXXXX" part of a filename (without extension)
int EquationOfStateV1::parse_eos_id_from_filename_stub(
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

static int parse_complicated_eos_params(
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
    processed_line = EOSUtils::trim_string(processed_line);

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
      key_part = EOSUtils::trim_string(processed_line.substr(0, eq_pos));
      value_part = EOSUtils::trim_string(processed_line.substr(eq_pos + 1));

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
        if (!EOSUtils::string_to_double_fortran_compat(token, val_converted)) {
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

// --- Shim Implementations ---
int EquationOfStateV1::call_air_eos_2000_shim(double rho, double T, double& P,
                                              double& E) {
  int istat;
  // Outputs for dPdT etc. are not provided by this analytic EOS, so pass
  // dummies or handle if API changes
  c_air_eos_2000(rho, T, &P, &E, &istat);
  return istat;
}

int EquationOfStateV1::call_carbon_eos_2001_shim(double rho, double T,
                                                 double& P, double& E) {
  int istat;
  c_carbon_eos_2001(rho, T, &P, &E, &istat);
  return istat;
}

EquationOfStateV1::EquationOfStateV1()
    : tfd_data_(std::make_unique<EOS_Internal::TFDMatrices>()) {  // Default TFD
                                                                  // version
  std::cout << "C++: EquationOfStateV1 instance created." << std::endl;

  // Initialize TFD data struct
  tfd_data_->initialized = false;
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
    std::cerr << "C++ Error: Failed to set complicated_eos_use_old_cold_term. "
                 "Fortran istat: "
              << istat << std::endl;
  }
  return istat;
}

int EquationOfStateV1::setUseTFDDataVer1(bool value) {
  tfd_use_ver1_setting_ = value;
  int istat;
  c_set_use_tfd_data_ver1(value, &istat);
  if (istat != 0) {
    std::cerr << "C++ Error: Failed to set use_tfd_data_ver1. Fortran istat: "
              << istat << std::endl;
  }
  return istat;
}

int EquationOfStateV1::loadTFDDataInternal(const std::string& hdf5_filepath) {
  tfd_data_->initialized = false;  // Reset
  std::cout << "C++ (loadTFDDataInternal): Attempting to load TFD data from "
               "HDF5 file: "
            << hdf5_filepath << std::endl;

  // Suppress HDF5 default error printing, we'll handle errors via return codes
  H5E_auto2_t old_func;
  void* old_client_data;
  H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);  // Turn off HDF5 auto error printing

  hid_t file_id = H5Fopen(hdf5_filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    std::cerr << "HDF5 Error: Could not open TFD HDF5 file: " << hdf5_filepath
              << std::endl;
    H5Eset_auto2(H5E_DEFAULT, old_func,
                 old_client_data);    // Restore error printing
    return EOS_ERROR_FILE_NOT_FOUND;  // Or a specific HDF5 error code
  }

  // Read file-level attributes
  std::string creation_date, tfd_version_attr;
  if (HDF5Utils::read_hdf5_string_attribute(file_id, "file_creation_date",
                                            creation_date)) {
    std::cout << "HDF5 Attr: File Creation Date = " << creation_date
              << std::endl;
  }
  if (HDF5Utils::read_hdf5_string_attribute(file_id, "tfd_model_version",
                                            tfd_version_attr)) {
    std::cout << "HDF5 Attr: TFD Model Version = " << tfd_version_attr
              << std::endl;
  }

  // Read model description dataset
  if (HDF5Utils::read_hdf5_scalar_vlen_string_dataset(
          file_id, "/model_description",
          tfd_data_->tfd_file_description)) {  // Assuming tfd_file_description
                                               // is a member
    std::cout << "HDF5 Dataset: Model Description loaded (length "
              << tfd_data_->tfd_file_description.length() << ")." << std::endl;
    // For brevity, don't print long descriptions here.
  } else {
    std::cout
        << "HDF5 Info: No '/model_description' dataset found or failed to read."
        << std::endl;
  }

  herr_t status;
  hid_t common_dims_group_id =
      H5Gopen2(file_id, "/common_dimensions", H5P_DEFAULT);
  if (common_dims_group_id < 0) {
    std::cerr << "HDF5 Error: Could not open group '/common_dimensions' in "
              << hdf5_filepath << std::endl;
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_FILE_PARSE;  // Or HDF5 specific
  }

  status = read_hdf5_scalar_int(common_dims_group_id, "N1", tfd_data_->N1);
  if (status < 0) {
    H5Gclose(common_dims_group_id);
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_FILE_PARSE;
  }
  status = read_hdf5_scalar_int(common_dims_group_id, "N2", tfd_data_->N2);
  if (status < 0) {
    H5Gclose(common_dims_group_id);
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_FILE_PARSE;
  }
  H5Gclose(common_dims_group_id);

  if (tfd_data_->N1 <= 0 || tfd_data_->N2 <= 0) {
    std::cerr << "HDF5 Error: Invalid dimensions N1=" << tfd_data_->N1
              << ", N2=" << tfd_data_->N2 << " read from " << hdf5_filepath
              << std::endl;
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_INVALID_DIMENSIONS;
  }
  std::cout << "HDF5: Read TFD dimensions: " << tfd_data_->N1 << "x"
            << tfd_data_->N2 << std::endl;

  // Determine which TFD version group to read from
  // std::string tfd_version_group = (tfd_use_ver1_setting_ ? "/tfd_version_1" :
  // "/tfd_version_2"); For simplicity, let's assume datasets are at root or a
  // fixed path for now, and the file name itself (tfd_ver1.h5 vs tfd_ver2.h5)
  // distinguishes versions. Or, if the file contains multiple versions: hid_t
  // version_group_id = H5Gopen2(file_id, tfd_version_group.c_str(),
  // H5P_DEFAULT); if (version_group_id < 0) { /* error */ } status =
  // read_hdf5_dataset_double(version_group_id, "MatrixA", ...);
  // H5Gclose(version_group_id);

  // Assuming MatrixA and MatrixB are top-level datasets in the HDF5 file for
  // now
  status = read_hdf5_dataset_double(file_id, "/matrix_A", tfd_data_->N1,
                                    tfd_data_->N2, tfd_data_->matrix_A);
  if (status < 0) {
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_TFD_LOAD_A;
  }
  std::cout << "HDF5: Successfully read Matrix A." << std::endl;

  status = read_hdf5_dataset_double(file_id, "/matrix_B", tfd_data_->N1,
                                    tfd_data_->N2, tfd_data_->matrix_B);
  if (status < 0) {
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_TFD_LOAD_B;
  }
  std::cout << "HDF5: Successfully read Matrix B." << std::endl;

  H5Fclose(file_id);
  H5Eset_auto2(H5E_DEFAULT, old_func,
               old_client_data);  // Restore HDF5 auto error printing

  tfd_data_->initialized = true;
  // tfd_data.loaded_tfd_filepath = hdf5_filepath; // Store which file was
  // loaded
  std::cout << "C++: TFD data loaded successfully from HDF5 file: "
            << hdf5_filepath << std::endl;
  return EOS_SUCCESS;
}

// --- Initialization ---
int EquationOfStateV1::initialize(const std::vector<int>& eos_id_list,
                                  const std::string& eos_data_dir) {
  // 0. Ensure Fortran side knows which TFD version we prefer.
  //    The C++ tfd_use_ver1_setting_ is the master; set it if user calls the
  //    setter. Then ensure Fortran is synced.
  int f_istat;
  c_set_use_tfd_data_ver1(tfd_use_ver1_setting_, &f_istat);
  if (f_istat != 0) {
    std::cerr << "Error setting Fortran TFD version preference." << std::endl;
    // Decide if this is a fatal error for initialize
  }

  // 1. Load TFD data
  std::string tfd_filename_to_load =
      (tfd_use_ver1_setting_ ? "tfd_ver1.h5" : "tfd_ver2.h5");
  std::string full_tfd_hdf5_path =
      eos_data_dir + "/" + tfd_filename_to_load;  // Basic join

  if (perform_signature_check_)
    check_file_signature(tfd_use_ver1_setting_ ? -1000 : -2000,
                         full_tfd_hdf5_path);

  int tfd_load_status = loadTFDDataInternal(full_tfd_hdf5_path);
  if (tfd_load_status != EOS_SUCCESS) {
    std::cerr << "Error initializing: Failed to load TFD data. Code: "
              << tfd_load_status << std::endl;
    return tfd_load_status;
  }

  // 2. Load data for each EOS ID
  // 2. Load data for each EOS ID in the list
  for (int current_eos_id : eos_id_list) {
    std::string mat_group_subpath;  // e.g., "mat100"
    std::string file_stub_rel_path = get_eos_relative_path_stub(
        current_eos_id, mat_group_subpath);  // e.g. "mat100/eos10000"

    std::string dat_file_rel_path = file_stub_rel_path + ".dat";
    std::string h5_file_rel_path = file_stub_rel_path + ".h5";

    std::string full_param_file_path;
    std::string chosen_relative_param_path;
    std::string file_extension_found;

    std::string full_dat_path =
        EOSUtils::simple_path_join(eos_data_dir, dat_file_rel_path);
    std::string full_h5_path =
        EOSUtils::simple_path_join(eos_data_dir, h5_file_rel_path);

    std::ifstream test_dat_file(full_dat_path);
    if (test_dat_file.is_open()) {
      test_dat_file.close();
      full_param_file_path = full_dat_path;
      chosen_relative_param_path = dat_file_rel_path;
      file_extension_found = "dat";
      std::cout << "Found material parameter file: " << full_param_file_path
                << std::endl;
    } else {
      std::ifstream test_h5_file(full_h5_path);
      if (test_h5_file
              .is_open()) {  // Basic existence check, HDF5 lib will truly open
        test_h5_file.close();
        full_param_file_path = full_h5_path;
        chosen_relative_param_path = h5_file_rel_path;
        file_extension_found = "h5";
        std::cout << "Found material parameter file: " << full_param_file_path
                  << std::endl;
      }
    }

    if (full_param_file_path.empty()) {
      // Special check for hardcoded/registered analytic IDs that don't need a
      // file This logic was in previous version of initialize, may need to be
      // adapted if (analytic_eos_registry_.count(current_eos_id)) { ... handle
      // ... continue; } For now, assume AnalyticEOS will also have a (possibly
      // minimal) .dat file
      std::cerr << "Error: No parameter file (.dat or .h5) found for EOS ID "
                << current_eos_id << " using stub " << file_stub_rel_path
                << std::endl;
      // return EOS_ERROR_FILE_NOT_FOUND; // Or collect errors
      continue;  // Skip this ID
    }

    // Determine model_type_str
    std::string model_type_str;
    if (file_extension_found == "h5") {
      // Read model_type_str from a specific dataset in the HDF5 file
      hid_t pfile_id =
          H5Fopen(full_param_file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if (pfile_id < 0) { /* ... error ... */
        continue;
      }
      // Assume helper for reading string dataset like "/model_info/type_name"
      // or "/model_type" For now, using a placeholder name, adjust to your HDF5
      // structure. static bool read_hdf5_scalar_vlen_string_dataset(hid_t
      // file_id, const char* dset_name, std::string& out_str)
      if (!HDF5Utils::read_hdf5_scalar_vlen_string_dataset(
              pfile_id, "/model_info/type_name",
              model_type_str)) {  // Use EOSUtils if moved
        std::cerr << "Error: Could not read model type from HDF5 param file: "
                  << full_param_file_path << std::endl;
        H5Fclose(pfile_id);
        continue;
      }
      H5Fclose(pfile_id);
    } else {  // .dat file
      std::ifstream param_desc_file(full_param_file_path);
      std::string line;
      bool found_model_type = false;
      while (std::getline(param_desc_file, line)) {
        line = EOSUtils::trim_string(line);
        std::string prefix = "#MODEL_TYPE:";
        if (line.rfind(prefix, 0) == 0) {
          model_type_str = EOSUtils::trim_string(line.substr(prefix.length()));
          found_model_type = true;
          break;
        }
      }
      param_desc_file.close();
      if (!found_model_type) {
        std::cerr << "Error: #MODEL_TYPE: not found in text param file: "
                  << full_param_file_path << std::endl;
        continue;
      }
    }
    if (model_type_str.empty()) {
      std::cerr << "Error: Failed to determine model type for "
                << chosen_relative_param_path << std::endl;
      continue;
    }
    std::transform(model_type_str.begin(), model_type_str.end(),
                   model_type_str.begin(), ::tolower);

    // Factory logic
    std::unique_ptr<MaterialEOS> material_ptr = nullptr;
    std::string material_name = file_stub_rel_path;  // Use stub as a base name

    // Ensure model_type_str is normalized (e.g. lowercase) before comparison
    if (model_type_str == "complicated_legacy_text") {
      material_ptr =
          std::make_unique<ComplicatedLegacyEOS>(current_eos_id, material_name);
    } else if (model_type_str ==
               "polynomial_v1_hdf5") {  // This string must match what's in HDF5
      material_ptr = std::make_unique<PolyEOS>(current_eos_id, material_name);
    } else if (model_type_str ==
               "analytic_air_2000") {  // Match types for AnalyticEOS
      material_ptr = std::make_unique<AnalyticEOS>(
          current_eos_id, AnalyticEOS::AnalyticForm::AIR_2000, material_name);
    } else if (model_type_str == "analytic_carbon_2001") {
      material_ptr = std::make_unique<AnalyticEOS>(
          current_eos_id, AnalyticEOS::AnalyticForm::CARBON_2001,
          material_name);
    } else {
      std::cerr << "Error: Unknown model type '" << model_type_str
                << "' for EOS ID " << current_eos_id << " from file "
                << chosen_relative_param_path << std::endl;
      continue;
    }

    if (material_ptr) {
      bool model_needs_tfd =
          (dynamic_cast<ComplicatedLegacyEOS*>(material_ptr.get()) != nullptr);
      int mat_init_stat = material_ptr->initialize(
          full_param_file_path, eos_data_dir,
          (model_needs_tfd ? tfd_data_.get() : nullptr));
      if (mat_init_stat == EOS_SUCCESS) {
        loaded_materials_[current_eos_id] = std::move(material_ptr);
        std::cout << "Successfully initialized EOS ID " << current_eos_id
                  << " (Type: " << model_type_str << ")" << std::endl;
      } else {
        std::cerr << "Error: Failed to initialize MaterialEOS for ID "
                  << current_eos_id << " from file "
                  << chosen_relative_param_path << ". Status: " << mat_init_stat
                  << std::endl;
        // return EOS_ERROR_MATERIAL_INIT_FAILED; // Or collect errors
      }
    }
  }  // end for eos_id_list

  return EOS_SUCCESS;
}

int EquationOfStateV1::check_eos_data_dir(
    const std::string& eos_data_dir, const std::vector<int>& eos_ids_to_check) {
  // Check TFD files
  std::string tfd_files_to_check[2] = {"tfd_ver1.h5", "tfd_ver2.h5"};

  for (const std::string& tfd_suffix : tfd_files_to_check) {
    std::string tfd_path_str = EOSUtils::simple_path_join(eos_data_dir, tfd_suffix);
    std::ifstream test_file(tfd_path_str);
    if (!test_file.is_open()) {
      std::cerr << "Check Error: TFD file not found or not accessible: "
                << tfd_path_str << std::endl;
      return EOS_ERROR_FILE_NOT_FOUND;
    }
    test_file.close();
  }

  // Check specific EOS ID files
  for (int id : eos_ids_to_check) {
    if (analytic_eos_registry_.count(id)) {
      std::cout << "Check Info: Analytic EOS ID " << id
                << " is registered, file existence check skipped/optional."
                << std::endl;
      continue;
    }
    std::string file_path_str = EOSUtils::get_eos_file_path(eos_data_dir, id);
    if (file_path_str.empty()) {
      std::cerr << "Check Error: Could not construct path for EOS ID " << id
                << std::endl;
      return EOS_ERROR_FILE_NOT_FOUND;  // Or a different error code
    }
    std::ifstream eos_file(file_path_str);
    if (!eos_file.is_open()) {
      std::cerr << "Check Error: EOS data file not found for ID " << id
                << " at " << file_path_str << std::endl;
      return EOS_ERROR_FILE_NOT_FOUND;
    }
    eos_file.close();
  }
  std::cout << "check_eos_data_dir: All checked files appear accessible."
            << std::endl;
  return EOS_SUCCESS;
}

// --- Computation ---
int EquationOfStateV1::compute(int eos_id, double rho, double T, double& P,
                               double& E, double& dPdT, double& dEdT,
                               double& dPdrho) {
  auto it = loaded_materials_.find(eos_id);
  if (it == loaded_materials_.end()) {
    std::cerr << "Error: EOS ID " << eos_id << " not initialized." << std::endl;
    return EOS_ERROR_UNKNOWN_EOS_ID;
  }
  P = 0.0;
  E = 0.0;
  dPdT = 0.0;
  dEdT = 0.0;
  dPdrho = 0.0;  // Initialize outputs
  return it->second->compute(rho, T, P, E, dPdT, dEdT, dPdrho);
}

void EquationOfStateV1::free_resources() {
  loaded_materials_
      .clear();  // unique_ptr handles deletion of MaterialEOS objects
  if (tfd_data_) {
    tfd_data_->clear();  // Clear data within TFDMatrices
  }
  std::cout << "C++: EquationOfStateV1 resources freed." << std::endl;
}

// Helper for packing data into buffer and advancing offset
template <typename T>
void pack_item(char*& current_ptr, const T& item) {
  std::memcpy(current_ptr, &item, sizeof(T));
  current_ptr += sizeof(T);
}

template <typename T>
void pack_vector_data(char*& current_ptr, const std::vector<T>& vec_data) {
  if (!vec_data.empty()) {
    std::memcpy(current_ptr, vec_data.data(), vec_data.size() * sizeof(T));
    current_ptr += vec_data.size() * sizeof(T);
  }
}

// Helper for unpacking data from buffer and advancing offset
template <typename T>
void unpack_item(const char*& current_ptr, T& item) {
  std::memcpy(&item, current_ptr, sizeof(T));
  current_ptr += sizeof(T);
}

template <typename T>
void unpack_vector_data(const char*& current_ptr, std::vector<T>& vec_data,
                        size_t num_elements) {
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

  required_size += sizeof(tfd_data_->N1);
  required_size += sizeof(tfd_data_->N2);
  if (tfd_data_->initialized) {
    required_size += tfd_data_->matrix_A.size() * sizeof(double);
    required_size += tfd_data_->matrix_B.size() * sizeof(double);
  } else {
    // If TFD not initialized, we still pack N1=0, N2=0. No matrix data.
  }

  required_size += sizeof(int);  // For num_loaded_materials_
  // for (const auto& pair : loaded_materials_) {
  //   const MaterialData& md = pair.second;
  //   required_size += sizeof(md.eos_id_key);
  //   required_size += sizeof(InternalEOSType);  // Pack the enum value
  //   required_size += sizeof(int);              // For num_params
  //   required_size += md.params.size() * sizeof(double);
  // }
  //
  // // --- Mode 1: Calculate size only ---
  // if (buffer == nullptr) {
  //   buffer_size = required_size;
  //   return EOS_SUCCESS;
  // }
  //
  // // --- Mode 2: Pack data if buffer is sufficient ---
  // if (buffer_size < required_size) {
  //   std::cerr << "Pack Error: Provided buffer too small. Need " << required_size
  //             << ", got " << buffer_size << std::endl;
  //   return -1;  // Or a specific error code for buffer too small
  // }
  //
  // char* current_ptr = buffer;
  //
  // pack_item(current_ptr, tfd_use_ver1_setting_);
  // pack_item(current_ptr, complicated_eos_use_old_cold_term_setting_);
  //
  // if (tfd_data_->initialized) {
  //   pack_item(current_ptr, tfd_data_->N1);
  //   pack_item(current_ptr, tfd_data_->N2);
  //   pack_vector_data(current_ptr, tfd_data_->matrix_A);
  //   pack_vector_data(current_ptr, tfd_data_->matrix_B);
  // } else {
  //   int zero_dim = 0;
  //   pack_item(current_ptr, zero_dim);  // N1 = 0
  //   pack_item(current_ptr, zero_dim);  // N2 = 0
  //                                      // No matrix data to pack
  // }
  //
  // int num_loaded = static_cast<int>(loaded_materials_.size());
  // pack_item(current_ptr, num_loaded);
  //
  // for (const auto& pair : loaded_materials_) {
  //   const MaterialData& md = pair.second;
  //   pack_item(current_ptr, md.eos_id_key);
  //   pack_item(current_ptr, md.internal_type);  // Packing enum directly
  //
  //   int num_params = static_cast<int>(md.params.size());
  //   pack_item(current_ptr, num_params);
  //   pack_vector_data(current_ptr, md.params);
  // }
  //
  // buffer_size = static_cast<int>(current_ptr - buffer);  // Actual size used
  // if (buffer_size != required_size) {
  //   // This should ideally not happen if logic is correct
  //   std::cerr << "Pack Error: Mismatch in calculated vs packed size. Packed: "
  //             << buffer_size << ", Required: " << required_size << std::endl;
  //   return -2;  // Packing logic error
  // }
  //
  // std::cout << "C++ (pack_data): Successfully packed " << buffer_size
  //           << " bytes." << std::endl;
  return EOS_SUCCESS;
}

int EquationOfStateV1::unpack_data(const char* buffer, int buffer_size) {
  // Clear existing state before unpacking
  free_resources();  // Important to reset state

  const char* current_ptr = buffer;
  const char* buffer_end = buffer + buffer_size;  // For boundary checks

  // Helper lambda for boundary check
  auto check_boundary = [&](size_t needed) -> bool {
    if (current_ptr + needed > buffer_end) {
      std::cerr << "Unpack Error: Buffer too small or data corrupted. "
                << "Attempting to read " << needed << " bytes with "
                << (buffer_end - current_ptr) << " bytes remaining."
                << std::endl;
      return false;
    }
    return true;
  };

  if (!check_boundary(sizeof(tfd_use_ver1_setting_))) return -1;
  unpack_item(current_ptr, tfd_use_ver1_setting_);

  if (!check_boundary(sizeof(complicated_eos_use_old_cold_term_setting_)))
    return -1;
  unpack_item(current_ptr, complicated_eos_use_old_cold_term_setting_);

  // Update Fortran side with unpacked control variables
  int f_istat;
  c_set_use_tfd_data_ver1(tfd_use_ver1_setting_,
                          &f_istat);  // Error check f_istat?
  c_set_complicated_eos_use_old_cold_term(
      complicated_eos_use_old_cold_term_setting_, &f_istat);  // Error check?

  std::cout << "C++ (unpack_data): Unpacked tfd_use_ver1_setting="
            << tfd_use_ver1_setting_
            << ", complicated_eos_use_old_cold_term_setting="
            << complicated_eos_use_old_cold_term_setting_ << std::endl;

  if (!check_boundary(sizeof(tfd_data_->N1))) return -1;
  unpack_item(current_ptr, tfd_data_->N1);
  if (!check_boundary(sizeof(tfd_data_->N2))) return -1;
  unpack_item(current_ptr, tfd_data_->N2);

  if (tfd_data_->N1 > 0 &&
      tfd_data_->N2 > 0) {  // Check if TFD data was actually packed
    size_t tfd_a_bytes =
        static_cast<size_t>(tfd_data_->N1) * tfd_data_->N2 * sizeof(double);
    if (!check_boundary(tfd_a_bytes)) return -1;
    unpack_vector_data(current_ptr, tfd_data_->matrix_A,
                       static_cast<size_t>(tfd_data_->N1) * tfd_data_->N2);

    size_t tfd_b_bytes =
        static_cast<size_t>(tfd_data_->N1) * tfd_data_->N2 * sizeof(double);
    if (!check_boundary(tfd_b_bytes)) return -1;
    unpack_vector_data(current_ptr, tfd_data_->matrix_B,
                       static_cast<size_t>(tfd_data_->N1) * tfd_data_->N2);
    tfd_data_->initialized = true;
    std::cout << "C++ (unpack_data): Unpacked TFD matrices (" << tfd_data_->N1
              << "x" << tfd_data_->N2 << ")." << std::endl;
  } else {
    tfd_data_->initialized = false;  // N1 or N2 was 0, so no TFD data
    std::cout << "C++ (unpack_data): TFD dimensions were zero, no TFD matrices "
                 "unpacked."
              << std::endl;
  }

  // int num_loaded_materials;
  // if (!check_boundary(sizeof(num_loaded_materials))) return -1;
  // unpack_item(current_ptr, num_loaded_materials);
  // std::cout << "C++ (unpack_data): Expecting " << num_loaded_materials
  //           << " materials." << std::endl;
  //
  // for (int i = 0; i < num_loaded_materials; ++i) {
  //   MaterialData md;
  //   if (!check_boundary(sizeof(md.eos_id_key))) return -1;
  //   unpack_item(current_ptr, md.eos_id_key);
  //
  //   if (!check_boundary(sizeof(md.internal_type))) return -1;
  //   unpack_item(current_ptr, md.internal_type);  // Unpacking enum directly
  //
  //   int num_params;
  //   if (!check_boundary(sizeof(num_params))) return -1;
  //   unpack_item(current_ptr, num_params);
  //
  //   if (num_params > 0) {
  //     size_t params_bytes = static_cast<size_t>(num_params) * sizeof(double);
  //     if (!check_boundary(params_bytes)) return -1;
  //     unpack_vector_data(current_ptr, md.params, num_params);
  //   }
  //
  //   // Re-establish analytic_func pointer based on internal_type and eos_id_key
  //   // This relies on analytic_eos_registry_ being populated (e.g. in
  //   // constructor) Or, find it in the registry by the unpacked md.eos_id_key.
  //   auto reg_it = analytic_eos_registry_.find(md.eos_id_key);
  //   if (reg_it != analytic_eos_registry_.end()) {
  //     // Cross-check if internal_type matches what's in registry for this ID
  //     if (reg_it->second.internal_type == md.internal_type) {
  //       md.analytic_func = reg_it->second.func_ptr;
  //     } else {
  //       std::cerr << "Unpack Warning: Mismatch in unpacked internal_type for "
  //                    "analytic EOS ID "
  //                 << md.eos_id_key << ". Using registry's type." << std::endl;
  //       md.internal_type = reg_it->second.internal_type;  // Prefer registry
  //       md.analytic_func = reg_it->second.func_ptr;
  //     }
  //   } else if (md.internal_type == InternalEOSType::ANALYTIC_AIR ||
  //              md.internal_type == InternalEOSType::ANALYTIC_CARBON) {
  //     // This case should ideally not happen if eos_id_key was used for registry
  //     // lookup
  //     std::cerr << "Unpack Error: Unpacked analytic EOS ID " << md.eos_id_key
  //               << " not found in registry, but type suggests it's analytic."
  //               << std::endl;
  //     return -1;  // Data integrity issue
  //   }
  //
  //   // Set type_from_file based on internal_type (approximate reverse)
  //   if (md.internal_type == InternalEOSType::ANALYTIC_AIR ||
  //       md.internal_type == InternalEOSType::ANALYTIC_CARBON) {
  //     md.type_from_file = EOSTypeFromFile::ANALYTIC;
  //   } else if (md.internal_type == InternalEOSType::COMPLICATED_V1) {
  //     md.type_from_file = EOSTypeFromFile::COMPLICATED;
  //   } else {
  //     md.type_from_file = EOSTypeFromFile::NOT_SET;  // Or UNKNOWN
  //   }
  //
  //   loaded_materials_[md.eos_id_key] = md;
  //   std::cout << "C++ (unpack_data): Unpacked material ID " << md.eos_id_key
  //             << ", type " << static_cast<int>(md.internal_type)
  //             << ", num_params " << num_params << std::endl;
  // }
  //
  // // Final check: did we consume the whole buffer?
  // if (current_ptr != buffer_end) {
  //   std::cerr
  //       << "Unpack Warning: " << (buffer_end - current_ptr)
  //       << " bytes remaining in buffer after unpack. Expected to consume all "
  //       << buffer_size << " bytes." << std::endl;
  //   // This might be an error depending on strictness.
  // }

  std::cout << "C++ (unpack_data): Successfully unpacked data." << std::endl;
  return EOS_SUCCESS;
}

bool EquationOfStateV1::get_perform_signature_check() const {
  return perform_signature_check_;
}

void EquationOfStateV1::set_perform_signature_check(bool enable) {
  perform_signature_check_ = enable;
}

void EquationOfStateV1::check_file_signature(
    const int eos_id, const std::string& full_filepath) const {
  if (!perform_signature_check_) {
    return;
  }

  auto it = official_file_checksums.find(eos_id);
  if (it == official_file_checksums.end()) {
    std::cout << "Signature Check Info: No official checksum registered for '"
              << eos_id << "'. Skipping check." << std::endl;
    return;
  }
  const std::string& official_md5 = it->second;

  std::string calculated_md5 =
      EOSCheckSumUtils::calculate_file_md5(full_filepath);

  if (calculated_md5.empty()) {
    std::cerr << "Signature Check Warning: Could not calculate MD5 for '"
              << full_filepath << "'. Checksum validation skipped."
              << std::endl;
  }

  if (calculated_md5 != official_md5) {
    std::cerr << "Signature Check Warning: MD5 for '" << full_filepath
              << " does not match the database!" << std::endl;
  }
}
