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
#include "eos_error_codes.h"
#include "materials/AnalyticEOS.h"
#include "materials/ComplicatedLegacyEOS.h"
#include "materials/MaterialEOS.h"  // Base class for EOS models
#include "materials/PolyEOS.h"
#include "materials/TFDMatrices.h"  // For TFDMatrices definition
#include "utils/checksum_utils.h"
#include "utils/hdf5_utils.h"
#include "utils/string_utils.h"

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

  status = HDF5Utils::read_hdf5_scalar_int(common_dims_group_id, "N1",
                                           tfd_data_->N1);
  if (status < 0) {
    H5Gclose(common_dims_group_id);
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_FILE_PARSE;
  }
  status = HDF5Utils::read_hdf5_scalar_int(common_dims_group_id, "N2",
                                           tfd_data_->N2);
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
  status = HDF5Utils::read_hdf5_dataset_double(
      file_id, "/matrix_A", tfd_data_->N1, tfd_data_->N2, tfd_data_->matrix_A);
  if (status < 0) {
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_TFD_LOAD_A;
  }
  std::cout << "HDF5: Successfully read Matrix A." << std::endl;

  status = HDF5Utils::read_hdf5_dataset_double(
      file_id, "/matrix_B", tfd_data_->N1, tfd_data_->N2, tfd_data_->matrix_B);
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
  free_resources();  // Clear previous state, also re-creates tfd_data_
                     // unique_ptr

  // 0. Sync Fortran control variables (call C-API setters)
  int f_istat_dummy;  // We might not strictly need to check status here if
                      // setters are simple
  c_set_use_tfd_data_ver1(tfd_use_ver1_setting_, &f_istat_dummy);
  c_set_complicated_eos_use_old_cold_term(
      complicated_eos_use_old_cold_term_setting_, &f_istat_dummy);

  // 1. Load TFD data (HDF5 based)
  std::string tfd_relative_filename =
      (tfd_use_ver1_setting_ ? "tfd_ver1.h5" : "tfd_ver2.h5");
  std::string full_tfd_filepath =
      EOStringUtils::simple_path_join(eos_data_dir, tfd_relative_filename);

  int tfd_load_stat = loadTFDDataInternal(full_tfd_filepath);
  if (tfd_load_stat != EOS_SUCCESS) {
    std::cerr << "Error (EquationOfStateV1::initialize): Failed to load TFD "
                 "data from "
              << full_tfd_filepath << ". Code: " << tfd_load_stat << std::endl;
    return tfd_load_stat;  // Critical failure
  }

  // 2. Load data for each EOS ID in the list
  for (int current_eos_id : eos_id_list) {
    std::string mat_group_subpath;  // e.g., "mat100"
    std::string file_stub_rel_path = EOStringUtils::get_eos_relative_path_stub(
        current_eos_id, mat_group_subpath);
    if (file_stub_rel_path.find("unknown") != std::string::npos ||
        file_stub_rel_path.empty()) {  // Error from get_eos_relative_path_stub
      std::cerr << "Error (EquationOfStateV1::initialize): Could not form "
                   "valid path stub for EOS ID "
                << current_eos_id << std::endl;
      continue;  // Skip this problematic ID
    }

    std::string dat_file_rel_path = file_stub_rel_path + ".dat";
    std::string h5_file_rel_path = file_stub_rel_path + ".h5";

    std::string full_param_file_path_found;
    std::string chosen_relative_param_path;
    std::string model_type_str;

    std::string full_dat_path =
        EOStringUtils::simple_path_join(eos_data_dir, dat_file_rel_path);
    std::string full_h5_path =
        EOStringUtils::simple_path_join(eos_data_dir, h5_file_rel_path);

    std::ifstream test_open_dat(full_dat_path);
    if (test_open_dat.is_open()) {
      test_open_dat.close();
      full_param_file_path_found = full_dat_path;
      chosen_relative_param_path = dat_file_rel_path;

      // Read #MODEL_TYPE from .dat file
      std::ifstream desc_file(full_param_file_path_found);
      std::string line;
      bool found_mt = false;
      while (std::getline(desc_file, line)) {
        line = EOStringUtils::trim_string(line);
        std::string prefix = "#MODEL_TYPE:";
        if (line.rfind(prefix, 0) == 0) {
          model_type_str =
              EOStringUtils::trim_string(line.substr(prefix.length()));
          found_mt = true;
          break;
        }
      }
      if (!found_mt) {
        std::cerr << "Error: #MODEL_TYPE not found in "
                  << full_param_file_path_found << std::endl;
        continue;  // Skip this file
      }
    } else {
      // Try .h5 file (basic existence check, HDF5 lib will truly open it later)
      std::ifstream test_open_h5(
          full_h5_path, std::ios::binary);  // Open H5 for existence check
      if (test_open_h5.is_open()) {
        test_open_h5.close();
        full_param_file_path_found = full_h5_path;
        chosen_relative_param_path = h5_file_rel_path;

        // Read model_type_str from HDF5 file (e.g., dataset
        // "/model_info/type_name") This requires opening the HDF5 file just for
        // this.
        hid_t pfile_id = H5Fopen(full_param_file_path_found.c_str(),
                                 H5F_ACC_RDONLY, H5P_DEFAULT);
        if (pfile_id < 0) {
          std::cerr << "Error: Could not open HDF5 param file for type read: "
                    << full_param_file_path_found << std::endl;
          continue;
        }
        // Ensure EOSUtils::read_hdf5_scalar_vlen_string_dataset is correctly
        // namespaced and available
        if (!HDF5Utils::read_hdf5_scalar_vlen_string_dataset(
                pfile_id, "/model_info/type_name", model_type_str)) {
          std::cerr << "Error: Could not read /model_info/type_name from HDF5 "
                       "param file: "
                    << full_param_file_path_found << std::endl;
          H5Fclose(pfile_id);
          continue;
        }
        H5Fclose(pfile_id);
      }
    }

    if (full_param_file_path_found.empty()) {
      std::cerr << "Error: No parameter file (.dat or .h5) found for EOS ID "
                << current_eos_id << " (stub: " << file_stub_rel_path << ")"
                << std::endl;
      continue;
    }
    if (model_type_str.empty()) {
      std::cerr << "Error: Could not determine model type for EOS ID "
                << current_eos_id << " from file " << chosen_relative_param_path
                << std::endl;
      continue;
    }
    std::transform(model_type_str.begin(), model_type_str.end(),
                   model_type_str.begin(), ::tolower);

    // Factory logic
    std::unique_ptr<MaterialEOS> material_ptr = nullptr;
    std::string material_name =
        chosen_relative_param_path;  // Use relative path as a base name

    if (model_type_str == "complicated_legacy_text") {
      material_ptr =
          std::make_unique<ComplicatedLegacyEOS>(current_eos_id, material_name);
    } else if (model_type_str == "polynomial_v1_hdf5") {
      material_ptr = std::make_unique<PolyEOS>(current_eos_id, material_name);
    } else if (model_type_str == "analytic_air_2000") {
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
      // Determine if this model type inherently needs TFD
      bool model_needs_tfd =
          (dynamic_cast<ComplicatedLegacyEOS*>(material_ptr.get()) != nullptr);

      const EOS_Internal::TFDMatrices* tfd_to_pass = nullptr;
      if (model_needs_tfd) {
        if (!tfd_data_ || !tfd_data_->initialized) {
          std::cerr << "Error: TFD data required by model type '"
                    << model_type_str << "' for EOS ID " << current_eos_id
                    << " but TFD is not initialized." << std::endl;
          continue;
        }
        tfd_to_pass = tfd_data_.get();
      }

      int mat_init_stat = material_ptr->initialize(full_param_file_path_found,
                                                   eos_data_dir, tfd_to_pass);
      if (mat_init_stat ==
          EOS_SUCCESS) {  // Assuming derived init returns EOS_SUCCESS
        loaded_materials_[current_eos_id] = std::move(material_ptr);
        std::cout << "Successfully initialized EOS ID " << current_eos_id
                  << " (Type: " << model_type_str << ")" << std::endl;
      } else {
        std::cerr << "Error: Failed to initialize MaterialEOS for ID "
                  << current_eos_id << ". File: " << chosen_relative_param_path
                  << ". Status: " << mat_init_stat << std::endl;
      }
    }
  }  // end for eos_id_list

  if (loaded_materials_.empty() && !eos_id_list.empty()) {
    std::cerr << "Warning: No EOS models were successfully loaded from the "
                 "provided list."
              << std::endl;
    // Depending on policy, this could be an error.
  }
  return EOS_SUCCESS;
}

int EquationOfStateV1::check_eos_data_dir(
    const std::string& eos_data_dir_root,
    const std::vector<int>& eos_ids_to_check) {
  // Check TFD HDF5 files
  std::string tfd_hdf5_files_to_check[2] = {"tfd_ver1.h5", "tfd_ver2.h5"};
  for (const std::string& tfd_suffix : tfd_hdf5_files_to_check) {
    std::string tfd_path_str =
        EOStringUtils::simple_path_join(eos_data_dir_root, tfd_suffix);
    std::ifstream test_file(tfd_path_str, std::ios::binary);
    if (!test_file.is_open()) { /* ... error ... */
      return EOS_ERROR_FILE_NOT_FOUND;
    }
    test_file.close();
  }

  // Check specific EOS ID files (probe for .dat then .h5)
  for (int id : eos_ids_to_check) {
    std::string mat_group_subpath;
    std::string file_stub_rel =
        EOStringUtils::get_eos_relative_path_stub(id, mat_group_subpath);
    if (file_stub_rel.find("unknown") != std::string::npos ||
        file_stub_rel.empty())
      continue;

    std::string dat_path_str = EOStringUtils::simple_path_join(
        eos_data_dir_root, file_stub_rel + ".dat");
    std::string h5_path_str = EOStringUtils::simple_path_join(
        eos_data_dir_root, file_stub_rel + ".h5");

    std::ifstream test_dat(dat_path_str);
    std::ifstream test_h5(h5_path_str, std::ios::binary);

    if (!test_dat.is_open() && !test_h5.is_open()) {
      std::cerr << "Check Error: EOS param file (.dat or .h5) not found for ID "
                << id << " (stub: " << file_stub_rel << ")" << std::endl;
      return EOS_ERROR_FILE_NOT_FOUND;
    }
    if (test_dat.is_open()) test_dat.close();
    if (test_h5.is_open()) test_h5.close();
  }
  std::cout << "check_eos_data_dir: All checked TFD files found and param "
               "files (.dat or .h5) exist."
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
  //   std::cerr << "Pack Error: Provided buffer too small. Need " <<
  //   required_size
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
  //   std::cerr << "Pack Error: Mismatch in calculated vs packed size. Packed:
  //   "
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
  //   // Re-establish analytic_func pointer based on internal_type and
  //   eos_id_key
  //   // This relies on analytic_eos_registry_ being populated (e.g. in
  //   // constructor) Or, find it in the registry by the unpacked
  //   md.eos_id_key. auto reg_it = analytic_eos_registry_.find(md.eos_id_key);
  //   if (reg_it != analytic_eos_registry_.end()) {
  //     // Cross-check if internal_type matches what's in registry for this ID
  //     if (reg_it->second.internal_type == md.internal_type) {
  //       md.analytic_func = reg_it->second.func_ptr;
  //     } else {
  //       std::cerr << "Unpack Warning: Mismatch in unpacked internal_type for
  //       "
  //                    "analytic EOS ID "
  //                 << md.eos_id_key << ". Using registry's type." <<
  //                 std::endl;
  //       md.internal_type = reg_it->second.internal_type;  // Prefer registry
  //       md.analytic_func = reg_it->second.func_ptr;
  //     }
  //   } else if (md.internal_type == InternalEOSType::ANALYTIC_AIR ||
  //              md.internal_type == InternalEOSType::ANALYTIC_CARBON) {
  //     // This case should ideally not happen if eos_id_key was used for
  //     registry
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
  //       << " bytes remaining in buffer after unpack. Expected to consume all
  //       "
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
