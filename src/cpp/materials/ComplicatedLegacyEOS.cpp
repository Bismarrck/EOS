#include "ComplicatedLegacyEOS.h"

#include "TFDMatrices.h"  // For EOS_Internal::TFDMatrices
#include "utils/string_utils.h"  // For EOSUtils::parse_complicated_eos_params if moved there
// or keep parse_complicated_eos_params as static helper here.
#include <fstream>
#include <iostream>
#include <sstream>  // For pack/unpack with stringstream

#include "eos_error_codes.h"
#include "utils/parser_utils.h"

// Forward declare the C-wrapper from fortran_c_wrappers.f90
// This should match the BIND(C, name='...') in Fortran
extern "C" {
void c_complicated_eos(const double* params, int num_params, double rho,
                       double T, const double* matrix_a_flat, int n1_a,
                       int n2_a, const double* matrix_b_flat, int n1_b,
                       int n2_b, double* P, double* E, double* dPdT,
                       double* dEdT, double* dPdrho, int* istat_out);
}

// Define the static member for parameter order
// This order MUST match what your Complicated_EOS Fortran routine expects.
const std::vector<std::string> ComplicatedLegacyEOS::s_param_order_ = {
    "key1", "key2", "key3",
    "key4", "key5", "key6"  // Populate with your actual keys in order
};

// If parse_complicated_eos_params is specific to this class, it can be a
// private static member. For now, assuming it's a free function in EOSUtils or
// defined here. Let's put a simplified version here for now (it was defined in
// EquationOfStateV1.cpp previously)
namespace EOSUtils { /* Assume parse_complicated_eos_params is available from
                        string_utils.h/cpp */
}

ComplicatedLegacyEOS::ComplicatedLegacyEOS(const int eos_id,
                                           const std::string& name)
    : MaterialEOS(eos_id, ModelCategory::COMPLICATED_TABLE_TFD, name),
      tfd_data_ptr_(nullptr) {}

int ComplicatedLegacyEOS::initialize(
    const std::string& material_param_filepath,
    const std::string& eos_data_dir_root,  // Not directly used if path is full
    const EOS_Internal::TFDMatrices* tfd_data) {
  std::cout << "Initializing ComplicatedLegacyEOS for ID " << eos_id_
            << " from " << material_param_filepath << std::endl;

  if (!tfd_data || !tfd_data->initialized) {
    std::cerr << "Error (ComplicatedLegacyEOS::initialize): TFD data is "
                 "required but not provided or not initialized."
              << std::endl;
    return EOS_ERROR_COMPLICATED_EOS_INIT_FAILED;  // Or a more specific error
  }
  tfd_data_ptr_ = tfd_data;

  std::ifstream param_file(material_param_filepath);
  if (!param_file.is_open()) {
    std::cerr << "Error (ComplicatedLegacyEOS::initialize): Could not open "
                 "parameter file: "
              << material_param_filepath << std::endl;
    return EOS_ERROR_COMPLICATED_EOS_INIT_FAILED;
  }

  // Skip #MODEL_TYPE line, as it was already read by EquationOfStateV1
  std::string line;
  bool model_type_skipped = false;
  std::streampos original_pos = param_file.tellg();
  while (std::getline(param_file, line)) {
    std::string temp_line = EOStringUtils::trim_string(line);
    if (temp_line.rfind("#MODEL_TYPE:", 0) == 0) {
      model_type_skipped = true;
      original_pos = param_file.tellg();  // Mark position after #MODEL_TYPE
      break;
    }
  }
  if (!model_type_skipped) {
    param_file.clear();  // Clear EOF flags if any
    param_file.seekg(
        0, std::ios::beg);  // Rewind to start if #MODEL_TYPE not found
    std::cout << "Warning (ComplicatedLegacyEOS::initialize): #MODEL_TYPE line "
                 "not found/skipped in "
              << material_param_filepath << ". Parsing from beginning."
              << std::endl;
  } else {
    param_file.seekg(original_pos);  // Seek to position after #MODEL_TYPE
  }

  // Use the robust key-value parser (ensure it's accessible, e.g., from
  // EOSUtils) It needs the ifstream, filepath for errors, the order, and output
  // vector.
  int parse_stat = EOSParserUtils::parse_complicated_eos_params(
      param_file, material_param_filepath, s_param_order_, params_);
  param_file.close();

  if (parse_stat != EOS_SUCCESS) {  // Assuming parse_complicated_eos_params
                                    // returns EOS_SUCCESS (0)
    std::cerr << "Error (ComplicatedLegacyEOS::initialize): Failed to parse "
                 "parameters from "
              << material_param_filepath << ". Status: " << parse_stat
              << std::endl;
    params_.clear();  // Ensure params are empty on failure
    return EOS_ERROR_COMPLICATED_EOS_INIT_FAILED;
  }

  std::cout << "ComplicatedLegacyEOS ID " << eos_id_ << " initialized with "
            << params_.size() << " parameters." << std::endl;
  return EOS_SUCCESS;
}

ComputeResult ComplicatedLegacyEOS::compute(double rho, double T) const {
  ComputeResult result;
  if (params_.empty()) {
    std::cerr << "Error (ComplicatedLegacyEOS::compute): Parameters not loaded "
                 "for EOS ID "
              << eos_id_ << std::endl;
    result.istat = EOS_ERROR_COMPLICATED_EOS_EMPTY_PARAMS;
    return result;
  }
  // if (!tfd_data_ptr_ || !tfd_data_ptr_->initialized) {
  //     std::cerr << "Error (ComplicatedLegacyEOS::compute): TFD data not
  //     available for EOS ID " << eos_id_ << std::endl; return
  //     MAT_EOS_COMPUTE_FAILED;
  // }

  c_complicated_eos(params_.data(), static_cast<int>(params_.size()), rho, T,
                    tfd_data_ptr_->matrix_A.data(), tfd_data_ptr_->N1,
                    tfd_data_ptr_->N2, tfd_data_ptr_->matrix_B.data(),
                    tfd_data_ptr_->N1,
                    tfd_data_ptr_->N2,  // Assuming N1/N2 same for A & B
                    &result.P, &result.E, &result.dPdT, &result.dEdT,
                    &result.dPdrho, &result.istat);

  return result;
}

int ComplicatedLegacyEOS::pack_parameters(std::ostream& os) const {
  size_t num_params = params_.size();
  os.write(reinterpret_cast<const char*>(&num_params), sizeof(num_params));
  if (!os.good()) return EOS_ERROR_PACK_FAILED;

  if (num_params > 0) {
    os.write(reinterpret_cast<const char*>(params_.data()),
             num_params * sizeof(double));
    if (!os.good()) return EOS_ERROR_PACK_FAILED;
  }
  return EOS_SUCCESS;
}

int ComplicatedLegacyEOS::unpack_parameters(std::istream& is) {
  params_.clear();
  size_t num_params = 0;
  is.read(reinterpret_cast<char*>(&num_params), sizeof(num_params));
  if (!is.good() && num_params > 0) {  // If expect params but read fails
    std::cerr
        << "Error (ComplicatedLegacyEOS::unpack): Failed to read num_params."
        << std::endl;
    return EOS_ERROR_UNPACK_FAILED;
  }
  if (num_params == 0 && !is.good() &&
      !is.eof()) {  // Read error but num_params could be 0
    std::cerr << "Error (ComplicatedLegacyEOS::unpack): Stream error after "
                 "reading num_params=0."
              << std::endl;
    return EOS_ERROR_UNPACK_FAILED;
  }

  if (num_params > 0) {
    try {
      params_.resize(num_params);  // Can throw std::bad_alloc
    } catch (const std::bad_alloc& e) {
      std::cerr << "Error (ComplicatedLegacyEOS::unpack): Failed to resize "
                   "params vector: "
                << e.what() << std::endl;
      return EOS_ERROR_UNPACK_FAILED;  // Or specific memory error
    }
    is.read(reinterpret_cast<char*>(params_.data()),
            num_params * sizeof(double));
    if (!is.good()) {
      std::cerr << "Error (ComplicatedLegacyEOS::unpack): Failed to read "
                << num_params << " parameters." << std::endl;
      params_.clear();  // Clear partially read data
      return EOS_ERROR_UNPACK_FAILED;
    }
  }
  // tfd_data_ptr_ will be set by EquationOfStateV1 after this call if needed
  return EOS_SUCCESS;
}
