#include "AnalyticEOS.h"

#include <iostream>
#include <sstream>
#include "eos_error_codes.h"

// Declare C-wrappers for analytic Fortran routines
extern "C" {
void c_air_eos_2000(double rho, double T, double* P, double* E, int* istat_out);
void c_carbon_eos_2001(double rho, double T, double* P, double* E,
                       int* istat_out);
// Add declarations for other analytic C-wrappers if you have them
}

AnalyticEOS::AnalyticEOS(int eos_id, AnalyticForm form, const std::string& name)
    : MaterialEOS(eos_id, ModelCategory::ANALYTIC, name), form_(form) {
  std::cout << "AnalyticEOS created for ID " << eos_id_ << " with form "
            << static_cast<int>(form_) << std::endl;
}

int AnalyticEOS::initialize(const std::string& material_param_filepath,
                            const std::string& eos_data_dir_root,
                            const EOS_Internal::TFDMatrices* tfd_data) {
  std::cout << "Initializing AnalyticEOS ID " << eos_id_
            << " (Form: " << static_cast<int>(form_) << ")" << std::endl;
  // material_param_filepath might be a .dat file with just #MODEL_TYPE:
  // ANALYTIC_AIR_2000 etc. Or it might be empty if EquationOfStateV1's factory
  // directly creates based on ID. For analytic models, usually no parameters
  // are loaded from a file for the instance itself. The 'form_' enum (or
  // similar mechanism) determines which hardcoded Fortran routine to call.
  if (tfd_data != nullptr) {
    std::cout << "Warning (AnalyticEOS::initialize): TFD data provided but not "
                 "used by AnalyticEOS."
              << std::endl;
  }
  // Optionally, open material_param_filepath to read descriptive metadata if
  // present, beyond just the #MODEL_TYPE line (which EquationOfStateV1 already
  // processed).
  return EOS_SUCCESS;
}

ComputeResult AnalyticEOS::compute(double rho, double T) const {
  // Initialize derivatives to 0 as these simple analytic forms might not
  // provide them
  ComputeResult result;

  switch (form_) {
    case AnalyticForm::AIR_2000:
      c_air_eos_2000(rho, T, &result.P, &result.E, &result.istat);
      break;
    case AnalyticForm::CARBON_2001:
      c_carbon_eos_2001(rho, T, &result.P, &result.E, &result.istat);
      break;
    // Add cases for other analytic forms
    case AnalyticForm::UNDEFINED:
    default:
      std::cerr << "Error (AnalyticEOS::compute): Undefined or unsupported "
                   "analytic form for EOS ID "
                << eos_id_ << std::endl;
      result.istat = EOS_ERROR_INVALID_EOS_TYPE;
      break;
  }
  return result;
}

int AnalyticEOS::pack_parameters(std::ostream& os) const {
  // Pack the 'form_' enum as its underlying integer type
  int form_val = static_cast<int>(form_);
  os.write(reinterpret_cast<const char*>(&form_val), sizeof(form_val));
  if (!os.good()) {
    std::cerr << "Error (AnalyticEOS::pack_parameters): Failed to write form "
                 "for EOS ID "
              << eos_id_ << std::endl;
    return EOS_ERROR_UNPACK_FAILED;
  }
  return EOS_SUCCESS;
}

int AnalyticEOS::unpack_parameters(std::istream& is) {
  int form_val = 0;
  is.read(reinterpret_cast<char*>(&form_val), sizeof(form_val));
  if (!is.good()) {
    std::cerr << "Error (AnalyticEOS::unpack_parameters): Failed to read form "
                 "for EOS ID "
              << eos_id_ << std::endl;
    form_ = AnalyticForm::UNDEFINED;  // Reset to a safe default
    return EOS_ERROR_UNPACK_FAILED;
  }

  // Validate the unpacked enum value
  if (form_val >= static_cast<int>(AnalyticForm::UNDEFINED) &&
      form_val <= static_cast<int>(
                      AnalyticForm::CARBON_2001)) {  // Assuming CARBON_2001 is
                                                     // the last valid form
    form_ = static_cast<AnalyticForm>(form_val);
  } else {
    std::cerr
        << "Error (AnalyticEOS::unpack_parameters): Invalid AnalyticForm value "
        << form_val << " unpacked for EOS ID " << eos_id_ << std::endl;
    form_ = AnalyticForm::UNDEFINED;
    return EOS_ERROR_UNPACK_FAILED;  // Or a specific error for invalid enum
  }
  return EOS_SUCCESS;
}
