#include "AnalyticEOS.h"
#include <iostream>
#include <sstream> // For pack/unpack

// Declare C-wrappers for analytic Fortran routines
extern "C" {
    void c_air_eos_2000(double rho, double T, double* P, double* E, int* istat_out);
    void c_carbon_eos_2001(double rho, double T, double* P, double* E, int* istat_out);
    // Add declarations for other analytic C-wrappers if you have them
}

AnalyticEOS::AnalyticEOS(int eos_id, AnalyticForm form, const std::string& name)
    : MaterialEOS(eos_id, ModelCategory::ANALYTIC, name), form_(form) {
    std::cout << "AnalyticEOS created for ID " << eos_id_ << " with form " << static_cast<int>(form_) << std::endl;
}

int AnalyticEOS::initialize(const std::string& material_param_filepath,
                           const std::string& eos_data_dir_root,
                           const EOS_Internal::TFDMatrices* tfd_data) {
    std::cout << "Initializing AnalyticEOS ID " << eos_id_ << " (Form: " << static_cast<int>(form_) << ")" << std::endl;
    // material_param_filepath might be a .dat file with just #MODEL_TYPE: ANALYTIC_AIR_2000 etc.
    // Or it might be empty if EquationOfStateV1's factory directly creates based on ID.
    // For analytic models, usually no parameters are loaded from a file for the instance itself.
    // The 'form_' enum (or similar mechanism) determines which hardcoded Fortran routine to call.
    if (tfd_data != nullptr) {
        std::cout << "Warning (AnalyticEOS::initialize): TFD data provided but not used by AnalyticEOS." << std::endl;
    }
    // Optionally, open material_param_filepath to read descriptive metadata if present,
    // beyond just the #MODEL_TYPE line (which EquationOfStateV1 already processed).
    return MAT_EOS_SUCCESS;
}

int AnalyticEOS::compute(double rho, double T,
                        double& P_out, double& E_out,
                        double& dPdT_out, double& dEdT_out, double& dPdrho_out) const {
    int istat_fortran = -1; // Default error

    // Initialize derivatives to 0 as these simple analytic forms might not provide them
    dPdT_out = 0.0;
    dEdT_out = 0.0;
    dPdrho_out = 0.0;

    switch (form_) {
        case AnalyticForm::AIR_2000:
            c_air_eos_2000(rho, T, &P_out, &E_out, &istat_fortran);
            break;
        case AnalyticForm::CARBON_2001:
            c_carbon_eos_2001(rho, T, &P_out, &E_out, &istat_fortran);
            break;
        // Add cases for other analytic forms
        case AnalyticForm::UNDEFINED:
        default:
            std::cerr << "Error (AnalyticEOS::compute): Undefined or unsupported analytic form for EOS ID " << eos_id_ << std::endl;
            P_out = 0.0; E_out = 0.0;
            return MAT_EOS_COMPUTE_FAILED; // Or specific error
    }
    return istat_fortran;
}

int AnalyticEOS::pack_parameters(std::ostream& os) const {
    // Analytic models might have no parameters to pack, or just their 'form' / 'type'
    int form_val = static_cast<int>(form_);
    os.write(reinterpret_cast<const char*>(&form_val), sizeof(form_val));
    return os.good() ? MAT_EOS_SUCCESS : MAT_EOS_PACK_FAILED;
}

int AnalyticEOS::unpack_parameters(std::istream& is) {
    int form_val = 0;
    is.read(reinterpret_cast<char*>(&form_val), sizeof(form_val));
    if (!is.good()) return MAT_EOS_UNPACK_FAILED;
    form_ = static_cast<AnalyticForm>(form_val);
    // Basic validation of unpacked form_
    if (form_ < AnalyticForm::UNDEFINED || form_ > AnalyticForm::CARBON_2001 /* adjust with max enum */) {
         std::cerr << "Warning (AnalyticEOS::unpack): Unpacked invalid AnalyticForm value: " << form_val << std::endl;
         form_ = AnalyticForm::UNDEFINED; // Reset to safe default
         return MAT_EOS_UNPACK_FAILED; // Or handle as error
    }
    return MAT_EOS_SUCCESS;
}
