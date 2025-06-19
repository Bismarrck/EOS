#ifndef ANALYTIC_EOS_H
#define ANALYTIC_EOS_H

#include "MaterialEOS.h"
#include <functional> // For std::function if used for dispatch

// Forward declare any specific Fortran C-wrapper functions if needed here
// extern "C" { ... }

class AnalyticEOS : public MaterialEOS {
public:
    enum class AnalyticForm {
        UNDEFINED,
        AIR_2000,    // Maps to air_eos_2000
        CARBON_2001  // Maps to carbon_eos_2001
        // Add more specific analytic forms here
    };

    AnalyticEOS(int eos_id, AnalyticForm form, const std::string& name = "");

    int initialize(const std::string& material_param_filepath,
                   const std::string& eos_data_dir_root,
                   const EOS_Internal::TFDMatrices* tfd_data) override;

    ComputeResult compute(double rho, double T) const override;

    int pack_parameters(std::ostream& os) const override;
    int unpack_parameters(std::istream& is) override;

private:
    AnalyticForm form_;
    // Could use a std::function to store pointer to the correct C++ shim or Fortran wrapper call
    // std::function<int(double, double, double&, double&, ...)> compute_func_;
    // For this example, we'll use a switch on form_ in the .cpp
};

#endif // ANALYTIC_EOS_H