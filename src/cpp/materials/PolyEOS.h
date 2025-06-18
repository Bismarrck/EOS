#ifndef POLY_EOS_H
#define POLY_EOS_H

#include "MaterialEOS.h"
#include <vector>
// #include <hdf5.h> // Only if HDF5 types are directly in header, prefer in .cpp

// Forward declare specific Fortran C-wrapper for poly_eos_compute
// extern "C" void c_poly_eos_compute(...);


class PolyEOS : public MaterialEOS {
public:
    PolyEOS(int eos_id, const std::string& name = "");

    int initialize(const std::string& material_param_filepath, // This will be an HDF5 file path
                   const std::string& eos_data_dir_root,
                   const EOS_Internal::TFDMatrices* tfd_data) override; // tfd_data will be nullptr

    int compute(double rho, double T,
                double& P_out, double& E_out,
                double& dPdT_out, double& dEdT_out, double& dPdrho_out) const override;

    int pack_parameters(std::ostream& os) const override;
    int unpack_parameters(std::istream& is) override;

private:
    // Members to store polynomial coefficients and other parameters loaded from HDF5
    // Example:
    std::vector<double> poly_coefficients_;
    double rho_min_validity_;
    double rho_max_validity_;
    // ... etc.
};

#endif // POLY_EOS_H