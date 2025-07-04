#include "PolyEOS.h"
#include <hdf5.h>
#include <iostream>
#include <sstream>  // For pack/unpack
#include "utils/hdf5_utils.h"
#include "eos_error_codes.h"

// Declare the C-wrapper for the Fortran poly_eos_compute routine
extern "C" {
void c_poly_eos_compute(const double* poly_coeffs, int num_coeffs, double rho,
                        double T, double* P, double* E, int* istat);
}

PolyEOS::PolyEOS(int eos_id, const std::string& name)
    : MaterialEOS(eos_id, ModelCategory::POLYNOMIAL_HDF5, name),
      rho_min_validity_(0.0),
      rho_max_validity_(0.0)  // Initialize members
{
  std::cout << "PolyEOS created for ID " << eos_id_ << std::endl;
}

int PolyEOS::initialize(
    const std::string& material_param_filepath_h5,
    const std::string& eos_data_dir_root,  // May not be needed if path is full
    const EOS_Internal::TFDMatrices* tfd_data) {  // tfd_data will be nullptr

  std::cout << "Initializing PolyEOS ID " << eos_id_
            << " from HDF5 file: " << material_param_filepath_h5 << std::endl;

  if (tfd_data != nullptr) {
    std::cout << "Warning (PolyEOS::initialize): TFD data provided but not "
                 "used by PolyEOS."
              << std::endl;
  }

  // Suppress HDF5 auto error printing for this scope
  H5E_auto2_t old_func;
  void* old_client_data;
  H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  hid_t file_id =
      H5Fopen(material_param_filepath_h5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    std::cerr
        << "Error (PolyEOS::initialize): Could not open HDF5 parameter file: "
        << material_param_filepath_h5 << std::endl;
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_POLY_EOS_INIT_FAILED;
  }

  // Read parameters specific to PolyEOS from the HDF5 file
  hid_t dset_id = H5Dopen2(file_id, "/polynomial/coefficients",
                           H5P_DEFAULT);  // Example path
  if (dset_id < 0) {
    std::cerr << "HDF5 Error (PolyEOS): Could not open dataset "
                 "'/polynomial/coefficients'."
              << std::endl;
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_POLY_EOS_INIT_FAILED;
  }
  hid_t dspace_id = H5Dget_space(dset_id);
  int rank = H5Sget_simple_extent_ndims(dspace_id);
  if (rank != 1) { /* error not 1D array */
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_POLY_EOS_INIT_FAILED;
  }

  hsize_t dims[1];
  H5Sget_simple_extent_dims(dspace_id, dims, NULL);
  poly_coefficients_.resize(dims[0]);
  herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, poly_coefficients_.data());

  H5Sclose(dspace_id);
  H5Dclose(dset_id);

  if (status < 0) {
    std::cerr
        << "HDF5 Error (PolyEOS): Failed to read '/polynomial/coefficients'."
        << std::endl;
    H5Fclose(file_id);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return EOS_ERROR_POLY_EOS_INIT_FAILED;
  }
  std::cout << "PolyEOS: Loaded " << poly_coefficients_.size()
            << " polynomial coefficients." << std::endl;

  // Example: Read validity range rho_min (scalar double)
  if (HDF5Utils::read_hdf5_scalar_double(file_id, "/validity/rho_min",
                                         rho_min_validity_) <
      0) {  // You'll need this helper
    std::cout << "Warning (PolyEOS): Could not read /validity/rho_min."
              << std::endl;
    // Set a default or handle as error
  }
  if (HDF5Utils::read_hdf5_scalar_double(file_id, "/validity/rho_max",
                                         rho_max_validity_) < 0) {
    std::cout << "Warning (PolyEOS): Could not read /validity/rho_max."
              << std::endl;
  }
  std::cout << "PolyEOS: rho_min_validity=" << rho_min_validity_
            << ", rho_max_validity=" << rho_max_validity_ << std::endl;

  // Read other necessary parameters from HDF5...

  H5Fclose(file_id);
  H5Eset_auto2(H5E_DEFAULT, old_func,
               old_client_data);  // Restore error printing

  std::cout << "PolyEOS ID " << eos_id_ << " initialized." << std::endl;
  return EOS_SUCCESS;
}

ComputeResult PolyEOS::compute(double rho, double T) const {
  ComputeResult result;
  if (poly_coefficients_.empty()) {
    std::cerr << "Error (PolyEOS::compute): Polynomial coefficients not loaded "
                 "for EOS ID "
              << eos_id_ << std::endl;
    result.istat = EOS_ERROR_POLY_EOS_EMPTY_PARAMS;
    return result;
  }
  c_poly_eos_compute(poly_coefficients_.data(),
                     static_cast<int>(poly_coefficients_.size()), rho, T,
                     &result.P, &result.E, &result.istat);
  return result;
}

int PolyEOS::pack_parameters(std::ostream& os) const {
  // 1. Pack poly_coefficients_
  size_t num_coeffs = poly_coefficients_.size();
  os.write(reinterpret_cast<const char*>(&num_coeffs), sizeof(num_coeffs));
  if (!os.good()) return EOS_ERROR_PACK_FAILED;
  if (num_coeffs > 0) {
    os.write(reinterpret_cast<const char*>(poly_coefficients_.data()),
             num_coeffs * sizeof(double));
    if (!os.good()) return EOS_ERROR_PACK_FAILED;
  }

  // 2. Pack validity range (example members)
  os.write(reinterpret_cast<const char*>(&rho_min_validity_),
           sizeof(rho_min_validity_));
  if (!os.good()) return EOS_ERROR_PACK_FAILED;
  os.write(reinterpret_cast<const char*>(&rho_max_validity_),
           sizeof(rho_max_validity_));
  if (!os.good()) return EOS_ERROR_PACK_FAILED;
  // Pack any other PolyEOS specific parameters here...

  return EOS_SUCCESS;
}

int PolyEOS::unpack_parameters(std::istream& is) {
  poly_coefficients_.clear();
  // 1. Unpack poly_coefficients_
  size_t num_coeffs = 0;
  is.read(reinterpret_cast<char*>(&num_coeffs), sizeof(num_coeffs));
  if (!is.good() && num_coeffs > 0) return EOS_ERROR_UNPACK_FAILED;
  if (num_coeffs == 0 && !is.good() && !is.eof()) return EOS_ERROR_UNPACK_FAILED;

  if (num_coeffs > 0) {
    try {
      poly_coefficients_.resize(num_coeffs);
    } catch (const std::bad_alloc& e) {
      std::cerr
          << "Error (PolyEOS::unpack): Failed to resize poly_coefficients_: "
          << e.what() << std::endl;
      return EOS_ERROR_UNPACK_FAILED;
    }
    is.read(reinterpret_cast<char*>(poly_coefficients_.data()),
            num_coeffs * sizeof(double));
    if (!is.good()) {
      poly_coefficients_.clear();
      return EOS_ERROR_UNPACK_FAILED;
    }
  }

  // 2. Unpack validity range
  is.read(reinterpret_cast<char*>(&rho_min_validity_),
          sizeof(rho_min_validity_));
  if (!is.good()) return EOS_ERROR_UNPACK_FAILED;
  is.read(reinterpret_cast<char*>(&rho_max_validity_),
          sizeof(rho_max_validity_));
  if (!is.good()) return EOS_ERROR_UNPACK_FAILED;

  // Unpack any other PolyEOS specific parameters here...

  return EOS_SUCCESS;
}
