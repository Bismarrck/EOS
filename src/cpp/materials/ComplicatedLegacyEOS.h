#ifndef COMPLICATED_LEGACY_EOS_H
#define COMPLICATED_LEGACY_EOS_H

#include <vector>

#include "MaterialEOS.h"

// Forward declare specific Fortran C-wrapper
// extern "C" void c_complicated_eos(...);

class ComplicatedLegacyEOS : public MaterialEOS {
 public:
  ComplicatedLegacyEOS(int eos_id, const std::string& name = "");

  int initialize(const std::string& material_param_filepath,
                 const std::string& eos_data_dir_root,
                 const EOS_Internal::TFDMatrices* tfd_data) override;

  ComputeResult compute(double rho, double T) const override;

  int pack_parameters(std::ostream& os) const override;
  int unpack_parameters(std::istream& is) override;

  void internal_set_tfd_pointer_after_unpack(
      const EOS_Internal::TFDMatrices* tfd_data) {
    tfd_data_ptr_ = tfd_data;
  }

 private:
  std::vector<double> params_;  // Loaded from text file
  const EOS_Internal::TFDMatrices* tfd_data_ptr_ =
      nullptr;  // Pointer to shared TFD data
  // Predefined order of keys for parsing this model's parameter file
  // Can be static const or initialized in constructor
  static const std::vector<std::string> s_param_order_;
};

#endif  // COMPLICATED_LEGACY_EOS_H