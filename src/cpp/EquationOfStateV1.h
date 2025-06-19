// src/cpp/EquationOfStateV1.h
#ifndef EQUATIONOFSTATEV1_H
#define EQUATIONOFSTATEV1_H

#include <hdf5.h>

#include <functional>
#include <map>
#include <memory>  // For std::unique_ptr
#include <string>
#include <vector>

// Forward declare MaterialEOS and TFDMatrices from their new location
class MaterialEOS;  // Assuming MaterialEOS.h will be included in .cpp
struct ComputeResult;
namespace EOS_Internal {
struct TFDMatrices;
}

// Forward declare Fortran-interfacing functions (actual definitions in .cpp)
// These are the C names given in BIND(C, name='...')
extern "C" {
void c_set_complicated_eos_use_old_cold_term(bool value, int* istat_out);
void c_set_use_tfd_data_ver1(bool value, int* istat_out);
// bool c_get_use_tfd_data_ver1(int* istat_out); // If getter is implemented
}

class EquationOfStateV1 {
 public:
  // EOSType to be read from files
  enum class EOSTypeFromFile { NOT_SET, ANALYTIC, COMPLICATED };

  // More specific internal types for dispatch
  enum class InternalEOSType {
    UNKNOWN,
    ANALYTIC_AIR,
    ANALYTIC_CARBON,
    COMPLICATED_V1  // Example if we have versions of complicated models
  };

  // C++ function signature for an analytic EOS compute call (after C++
  // shimming) Takes (rho, T, P_out, E_out, istat_out)
  using AnalyticEOSFunc = std::function<int(double, double, double&, double&)>;

  // C++ function signature for a complicated EOS compute call
  // Takes (params_vec, rho, T, P_out, E_out, dPdT_out, dEdT_out, dPdrho_out,
  // istat_out)
  using ComplicatedEOSFunc =
      std::function<int(const std::vector<double>&, double, double, double&,
                        double&, double&, double&, double&)>;

  struct MaterialData {
    int eos_id_key;  // The ID used to register/find this material
    EOSTypeFromFile type_from_file;  // Type read from #TYPE line
    InternalEOSType internal_type;   // Specific internal type for dispatch

    std::vector<double> params;  // For COMPLICATED type
    AnalyticEOSFunc
        analytic_func;  // For ANALYTIC type (points to a C++ lambda/shim)

    MaterialData()
        : eos_id_key(0),
          type_from_file(EOSTypeFromFile::NOT_SET),
          internal_type(InternalEOSType::UNKNOWN),
          analytic_func(nullptr) {}
  };

  EquationOfStateV1();
  ~EquationOfStateV1();

  // --- Initialization and Configuration ---
  int initialize(const std::vector<int>& eos_id_list,
                 const std::string& eos_data_dir);
  int check_eos_data_dir(
      const std::string& eos_data_dir,
      const std::vector<int>& eos_ids_to_check);  // Basic check

  // Control Variable Setters (from Phase 1)
  int setComplicatedEOSUseOldColdTerm(bool value);
  int setUseTFDDataVer1(
      bool value);  // This will influence TFD file choice in initialize

  void set_perform_signature_check(bool enable);
  bool get_perform_signature_check() const;

  // Main TFD loading function
  int loadTFDDataInternal(
      const std::string& hdf5_filepath);  // Takes base dir for TFD files

  // --- Computation ---
  ComputeResult compute(int eos_id, double rho, double T);

  void free_resources();  // Renamed from free() to avoid conflict if any

  // --- MPI Data Packing/Unpacking ---
  // Calculates required buffer size, or packs data if buffer is provided.
  int pack_data(char*& buffer, int& buffer_size)
      const;  // Made const as it shouldn't change state

  // Unpacks data from buffer and initializes the object.
  int unpack_data(const char* buffer, int buffer_size);

 private:
  // --- Private Helper Methods ---

  // --- Members ---
  std::unique_ptr<EOS_Internal::TFDMatrices> tfd_data_;  // Owns the TFD data
  std::map<int, std::unique_ptr<MaterialEOS>>
      loaded_materials_;  // Keyed by integer eos_id

  // Registry for known analytic EOS types and their corresponding C++ shim
  // functions and internal types Key: The eos_id as it would appear in
  // eosXXXXX.dat (e.g. 10000 for a complicated, or a specific ID for analytic)
  // Value: Pair of (InternalEOSType, C++ shim function pointer)
  struct AnalyticRegistryEntry {
    InternalEOSType internal_type;
    AnalyticEOSFunc func_ptr;
  };
  std::map<int, AnalyticRegistryEntry> analytic_eos_registry_;

  bool tfd_use_ver1_setting_ =
      true;  // C++ mirror of the Fortran control, default true
  bool complicated_eos_use_old_cold_term_setting_ = false;

  // Signature check
  bool perform_signature_check_ = true;
  void check_file_signature(int eos_id, const std::string& full_filepath) const;
};

#endif  // EQUATIONOFSTATEV1_H