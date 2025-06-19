#ifndef MATERIAL_EOS_H
#define MATERIAL_EOS_H

#include <iosfwd>  // For std::ostream&, std::istream& forward declarations
#include <string>

// Forward declaration for TFDMatrices to avoid circular include if MaterialEOS
// methods took it by value/ref However, since initialize takes it by pointer, a
// forward declaration is fine.
namespace EOS_Internal {
struct TFDMatrices;
}


struct ComputeResult {
  double P = 0.0;
  double E = 0.0;
  double dPdT = 0.0;
  double dEdT = 0.0;
  double dPdrho = 0.0;
  // Add other common results like sound_speed, Cv, Cp if needed later
  int istat = -1; // Status of the computation (0 for success)

  explicit ComputeResult(const int status = -1) : istat(status) {} // Default constructor
};

class MaterialEOS {
 public:
  // Enum to distinguish the high-level model category
  // This might be redundant if the factory in EquationOfStateV1 uses strings,
  // but can be useful internally for the MaterialEOS objects themselves.
  enum class ModelCategory {
    UNSET,
    ANALYTIC,
    COMPLICATED_TABLE_TFD,  // Uses TFD data
    POLYNOMIAL_HDF5,        // Reads parameters from its own HDF5
    // Add more categories as new fundamental model types are introduced
  };

  MaterialEOS(int eos_id, ModelCategory category, std::string name = "")
      : eos_id_(eos_id), model_category_(category), name_(std::move(name)) {}

  virtual ~MaterialEOS() = default;

  // Initialization method for the specific material instance.
  // material_param_filepath: Full path to the specific parameter file for this
  // material (e.g., eosXXXXX.dat or eosYYYYY.h5).
  // eos_data_dir_root: Root path of the global EOS data directory (for
  // resolving relative paths if any). tfd_data: Pointer to the globally loaded
  // TFDMatrices (can be nullptr if this model doesn't use TFD).
  virtual int initialize(const std::string& material_param_filepath,
                         const std::string& eos_data_dir_root,
                         const EOS_Internal::TFDMatrices* tfd_data) = 0;

  // Core computation method.
  virtual ComputeResult compute(double rho, double T) const = 0;

  // Methods for MPI packing/unpacking of *this material's specific parameters*.
  // TFD data is handled by EquationOfStateV1.
  // The stream can be a custom buffer wrapper or directly std::ostream/istream
  // on a stringstream.
  virtual int pack_parameters(std::ostream& os) const = 0;
  virtual int unpack_parameters(
      std::istream& is) = 0;  // Non-const as it modifies the object

  // Accessors
  int get_eos_id() const { return eos_id_; }
  ModelCategory get_model_category() const { return model_category_; }
  const std::string& get_name() const { return name_; }
  void set_name(const std::string& name) { name_ = name; }

 protected:
  int eos_id_;
  ModelCategory model_category_;
  std::string name_;
};

#endif  // MATERIAL_EOS_H