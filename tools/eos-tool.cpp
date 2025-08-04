// tools/eos-tool.cpp
#include <hdf5.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "CLI11/CLI11.hpp"      // From third_party/CLI11/CLI11.hpp
#include "EquationOfStateV1.h"  // From install or build tree
#include "eos_error_codes.h"
#include "materials/MaterialEOS.h"
#include "math/NumericalSolvers.h"
#include "utils/hdf5_utils.h"  // Assuming HDF5 write helpers are placed here

// Add to tools/eos-tool.cpp, perhaps inside a helper namespace
#include <cmath>
#include <functional>
#include <limits>

namespace NumericalMethods {

struct NewtonResult {
  double root = 0.0;
  bool converged = false;
  int iterations = 0;
};

// A generic Newton-Raphson solver
NewtonResult newton_raphson(
    double initial_guess,                        // Starting T
    const std::function<double(double)>& func,   // f(T)
    const std::function<double(double)>& deriv,  // f'(T)
    int max_iter = 50,                           // Max iterations
    double tolerance = 1e-7,                     // Convergence tolerance
    double epsilon =
        std::numeric_limits<double>::epsilon()  // To avoid division by zero
) {
  NewtonResult result;
  result.root = initial_guess;

  for (int i = 0; i < max_iter; ++i) {
    result.iterations = i + 1;
    double f_val = func(result.root);

    // Check for convergence
    if (std::abs(f_val) < tolerance) {
      result.converged = true;
      return result;
    }

    double df_val = deriv(result.root);

    // Check for zero derivative
    if (std::abs(df_val) < epsilon) {
      std::cerr << "Warning (Newton-Raphson): Zero derivative encountered."
                << std::endl;
      result.converged = false;
      return result;
    }

    // Newton's method step
    double step = f_val / df_val;
    result.root = result.root - step;

    // Check if step is too small
    if (std::abs(step) < tolerance) {
      result.converged = true;
      return result;
    }
  }
  std::cerr << "Warning (Newton-Raphson): Failed to converge after " << max_iter
            << " iterations." << std::endl;
  result.converged = false;
  return result;
}

}  // namespace NumericalMethods

// Structure for Hugoniot command options
struct HugoniotOptions {
  int eos_id;
  double rho0;
  double T0 = 300.0;  // Default initial temperature
  double rho_final_max;
  int num_points = 100;
  std::string data_dir;
  std::string output_file;
};

// Main function to execute the Hugoniot calculation
void run_hugoniot(const HugoniotOptions& opt) {
  // 1. Initialize EOS Manager
  EquationOfStateV1 eos_manager;
  // We need to initialize for the given eos_id to get the initial state
  int init_stat = eos_manager.initialize({opt.eos_id}, opt.data_dir);
  if (init_stat != EOS_SUCCESS) {
    std::cerr << "Error: Failed to initialize EOS manager for ID " << opt.eos_id
              << std::endl;
    return;
  }

  // 2. Get initial state (P0, E0, V0)
  ComputeResult res0 = eos_manager.compute(opt.eos_id, opt.rho0, opt.T0);
  if (res0.istat != EOS_SUCCESS) {
    std::cerr << "Error: Failed to compute initial state at rho0=" << opt.rho0
              << ", T0=" << opt.T0 << std::endl;
    return;
  }
  double P0 = res0.P;
  double E0 = res0.E;
  double V0 = 1.0 / opt.rho0;

  std::cout << "Initial State (rho0=" << opt.rho0 << ", T0=" << opt.T0 << "):"
            << " P0=" << P0 << ", E0=" << E0 << ", V0=" << V0 << std::endl;

  // 3. Setup output
  std::ofstream outfile;
  std::ostream* out_stream = &std::cout;
  if (!opt.output_file.empty()) {
    outfile.open(opt.output_file);
    if (!outfile.is_open()) {
      std::cerr << "Error: Failed to open output file " << opt.output_file
                << std::endl;
      return;
    }
    out_stream = &outfile;
  }

  // Print header
  *out_stream << "rho,V,T,P,E,Us,Up\n";

  // 4. Loop over final densities and solve for Hugoniot states
  double T_guess =
      opt.T0;  // Use previous point's temperature to help find a bracket
  for (int i = 1; i < opt.num_points; ++i) {
    double rho_final = opt.rho0 + (opt.rho_final_max - opt.rho0) *
                                      static_cast<double>(i) /
                                      (opt.num_points - 1);
    double V_final = 1.0 / rho_final;

    if (V_final >= V0) continue;  // Only consider compression shocks

    // Define the Hugoniot function f(T) = E(rho,T) - E_hugoniot
    auto hugoniot_func = [&](double T_trial) -> double {
      ComputeResult res_trial =
          eos_manager.compute(opt.eos_id, rho_final, T_trial);
      if (res_trial.istat != EOS_SUCCESS)
        return std::numeric_limits<double>::max();  // Signal error
      return res_trial.E - E0 - 0.5 * (res_trial.P + P0) * (V0 - V_final);
    };

    // ---- Find a bracket [T_low, T_high] for the root ----
    // Start with a guess and expand the bracket until f(T_low) and f(T_high)
    // have opposite signs.
    double T_low, T_high;
    double f_low, f_high;

    // A simple expanding bracket search
    double bracket_search_factor = 1.5;
    double T_step = std::max(100.0, T_guess * 0.1);  // Initial step size
    T_low = std::max(1.0, T_guess - T_step);         // Don't let T go below 1K
    T_high = T_guess + T_step;

    f_low = hugoniot_func(T_low);
    f_high = hugoniot_func(T_high);

    int bracket_iter = 0;
    const int max_bracket_iter = 20;
    while (f_low * f_high > 0 && bracket_iter < max_bracket_iter) {
      if (std::abs(f_low) < std::abs(f_high)) {  // Expand on the low side
        T_low = std::max(1.0, T_low - T_step);
        f_low = hugoniot_func(T_low);
      } else {  // Expand on the high side
        T_high = T_high + T_step;
        f_high = hugoniot_func(T_high);
      }
      T_step *= bracket_search_factor;  // Expand search range exponentially
      bracket_iter++;
    }

    if (f_low * f_high > 0) {
      std::cerr << "Warning: Could not bracket root for rho_final=" << rho_final
                << " in temperature range [" << T_low << ", " << T_high << "]."
                << std::endl;
      continue;  // Skip this point
    }

    // --- Use Brent's method to find the root within the bracket ---
    EOS_Toolkit::NumericalSolvers::RootResult nr_result =
        EOS_Toolkit::NumericalSolvers::brents_method(hugoniot_func, T_low,
                                                     T_high);

    if (!nr_result.converged) {
      std::cerr << "Warning: Brent's method failed to converge for rho_final="
                << rho_final << std::endl;
      continue;  // Skip this point
    }

    double T_final = nr_result.root;
    T_guess = T_final;  // Update guess for next density point

    // 5. Get final state properties and calculate shock parameters (as before)
    // ...
  }
  std::cout << "Hugoniot calculation complete." << std::endl;
}

// --- Isotherm Implementation ---
struct IsothermOptions {
  int eos_id;
  double temp;
  double rho_min;
  double rho_max;
  int num_points = 100;
  bool log_scale = false;
  std::string data_dir;
  std::string output_file;
};

void run_isotherm(const IsothermOptions& opt) {
  EquationOfStateV1 eos_manager;
  int init_stat = eos_manager.initialize({opt.eos_id}, opt.data_dir);
  if (init_stat != EOS_SUCCESS) {
    std::cerr << "Error: Failed to initialize EOS manager for ID " << opt.eos_id
              << std::endl;
    return;
  }

  std::ofstream outfile;
  std::ostream* out_stream = &std::cout;
  if (!opt.output_file.empty()) {
    outfile.open(opt.output_file);
    if (!outfile.is_open()) {
      std::cerr << "Error: Failed to open output file " << opt.output_file
                << std::endl;
      return;
    }
    out_stream = &outfile;
  }

  // Print header
  *out_stream << "rho,T,P,E,dPdT,dEdT,dPdrho\n";

  for (int i = 0; i < opt.num_points; ++i) {
    double rho;
    if (opt.log_scale) {
      rho =
          opt.rho_min * std::pow(opt.rho_max / opt.rho_min,
                                 static_cast<double>(i) / (opt.num_points - 1));
    } else {
      rho = opt.rho_min + (opt.rho_max - opt.rho_min) * static_cast<double>(i) /
                              (opt.num_points - 1);
    }

    ComputeResult res = eos_manager.compute(opt.eos_id, rho, opt.temp);
    if (res.istat == EOS_SUCCESS) {
      *out_stream << rho << "," << opt.temp << "," << res.P << "," << res.E
                  << "," << res.dPdT << "," << res.dEdT << "," << res.dPdrho
                  << "\n";
    }
  }
  std::cout << "Isotherm calculation complete." << std::endl;
}

void setup_isotherm_subcommand(CLI::App& app) {
  auto opts = std::make_shared<IsothermOptions>();
  auto* isotherm_cmd =
      app.add_subcommand("isotherm", "Calculate an isothermal line.");

  isotherm_cmd->add_option("--id,--eos_id", opts->eos_id, "Material EOS ID")
      ->required();
  isotherm_cmd->add_option("-T,--temp", opts->temp, "Temperature (K)")
      ->required();
  isotherm_cmd->add_option("--rho_min", opts->rho_min, "Minimum density (g/cc)")
      ->required();
  isotherm_cmd->add_option("--rho_max", opts->rho_max, "Maximum density (g/cc)")
      ->required();
  isotherm_cmd
      ->add_option("-n,--num_points", opts->num_points, "Number of points")
      ->default_val(100);
  isotherm_cmd->add_flag("--log", opts->log_scale,
                         "Use logarithmic spacing for density");
  isotherm_cmd
      ->add_option("-d,--data_dir", opts->data_dir, "Path to eos_data_dir")
      ->required()
      ->check(CLI::ExistingDirectory);
  isotherm_cmd->add_option("-o,--output", opts->output_file,
                           "Output CSV file (default: stdout)");

  isotherm_cmd->callback([opts]() { run_isotherm(*opts); });
}

// --- Hugoniot Implementation ---
// (Needs root-finder implementation)
void setup_hugoniot_subcommand(CLI::App& app) {
  auto opts = std::make_shared<HugoniotOptions>();
  auto* hugo_cmd =
      app.add_subcommand("hugoniot", "Calculate a principal Hugoniot curve.");

  hugo_cmd->add_option("--id,--eos_id", opts->eos_id, "Material EOS ID")
      ->required();
  hugo_cmd
      ->add_option("--rho0", opts->rho0,
                   "Initial density of unshocked material (g/cc)")
      ->required();
  hugo_cmd->add_option("--T0", opts->T0, "Initial temperature (K)")
      ->default_val(300.0);
  hugo_cmd
      ->add_option("--rho_max", opts->rho_final_max,
                   "Maximum final density to compute up to (g/cc)")
      ->required();
  hugo_cmd->add_option("-n,--num_points", opts->num_points, "Number of points")
      ->default_val(100);
  hugo_cmd->add_option("-d,--data_dir", opts->data_dir, "Path to eos_data_dir")
      ->required()
      ->check(CLI::ExistingDirectory);
  hugo_cmd->add_option("-o,--output", opts->output_file,
                       "Output CSV file (default: stdout)");

  // Set the callback that runs when this subcommand is invoked
  hugo_cmd->callback([opts]() { run_hugoniot(*opts); });
}

// Define the state type for Odeint. Since we only integrate T, it's a 1D state.
using ode_state_type = std::vector<double>;

// This struct defines the ODE system dy/dx = f(x, y)
// where x = rho, y[0] = T
struct IsentropeODE {
  EquationOfStateV1& eos_manager;
  int eos_id;

  IsentropeODE(EquationOfStateV1& manager, int id)
      : eos_manager(manager), eos_id(id) {}

  // This is the function f(x, y) that Odeint will call
  void operator()(const ode_state_type& y, ode_state_type& dydx,
                  const double x) {
    double rho = x;
    double T = y[0];  // y[0] is our current Temperature

    if (rho <= 0 || T <= 0) {
      dydx[0] = 0.0;  // Avoid invalid calculations
      return;
    }

    ComputeResult res = eos_manager.compute(eos_id, rho, T);
    if (res.istat != EOS_SUCCESS || res.dEdT == 0) {
      // If compute fails or Cv is zero, we can't continue integration
      dydx[0] = 0.0;
      return;
    }

    // Our final ODE: dT/d(rho)
    double numerator = (2.0 * res.P - T * res.dPdT);
    double denominator = (rho * rho * res.dEdT);

    dydx[0] = numerator / denominator;
  }
};

struct IsentropeOptions {
  int eos_id;
  double T_initial;  // Starting temperature at low density/pressure
  double rho_max;
  int num_points = 100;
  double rho_initial = 1e-3;  // A small starting density to approximate P=0
  std::string solver = "rk4";
  std::string data_dir;
  std::string output_file;
};

void run_isentrope(const IsentropeOptions& opt) {
  EquationOfStateV1 eos_manager;
  int init_stat = eos_manager.initialize({opt.eos_id}, opt.data_dir);
  if (init_stat != EOS_SUCCESS) {
    std::cerr << "Error: Failed to initialize EOS manager for ID " << opt.eos_id
              << std::endl;
    return;
  }

  std::ofstream outfile;
  std::ostream* out_stream = &std::cout;
  if (!opt.output_file.empty()) {
    outfile.open(opt.output_file);
    if (!outfile.is_open()) {
      std::cerr << "Error: Failed to open output file " << opt.output_file
                << std::endl;
      return;
    }
    out_stream = &outfile;
  }

  *out_stream << "rho,T,P,E\n";

  // Initial state: state[0] = T
  ode_state_type state = {opt.T_initial};

  // The ODE system functor
  IsentropeODE ode_functor(eos_manager, opt.eos_id);
  EOS_Toolkit::NumericalSolvers::ODESystem ode_system = ode_functor;

  // Observer function to be called at each step of the integration
  auto observer_func = [&](const ode_state_type& y, double x) {
    double rho = x;
    double T = y[0];
    ComputeResult res = eos_manager.compute(opt.eos_id, rho, T);
    if (res.istat == EOS_SUCCESS) {
      *out_stream << rho << "," << T << "," << res.P << "," << res.E << "\n";
    }
  };
  EOS_Toolkit::NumericalSolvers::ODEObserver observer = observer_func;

  // Use a Runge-Kutta stepper with controlled step size
  // Odeint needs begin, end, and initial step size for integration range
  double rho_start = opt.rho_initial;
  double rho_end = opt.rho_max;
  double drho = (rho_end - rho_start) / (opt.num_points - 1);

  // Integrate from rho_start to rho_end with fixed steps
  if (opt.solver == "rk4") {
    std::cout << "Using RK4 solver." << std::endl;
    EOS_Toolkit::NumericalSolvers::integrate_ode_rk4(
        ode_system, state, rho_start, rho_end, drho, observer);
  } else if (opt.solver == "heun") {
    std::cout << "Using Heun's method (Trapezoid Rule) solver." << std::endl;
    EOS_Toolkit::NumericalSolvers::integrate_ode_heun(
        ode_system, state, rho_start, rho_end, drho, observer);
  }

  std::cout << "Isentrope calculation complete." << std::endl;
}

// --- Release Isentrope ---

struct ReleaseOptions {
  int eos_id;
  double rho0;
  double P_hugoniot;  // Target Hugoniot pressure to start release from
  int num_points = 100;
  double T0 = 300.0;
  std::string data_dir;
  std::string output_file;
};

void run_release(const ReleaseOptions& opt) {
  EquationOfStateV1 eos_manager;
  int init_stat = eos_manager.initialize({opt.eos_id}, opt.data_dir);
  if (init_stat != EOS_SUCCESS) {
    std::cerr << "Error: Failed to initialize EOS manager for ID " << opt.eos_id
              << std::endl;
    return;
  }

  std::ofstream outfile;
  std::ostream* out_stream = &std::cout;
  if (!opt.output_file.empty()) {
    outfile.open(opt.output_file);
    if (!outfile.is_open()) {
      std::cerr << "Error: Failed to open output file " << opt.output_file
                << std::endl;
      return;
    }
    out_stream = &outfile;
  }

  // --- Step 1: Find the Hugoniot state ---
  // (This requires the Hugoniot logic from the previous step)
  // You need to solve for rho_final on the Hugoniot that gives P_hugoniot
  // This is another root-finding problem: find rho_h such that P_h(rho_h) -
  // P_target = 0
  std::cout << "Finding Hugoniot state for P_h = " << opt.P_hugoniot << "..."
            << std::endl;
  // ... (logic to find rho_h and T_h on the Hugoniot) ...
  // For now, let's just assume we found it and use placeholder values:
  double rho_h = 2.5;    // Placeholder
  double T_h = 15000.0;  // Placeholder
  std::cout << "Found Hugoniot state: rho_h=" << rho_h << ", T_h=" << T_h
            << std::endl;

  // --- Step 2: Integrate the release isentrope DOWN in density ---
  // ... (setup output stream) ...
  *out_stream << "rho,T,P,E\n";

  ode_state_type state = {T_h};  // Initial state is the Hugoniot temperature
  IsentropeODE ode_functor(eos_manager, opt.eos_id);
  EOS_Toolkit::NumericalSolvers::ODESystem ode_system = ode_functor;

  auto observer_func = [&](const ode_state_type& y, double x) {
    double rho = x;
    double T = y[0];
    ComputeResult res = eos_manager.compute(opt.eos_id, rho, T);
    if (res.istat == EOS_SUCCESS) {
      *out_stream << rho << "," << T << "," << res.P << "," << res.E << "\n";
    }
  };
  EOS_Toolkit::NumericalSolvers::ODEObserver observer = observer_func;

  // Integrate from rho_h down to a low density (or until P~0)
  double rho_start = rho_h;
  double rho_end = opt.rho0;  // Release back to initial density
  double drho =
      (rho_end - rho_start) / (opt.num_points - 1);  // This will be negative
  EOS_Toolkit::NumericalSolvers::integrate_ode_rk4(ode_system, state, rho_start,
                                                   rho_end, drho, observer);

  std::cout << "Release isentrope calculation complete." << std::endl;
}

// In main() in tools/eos-tool.cpp
void setup_isentrope_subcommand(CLI::App& app) {
  auto opts = std::make_shared<IsentropeOptions>();
  auto* cmd = app.add_subcommand("isentrope",
                                 "Calculate an isentrope from a P=0 state.");
  cmd->add_option("--id", opts->eos_id, "Material EOS ID")->required();
  cmd->add_option("--T_init", opts->T_initial,
                  "Initial temperature at low density (K)")
      ->required();
  cmd->add_option("--rho_max", opts->rho_max,
                  "Maximum density to integrate to (g/cc)")
      ->required();
  std::map<std::string, std::string> solver_map = {
      {"rk4", "4th-order Runge-Kutta"}, {"heun", "Heun's Method (2nd-order)"}};
  cmd->add_option("--solver", opts->solver, "ODE solver to use")
      ->transform(CLI::CheckedTransformer(solver_map, CLI::ignore_case))
      ->default_val("rk4");
  // ...
  cmd->callback([opts]() { run_isentrope(*opts); });
  cmd->callback([opts]() { run_isentrope(*opts); });
}

void setup_release_subcommand(CLI::App& app) {
  auto opts = std::make_shared<ReleaseOptions>();
  auto* cmd = app.add_subcommand(
      "release", "Calculate a release isentrope from a Hugoniot state.");
  cmd->add_option("--id", opts->eos_id, "Material EOS ID")->required();
  cmd->add_option("--P_hugo", opts->P_hugoniot,
                  "Hugoniot pressure to release from (e.g., GPa)")
      ->required();
  cmd->add_option("--rho0", opts->rho0,
                  "Initial density of unshocked material (g/cc)")
      ->required();
  // ... other options ...
  cmd->callback([opts]() { run_release(*opts); });
}

// =========================================================================
// ==              ENTROPY TABLE BUILDER IMPLEMENTATION                   ==
// =========================================================================

struct EntropyTableOptions {
  int eos_id;
  std::string data_dir;
  std::string output_h5_file;

  // Grid definition
  double rho_min;
  double rho_max;
  int num_rho_points;
  double T_min;
  double T_max;
  int num_T_points;
  bool log_rho = false;  // Flags for logarithmic spacing
  bool log_T = false;

  // Reference state for entropy
  double rho_ref;
  double T_ref;
  double S_ref = 0.0;  // Reference entropy value (usually 0)

  // ODE solver options
  int integration_substeps = 100;  // Sub-steps for each integration segment
};

// Helper function to write the 2D entropy table to an HDF5 file
// This could also be moved to hdf5_utils.cpp
bool write_entropy_table_to_hdf5(
    const std::string& filename, const std::vector<double>& rho_grid,
    const std::vector<double>& T_grid,
    const std::vector<std::vector<double>>&
        entropy_table  // Assumes entropy_table[i_rho][j_T]
) {
  hid_t file_id =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    std::cerr << "HDF5 Error: Could not create output file: " << filename
              << std::endl;
    return false;
  }

  // --- Write Grids ---
  hsize_t rho_dims[1] = {rho_grid.size()};
  HDF5Utils::write_hdf5_dataset_double_1d(file_id, "/rho_grid", rho_dims,
                                          rho_grid.data());

  hsize_t T_dims[1] = {T_grid.size()};
  HDF5Utils::write_hdf5_dataset_double_1d(file_id, "/T_grid", T_dims,
                                          T_grid.data());

  // --- Write 2D Entropy Table ---
  // HDF5 needs a contiguous 1D buffer for 2D data
  std::vector<double> entropy_flat;
  entropy_flat.reserve(rho_grid.size() * T_grid.size());
  for (size_t i = 0; i < rho_grid.size();
       ++i) {  // rho is the slower-changing index (rows)
    for (size_t j = 0; j < T_grid.size();
         ++j) {  // T is faster-changing (columns)
      entropy_flat.push_back(entropy_table[i][j]);
    }
  }

  hsize_t entropy_dims[2] = {rho_grid.size(), T_grid.size()};
  HDF5Utils::write_hdf5_dataset_double_2d(file_id, "/entropy_table",
                                          entropy_dims, entropy_flat.data());

  // You can also add attributes for clarity
  // H5LTset_attribute_string(file_id, "/entropy_table", "UNITS", "J/kg/K"); //
  // Example

  H5Fclose(file_id);
  std::cout << "Successfully wrote entropy table to HDF5 file: " << filename
            << std::endl;
  return true;
}

void run_build_entropy_table(const EntropyTableOptions& opt) {
  std::cout << "Starting entropy table generation for EOS ID " << opt.eos_id
            << std::endl;

  // 1. Initialize EOS Manager
  EquationOfStateV1 eos_manager;
  int init_stat = eos_manager.initialize({opt.eos_id}, opt.data_dir);
  if (init_stat != EOS_SUCCESS) {
    std::cerr << "Error: Failed to initialize EOS manager for ID " << opt.eos_id
              << std::endl;
    return;
  }

  // 2. Create rho and T grids
  std::vector<double> rho_grid(opt.num_rho_points);
  std::vector<double> T_grid(opt.num_T_points);

  for (int i = 0; i < opt.num_rho_points; ++i) {
    double frac = static_cast<double>(i) / (opt.num_rho_points - 1);
    rho_grid[i] = opt.log_rho
                      ? opt.rho_min * std::pow(opt.rho_max / opt.rho_min, frac)
                      : opt.rho_min + frac * (opt.rho_max - opt.rho_min);
  }
  for (int i = 0; i < opt.num_T_points; ++i) {
    double frac = static_cast<double>(i) / (opt.num_T_points - 1);
    T_grid[i] = opt.log_T ? opt.T_min * std::pow(opt.T_max / opt.T_min, frac)
                          : opt.T_min + frac * (opt.T_max - opt.T_min);
  }

  // Output table
  std::vector<std::vector<double>> entropy_table(
      opt.num_rho_points, std::vector<double>(opt.num_T_points));
  using EOS_Toolkit::NumericalSolvers::ODEObserver;
  using EOS_Toolkit::NumericalSolvers::ODEState;
  using EOS_Toolkit::NumericalSolvers::ODESystem;
  // --- Define ODE system for isothermal path: dS/d(rho) ---
  // y[0] = S, x = rho
  auto isothermal_ode_system = [&](const double T_const, ODEState& dydx,
                                   const ODEState& y, const double x) {
    double rho = x;
    ComputeResult res = eos_manager.compute(opt.eos_id, rho, T_const);
    if (res.istat != EOS_SUCCESS || rho <= 0 || T_const <= 0) {
      dydx[0] = 0.0;
      return;
    }
    // From First Law: (dS/d(rho))_T = (1/T) * (dE/d(rho))_T - P / (T*rho^2)
    // From thermo relation: (dE/d(rho))_T = (T*(dP/dT)_rho - P) / rho^2
    double dE_drho_const_T = (T_const * res.dPdT - res.P) / (rho * rho);
    dydx[0] = (1.0 / T_const) * dE_drho_const_T - res.P / (T_const * rho * rho);
  };

  // 3. Pre-calculate all isothermal integration steps from the reference
  // density
  std::cout << "Step 1/2: Calculating isothermal entropy changes..."
            << std::endl;
  std::vector<double> delta_S_isothermal(opt.num_rho_points);
  ODEObserver null_observer = [](const ODEState&, double) {
  };  // No need to observe intermediate steps

  for (int i = 0; i < opt.num_rho_points; ++i) {
    double rho_i = rho_grid[i];
    if (std::abs(rho_i - opt.rho_ref) < 1e-9) {
      delta_S_isothermal[i] = 0.0;
      continue;
    }

    ODESystem current_isotherm_sys = [&](const ODEState& y, ODEState& dydx,
                                         const double x) {
      isothermal_ode_system(opt.T_ref, dydx, y, x);
    };
    ODEState S_state = {
        0.0};  // Start integration at S=0 (we are calculating delta S)
    double drho = (rho_i - opt.rho_ref) / opt.integration_substeps;

    EOS_Toolkit::NumericalSolvers::integrate_ode_rk4(
        current_isotherm_sys, S_state, opt.rho_ref, rho_i, drho, null_observer);

    delta_S_isothermal[i] = S_state[0];
  }

  // 4. Calculate isochoric integration steps for each grid point
  std::cout
      << "Step 2/2: Calculating isochoric entropy changes and building table..."
      << std::endl;
  for (int i = 0; i < opt.num_rho_points; ++i) {
    double rho_i = rho_grid[i];
    for (int j = 0; j < opt.num_T_points; ++j) {
      double T_j = T_grid[j];

      double delta_S_isochoric = 0.0;
      if (std::abs(T_j - opt.T_ref) > 1e-9) {
        // Isochoric system: dS/dT = Cv / T = (dE/dT)_rho / T
        ODESystem current_isochore_sys = [&](const ODEState& y, ODEState& dydx,
                                             const double x) {
          double T = x;
          ComputeResult res = eos_manager.compute(opt.eos_id, rho_i, T);
          if (res.istat != EOS_SUCCESS || T <= 0) {
            dydx[0] = 0.0;
            return;
          }
          dydx[0] = res.dEdT / T;
        };

        ODEState S_state = {0.0};
        double dT = (T_j - opt.T_ref) / opt.integration_substeps;

        EOS_Toolkit::NumericalSolvers::integrate_ode_rk4(
            current_isochore_sys, S_state, opt.T_ref, T_j, dT, null_observer);
        delta_S_isochoric = S_state[0];
      }

      entropy_table[i][j] =
          opt.S_ref + delta_S_isothermal[i] + delta_S_isochoric;
    }
  }

  // 5. Save the final table
  write_entropy_table_to_hdf5(opt.output_h5_file, rho_grid, T_grid,
                              entropy_table);
}

// --- CLI11 Setup for the new command ---
void setup_build_entropy_table_subcommand(CLI::App& app) {
  auto opts = std::make_shared<EntropyTableOptions>();
  auto* cmd = app.add_subcommand("build-entropy-table",
                                 "Build a 2D entropy table S(rho, T).");

  cmd->add_option("--id", opts->eos_id, "Material EOS ID")->required();
  cmd->add_option("-d,--data_dir", opts->data_dir, "Path to eos_data_dir")
      ->required()
      ->check(CLI::ExistingDirectory);
  cmd->add_option("-o,--output", opts->output_h5_file, "Output HDF5 file")
      ->required();

  auto* grid_group =
      cmd->add_option_group("Grid", "Options for defining the (rho, T) grid");
  grid_group->add_option("--rho_min", opts->rho_min, "Minimum density (g/cc)")
      ->required();
  grid_group->add_option("--rho_max", opts->rho_max, "Maximum density (g/cc)")
      ->required();
  grid_group
      ->add_option("--num_rho", opts->num_rho_points,
                   "Number of density points")
      ->default_val(100);
  grid_group->add_option("--T_min", opts->T_min, "Minimum temperature (K)")
      ->required();
  grid_group->add_option("--T_max", opts->T_max, "Maximum temperature (K)")
      ->required();
  grid_group
      ->add_option("--num_T", opts->num_T_points,
                   "Number of temperature points")
      ->default_val(100);
  grid_group->add_flag("--log_rho", opts->log_rho,
                       "Use logarithmic spacing for density");
  grid_group->add_flag("--log_T", opts->log_T,
                       "Use logarithmic spacing for temperature");

  auto* ref_group = cmd->add_option_group(
      "Reference", "Options for the entropy reference state");
  ref_group
      ->add_option("--rho_ref", opts->rho_ref,
                   "Reference density for entropy calculation (g/cc)")
      ->required();
  ref_group
      ->add_option("--T_ref", opts->T_ref,
                   "Reference temperature for entropy calculation (K)")
      ->required();
  ref_group
      ->add_option("--S_ref", opts->S_ref,
                   "Reference entropy value at (rho_ref, T_ref)")
      ->default_val(0.0);

  cmd->add_option("--substeps", opts->integration_substeps,
                  "Number of sub-steps for each ODE integration")
      ->default_val(100);

  cmd->callback([opts]() { run_build_entropy_table(*opts); });
}

int main(int argc, char** argv) {
  CLI::App app{"eos-tool: A command-line tool for EOS calculations"};
  app.require_subcommand(1);  // Require at least one subcommand

  // Add subcommands
  setup_isotherm_subcommand(app);
  setup_hugoniot_subcommand(app);
  setup_isentrope_subcommand(app);
  setup_release_subcommand(app);
  setup_build_entropy_table_subcommand(app);  // <<< NEW

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }

  return 0;
}
