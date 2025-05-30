// src/cpp/EquationOfStateV1.h
#ifndef EQUATIONOFSTATEV1_H
#define EQUATIONOFSTATEV1_H

#include <vector>
#include <string>
#include <map>

// Forward declare Fortran-interfacing functions (actual definitions in .cpp)
// These are the C names given in BIND(C, name='...')
extern "C" {
    void c_air_eos_2000(double rho, double T, double* P, double* E, int* istat_out);
    void c_carbon_eos_2001(double rho, double T, double* P, double* E, int* istat_out);

    void c_set_complicated_eos_use_old_cold_term(bool value, int* istat_out);
    void c_set_use_tfd_data_ver1(bool value, int* istat_out);
    // bool c_get_use_tfd_data_ver1(int* istat_out); // If getter is implemented
}

extern "C" {
    // Assuming fixed dimensions 100x50 for now in the C wrapper arguments
    // If dimensions can vary and need to be passed, the signature would change.
    void c_tfd_eos_compute(double rho, double T,
                           const double* matrix_a_flat, const double* matrix_b_flat, // Pass as flat arrays
                           int dim1_a, int dim2_a, // Dimensions of A
                           int dim1_b, int dim2_b, // Dimensions of B
                           double* tfd_result_x, double* tfd_result_y,
                           int* istat_out); // Good practice to always have an istat
}

// Structure to hold TFD matrix data in C++
struct TFDMatrices {
    std::vector<double> matrix_A; // Stored flat
    std::vector<double> matrix_B; // Stored flat
    int N1_A = 0, N2_A = 0;       // Dimensions for A
    int N1_B = 0, N2_B = 0;       // Dimensions for B
    bool initialized = false;

    // Helper to access elements as if 2D (for C++ side if needed, Fortran gets flat)
    // Example for matrix_A: A(row, col) = matrix_A[row * N2_A + col]
};

// New extern "C" for Complicated EOS
extern "C" {
    void c_complicated_eos(const double* params, int num_params, double rho, double T,
                           const double* matrix_a_flat, int n1_a, int n2_a,
                           const double* matrix_b_flat, int n1_b, int n2_b,
                           double* P, double* E, double* dPdT, double* dEdT, double* dPdrho,
                           int* istat_out);
}


class EquationOfStateV1 {
public:
    EquationOfStateV1();
    ~EquationOfStateV1();

    // Enum for identifying EOS types (will be used more in Phase 4)
    enum class EOSType {
        UNKNOWN,
        ANALYTIC_AIR_2000,      // Example specific analytic type
        ANALYTIC_CARBON_2001,   // Example specific analytic type
        COMPLICATED_GENERIC     // Placeholder for complicated type
    };

    // Placeholder for material data (will be expanded)
    struct MaterialInfo {
        EOSType type;
        int eos_id; // The unique ID like 10000, 10001, or a special ID for analytic
        // For analytic, we might just use the EOSType to dispatch
    };

    // --- TFD Data Management ---
    // Loads TFD data from specified files.
    // For Phase 2, we might load one set. Later, initialize might choose.
    int loadTFDData(const std::string& file_A_path, const std::string& file_B_path);

    // --- TFD Computation Call (for testing Phase 2 directly) ---
    // In reality, Complicated_EOS would trigger this.
    int computeTFD(double rho, double T, double& result_x, double& result_y);


    // --- Control Variable Management ---
    int setComplicatedEOSUseOldColdTerm(bool value);
    int setUseTFDDataVer1(bool value);
    // bool getUseTFDDataVer1(int* istat_out); // If getter is needed

    // --- EOS Computation (very basic for Phase 1) ---
    // This will evolve significantly into the main compute method
    int computeAirEOS(double rho, double T, double& P, double& E);
    int computeCarbonEOS(double rho, double T, double& P, double& E);

    // --- Complicated EOS Computation Call (for testing Phase 3 directly) ---
    // This will be part of the main `compute(eos_id, ...)` method in Phase 4.
    int computeComplicatedEOS(const std::vector<double>& params, // Material-specific params
                              double rho, double T,
                              double& P, double& E, double& dPdT, double& dEdT, double& dPdrho);

    // --- Methods to be implemented in later phases ---
    // int initialize(const std::vector<int>& eos_id_list, const std::string& eos_data_dir);
    // int check_eos_data_dir(/* ... */);
    // int pack(char *&buffer, int &buffer_size);
    // int unpack(const char *buffer, int buffer_size);
    // int compute(int eos_id, double rho, double T, double& P, double& E, ...); // Main compute
    // void free();

private:

    // --- TFD ---
    TFDMatrices tfd_data; // Member to store loaded TFD matrices

    // Helper function to read a single TFD matrix from a file
    // Returns 0 on success, error code on failure.
    static int readTFDMatrixFromFile(const std::string& file_path,
                                     std::vector<double>& matrix_data,
                                     int& N1, int& N2);

    // Example: mapping an internal C++ EOSType to a function pointer or lambda
    // for dispatching compute calls. (More for Phase 4)
    // using AnalyticEOSFunc = int(EquationOfStateV1::*)(double, double, double&, double&);
    // std::map<EOSType, AnalyticEOSFunc> analytic_dispatch_table;
};

#endif // EQUATIONOFSTATEV1_H