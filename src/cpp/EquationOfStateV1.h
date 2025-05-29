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


    // --- Control Variable Management ---
    int setComplicatedEOSUseOldColdTerm(bool value);
    int setUseTFDDataVer1(bool value);
    // bool getUseTFDDataVer1(int* istat_out); // If getter is needed

    // --- EOS Computation (very basic for Phase 1) ---
    // This will evolve significantly into the main compute method
    int computeAirEOS(double rho, double T, double& P, double& E);
    int computeCarbonEOS(double rho, double T, double& P, double& E);

    // --- Methods to be implemented in later phases ---
    // int initialize(const std::vector<int>& eos_id_list, const std::string& eos_data_dir, bool use_tfd_version_1);
    // int check_eos_data_dir(/* ... */);
    // int pack(char *&buffer, int &buffer_size);
    // int unpack(const char *buffer, int buffer_size);
    // int compute(int eos_id, double rho, double T, double& P, double& E, ...); // Main compute
    // void free();

private:
    // For Phase 1, not much state is managed yet.
    // This will grow to include TFD data, material params, etc.

    // Example: mapping an internal C++ EOSType to a function pointer or lambda
    // for dispatching compute calls. (More for Phase 4)
    // using AnalyticEOSFunc = int(EquationOfStateV1::*)(double, double, double&, double&);
    // std::map<EOSType, AnalyticEOSFunc> analytic_dispatch_table;
};

#endif // EQUATIONOFSTATEV1_H