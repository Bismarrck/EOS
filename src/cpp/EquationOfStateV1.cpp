// src/cpp/EquationOfStateV1.cpp
#include "EquationOfStateV1.h"
#include <iostream> // For debug prints

EquationOfStateV1::EquationOfStateV1() {
    std::cout << "C++: EquationOfStateV1 instance created." << std::endl;
    // Initialize dispatch table if used for analytic functions later
    // analytic_dispatch_table[EOSType::ANALYTIC_AIR_2000] = &EquationOfStateV1::computeAirEOS_internal_wrapper;
}

EquationOfStateV1::~EquationOfStateV1() {
    std::cout << "C++: EquationOfStateV1 instance destroyed." << std::endl;
}

// --- Control Variable Management ---
int EquationOfStateV1::setComplicatedEOSUseOldColdTerm(bool value) {
    int istat;
    c_set_complicated_eos_use_old_cold_term(value, &istat);
    if (istat != 0) {
        std::cerr << "C++ Error: Failed to set complicated_eos_use_old_cold_term. Fortran istat: " << istat << std::endl;
    }
    return istat;
}

int EquationOfStateV1::setUseTFDDataVer1(bool value) {
    int istat;
    c_set_use_tfd_data_ver1(value, &istat);
    if (istat != 0) {
        std::cerr << "C++ Error: Failed to set use_tfd_data_ver1. Fortran istat: " << istat << std::endl;
    }
    return istat;
}

// --- EOS Computation (very basic for Phase 1) ---
int EquationOfStateV1::computeAirEOS(double rho, double T, double& P, double& E) {
    int istat;
    c_air_eos_2000(rho, T, &P, &E, &istat); // P, E passed by pointer
    if (istat != 0) {
        std::cerr << "C++ Error: c_air_eos_2000 failed. Fortran istat: " << istat << std::endl;
    }
    return istat;
}

int EquationOfStateV1::computeCarbonEOS(double rho, double T, double& P, double& E) {
    int istat;
    c_carbon_eos_2001(rho, T, &P, &E, &istat); // P, E passed by pointer
    if (istat != 0) {
        std::cerr << "C++ Error: c_carbon_eos_2001 failed. Fortran istat: " << istat << std::endl;
        // Potentially map Fortran istat to C++ exceptions or more detailed error codes here
    }
    return istat;
}