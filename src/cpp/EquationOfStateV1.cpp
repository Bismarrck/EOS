// src/cpp/EquationOfStateV1.cpp
#include "EquationOfStateV1.h"
#include <iostream>     // For debug prints
#include <fstream>      // for std::ifstream
#include <sstream>      // For std::istringstream
#include <iomanip>      // For std::fixed, std::setprecision in main, can be useful here too


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

// --- TFD Data Management ---

// Helper function to read a single TFD matrix
int EquationOfStateV1::readTFDMatrixFromFile(const std::string& file_path,
                                             std::vector<double>& matrix_data,
                                             int& N1, int& N2) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open TFD file: " << file_path << std::endl;
        return 1; // Error code for file not found/open
    }

    std::string line;
    // Read dimensions
    if (!std::getline(infile, line)) {
        std::cerr << "Error: Could not read dimensions line from " << file_path << std::endl;
        return 2; // Error code for read error
    }
    std::istringstream dim_ss(line);
    if (!(dim_ss >> N1 >> N2)) {
        std::cerr << "Error: Could not parse dimensions N1, N2 from " << file_path << std::endl;
        return 3; // Error code for parse error
    }

    if (N1 <= 0 || N2 <= 0) {
        std::cerr << "Error: Invalid dimensions N1=" << N1 << ", N2=" << N2 << " in " << file_path << std::endl;
        return 4; // Error code for invalid dimensions
    }

    matrix_data.clear();
    matrix_data.reserve(static_cast<size_t>(N1) * N2); // C++20: use static_cast<size_t> for N1*N2

    double value;
    while (std::getline(infile, line)) {
        std::istringstream val_ss(line);
        while (val_ss >> value) {
            matrix_data.push_back(value);
        }
    }

    if (matrix_data.size() != static_cast<size_t>(N1 * N2)) {
        std::cerr << "Error: Read " << matrix_data.size() << " values, but expected "
                  << (N1 * N2) << " for dimensions " << N1 << "x" << N2
                  << " in file " << file_path << std::endl;
        return 5; // Error code for data size mismatch
    }

    std::cout << "Successfully read " << N1 << "x" << N2 << " matrix from " << file_path << std::endl;
    return 0; // Success
}

int EquationOfStateV1::loadTFDData(const std::string& file_A_path, const std::string& file_B_path) {
    int istat_A = readTFDMatrixFromFile(file_A_path, tfd_data.matrix_A, tfd_data.N1_A, tfd_data.N2_A);
    if (istat_A != 0) {
        tfd_data.initialized = false;
        return 10 + istat_A; // Error loading matrix A
    }

    int istat_B = readTFDMatrixFromFile(file_B_path, tfd_data.matrix_B, tfd_data.N1_B, tfd_data.N2_B);
    if (istat_B != 0) {
        tfd_data.initialized = false;
        return 20 + istat_B; // Error loading matrix B
    }

    // Optional: Add checks if N1_A, N2_A, N1_B, N2_B match expected (e.g., 100, 50)
    // For now, we assume the files provide dimensions that Fortran expects or can handle.
    // The c_tfd_eos_compute wrapper will pass these dimensions.

    tfd_data.initialized = true;
    std::cout << "C++: TFD data loaded and initialized." << std::endl;
    return 0; // Success
}

// --- TFD Computation Call ---
int EquationOfStateV1::computeTFD(double rho, double T, double& result_x, double& result_y) {
    if (!tfd_data.initialized) {
        std::cerr << "C++ Error: TFD data not initialized before calling computeTFD." << std::endl;
        return -1; // Or some other error code
    }

    // Assuming the Fortran routine expects 100x50 or can handle variable sizes
    // passed via dim1_a, dim2_a etc.
    // If Fortran is strictly 100x50, add checks here:
    // if (tfd_data.N1_A != 100 || tfd_data.N2_A != 50 || /* similar for B */) {
    //     std::cerr << "C++ Error: Loaded TFD matrix dimensions mismatch expected 100x50." << std::endl;
    //     return -2;
    // }

    int istat_fortran;
    c_tfd_eos_compute(rho, T,
                      tfd_data.matrix_A.data(), tfd_data.matrix_B.data(),
                      tfd_data.N1_A, tfd_data.N2_A,
                      tfd_data.N1_B, tfd_data.N2_B,
                      &result_x, &result_y, &istat_fortran);

    if (istat_fortran != 0) {
        std::cerr << "C++ Error: c_tfd_eos_compute failed. Fortran istat: " << istat_fortran << std::endl;
    }
    return istat_fortran;
}
