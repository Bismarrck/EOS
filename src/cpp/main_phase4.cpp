// src/cpp/main_phase4.cpp
#include "EquationOfStateV1.h" // Assuming EquationOfStateV1.h is in the same dir or include path set by CMake
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

// Helper to print results
void print_eos_results(int eos_id, int istat, double P, double E, double dPdT, double dEdT, double dPdrho) {
    std::cout << std::fixed << std::setprecision(5);
    if (istat == EOS_SUCCESS) {
        std::cout << "  EOS ID " << eos_id << " SUCCESS:" << std::endl;
        std::cout << "    P      = " << P << std::endl;
        std::cout << "    E      = " << E << std::endl;
        std::cout << "    dPdT   = " << dPdT << std::endl;
        std::cout << "    dEdT   = " << dEdT << std::endl;
        std::cout << "    dPdrho = " << dPdrho << std::endl;
    } else {
        std::cout << "  EOS ID " << eos_id << " FAILED. istat=" << istat << std::endl;
    }
}


int main() {
    std::cout << "C++: Phase 4 Test Program" << std::endl;
    std::cout << "-------------------------" << std::endl;

    EquationOfStateV1 eos_manager;

    // Define eos_data_dir. Adjust path as necessary.
    // For testing, often easiest to run executable from project root, or copy data.
    std::string eos_data_path = "../../eos_data_dir"; // If build is in Project/build
    // If running from build dir: std::string eos_data_path = "../eos_data_dir";
    // Or make sure your dummy "eos_data_dir" from project root is copied to the executable's CWD.
    // For simplicity, I will assume that when the program runs, "eos_data_dir" is a subdirectory
    // of the current working directory. You may need to adjust this path.
    // If executable is in "EquationOfStateProject/build/", then eos_data_path = "../eos_data_dir/"
    // For now, let's assume we run it from "EquationOfStateProject/" and data is in "eos_data_dir/"
    // This usually means copying files or configuring CMake to run from a specific dir.
    // Simplest for now: current dir is `build`, data is `../eos_data_dir`
    // If your executable is in `build/`, try:
    // std::string eos_data_path = "../eos_data_dir";
    // If your executable is in `build/src/cpp/` try:
    // std::string eos_data_path = "../../../eos_data_dir";
    // Let's set a placeholder that often works if the CWD is the build directory:
    eos_data_path = "../eos_data_dir"; // Common if executable is in build/


    // --- Test check_eos_data_dir ---
    std::cout << "\n--- Checking EOS Data Directory ---" << std::endl;
    // For check_eos_data_dir, list IDs you expect files for (excluding registered analytic if files are optional)
    std::vector<int> ids_to_check_files_for = {10000}; // Check the complicated one
    int check_stat = eos_manager.check_eos_data_dir(eos_data_path, ids_to_check_files_for);
    if (check_stat != EOS_SUCCESS) {
        std::cerr << "EOS Data Directory check failed. Aborting." << std::endl;
        return 1;
    }
    std::cout << "EOS Data Directory check passed." << std::endl;


    // --- Test Initialization ---
    std::cout << "\n--- Initializing EOS Manager (TFD v1 by default) ---" << std::endl;
    std::vector<int> eos_to_load = {
        1,       // Analytic Air (registered, file mat000/eos00001.dat may exist with #TYPE)
        2,       // Analytic Carbon (registered, file mat000/eos00002.dat may exist with #TYPE)
        10000    // Complicated (file mat100/eos10000.dat with #TYPE and params)
    };
    int init_stat = eos_manager.initialize(eos_to_load, eos_data_path);
    if (init_stat != EOS_SUCCESS) {
        std::cerr << "EOS Manager initialization failed. Error code: " << init_stat << std::endl;
        return 1;
    }
    std::cout << "EOS Manager initialized successfully." << std::endl;

    // --- Test Computation ---
    std::cout << "\n--- Testing Computations ---" << std::endl;
    double rho, T, P, E, dPdT, dEdT, dPdrho;
    int comp_stat;

    // Test Analytic Air (ID 1)
    rho = 1.2; T = 300.0;
    std::cout << "Computing for EOS ID 1 (Air Analytic) with rho=" << rho << ", T=" << T << std::endl;
    comp_stat = eos_manager.compute(1, rho, T, P, E, dPdT, dEdT, dPdrho);
    print_eos_results(1, comp_stat, P, E, dPdT, dEdT, dPdrho);

    // Test Analytic Carbon (ID 2) - success
    rho = 2.0; T = 1000.0;
    std::cout << "Computing for EOS ID 2 (Carbon Analytic) with rho=" << rho << ", T=" << T << std::endl;
    comp_stat = eos_manager.compute(2, rho, T, P, E, dPdT, dEdT, dPdrho);
    print_eos_results(2, comp_stat, P, E, dPdT, dEdT, dPdrho);

    // Test Analytic Carbon (ID 2) - error case (if applicable from its Fortran logic)
    rho = -1.0; T = 500.0; // Expecting error from carbon_eos_2001
    std::cout << "Computing for EOS ID 2 (Carbon Analytic) with rho=" << rho << ", T=" << T << " (expecting error)" << std::endl;
    comp_stat = eos_manager.compute(2, rho, T, P, E, dPdT, dEdT, dPdrho);
    print_eos_results(2, comp_stat, P, E, dPdT, dEdT, dPdrho);


    // Test Complicated EOS (ID 10000)
    rho = 1.8; T = 1200.0;
    std::cout << "Computing for EOS ID 10000 (Complicated) with rho=" << rho << ", T=" << T << std::endl;
    comp_stat = eos_manager.compute(10000, rho, T, P, E, dPdT, dEdT, dPdrho);
    print_eos_results(10000, comp_stat, P, E, dPdT, dEdT, dPdrho);


    // --- Test with TFD v2 ---
    std::cout << "\n--- Re-initializing with TFD v2 ---" << std::endl;
    EquationOfStateV1 eos_manager_v2;
    eos_manager_v2.setUseTFDDataVer1(false); // Request TFD v2
    init_stat = eos_manager_v2.initialize(eos_to_load, eos_data_path);
    if (init_stat != EOS_SUCCESS) {
        std::cerr << "EOS Manager (v2) initialization failed. Error code: " << init_stat << std::endl;
        return 1;
    }
    std::cout << "EOS Manager (v2) initialized successfully." << std::endl;
    std::cout << "Computing for EOS ID 10000 (Complicated with TFD v2) with rho=" << rho << ", T=" << T << std::endl;
    comp_stat = eos_manager_v2.compute(10000, rho, T, P, E, dPdT, dEdT, dPdrho);
    print_eos_results(10000, comp_stat, P, E, dPdT, dEdT, dPdrho);


    std::cout << "\n-------------------------" << std::endl;
    std::cout << "C++: Phase 4 Test Complete." << std::endl;
    return 0;
}
