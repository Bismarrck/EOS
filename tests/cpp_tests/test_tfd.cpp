// src/cpp/main_phase3.cpp
#include "../../src/cpp/EquationOfStateV1.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

int main() {
    std::cout << "C++: Phase 3 Test Program" << std::endl;
    std::cout << "-------------------------" << std::endl;

    EquationOfStateV1 eos_instance;

    // --- Setup: Load TFD Data (from Phase 2) ---
    // Ensure TFD files are accessible (e.g., in build dir or correct relative path)
    std::string tfd_file_A = "../eos_data_dir/tfd_matA_v1.dat"; // Use the 2x3 dummy files
    std::string tfd_file_B = "../eos_data_dir/tfd_matB_v1.dat";
    std::cout << "\nLoading TFD Data..." << std::endl;
    int load_istat = eos_instance.loadTFDData(tfd_file_A, tfd_file_B);
    if (load_istat != 0) {
        std::cerr << "  Failed to load TFD data. Error: " << load_istat << ". Aborting test." << std::endl;
        return 1;
    }
    std::cout << "  TFD data loaded successfully." << std::endl;

    // --- Test Complicated EOS ---
    std::cout << "\nTesting Complicated EOS computation:" << std::endl;
    // Dummy parameters for a "material"
    std::vector<double> material_params = {10.0, 20.0, 0.5, 0.7}; // Example: 4 parameters

    double rho = 1.8;
    double T = 1200.0;
    double P, E, dPdT, dEdT, dPdrho;
    int eos_istat;

    std::cout << "Calling Complicated EOS with rho=" << rho << ", T=" << T
              << ", params[0]=" << material_params[0] << std::endl;

    eos_istat = eos_instance.computeComplicatedEOS(material_params, rho, T, P, E, dPdT, dEdT, dPdrho);

    std::cout << std::fixed << std::setprecision(5);
    if (eos_istat == 0) {
        std::cout << "  Complicated EOS success:" << std::endl;
        std::cout << "    P      = " << P << std::endl;
        std::cout << "    E      = " << E << std::endl;
        std::cout << "    dPdT   = " << dPdT << std::endl;
        std::cout << "    dEdT   = " << dEdT << std::endl;
        std::cout << "    dPdrho = " << dPdrho << std::endl;
    } else {
        std::cout << "  Complicated EOS failed. istat=" << eos_istat << std::endl;
    }

    // Test case with insufficient parameters for Complicated_EOS stub
    std::cout << "\nTesting Complicated EOS with insufficient params:" << std::endl;
    std::vector<double> short_params = {5.0}; // Only 1 param, stub expects at least 2
     eos_istat = eos_instance.computeComplicatedEOS(short_params, rho, T, P, E, dPdT, dEdT, dPdrho);
    if (eos_istat != 0) {
        std::cout << "  Complicated EOS failed as expected with insufficient params. istat=" << eos_istat << std::endl;
    } else {
        std::cout << "  Complicated EOS succeeded unexpectedly with insufficient params." << std::endl;
    }


    std::cout << "\n-------------------------" << std::endl;
    std::cout << "C++: Phase 3 Test Complete." << std::endl;
    return 0;
}
