// src/cpp/main_phase1.cpp
#include "EquationOfStateV1.h" // Assuming EquationOfStateV1.h is in the same dir or include path set by CMake
#include <iostream>
#include <iomanip> // For std::fixed, std::setprecision

// Keep Phase 0 test functions for now, just for linkage check if desired
extern "C" {
    void greet_from_fortran(const char* message);
    int add_numbers_fortran(int a, int b, int* istat_out);
}


int main() {
    std::cout << "C++: Phase 1 Test Program" << std::endl;
    std::cout << "-------------------------" << std::endl;

    std::cout << "\n--- Testing EquationOfStateV1 ---" << std::endl;
    EquationOfStateV1 eos_instance;

    // Test setting control variables
    std::cout << "\nSetting control variables:" << std::endl;
    int istat_ctrl;
    istat_ctrl = eos_instance.setComplicatedEOSUseOldColdTerm(true);
    std::cout << "C++: setComplicatedEOSUseOldColdTerm called. C++ istat: " << istat_ctrl << std::endl;
    istat_ctrl = eos_instance.setUseTFDDataVer1(false);
    std::cout << "C++: setUseTFDDataVer1 called. C++ istat: " << istat_ctrl << std::endl;
    istat_ctrl = eos_instance.setUseTFDDataVer1(true); // Set back
    std::cout << "C++: setUseTFDDataVer1 called. C++ istat: " << istat_ctrl << std::endl;


    // Test EOS computations
    std::cout << "\nTesting EOS computations:" << std::endl;
    double rho, T, P, E;
    int eos_istat;

    // Air EOS
    rho = 1.2; T = 300.0;
    std::cout << "Calling Air EOS with rho=" << rho << ", T=" << T << std::endl;
    eos_istat = eos_instance.computeAirEOS(rho, T, P, E);
    std::cout << std::fixed << std::setprecision(5);
    if (eos_istat == 0) {
        std::cout << "  Air EOS success: P=" << P << ", E=" << E << std::endl;
    } else {
        std::cout << "  Air EOS failed. istat=" << eos_istat << std::endl;
    }

    // Carbon EOS - success case
    rho = 2.0; T = 1000.0;
    std::cout << "Calling Carbon EOS with rho=" << rho << ", T=" << T << std::endl;
    eos_istat = eos_instance.computeCarbonEOS(rho, T, P, E);
    if (eos_istat == 0) {
        std::cout << "  Carbon EOS success: P=" << P << ", E=" << E << std::endl;
    } else {
        std::cout << "  Carbon EOS failed. istat=" << eos_istat << std::endl;
    }

    // Carbon EOS - error case (e.g., negative density)
    rho = -1.0; T = 500.0;
    std::cout << "Calling Carbon EOS with rho=" << rho << ", T=" << T << " (expecting error)" << std::endl;
    eos_istat = eos_instance.computeCarbonEOS(rho, T, P, E);
    if (eos_istat == 0) {
        std::cout << "  Carbon EOS success (unexpected): P=" << P << ", E=" << E << std::endl;
    } else {
        std::cout << "  Carbon EOS failed as expected. istat=" << eos_istat << std::endl;
    }

    std::cout << "\n-------------------------" << std::endl;
    std::cout << "C++: Phase 1 Test Complete." << std::endl;
    return 0;
}
