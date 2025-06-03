// src/cpp/main_phase5.cpp
#include "EquationOfStateV1.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstring> // For std::memcmp

// Helper from main_phase4.cpp
void print_eos_results(int eos_id, int istat, double P, double E, double dPdT, double dEdT, double dPdrho) {
    std::cout << std::fixed << std::setprecision(5);
    if (istat == EOS_SUCCESS) {
        std::cout << "  EOS ID " << eos_id << " SUCCESS:" << std::endl;
        std::cout << "    P      = " << P << std::endl; /* ... same as before ... */ }
    else { /* ... */ }
}


int main() {
    std::cout << "C++: Phase 5 Test Program (Pack/Unpack)" << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    EquationOfStateV1 eos_manager_rank0; // Simulates rank 0
    std::string eos_data_path = "../eos_data_dir"; // Adjust as needed

    // --- Initialize on Rank 0 ---
    std::cout << "\n--- Initializing EOS Manager on Rank 0 ---" << std::endl;
    std::vector<int> eos_to_load = {1, 2, 10000};
    eos_manager_rank0.setComplicatedEOSUseOldColdTerm(true); // Set a non-default control
    eos_manager_rank0.setUseTFDDataVer1(true);              // Explicitly use TFD v1

    int init_stat = eos_manager_rank0.initialize(eos_to_load, eos_data_path);
    if (init_stat != EOS_SUCCESS) {
        std::cerr << "Rank 0 EOS Manager initialization failed. Error code: " << init_stat << std::endl;
        return 1;
    }
    std::cout << "Rank 0 EOS Manager initialized successfully." << std::endl;

    // --- Pack Data on Rank 0 ---
    std::cout << "\n--- Packing data on Rank 0 ---" << std::endl;
    char* packed_buffer = nullptr;
    int buffer_size = 0;

    // First call: get size
    int pack_stat = eos_manager_rank0.pack_data(packed_buffer, buffer_size);
    if (pack_stat != EOS_SUCCESS || buffer_size == 0) {
        std::cerr << "Rank 0: Failed to get packed data size. Stat: " << pack_stat << std::endl;
        return 1;
    }
    std::cout << "Rank 0: Calculated packed data size: " << buffer_size << " bytes." << std::endl;

    // Allocate buffer
    packed_buffer = new char[buffer_size];

    // Second call: pack data
    int actual_packed_size = buffer_size; // Pass the allocated size
    pack_stat = eos_manager_rank0.pack_data(packed_buffer, actual_packed_size);
    if (pack_stat != EOS_SUCCESS) {
        std::cerr << "Rank 0: Failed to pack data. Stat: " << pack_stat << std::endl;
        delete[] packed_buffer;
        return 1;
    }
    if (actual_packed_size != buffer_size) {
         std::cerr << "Rank 0: Mismatch between calculated size and actual packed size!" << std::endl;
         // This might indicate an issue if actual_packed_size is less than buffer_size after packing
    }
    std::cout << "Rank 0: Data packed successfully. Actual size used: " << actual_packed_size << " bytes." << std::endl;


    // --- Simulate MPI Broadcast: Unpack on Rank > 0 ---
    std::cout << "\n--- Unpacking data on Rank > 0 (Simulated) ---" << std::endl;
    EquationOfStateV1 eos_manager_rankN; // Simulates rank > 0

    int unpack_stat = eos_manager_rankN.unpack_data(packed_buffer, actual_packed_size);
    if (unpack_stat != EOS_SUCCESS) {
        std::cerr << "Rank N: Failed to unpack data. Stat: " << unpack_stat << std::endl;
        delete[] packed_buffer;
        return 1;
    }
    std::cout << "Rank N: Data unpacked successfully." << std::endl;

    // --- Verify by Computing on Rank N and comparing with Rank 0 (Optional but good) ---
    std::cout << "\n--- Verifying by computing on Rank N ---" << std::endl;
    double rho = 1.8, T = 1200.0;
    double P_r0, E_r0, dPdT_r0, dEdT_r0, dPdrho_r0;
    double P_rN, E_rN, dPdT_rN, dEdT_rN, dPdrho_rN;
    int comp_stat_r0, comp_stat_rN;

    // Complicated EOS on Rank 0
    comp_stat_r0 = eos_manager_rank0.compute(10000, rho, T, P_r0, E_r0, dPdT_r0, dEdT_r0, dPdrho_r0);
    std::cout << "Rank 0 computation for ID 10000:" << std::endl;
    print_eos_results(10000, comp_stat_r0, P_r0, E_r0, dPdT_r0, dEdT_r0, dPdrho_r0);

    // Complicated EOS on Rank N
    comp_stat_rN = eos_manager_rankN.compute(10000, rho, T, P_rN, E_rN, dPdT_rN, dEdT_rN, dPdrho_rN);
    std::cout << "Rank N computation for ID 10000:" << std::endl;
    print_eos_results(10000, comp_stat_rN, P_rN, E_rN, dPdT_rN, dEdT_rN, dPdrho_rN);

    if (comp_stat_r0 == EOS_SUCCESS && comp_stat_rN == EOS_SUCCESS &&
        std::abs(P_r0 - P_rN) < 1e-9 && std::abs(E_r0 - E_rN) < 1e-9) { // Basic check
        std::cout << "Verification: Computations on Rank 0 and Rank N match for ID 10000." << std::endl;
    } else {
        std::cout << "Verification ERROR: Computations on Rank 0 and Rank N DO NOT match for ID 10000." << std::endl;
    }

    // Test an analytic EOS too
    rho = 1.2; T = 300.0;
    comp_stat_r0 = eos_manager_rank0.compute(1, rho, T, P_r0, E_r0, dPdT_r0, dEdT_r0, dPdrho_r0);
    comp_stat_rN = eos_manager_rankN.compute(1, rho, T, P_rN, E_rN, dPdT_rN, dEdT_rN, dPdrho_rN);
     if (comp_stat_r0 == EOS_SUCCESS && comp_stat_rN == EOS_SUCCESS &&
        std::abs(P_r0 - P_rN) < 1e-9 && std::abs(E_r0 - E_rN) < 1e-9) {
        std::cout << "Verification: Computations on Rank 0 and Rank N match for ID 1." << std::endl;
    } else {
        std::cout << "Verification ERROR: Computations on Rank 0 and Rank N DO NOT match for ID 1." << std::endl;
    }


    delete[] packed_buffer;
    std::cout << "\n---------------------------------------" << std::endl;
    std::cout << "C++: Phase 5 Test Complete." << std::endl;
    return 0;
}
