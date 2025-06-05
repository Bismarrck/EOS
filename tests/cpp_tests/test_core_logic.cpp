#include "../../third_party/doctest/doctest.h"
#include "../../src/cpp/EquationOfStateV1.h"
#include <iostream>
#include <vector>
#include <string>


TEST_CASE("EquationOfStateV1 Core Logic Tests") {
    std::cout << "C++: Phase 4 Test Program" << std::endl;
    std::cout << "-------------------------" << std::endl;

    EquationOfStateV1 eos_manager;
    std::string eos_data_path = "../../../eos_data_dir";

    // --- Test check_eos_data_dir ---
    SUBCASE("Check EOS Data Directory") {
        INFO("Checking EOS Data Directory: " << eos_data_path);
        std::vector<int> ids_to_check_files_for = {10000};
        int check_stat = eos_manager.check_eos_data_dir(eos_data_path, ids_to_check_files_for);
        REQUIRE(check_stat == EOS_SUCCESS); // REQUIRE will stop test case on failure
    }

    // --- Test Initialization ---
    SUBCASE("Initialize EOS Manager (TFD v1)") {
        INFO("Initializing EOS Manager (TFD v1 by default)");
        std::vector<int> eos_to_load = {1, 2, 10000};
        int init_stat = eos_manager.initialize(eos_to_load, eos_data_path);
        REQUIRE(init_stat == EOS_SUCCESS);

        // --- Test Computation after init ---
        double rho, T, P, E, dPdT, dEdT, dPdrho;
        int comp_stat;

        // Test Analytic Air (ID 1)
        rho = 1.2; T = 300.0;
        INFO("Computing for EOS ID 1 (Air Analytic) with rho=" << rho << ", T=" << T);
        comp_stat = eos_manager.compute(1, rho, T, P, E, dPdT, dEdT, dPdrho);
        CHECK(comp_stat == EOS_SUCCESS);
        CHECK(P == doctest::Approx(360.0)); // Example check for expected value
        CHECK(E == doctest::Approx(540.0));

        // Test Analytic Carbon (ID 2) - success
        rho = 2.0; T = 1000.0;
        INFO("Computing for EOS ID 2 (Carbon Analytic) with rho=" << rho << ", T=" << T);
        comp_stat = eos_manager.compute(2, rho, T, P, E, dPdT, dEdT, dPdrho);
        CHECK(comp_stat == EOS_SUCCESS);
        CHECK(P == doctest::Approx(4000.0));
        CHECK(E == doctest::Approx(5000.0));


        // Test Analytic Carbon (ID 2) - error case
        rho = -1.0; T = 500.0;
        INFO("Computing for EOS ID 2 (Carbon Analytic) with rho=" << rho << ", T=" << T << " (expecting error)");
        comp_stat = eos_manager.compute(2, rho, T, P, E, dPdT, dEdT, dPdrho);
        CHECK(comp_stat != EOS_SUCCESS); // Or check for specific error code (e.g., 101 from Fortran)
        CHECK(comp_stat == 101); // If 101 is the specific Fortran error

        // Test Complicated EOS (ID 10000)
        rho = 1.8; T = 1200.0;
        INFO("Computing for EOS ID 10000 (Complicated) with rho=" << rho << ", T=" << T);
        comp_stat = eos_manager.compute(10000, rho, T, P, E, dPdT, dEdT, dPdrho);
        CHECK(comp_stat == EOS_SUCCESS);
        // Add CHECKs for P, E, dPdT etc. based on your dummy TFD data and logic
        // e.g. CHECK(P == doctest::Approx(19.98000));
    }

    SUBCASE("Initialize and Compute with TFD v2") {
        EquationOfStateV1 eos_manager_v2; // New instance for this subcase
        eos_manager_v2.setUseTFDDataVer1(false);
        INFO("Initializing with TFD v2");
        std::vector<int> eos_to_load = {10000}; // Just load one for this test
        int init_stat = eos_manager_v2.initialize(eos_to_load, eos_data_path);
        REQUIRE(init_stat == EOS_SUCCESS);

        double rho = 1.8, T = 1200.0;
        double P, E, dPdT, dEdT, dPdrho;
        int comp_stat;
        INFO("Computing for EOS ID 10000 (Complicated with TFD v2) with rho=" << rho << ", T=" << T);
        comp_stat = eos_manager_v2.compute(10000, rho, T, P, E, dPdT, dEdT, dPdrho);
        CHECK(comp_stat == EOS_SUCCESS);
        // Add CHECKs for P, E for TFD v2. Values will be different.
        // Example if tfd_ver2.dat had different first elements:
        // P_v2 = params_in(1) * rho + tfd_v2_result_x
        // CHECK(P == doctest::Approx(expected_P_with_TFDv2));
    }
}
