#include "../../third_party/doctest/doctest.h"
#include "../../src/cpp/EquationOfStateV1.h"
#include <vector>
#include <string>


TEST_CASE("EquationOfStateV1 Pack/Unpack Tests") {
    EquationOfStateV1 eos_manager_rank0;
    std::string eos_data_path = "../../../eos_data_dir"; // Adjust path

    INFO("Initializing EOS Manager on Rank 0 for Pack/Unpack Test");
    std::vector<int> eos_to_load = {1, 2, 10000};
    eos_manager_rank0.setComplicatedEOSUseOldColdTerm(true);
    eos_manager_rank0.setUseTFDDataVer1(true);

    int init_stat = eos_manager_rank0.initialize(eos_to_load, eos_data_path);
    REQUIRE(init_stat == EOS_SUCCESS);

    // --- Pack Data ---
    char* packed_buffer = nullptr;
    int buffer_size = 0;

    INFO("Getting packed data size");
    int pack_stat = eos_manager_rank0.pack_data(packed_buffer, buffer_size);
    REQUIRE(pack_stat == EOS_SUCCESS);
    REQUIRE(buffer_size > 0);
    packed_buffer = new char[buffer_size];

    INFO("Packing data");
    int actual_packed_size = buffer_size;
    pack_stat = eos_manager_rank0.pack_data(packed_buffer, actual_packed_size);
    REQUIRE(pack_stat == EOS_SUCCESS);
    CHECK(actual_packed_size == buffer_size);

    // --- Unpack Data ---
    EquationOfStateV1 eos_manager_rankN;
    INFO("Unpacking data");
    int unpack_stat = eos_manager_rankN.unpack_data(packed_buffer, actual_packed_size);
    REQUIRE(unpack_stat == EOS_SUCCESS);

    // --- Verification ---
    SUBCASE("Verify by computing after unpack") {
        double rho = 1.8, T = 1200.0;
        double P_r0, E_r0, dPdT_r0, dEdT_r0, dPdrho_r0;
        double P_rN, E_rN, dPdT_rN, dEdT_rN, dPdrho_rN;
        int comp_stat_r0, comp_stat_rN;

        // Complicated EOS
        INFO("Comparing Complicated EOS (ID 10000) results");
        comp_stat_r0 = eos_manager_rank0.compute(10000, rho, T, P_r0, E_r0, dPdT_r0, dEdT_r0, dPdrho_r0);
        comp_stat_rN = eos_manager_rankN.compute(10000, rho, T, P_rN, E_rN, dPdT_rN, dEdT_rN, dPdrho_rN);
        REQUIRE(comp_stat_r0 == EOS_SUCCESS);
        REQUIRE(comp_stat_rN == EOS_SUCCESS);
        CHECK(P_r0 == doctest::Approx(P_rN));
        CHECK(E_r0 == doctest::Approx(E_rN));
        CHECK(dPdT_r0 == doctest::Approx(dPdT_rN));
        CHECK(dEdT_r0 == doctest::Approx(dEdT_rN));
        CHECK(dPdrho_r0 == doctest::Approx(dPdrho_rN));

        // Analytic EOS
        rho = 1.2; T = 300.0;
        INFO("Comparing Analytic EOS (ID 1) results");
        comp_stat_r0 = eos_manager_rank0.compute(1, rho, T, P_r0, E_r0, dPdT_r0, dEdT_r0, dPdrho_r0);
        comp_stat_rN = eos_manager_rankN.compute(1, rho, T, P_rN, E_rN, dPdT_rN, dEdT_rN, dPdrho_rN);
        REQUIRE(comp_stat_r0 == EOS_SUCCESS);
        REQUIRE(comp_stat_rN == EOS_SUCCESS);
        CHECK(P_r0 == doctest::Approx(P_rN));
        CHECK(E_r0 == doctest::Approx(E_rN));
    }
    delete[] packed_buffer;
}
