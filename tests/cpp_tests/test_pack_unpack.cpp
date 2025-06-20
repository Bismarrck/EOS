#include <string>
#include <vector>

#include "../../src/cpp/EquationOfStateV1.h"
#include "../../src/cpp/eos_error_codes.h"
#include "../../third_party/doctest/doctest.h"
#include "materials/MaterialEOS.h"

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
        ComputeResult result_r0, result_rN;

        // Complicated EOS
        INFO("Comparing Complicated EOS (ID 10000) results");
        result_r0 = eos_manager_rank0.compute(10000, rho, T);
        result_rN = eos_manager_rankN.compute(10000, rho, T);
        REQUIRE(result_r0.istat == EOS_SUCCESS);
        REQUIRE(result_rN.istat == EOS_SUCCESS);
        CHECK(result_r0.P == doctest::Approx(result_rN.P));
        CHECK(result_r0.E == doctest::Approx(result_rN.E));
        CHECK(result_r0.dPdT == doctest::Approx(result_rN.dPdT));
        CHECK(result_r0.dEdT == doctest::Approx(result_rN.dEdT));
        CHECK(result_r0.dPdrho == doctest::Approx(result_rN.dPdrho));

        // Analytic EOS
        rho = 1.2; T = 300.0;
        INFO("Comparing Analytic EOS (ID 1) results");
        result_r0 = eos_manager_rank0.compute(1, rho, T);
        result_rN = eos_manager_rankN.compute(1, rho, T);
        REQUIRE(result_r0.istat == EOS_SUCCESS);
        REQUIRE(result_rN.istat == EOS_SUCCESS);
        CHECK(result_r0.P == doctest::Approx(result_rN.P));
        CHECK(result_r0.E == doctest::Approx(result_rN.E));
    }
    delete[] packed_buffer;
}
