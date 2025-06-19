#include <algorithm>
#include <fstream>
#include <iomanip>  // For std::stod in a robust way if needed
#include <sstream>
#include <string>
#include <vector>

#include "../../third_party/doctest/doctest.h"
#include "EquationOfStateV1.h"
#include "eos_error_codes.h"
#include "materials/MaterialEOS.h"
#include "utils/string_utils.h"

struct NumericalTestCase {
  int eos_id = -1;
  double rho = 0.0;
  double T = 0.0;
  double ref_P = 0.0;
  double ref_E = 0.0;
  double ref_dPdT = 0.0;
  double ref_dEdT = 0.0;
  double ref_dPdrho = 0.0;
  std::string description;
  int line_num = -1;  // For error reporting

  // Configuration for this test
  bool use_tfd_v1_for_ref = true;
};

// Function to load test cases from a CSV file
// (Simplified CSV parsing, assumes no commas within fields, handles '#'
// comments)
std::vector<NumericalTestCase> load_numerical_test_cases(
    const std::string& filepath) {
  std::vector<NumericalTestCase> cases;
  std::ifstream file(filepath.c_str());
  std::string line;
  int current_line_num = 0;

  if (!file.is_open()) {
    FAIL("Failed to open reference data file: " << filepath);  // doctest FAIL
    return cases;                                              // Return empty
  }

  // Skip header line if present (optional, can be adapted)
  std::getline(file, line);
  current_line_num++;

  while (std::getline(file, line)) {
    current_line_num++;
    std::string processed_line = line;
    size_t comment_pos = processed_line.find('#');
    if (comment_pos != std::string::npos) {
      processed_line = processed_line.substr(0, comment_pos);
    }
    processed_line = EOStringUtils::trim_string(
        processed_line);  // Use trim_string from EquationOfStateV1.cpp context

    if (processed_line.empty()) continue;

    std::stringstream ss(processed_line);
    std::string field;
    NumericalTestCase tc;
    tc.line_num = current_line_num;

    try {
      // eos_id
      if (!std::getline(ss, field, ','))
        throw std::runtime_error("Missing eos_id");
      tc.eos_id = std::stoi(EOStringUtils::trim_string(field));
      // rho
      if (!std::getline(ss, field, ','))
        throw std::runtime_error("Missing rho");
      EOStringUtils::string_to_double_fortran_compat(
          EOStringUtils::trim_string(field),
          tc.rho);  // Use our D-exponent aware converter
      // T
      if (!std::getline(ss, field, ',')) throw std::runtime_error("Missing T");
      EOStringUtils::string_to_double_fortran_compat(
          EOStringUtils::trim_string(field), tc.T);
      // ref_P
      if (!std::getline(ss, field, ','))
        throw std::runtime_error("Missing ref_P");
      EOStringUtils::string_to_double_fortran_compat(
          EOStringUtils::trim_string(field), tc.ref_P);
      // ref_E
      if (!std::getline(ss, field, ','))
        throw std::runtime_error("Missing ref_E");
      EOStringUtils::string_to_double_fortran_compat(
          EOStringUtils::trim_string(field), tc.ref_E);
      // ref_dPdT
      if (!std::getline(ss, field, ','))
        throw std::runtime_error("Missing ref_dPdT");
      EOStringUtils::string_to_double_fortran_compat(
          EOStringUtils::trim_string(field), tc.ref_dPdT);
      // ref_dEdT
      if (!std::getline(ss, field, ','))
        throw std::runtime_error("Missing ref_dEdT");
      EOStringUtils::string_to_double_fortran_compat(
          EOStringUtils::trim_string(field), tc.ref_dEdT);
      // ref_dPdrho
      if (!std::getline(ss, field, ','))
        throw std::runtime_error("Missing ref_dPdrho");
      EOStringUtils::string_to_double_fortran_compat(
          EOStringUtils::trim_string(field), tc.ref_dPdrho);
      // description (optional, rest of the line)
      if (std::getline(ss, field, ',')) {  // if there's more after dPdrho
        tc.description = EOStringUtils::trim_string(field);
      } else {
        tc.description = "N/A";
      }
      if (std::getline(ss, field, ',')) {
        std::string tfd_v1_str = EOStringUtils::trim_string(field);
        std::transform(tfd_v1_str.begin(), tfd_v1_str.end(), tfd_v1_str.begin(),
                       ::tolower);
        if (tfd_v1_str == "false" || tfd_v1_str == "0")
          tc.use_tfd_v1_for_ref = false;
        else if (tfd_v1_str == "true" || tfd_v1_str == "1" ||
                 tfd_v1_str.empty())
          tc.use_tfd_v1_for_ref = true;
        else
          INFO("Warning: Unrecognized value for use_tfd_v1 '"
               << field << "', defaulting to true.");
      }

      cases.push_back(tc);
    } catch (const std::exception& e) {
      FAIL_CHECK("Error parsing reference data file "
                 << filepath << " at line " << current_line_num << ": "
                 << e.what() << " on content: \"" << processed_line << "\"");
    }
  }
  return cases;
}

TEST_CASE("Numerical Validation from Reference Data") {
  INFO("Starting Numerical Validation Test Case");

  // Path to reference data. Adjust relative to CTest execution directory.
  std::string ref_data_filepath =
      "../../../tests/test_data/reference_eos_data.csv";
  std::vector<NumericalTestCase> test_cases =
      load_numerical_test_cases(ref_data_filepath);
  REQUIRE_FALSE(test_cases.empty());  // Fail if no test cases were loaded

  // Group test cases by their required TFD version setting
  std::map<bool, std::vector<NumericalTestCase>> grouped_cases;
  for (const auto& tc : test_cases) {
    grouped_cases[tc.use_tfd_v1_for_ref].push_back(tc);
  }

  // EOS Manager setup (assuming eos_data_dir is needed for initialize)
  // Path to the actual eos_data_dir used by EquationOfStateV1
  std::string eos_data_path = "../../../eos_data_dir";

  EquationOfStateV1 eos_manager;

  // We need to initialize the eos_manager with all unique eos_ids found in
  // test_cases.
  std::vector<int> all_eos_ids_to_init;
  std::map<int, bool> unique_ids;  // Using map for uniqueness
  for (const auto& tc : test_cases) {
    unique_ids[tc.eos_id] = true;
  }
  all_eos_ids_to_init.reserve(unique_ids.size());
  for (const auto& pair : unique_ids) {
    all_eos_ids_to_init.push_back(pair.first);
  }

  INFO(
      "Initializing EOS manager for numerical tests with EOS IDs: (listing "
      "IDs)");
  for (int id : all_eos_ids_to_init) INFO("  ID: " << id);  // Can be verbose

  int init_stat = eos_manager.initialize(all_eos_ids_to_init, eos_data_path);
  REQUIRE_MESSAGE(init_stat == EOS_SUCCESS,
                  "EOS Manager initialization failed with code: " << init_stat);

  // Define a tolerance for floating point comparisons
  constexpr double rel_tolerance = 1e-7;  // Relative tolerance
  constexpr double abs_tolerance =
      1e-9;  // Absolute tolerance for values near zero

  for (auto const& group : grouped_cases) {
    bool use_v1 = group.first;
    auto case_group = group.second;

    MESSAGE("Testing with TFD Configuration: use_tfd_v1 = "
            << (use_v1 ? "true" : "false"));
    if (case_group.empty()) continue;

    EquationOfStateV1 eos_lib;
    eos_lib.setUseTFDDataVer1(use_v1);  // Set TFD version for this group

    std::vector<int> ids_for_this_group;
    std::map<int, bool> unique_ids_in_group;
    for (const auto& tc : case_group) {
      unique_ids_in_group[tc.eos_id] = true;
    }
    ids_for_this_group.reserve(unique_ids_in_group.size());
    for (const auto& pair : unique_ids_in_group) {
      ids_for_this_group.push_back(pair.first);
    }

    INFO("Initializing EOS manager for TFD config use_v1="
         << use_v1 << " with EOS IDs (count: " << ids_for_this_group.size()
         << ")");

    init_stat = eos_lib.initialize(ids_for_this_group, eos_data_path);
    REQUIRE_MESSAGE(init_stat == EOS_SUCCESS,
                    "EOS Manager initialization failed for TFD config use_v1="
                        << use_v1 << " with code: " << init_stat);

    for (const auto& tc : case_group) {
      // Create a unique subcase name
      std::string subcase_name = tc.description;
      if (subcase_name.empty() || subcase_name == "N/A") {
        subcase_name = "EOSID_" + std::to_string(tc.eos_id) + "_rho_" +
                       std::to_string(tc.rho) +  // Potentially long names
                       "_T_" + std::to_string(tc.T);
      }
      subcase_name += (use_v1 ? "_TFDv1" : "_TFDv2");

      SUBCASE(subcase_name.c_str()) {
        INFO("Test Case: " << (tc.description.empty()
                                   ? "Line " + std::to_string(tc.line_num)
                                   : tc.description));
        INFO("  Input: eos_id=" << tc.eos_id << ", rho=" << tc.rho
                                << ", T=" << tc.T);
        INFO("  Reference: P=" << tc.ref_P
                               << ", E=" << tc.ref_E /* << other refs */);

        ComputeResult result = eos_lib.compute(tc.eos_id, tc.rho, tc.T);
        REQUIRE(result.istat == EOS_SUCCESS);

        // Compare calculated values with reference values
        // Using doctest::Approx with relative and absolute epsilon might be
        // better.
        auto check_approx_equal = [&](double val_calc, double val_ref,
                                      const std::string& name) {
          if (std::abs(val_ref) > abs_tolerance * 1000) {
            // If ref value is not too small, relative makes sense
            REQUIRE_MESSAGE(
                std::abs(val_calc - val_ref) <=
                    rel_tolerance * std::abs(val_ref),
                name << ": Calc=" << val_calc << ", Ref=" << val_ref
                     << ", RelDiff="
                     << (std::abs(val_ref) > 0
                             ? std::abs(val_calc - val_ref) / std::abs(val_ref)
                             : 0.0));
          } else {
            // For small ref values, absolute tolerance is more robust
            REQUIRE_MESSAGE(
                std::abs(val_calc - val_ref) <= abs_tolerance,
                name << ": Calc=" << val_calc << ", Ref=" << val_ref
                     << ", AbsDiff=" << std::abs(val_calc - val_ref));
          }
        };

        check_approx_equal(result.P, tc.ref_P, "P");
        check_approx_equal(result.E, tc.ref_E, "E");
        check_approx_equal(result.dPdT, tc.ref_dPdT, "dPdT");
        check_approx_equal(result.dEdT, tc.ref_dEdT, "dEdT");
        check_approx_equal(result.dPdrho, tc.ref_dPdrho, "dPdrho");
      }
    }
  }
}
