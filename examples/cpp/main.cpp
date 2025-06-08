// examples/cpp_example/main.cpp
#if defined(BUILT_WITHIN_MAIN_PROJECT)
  #include "EquationOfStateV1.h" // Found via target_include_directories above
  #include "utils/string_utils.h"
#else
  #include <EOSCore/EquationOfStateV1.h> // For installed package
  #include <EOSCore/string_utils.h>
#endif
#include <string_utils.h>    // If string_utils.h is also public/installed
#include <iostream>
#include <vector>
#include <string>
#include <iomanip> // For output formatting

// A simple helper to print results, can be more elaborate
void print_results(const std::string& title, int istat, double P, double E) {
    std::cout << title << " (istat: " << istat << ")" << std::endl;
    if (istat == 0) { // Assuming 0 is EOS_SUCCESS
        std::cout << "  P = " << std::fixed << std::setprecision(5) << P
                  << ", E = " << E << std::endl;
    } else {
        std::cout << "  Computation failed." << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "C++ EOS Library Example" << std::endl;

    std::string data_dir_path;
    if (argc > 1) {
        data_dir_path = argv[1];
    } else {
        // Default path relative to where example might be run from build tree
        // For CI, we'll pass this path carefully.
        data_dir_path = "../../../eos_data_dir"; // Adjust if needed
        std::cout << "Usage: CppEOSDemo <path_to_eos_data_dir>" << std::endl;
        std::cout << "No path provided, using default: " << data_dir_path << std::endl;
    }


    EquationOfStateV1 eos_manager; // Using namespace if header uses it

    std::vector<int> eos_ids_to_load = {1, 10000}; // Load analytic air and one complicated

    int istat = eos_manager.initialize(eos_ids_to_load, data_dir_path);
    if (istat != 0) { // Assuming 0 is EOS_SUCCESS
        std::cerr << "Failed to initialize EOS manager. Error: " << istat << std::endl;
        return 1;
    }
    std::cout << "EOS manager initialized." << std::endl;

    double rho, T, P, E, dPdT, dEdT, dPdrho;

    // Test Analytic Air (ID 1)
    rho = 1.2; T = 300.0;
    istat = eos_manager.compute(1, rho, T, P, E, dPdT, dEdT, dPdrho);
    print_results("Analytic Air (ID 1)", istat, P, E);

    // Test Complicated (ID 10000)
    rho = 1.8; T = 1200.0;
    istat = eos_manager.compute(10000, rho, T, P, E, dPdT, dEdT, dPdrho);
    print_results("Complicated EOS (ID 10000)", istat, P, E);

    eos_manager.free_resources();
    std::cout << "Example finished." << std::endl;
    return 0;
}
