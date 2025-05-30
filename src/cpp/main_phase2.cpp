// src/cpp/main_phase2.cpp
#include "EquationOfStateV1.h"
#include <iostream>
#include <iomanip>
#include <string>   // For std::string paths
#include <vector>   // For file paths

int main() {
    std::cout << "C++: Phase 2 Test Program" << std::endl;
    std::cout << "-------------------------" << std::endl;

    EquationOfStateV1 eos_instance;

    // Test setting control variables (from Phase 1)
    std::cout << "\nSetting control variables:" << std::endl;
    eos_instance.setComplicatedEOSUseOldColdTerm(true);
    eos_instance.setUseTFDDataVer1(false);

    // Define paths to TFD data files
    // IMPORTANT: Adjust these paths relative to where your executable will run
    // or copy eos_data_dir to your build/executable directory.
    // For CMake, if eos_data_dir is at project root:
    // std::string base_path = "../../eos_data_dir/"; // If running from build/src/cpp/
    // Or more robustly, pass path via command line or configure it.
    // For simplicity in this test, assume files are accessible.
    std::string tfd_file_A = "../eos_data_dir/tfd_matA_v1.dat"; // Expect this in current working dir of exe
    std::string tfd_file_B = "../eos_data_dir/tfd_matB_v1.dat"; // Or provide full/relative paths

    // Path considerations:
    // If eos_data_dir is at project root:
    // One way is to copy data files to build directory using CMake's add_custom_command or file(COPY...)
    // For now, manually ensure the files are where the executable can find them,
    // e.g., by running the executable from the directory containing them, or using absolute paths.
    // Let's assume you'll run the executable from a directory where these files are directly.
    // Or for a quick test, place them in your build directory alongside the executable.

    // You might need to adjust these paths based on your build setup and execution CWD.
    // If your project root has eos_data_dir:
    // std::string project_root_data_path = "../../eos_data_dir/"; // if exe is in build/src/cpp
    // std::string tfd_file_A = project_root_data_path + "tfd_matA_v1.dat";
    // std::string tfd_file_B = project_root_data_path + "tfd_matB_v1.dat";
    // For this example, let's assume files are in the same directory as the executable
    // (e.g. copy them to your build directory)


    std::cout << "\nLoading TFD Data:" << std::endl;
    int load_istat = eos_instance.loadTFDData(tfd_file_A, tfd_file_B);

    if (load_istat != 0) {
        std::cerr << "  Failed to load TFD data. Error code: " << load_istat << std::endl;
        std::cerr << "  Ensure '" << tfd_file_A << "' and '" << tfd_file_B
                  << "' exist and are readable relative to the executable path." << std::endl;
    } else {
        std::cout << "  TFD data loaded successfully." << std::endl;

        // Test TFD computation
        std::cout << "\nTesting TFD computation:" << std::endl;
        double rho = 1.5, T = 500.0;
        double res_x, res_y;
        int tfd_istat = eos_instance.computeTFD(rho, T, res_x, res_y);

        std::cout << std::fixed << std::setprecision(5);
        if (tfd_istat == 0) {
            std::cout << "  TFD compute success: rho=" << rho << ", T=" << T
                      << " -> result_X=" << res_x << ", result_Y=" << res_y << std::endl;
        } else {
            std::cout << "  TFD compute failed. istat=" << tfd_istat << std::endl;
        }
    }


    std::cout << "\n-------------------------" << std::endl;
    std::cout << "C++: Phase 2 Test Complete." << std::endl;
    return 0;
}
