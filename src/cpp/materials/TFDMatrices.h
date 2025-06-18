#ifndef TFD_MATRICES_H
#define TFD_MATRICES_H

#include <vector>
#include <string> // For storing descriptions later if needed

namespace EOS_Internal { // Or your chosen internal namespace

    struct TFDMatrices {
        std::vector<double> matrix_A;
        std::vector<double> matrix_B;
        int N1 = 0; // Common dimension 1
        int N2 = 0; // Common dimension 2
        bool initialized = false;

        // Optional: Metadata loaded from HDF5
        std::string tfd_file_description;
        std::string matrix_a_description;
        std::string matrix_b_description;
        // std::string loaded_tfd_filepath; // To know which file was source

        void clear() {
            matrix_A.clear();
            matrix_B.clear();
            N1 = 0;
            N2 = 0;
            initialized = false;
            tfd_file_description.clear();
            matrix_a_description.clear();
            matrix_b_description.clear();
        }
    };

} // namespace EOS_Internal
#endif // TFD_MATRICES_H
