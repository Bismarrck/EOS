// src/cpp/EquationOfStateV1.cpp
#include "EquationOfStateV1.h"
#include <iostream>     // For debug prints
#include <fstream>      // for std::ifstream
#include <sstream>      // For std::istringstream
#include <iomanip>      // For std::fixed, std::setprecision in main, can be useful here too
#include <algorithm> // For std::transform for string tolower


// Helper for basic path joining (add to EquationOfStateV1.cpp or a common utils file)
// WARNING: This is a very basic path joiner, doesn't handle all edge cases
// like multiple separators or ensuring correct separator for OS.
// For robust cross-platform paths without C++17, a small utility or careful manual
// construction is needed.
static std::string simple_path_join(const std::string& p1, const std::string& p2) {
    char sep = '/'; // Common separator
    std::string tmp = p1;

#ifdef _WIN32 // Or other Windows checks
    sep = '\\';
#endif

    if (!p1.empty() && p1.back() != sep && !p2.empty() && p2.front() != sep) {
        tmp += sep;
    } else if (!p1.empty() && p1.back() == sep && !p2.empty() && p2.front() == sep) {
        // Remove one separator if both have it
        tmp.pop_back();
    }
    tmp += p2;
    return tmp;
}

// Update get_eos_file_path in EquationOfStateV1.cpp
static std::string get_eos_file_path(const std::string& base_dir, int eos_id_full) {
    std::ostringstream eos_id_ss;
    eos_id_ss << std::setw(5) << std::setfill('0') << eos_id_full;
    std::string eos_id_str = eos_id_ss.str();

    if (eos_id_str.length() < 3) return "";

    std::string material_id_str = eos_id_str.substr(0, 3);

    std::string path = base_dir;
    path = simple_path_join(path, "mat" + material_id_str);
    path = simple_path_join(path, "eos" + eos_id_str + ".dat");
    return path;
}

// --- Shim Implementations ---
int EquationOfStateV1::call_air_eos_2000_shim(double rho, double T, double& P, double& E) {
    int istat;
    // Outputs for dPdT etc. are not provided by this analytic EOS, so pass dummies or handle if API changes
    c_air_eos_2000(rho, T, &P, &E, &istat);
    return istat;
}

int EquationOfStateV1::call_carbon_eos_2001_shim(double rho, double T, double& P, double& E) {
    int istat;
    c_carbon_eos_2001(rho, T, &P, &E, &istat);
    return istat;
}

EquationOfStateV1::EquationOfStateV1() : tfd_use_ver1_setting_(true) { // Default TFD version
    std::cout << "C++: EquationOfStateV1 instance created." << std::endl;

    // Populate the analytic EOS registry
    // The key is the `eos_id` that would be found in the filename (e.g., from materials_list)
    // For this example, let's say analytic EOS "air" has ID 1, "carbon" has ID 2
    // These IDs need to be consistent with how you identify them in your eos_data_dir/#TYPE lines
    // or in the `eos_id_list` passed to initialize.
    // If analytic EOS files also follow the matXXX/eosXXXXX.dat pattern, then those full XXXXX IDs are used.
    // For simplicity, let's assume special IDs for these hardcoded analytic types for now.
    // If eos10000.dat said #TYPE: ANALYTIC and was for air_eos_2000:
    // analytic_eos_registry_[10000] = {InternalEOSType::ANALYTIC_AIR, &EquationOfStateV1::call_air_eos_2000_shim};
    // Or, if you have designated IDs for them:
    analytic_eos_registry_[1] = {InternalEOSType::ANALYTIC_AIR, &EquationOfStateV1::call_air_eos_2000_shim}; // User must provide '1' in eos_id_list
    analytic_eos_registry_[2] = {InternalEOSType::ANALYTIC_CARBON, &EquationOfStateV1::call_carbon_eos_2001_shim}; // User must provide '2'

    // Initialize TFD data struct
    tfd_data.initialized = false;
}

EquationOfStateV1::~EquationOfStateV1() {
    std::cout << "C++: EquationOfStateV1 instance destroyed." << std::endl;
    free_resources();
}

// --- Control Variable Management ---
int EquationOfStateV1::setComplicatedEOSUseOldColdTerm(bool value) {
    int istat;
    c_set_complicated_eos_use_old_cold_term(value, &istat);
    if (istat != 0) {
        std::cerr << "C++ Error: Failed to set complicated_eos_use_old_cold_term. Fortran istat: " << istat << std::endl;
    }
    return istat;
}

int EquationOfStateV1::setUseTFDDataVer1(bool value) {
    int istat;
    c_set_use_tfd_data_ver1(value, &istat);
    if (istat != 0) {
        std::cerr << "C++ Error: Failed to set use_tfd_data_ver1. Fortran istat: " << istat << std::endl;
    }
    return istat;
}

// New helper function in EquationOfStateV1.cpp
int EquationOfStateV1::read_matrix_values_from_stream(std::ifstream& infile, int N1, int N2,
                                                    std::vector<double>& matrix_data, const std::string& matrix_name_for_error) {
    matrix_data.clear();
    matrix_data.reserve(static_cast<size_t>(N1) * N2);
    double value;
    std::string line;

    size_t values_to_read = static_cast<size_t>(N1) * N2;
    size_t values_read = 0;

    // Read values, handling multiple values per line
    while (values_read < values_to_read && std::getline(infile, line)) {
        std::istringstream val_ss(line);
        while (val_ss >> value) {
            if (values_read < values_to_read) {
                matrix_data.push_back(value);
                values_read++;
            } else {
                std::cerr << "Warning: Extra data found after reading " << matrix_name_for_error << std::endl;
                // break from inner and outer loop or just return success if this is acceptable
                goto end_read_loop; // Ugly, but avoids nested break flags
            }
        }
        if (val_ss.fail() && !val_ss.eof()) { // Check for non-double data on a line
             std::cerr << "Error: Non-double value encountered while reading " << matrix_name_for_error << std::endl;
             return EOS_ERROR_FILE_PARSE;
        }
    }
end_read_loop:;

    if (values_read != values_to_read) {
        std::cerr << "Error: Read " << values_read << " values for " << matrix_name_for_error
                  << ", but expected " << values_to_read
                  << " for dimensions " << N1 << "x" << N2 << std::endl;
        return EOS_ERROR_FILE_PARSE; // Data size mismatch
    }
    return EOS_SUCCESS;
}


int EquationOfStateV1::loadTFDDataInternal(const std::string& tfd_base_dir) {
    std::string tfd_file_name = (tfd_use_ver1_setting_ ? "tfd_ver1.dat" : "tfd_ver2.dat");
    std::string full_tfd_path = simple_path_join(tfd_base_dir, tfd_file_name);

    std::ifstream infile(full_tfd_path);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open TFD file: " << full_tfd_path << std::endl;
        tfd_data.initialized = false;
        return EOS_ERROR_FILE_NOT_FOUND;
    }

    std::string line;
    // Read common dimensions N1, N2
    if (!std::getline(infile, line)) {
        std::cerr << "Error: Could not read dimensions line from " << full_tfd_path << std::endl;
        tfd_data.initialized = false;
        return EOS_ERROR_FILE_PARSE;
    }
    std::istringstream dim_ss(line);
    if (!(dim_ss >> tfd_data.N1 >> tfd_data.N2)) {
        std::cerr << "Error: Could not parse dimensions N1, N2 from " << full_tfd_path << std::endl;
        tfd_data.initialized = false;
        return EOS_ERROR_FILE_PARSE;
    }

    if (tfd_data.N1 <= 0 || tfd_data.N2 <= 0) {
        std::cerr << "Error: Invalid dimensions N1=" << tfd_data.N1 << ", N2=" << tfd_data.N2
                  << " in " << full_tfd_path << std::endl;
        tfd_data.initialized = false;
        return EOS_ERROR_INVALID_DIMENSIONS;
    }
    std::cout << "TFD dimensions from " << full_tfd_path << ": " << tfd_data.N1 << "x" << tfd_data.N2 << std::endl;

    // Read Matrix A
    int istat_A = read_matrix_values_from_stream(infile, tfd_data.N1, tfd_data.N2, tfd_data.matrix_A, "Matrix A");
    if (istat_A != EOS_SUCCESS) {
        tfd_data.initialized = false;
        return EOS_ERROR_TFD_LOAD_A + istat_A;
    }
    std::cout << "Successfully read Matrix A data." << std::endl;

    // Read Matrix B
    int istat_B = read_matrix_values_from_stream(infile, tfd_data.N1, tfd_data.N2, tfd_data.matrix_B, "Matrix B");
    if (istat_B != EOS_SUCCESS) {
        tfd_data.initialized = false;
        return EOS_ERROR_TFD_LOAD_B + istat_B;
    }
    std::cout << "Successfully read Matrix B data." << std::endl;

    // If there were more matrices, continue reading them here.

    tfd_data.initialized = true;
    std::cout << "C++: TFD data (" << (tfd_use_ver1_setting_ ? "ver1" : "ver2") << ") loaded successfully from " << full_tfd_path << std::endl;
    return EOS_SUCCESS;
}

// --- Initialization ---
int EquationOfStateV1::initialize(const std::vector<int>& eos_id_list, const std::string& eos_data_dir) {
    // 0. Ensure Fortran side knows which TFD version we prefer.
    //    The C++ tfd_use_ver1_setting_ is the master; set it if user calls the setter.
    //    Then ensure Fortran is synced.
    int f_istat;
    c_set_use_tfd_data_ver1(tfd_use_ver1_setting_, &f_istat);
    if (f_istat != 0) {
        std::cerr << "Error setting Fortran TFD version preference." << std::endl;
        // Decide if this is a fatal error for initialize
    }

    // 1. Load TFD data
    int tfd_load_status = loadTFDDataInternal(eos_data_dir); // Pass base data dir for TFD files
    if (tfd_load_status != EOS_SUCCESS) {
        std::cerr << "Error initializing: Failed to load TFD data. Code: " << tfd_load_status << std::endl;
        return tfd_load_status;
    }

    // 2. Load data for each EOS ID
    for (int id : eos_id_list) {
        std::string file_path_str = get_eos_file_path(eos_data_dir, id);
        if (file_path_str.empty()) {
            std::cerr << "Error: Could not construct path for EOS ID " << id << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND; // Or collect errors and report all
        }

        std::ifstream eos_file(file_path_str);
        if (!eos_file.is_open()) {
            // Special check: is this an ID for a known, hardcoded analytic EOS?
            auto it_analytic_reg = analytic_eos_registry_.find(id);
            if (it_analytic_reg != analytic_eos_registry_.end()) {
                MaterialData md;
                md.eos_id_key = id;
                md.type_from_file = EOSTypeFromFile::ANALYTIC; // By virtue of being in registry
                md.internal_type = it_analytic_reg->second.internal_type;
                md.analytic_func = it_analytic_reg->second.func_ptr;
                loaded_materials_[id] = md;
                std::cout << "Initialized hardcoded Analytic EOS ID " << id << " (" << file_path_str << " not strictly needed but implies registration)." << std::endl;
                continue; // Successfully processed this analytic ID
            }
            std::cerr << "Error: Could not open EOS data file: " << file_path_str << " for ID " << id << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND;
        }

        std::string line;
        std::string type_str_from_file;
        if (std::getline(eos_file, line)) {
            std::string prefix = "#TYPE:";
            size_t prefix_pos = line.find(prefix);
            if (prefix_pos != std::string::npos) {
                type_str_from_file = line.substr(prefix_pos + prefix.length());
                // Trim whitespace
                type_str_from_file.erase(0, type_str_from_file.find_first_not_of(" \t\n\r\f\v"));
                type_str_from_file.erase(type_str_from_file.find_last_not_of(" \t\n\r\f\v") + 1);
                std::transform(type_str_from_file.begin(), type_str_from_file.end(), type_str_from_file.begin(), ::tolower);
            } else {
                std::cerr << "Error: Missing or malformed #TYPE: line in " << file_path_str << std::endl;
                return EOS_ERROR_FILE_PARSE;
            }
        } else {
            std::cerr << "Error: Could not read #TYPE: line from " << file_path_str << std::endl;
            return EOS_ERROR_FILE_PARSE;
        }

        MaterialData md;
        md.eos_id_key = id;

        if (type_str_from_file == "analytic") {
            md.type_from_file = EOSTypeFromFile::ANALYTIC;
            auto it_analytic_reg = analytic_eos_registry_.find(id);
            if (it_analytic_reg != analytic_eos_registry_.end()) {
                md.internal_type = it_analytic_reg->second.internal_type;
                md.analytic_func = it_analytic_reg->second.func_ptr;
            } else {
                std::cerr << "Error: Analytic EOS ID " << id << " from file " << file_path_str
                          << " is not registered in C++ analytic_eos_registry_." << std::endl;
                return EOS_ERROR_UNKNOWN_EOS_ID;
            }
            // For analytic, the .dat file might contain no other parameters.
        } else if (type_str_from_file == "complicated") {
            md.type_from_file = EOSTypeFromFile::COMPLICATED;
            md.internal_type = InternalEOSType::COMPLICATED_V1; // Assuming one type for now

            // Read parameters for complicated EOS
            double param_val;
            while(eos_file >> param_val) { // Simple sequential read of doubles
                md.params.push_back(param_val);
            }
            if (eos_file.eof()){ /* good */ }
            else if (eos_file.fail() && !eos_file.eof()){ // check if fail was not due to eof
                 std::cerr << "Error: Non-double value encountered while reading params from " << file_path_str << std::endl;
                 return EOS_ERROR_FILE_PARSE;
            }

            if (md.params.empty()) {
                std::cerr << "Warning: No parameters found for complicated EOS ID " << id << " in " << file_path_str << std::endl;
                // Potentially return an error if params are mandatory
            }
             std::cout << "Read " << md.params.size() << " params for complicated EOS ID " << id << std::endl;

        } else {
            std::cerr << "Error: Unknown #TYPE: '" << type_str_from_file << "' in " << file_path_str << std::endl;
            return EOS_ERROR_INVALID_EOS_TYPE;
        }
        loaded_materials_[id] = md;
        std::cout << "Successfully loaded EOS ID " << id << " from " << file_path_str << std::endl;
    }
    return EOS_SUCCESS;
}

int EquationOfStateV1::check_eos_data_dir(const std::string& eos_data_dir, const std::vector<int>& eos_ids_to_check) {
    // Check TFD files
    std::string tfd_files_to_check[2] = { "tfd_ver1.dat", "tfd_ver2.dat" };

    for (const std::string& tfd_suffix : tfd_files_to_check) {
        std::string tfd_path_str = simple_path_join(eos_data_dir, tfd_suffix);
        std::ifstream test_file(tfd_path_str);
        if (!test_file.is_open()) {
            std::cerr << "Check Error: TFD file not found or not accessible: " << tfd_path_str << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND;
        }
        test_file.close();
    }

    // Check specific EOS ID files
    for (int id : eos_ids_to_check) {
        if (analytic_eos_registry_.count(id)) {
            std::cout << "Check Info: Analytic EOS ID " << id << " is registered, file existence check skipped/optional." << std::endl;
            continue;
        }
        std::string file_path_str = get_eos_file_path(eos_data_dir, id);
        if (file_path_str.empty()){
            std::cerr << "Check Error: Could not construct path for EOS ID " << id << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND; // Or a different error code
        }
        std::ifstream eos_file(file_path_str);
        if (!eos_file.is_open()) {
            std::cerr << "Check Error: EOS data file not found for ID " << id << " at " << file_path_str << std::endl;
            return EOS_ERROR_FILE_NOT_FOUND;
        }
        eos_file.close();
    }
    std::cout << "check_eos_data_dir: All checked files appear accessible." << std::endl;
    return EOS_SUCCESS;
}

// --- Computation ---
int EquationOfStateV1::compute(int eos_id, double rho, double T,
                             double& P, double& E,
                             double& dPdT, double& dEdT, double& dPdrho) {
    auto it = loaded_materials_.find(eos_id);
    if (it == loaded_materials_.end()) {
        std::cerr << "Error: EOS ID " << eos_id << " not initialized." << std::endl;
        return EOS_ERROR_UNKNOWN_EOS_ID;
    }

    const MaterialData& md = it->second;
    P = 0.0; E = 0.0; dPdT = 0.0; dEdT = 0.0; dPdrho = 0.0; // Initialize outputs

    int istat = EOS_SUCCESS;

    switch (md.internal_type) {
        case InternalEOSType::ANALYTIC_AIR:
        case InternalEOSType::ANALYTIC_CARBON:
            if (md.analytic_func) {
                // Current analytic shims only set P, E. Derivatives would be zero.
                // If analytic functions also give derivatives, the AnalyticEOSFunc signature
                // and shims need to be updated.
                istat = md.analytic_func(rho, T, P, E);
                // dPdT, dEdT, dPdrho remain 0 for these specific analytic stubs
            } else {
                std::cerr << "Error: No analytic function registered for EOS ID " << eos_id << std::endl;
                istat = EOS_ERROR_ANALYTIC_DISPATCH;
            }
            break;

        case InternalEOSType::COMPLICATED_V1:
            if (!tfd_data.initialized) {
                std::cerr << "Error: TFD data not initialized for complicated EOS ID " << eos_id << std::endl;
                return EOS_ERROR_TFD_NOT_INIT;
            }
            if (md.params.empty() && md.type_from_file == EOSTypeFromFile::COMPLICATED) {
                 // This check depends on whether complicated EOS can ever have zero params
                std::cerr << "Warning: Complicated EOS ID " << eos_id << " has no parameters loaded." << std::endl;
                // return EOS_ERROR_PARAMS_NOT_LOADED; // Or proceed if allowed
            }
            // Call the C wrapper from Phase 3 directly for now
            c_complicated_eos(md.params.data(), static_cast<int>(md.params.size()), rho, T,
                              tfd_data.matrix_A.data(), tfd_data.N1, tfd_data.N2,
                              tfd_data.matrix_B.data(), tfd_data.N1, tfd_data.N2,
                              &P, &E, &dPdT, &dEdT, &dPdrho,
                              &istat);
            break;

        default:
            std::cerr << "Error: Unknown internal EOS type for ID " << eos_id << std::endl;
            istat = EOS_ERROR_UNKNOWN_EOS_ID;
            break;
    }

    if (istat != 0) {
        // Fortran istat might be directly returned, or mapped to C++ error codes
        std::cerr << "C++ compute: Fortran calculation for EOS ID " << eos_id << " failed with istat: " << istat << std::endl;
    }
    return istat; // Return Fortran's istat or a mapped C++ error code
}


void EquationOfStateV1::free_resources() {
    loaded_materials_.clear();
    tfd_data.matrix_A.clear();
    tfd_data.matrix_B.clear();
    tfd_data.initialized = false;
    std::cout << "C++: EquationOfStateV1 resources freed." << std::endl;
}
