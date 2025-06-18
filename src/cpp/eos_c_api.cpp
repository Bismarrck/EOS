// src/cpp/eos_c_api.cpp
#include "EquationOfStateV1.h" // Your main C++ class header
#include "eos_error_codes.h"
#include <string>
#include <vector>
#include <memory>   // For std::unique_ptr
#include <iostream> // For C-API error messages/debug
#include <algorithm> // For std::transform on control_name

// --- Global Instance Management ---
// This unique_ptr will manage the lifetime of the single global EOS instance.
static std::unique_ptr<EquationOfStateV1> global_eos_instance = nullptr;

// Helper to get the instance or report an error if not initialized
// Returns nullptr on error, prints to cerr.
static EquationOfStateV1* get_eos_instance_for_api(const char* func_name) {
    if (!global_eos_instance) {
        std::cerr << "C-API Error (" << (func_name ? func_name : "unknown func")
                  << "): EOS instance not initialized. Call eos_c_initialize_instance first." << std::endl;
        return nullptr;
    }
    return global_eos_instance.get();
}


extern "C" {

// --- Instance Management ---

int eos_c_initialize_instance(const char* data_dir_c_str, int data_dir_len,
                              const int* eos_ids_arr, int num_eos_ids) {
    if (!data_dir_c_str || data_dir_len <= 0) {
        std::cerr << "C-API Error (eos_c_initialize_instance): Invalid data_dir_c_str or data_dir_len." << std::endl;
        return -1; // Invalid argument
    }
    // num_eos_ids can be 0 for an empty list, eos_ids_arr can be null then.
    if (num_eos_ids < 0 || (num_eos_ids > 0 && !eos_ids_arr)) {
        std::cerr << "C-API Error (eos_c_initialize_instance): Invalid eos_ids_arr or num_eos_ids." << std::endl;
        return -2; // Invalid argument
    }

    try {
        std::string data_dir(data_dir_c_str, static_cast<size_t>(data_dir_len));
        std::vector<int> eos_ids_vec;
        if (num_eos_ids > 0) {
            eos_ids_vec.assign(eos_ids_arr, eos_ids_arr + num_eos_ids);
        }

        // Create or re-create the instance
        global_eos_instance = std::make_unique<EquationOfStateV1>();
        if (!global_eos_instance) { // Should not happen with make_unique unless out of memory
             std::cerr << "C-API Error (eos_c_initialize_instance): Failed to allocate EquationOfStateV1 instance." << std::endl;
            return -3; // Allocation failure
        }
        std::cout << "C-API (eos_c_initialize_instance): Instance created. Calling C++ initialize..." << std::endl;
        return global_eos_instance->initialize(eos_ids_vec, data_dir);
    } catch (const std::bad_alloc& ba) {
        std::cerr << "C-API Exception (eos_c_initialize_instance): Memory allocation failed: " << ba.what() << std::endl;
        global_eos_instance.reset(); // Ensure no dangling partial state
        return -98; // Memory allocation exception
    } catch (const std::exception& e) {
        std::cerr << "C-API Exception (eos_c_initialize_instance): " << e.what() << std::endl;
        global_eos_instance.reset();
        return -99; // Generic C++ exception
    } catch (...) {
        std::cerr << "C-API Unknown exception in eos_c_initialize_instance." << std::endl;
        global_eos_instance.reset();
        return -100; // Unknown exception
    }
}

void eos_c_free_instance() {
    if (global_eos_instance) {
        std::cout << "C-API (eos_c_free_instance): Freeing EOS instance." << std::endl;
        // global_eos_instance->free_resources(); // Call explicit C++ free if it does more than destructor
        global_eos_instance.reset(); // Destroys the managed object
    } else {
        std::cout << "C-API (eos_c_free_instance): Instance was already null." << std::endl;
    }
}

// --- Control Variables ---

int eos_c_set_bool_control(const char* control_name_c_str, int name_len, bool value) {
    EquationOfStateV1* instance = get_eos_instance_for_api("eos_c_set_bool_control");
    if (!instance) return EOS_ERROR_UNKNOWN_EOS_ID -10; // Re-using an error code, or define C_API specific ones

    if (!control_name_c_str || name_len <= 0) {
        std::cerr << "C-API Error (eos_c_set_bool_control): Invalid control_name_c_str or name_len." << std::endl;
        return -1;
    }
    std::string control_name(control_name_c_str, static_cast<size_t>(name_len));
    // Optional: convert control_name to lowercase for case-insensitive matching
    std::transform(control_name.begin(), control_name.end(), control_name.begin(), ::tolower);


    std::cout << "C-API (eos_c_set_bool_control): Setting '" << control_name << "' to " << (value ? "true" : "false") << std::endl;
    if (control_name == "use_tfd_data_ver1") {
        return instance->setUseTFDDataVer1(value);
    } else if (control_name == "complicated_eos_use_old_cold_term") {
        return instance->setComplicatedEOSUseOldColdTerm(value);
    }
    // Add more boolean controls here if needed
    std::cerr << "C-API Error (eos_c_set_bool_control): Unknown boolean control name: " << control_name << std::endl;
    return -2; // Unknown control name
}

// --- Computation ---

int eos_c_compute(int eos_id, double rho, double T,
                  double* P_out, double* E_out,
                  double* dPdT_out, double* dEdT_out, double* dPdrho_out) {
    EquationOfStateV1* instance = get_eos_instance_for_api("eos_c_compute");
    if (!instance) return EOS_ERROR_UNKNOWN_EOS_ID -10;

    if (!P_out || !E_out || !dPdT_out || !dEdT_out || !dPdrho_out) {
        std::cerr << "C-API Error (eos_c_compute): One or more output pointers are null." << std::endl;
        return -1; // Invalid output pointers
    }
    // std::cout << "C-API (eos_c_compute): Calling C++ compute for ID " << eos_id << std::endl; // Can be verbose
    return instance->compute(eos_id, rho, T, *P_out, *E_out, *dPdT_out, *dEdT_out, *dPdrho_out);
}

// --- Pack/Unpack ---

int eos_c_pack_data_get_size(int* required_buffer_size_out) {
    EquationOfStateV1* instance = get_eos_instance_for_api("eos_c_pack_data_get_size");
    if (!instance) return EOS_ERROR_UNKNOWN_EOS_ID -10;

    if (!required_buffer_size_out) {
         std::cerr << "C-API Error (eos_c_pack_data_get_size): Null output pointer for size." << std::endl;
        return -1;
    }

    char* temp_buf_ptr = nullptr; // Pass nullptr to C++ pack_data to get size
    int size = 0;
    int stat = instance->pack_data(temp_buf_ptr, size); // C++ pack_data updates 'size'
    if (stat == EOS_SUCCESS) {
        *required_buffer_size_out = size;
    } else {
        *required_buffer_size_out = 0; // Ensure it's zero on error
    }
    // std::cout << "C-API (eos_c_pack_data_get_size): C++ pack_data returned size " << size << " with status " << stat << std::endl;
    return stat;
}

int eos_c_pack_data_to_buffer(char* buffer_inout, int* buffer_size_inout) {
    EquationOfStateV1* instance = get_eos_instance_for_api("eos_c_pack_data_to_buffer");
    if (!instance) return EOS_ERROR_UNKNOWN_EOS_ID -10;

    if (!buffer_inout || !buffer_size_inout || *buffer_size_inout <= 0) {
        std::cerr << "C-API Error (eos_c_pack_data_to_buffer): Invalid buffer or size_inout pointer, or non-positive capacity." << std::endl;
        return -1;
    }
    // std::cout << "C-API (eos_c_pack_data_to_buffer): Calling C++ pack_data with buffer capacity " << *buffer_size_inout << std::endl;
    // The C++ pack_data(char*& buffer, int& buffer_size) expects:
    // - `buffer` (1st arg): if non-null, it's the pre-allocated buffer to pack into.
    // - `buffer_size` (2nd arg): on input, capacity of buffer; on output, actual bytes written.
    // The char*& for the first argument in C++ is to allow the C++ function to easily advance an internal pointer.
    // When called from C, a char* is passed. This is fine.
    return instance->pack_data(buffer_inout, *buffer_size_inout);
}

int eos_c_unpack_data(const char* buffer_in, int buffer_size_in) {
    if (!buffer_in || buffer_size_in <= 0) {
        std::cerr << "C-API Error (eos_c_unpack_data): Invalid buffer_in or buffer_size_in." << std::endl;
        return -1;
    }

    // Unpack can re-initialize, so it might create the instance.
    if (!global_eos_instance) {
        std::cout << "C-API (eos_c_unpack_data): Instance is null, creating new one for unpack." << std::endl;
        global_eos_instance = std::make_unique<EquationOfStateV1>();
        if (!global_eos_instance) {
             std::cerr << "C-API Error (eos_c_unpack_data): Failed to allocate EquationOfStateV1 instance for unpack." << std::endl;
            return -3; // Allocation failure
        }
    }
    EquationOfStateV1* instance = global_eos_instance.get(); // Should exist now
    std::cout << "C-API (eos_c_unpack_data): Calling C++ unpack_data with buffer size " << buffer_size_in << std::endl;
    return instance->unpack_data(buffer_in, buffer_size_in);
}

// --- Data Check ---

int eos_c_check_eos_data_dir(const char* data_dir_c_str, int data_dir_len,
                             const int* eos_ids_to_check_arr, int num_ids_to_check) {
    // For check_eos_data_dir, an instance might not be strictly necessary if it only checks files.
    // However, our current C++ check_eos_data_dir uses analytic_eos_registry_ from an instance.
    // So, either create a temporary instance or require initialization first.
    // Let's require initialization first, like other operations.
    EquationOfStateV1* instance = get_eos_instance_for_api("eos_c_check_eos_data_dir");
    if (!instance) return EOS_ERROR_UNKNOWN_EOS_ID -10;


    if (!data_dir_c_str || data_dir_len <= 0) {
         std::cerr << "C-API Error (eos_c_check_eos_data_dir): Invalid data_dir_c_str or data_dir_len." << std::endl;
        return -1;
    }
    if (num_ids_to_check < 0 || (num_ids_to_check > 0 && !eos_ids_to_check_arr)) {
         std::cerr << "C-API Error (eos_c_check_eos_data_dir): Invalid eos_ids_to_check_arr or num_ids_to_check." << std::endl;
        return -2;
    }

    std::string data_dir(data_dir_c_str, static_cast<size_t>(data_dir_len));
    std::vector<int> ids_to_check_vec;
    if (num_ids_to_check > 0) {
        ids_to_check_vec.assign(eos_ids_to_check_arr, eos_ids_to_check_arr + num_ids_to_check);
    }
    // std::cout << "C-API (eos_c_check_eos_data_dir): Calling C++ check_eos_data_dir." << std::endl;
    return instance->check_eos_data_dir(data_dir, ids_to_check_vec);
}

} // extern "C"