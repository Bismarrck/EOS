// src/cpp/python_bindings.cpp
#include "EquationOfStateV1.h" // Your C++ class
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>       // For std::vector, std::map conversions
#include <pybind11/iostream.h>  // For redirecting C++ cout/cerr to Python
#include <stdexcept>            // For throwing Python-catchable exceptions
#include <string>               // For std::to_string

namespace py = pybind11;

// Helper to convert C++ error codes to Python exceptions if desired
// Or Python side can just check the integer status.
// For now, C++ methods return int status; Python methods will also return this.
// A more Pythonic way is to throw exceptions on error.

PYBIND11_MODULE(PyEOS, m) { // Module name prefixed with _ is common for C extensions
    m.doc() = "PyEOS: Python bindings for the EquationOfStateV1 C++ library";

    // Expose error constants (optional, but can be useful for Python side)
    m.attr("EOS_SUCCESS") = py::int_(EOS_SUCCESS);
    m.attr("EOS_ERROR_FILE_NOT_FOUND") = py::int_(EOS_ERROR_FILE_NOT_FOUND);
    // ... add other relevant error codes you defined in EquationOfStateV1.h

    py::class_<EquationOfStateV1>(m, "EquationOfState") // Python class name
        .def(py::init<>()) // Default constructor
        .def("initialize", &EquationOfStateV1::initialize,
             py::arg("eos_id_list"), py::arg("eos_data_dir"),
             "Initializes the EOS manager with a list of EOS IDs and the data directory.")
        .def("free_resources", &EquationOfStateV1::free_resources,
             "Frees loaded EOS data and TFD matrices.")
        .def("check_eos_data_dir", &EquationOfStateV1::check_eos_data_dir,
             py::arg("eos_data_dir"), py::arg("eos_ids_to_check"),
             "Checks if the EOS data directory contains necessary files.")
        .def("set_use_tfd_data_ver1", &EquationOfStateV1::setUseTFDDataVer1,
             py::arg("value"),
             "Sets the preference for TFD data version (true for v1, false for v2).")
        .def("set_complicated_eos_use_old_cold_term", &EquationOfStateV1::setComplicatedEOSUseOldColdTerm,
             py::arg("value"),
             "Sets the control variable for using old cold term in complicated EOS.")
        .def("compute",
             [](EquationOfStateV1 &self, int eos_id, double rho, double T) {
                 double P, E, dPdT, dEdT, dPdrho;
                 // Redirect C++ streams to Python's stdout/stderr for the duration of this call
                 py::scoped_ostream_redirect stream_redir_out(std::cout, py::module_::import("sys").attr("stdout"));
                 py::scoped_ostream_redirect stream_redir_err(std::cerr, py::module_::import("sys").attr("stderr"));

                 int istat = self.compute(eos_id, rho, T, P, E, dPdT, dEdT, dPdrho);
                 // In Python, it's common to return a tuple or a custom result object.
                 // And to raise an exception on failure instead of returning status code.
                 if (istat != EOS_SUCCESS) {
                     throw std::runtime_error("EOS computation failed for ID " + std::to_string(eos_id) +
                                              " with istat: " + std::to_string(istat));
                 }
                 // Return a dictionary for named results, or a tuple
                 py::dict results;
                 results["P"] = P;
                 results["E"] = E;
                 results["dPdT"] = dPdT;
                 results["dEdT"] = dEdT;
                 results["dPdrho"] = dPdrho;
                 results["istat"] = istat; // Still include istat for completeness
                 return results;
             },
             py::arg("eos_id"), py::arg("rho"), py::arg("T"),
             "Computes EOS properties. Returns a dictionary with P, E, dPdT, dEdT, dPdrho, istat.")
        .def("pack_data",
             [](EquationOfStateV1 &self) {
                 // Redirect C++ streams
                 py::scoped_ostream_redirect stream_redir_out(std::cout, py::module_::import("sys").attr("stdout"));
                 py::scoped_ostream_redirect stream_redir_err(std::cerr, py::module_::import("sys").attr("stderr"));

                 char* buf_ptr = nullptr;
                 int size = 0;
                 int stat_get_size = self.pack_data(buf_ptr, size); // Get size
                 if (stat_get_size != EOS_SUCCESS || size == 0) {
                     // If size is 0, it could be valid (no data to pack), or an error.
                     // pack_data should ideally distinguish this. For now, assume error if stat isn't success.
                     if (stat_get_size != EOS_SUCCESS) {
                        throw std::runtime_error("pack_data: failed to get buffer size, istat: " + std::to_string(stat_get_size));
                     }
                     return py::bytes(""); // Return empty bytes if size is 0 and successful
                 }

                 std::vector<char> buffer_vec(size);
                 buf_ptr = buffer_vec.data(); // pybind11 examples show this pattern
                                             // but EquationOfStateV1::pack_data expects char*&.
                                             // The C-API used char*, let's assume pack_data is:
                                             // int pack_data(char* buffer_to_fill, int& size_capacity_in_actual_out)
                                             // If pack_data remains char*&, this needs care.
                                             // My current C++ pack_data `int pack_data(char*& buffer, int& buffer_size)`
                                             // in mode 2 (buffer not null) expects `buffer` to be a valid pointer
                                             // and `buffer_size` to be its capacity. It writes into it and updates
                                             // `buffer_size` to actual written. The `char*&` for `buffer` is for
                                             // internal pointer advancement in the implementation of pack_data, not
                                             // for reallocating or returning a new buffer pointer.
                                             // So passing `buffer_vec.data()` should be fine.

                 int actual_size = size; // Input: capacity
                 int stat_pack = self.pack_data(buf_ptr, actual_size); // Output: actual_size
                 if (stat_pack != EOS_SUCCESS) {
                     throw std::runtime_error("pack_data: failed to pack data, istat: " + std::to_string(stat_pack));
                 }
                 return py::bytes(buffer_vec.data(), actual_size);
             },
             "Packs the EOS manager's state into a Python bytes object.")
        .def("unpack_data",
             [](EquationOfStateV1 &self, py::bytes buffer_py) {
                 // Redirect C++ streams
                 py::scoped_ostream_redirect stream_redir_out(std::cout, py::module_::import("sys").attr("stdout"));
                 py::scoped_ostream_redirect stream_redir_err(std::cerr, py::module_::import("sys").attr("stderr"));

                 std::string buffer_str(buffer_py); // py::bytes to std::string (implicitly gets char* and len)
                 int istat = self.unpack_data(buffer_str.data(), buffer_str.length());
                 if (istat != EOS_SUCCESS) {
                     throw std::runtime_error("unpack_data: failed with istat: " + std::to_string(istat));
                 }
                 return istat; // Or return None on success if exceptions are primary error mechanism
             },
             py::arg("buffer"),
             "Unpacks data from a Python bytes object to initialize the EOS manager.");

    // Add other methods or enums as needed.
    // py::enum_<EquationOfStateV1::EOSTypeFromFile>(m, "EOSTypeFromFile")
    //     .value("NOT_SET", EquationOfStateV1::EOSTypeFromFile::NOT_SET)
    //     .value("ANALYTIC", EquationOfStateV1::EOSTypeFromFile::ANALYTIC)
    //     .value("COMPLICATED", EquationOfStateV1::EOSTypeFromFile::COMPLICATED)
    //     .export_values();
}
