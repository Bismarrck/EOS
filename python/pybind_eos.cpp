// python/pybind_eos.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>      // For std::vector, std::map automatic conversion
#include <pybind11/functional.h> // For std::function if needed for callbacks

#include "../src/cpp/EquationOfStateV1.h"

namespace py = pybind11;

// Helper to convert C-style istat to Python exception
// You can create custom Python exceptions if desired.
void check_eos_status(int istat, const std::string& operation_name = "EOS operation") {
    if (istat != EOS_SUCCESS) { // Assuming EOS_SUCCESS = 0
        // You might want to map specific istat values to different Python exception types
        // or include more details from the C++ error codes.
        throw std::runtime_error(operation_name + " failed with status code: " + std::to_string(istat));
    }
}

PYBIND11_MODULE(PyEOS, m) { // Module name is PyEOS
    m.doc() = "Python bindings for the Equation of State (EOS) library";

    // Bind the EquationOfStateV1 class
    py::class_<EquationOfStateV1>(m, "EquationOfState") // Python class name will be "EquationOfState"
        .def(py::init<>(), "Default constructor")

        // Bind methods from EquationOfStateV1
        .def("initialize", [](EquationOfStateV1 &self, const std::vector<int>& eos_id_list, const std::string& eos_data_dir) {
            int istat = self.initialize(eos_id_list, eos_data_dir);
            check_eos_status(istat, "initialize");
        }, py::arg("eos_id_list"), py::arg("eos_data_dir"),
           "Initializes the EOS manager with a list of EOS IDs and the data directory path.")

        .def("set_use_tfd_data_ver1", [](EquationOfStateV1 &self, bool value) {
            int istat = self.setUseTFDDataVer1(value);
            check_eos_status(istat, "set_use_tfd_data_ver1");
        }, py::arg("value"),
           "Sets the preference for TFD data version (true for v1, false for v2).")

        .def("set_complicated_eos_use_old_cold_term", [](EquationOfStateV1 &self, bool value) {
            int istat = self.setComplicatedEOSUseOldColdTerm(value);
            check_eos_status(istat, "set_complicated_eos_use_old_cold_term");
        }, py::arg("value"),
           "Sets the control variable for using the old cold term in complicated EOS.")

        .def("compute", [](EquationOfStateV1 &self, int eos_id, double rho, double T) {
            double P, E, dPdT, dEdT, dPdrho;
            int istat = self.compute(eos_id, rho, T, P, E, dPdT, dEdT, dPdrho);
            check_eos_status(istat, "compute for EOS ID " + std::to_string(eos_id));
            // Return results as a Python tuple
            return std::make_tuple(P, E, dPdT, dEdT, dPdrho);
        }, py::arg("eos_id"), py::arg("rho"), py::arg("T"),
           "Computes EOS properties. Returns (P, E, dPdT, dEdT, dPdrho).")

        .def("pack_data", [](EquationOfStateV1 &self) {
            char* buf_ptr = nullptr;
            int size = 0;
            int stat_get_size = self.pack_data(buf_ptr, size); // Get size
            check_eos_status(stat_get_size, "pack_data (get_size)");

            if (size == 0) return py::bytes(""); // Return empty Python bytes object

            std::vector<char> buffer_vec(size); // Use std::vector as intermediate buffer
            char* actual_buf_ptr = buffer_vec.data();
            // The C++ pack_data expects char*& for its first argument in the internal logic,
            // but when buffer is not nullptr, it means it's an output buffer.
            // We pass the pointer to the data in our vector.
            int stat_pack = self.pack_data(actual_buf_ptr, size); // size is in/out for actual bytes
            check_eos_status(stat_pack, "pack_data (to_buffer)");

            return py::bytes(buffer_vec.data(), size); // Create Python bytes from the vector's data
        }, "Packs the EOS manager's state into a Python bytes object.")

        .def("unpack_data", [](EquationOfStateV1 &self, py::bytes buffer_py) {
            std::string buffer_str(buffer_py); // py::bytes to std::string (implicitly gets char* and len)
            int istat = self.unpack_data(buffer_str.data(), buffer_str.length());
            check_eos_status(istat, "unpack_data");
        }, py::arg("buffer"),
           "Unpacks data from a Python bytes object and initializes the EOS manager.")

        .def("check_eos_data_dir", [](EquationOfStateV1 &self, const std::string& eos_data_dir, const std::vector<int>& eos_ids_to_check) {
            int istat = self.check_eos_data_dir(eos_data_dir, eos_ids_to_check);
            check_eos_status(istat, "check_eos_data_dir");
        }, py::arg("eos_data_dir"), py::arg("eos_ids_to_check"),
           "Checks if the EOS data directory contains necessary files.")

        .def("free_resources", &EquationOfStateV1::free_resources,
             "Frees resources held by the EOS manager.");

    // You can also bind global constants like EOS_SUCCESS if needed
    m.attr("EOS_SUCCESS") = py::int_(EOS_SUCCESS);
    // Add other error codes if you want them accessible from Python
    m.attr("EOS_ERROR_FILE_NOT_FOUND") = py::int_(EOS_ERROR_FILE_NOT_FOUND);


    // If you had enums in EquationOfStateV1.h visible to pybind11, you'd bind them:
    // py::enum_<EquationOfStateV1::InternalEOSType>(m, "InternalEOSType")
    //     .value("UNKNOWN", EquationOfStateV1::InternalEOSType::UNKNOWN)
    //     // ... other values
    //     .export_values();
}