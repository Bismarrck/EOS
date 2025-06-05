# python/CMakeLists.txt
if(NOT TARGET EquationOfStateV1_lib OR NOT TARGET eos_fortran_core_obj)
    message(FATAL_ERROR "EquationOfStateV1_lib or eos_fortran_core_obj target not found. Ensure src/ is built first.")
endif()

# Find pybind11 using its own tools (assuming third_party/pybind11/tools is present)
# Adjust path to pybind11 if it's elsewhere relative to this CMakeLists.txt
# The THIRDPARTY_DIR is defined in the root CMakeLists.txt
include(${THIRDPARTY_DIR}/pybind11/tools/pybind11Tools.cmake)
# Alternatively, if you have pybind11 installed or use find_package:
# find_package(pybind11 REQUIRED)

# Define the Python extension module
# The output file will be <name>.so or <name>.pyd
pybind11_add_module(PyEOS SHARED
        pybind_eos.cpp
)

# Link the Python module against our C++/Fortran libraries
target_link_libraries(PyEOS PRIVATE
        EquationOfStateV1_lib   # C++ static library (includes C-API if you were using it)
        eos_fortran_core_obj    # Fortran object library
)

# Ensure pybind_eos.cpp can find EquationOfStateV1.h
# This should be handled by EquationOfStateV1_lib's PUBLIC include directory.
# And pybind11 headers should be found via pybind11_add_module or global include path.

# Optional: Installation rules for the Python module
# This usually involves installing to Python's site-packages or a specific local directory.
# pybind11_install_module(PyEOS DESTINATION my_python_packages) # Example
# For now, we'll build it in the build tree.