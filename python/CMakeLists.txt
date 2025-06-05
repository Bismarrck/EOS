# EquationOfStateProject/python/CMakeLists.txt
if(NOT TARGET EquationOfStateV1_lib OR NOT TARGET eos_fortran_core_obj)
    message(FATAL_ERROR "EquationOfStateV1_lib or eos_fortran_core_obj target not found. Build src/ first.")
endif()

# Find pybind11 - Option 1: If pybind11 was added as a proper subdirectory to root CMake
# find_package(pybind11 CONFIG REQUIRED) # This requires pybind11 to have been built/installed or properly subdirected

# Find pybind11 - Option 2: Manually include its tools and headers (more direct for vendored headers)
# Adjust path to where you placed pybind11's tools directory
include(${CMAKE_SOURCE_DIR}/third_party/pybind11/tools/pybind11Tools.cmake)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/third_party/pybind11/include)

# Define the Python extension module
# Module name should match PYBIND11_MODULE macro in python_bindings.cpp (PyEOS)
# Sources are relative to THIS CMakeLists.txt file.
pybind11_add_module(PyEOS
        SHARED # Or MODULE
        ${CMAKE_SOURCE_DIR}/src/cpp/python_bindings.cpp # Path to your binding code
)

# Link against your C++ static library and the Fortran object library
target_link_libraries(PyEOS PRIVATE
        EquationOfStateV1_lib
        eos_fortran_core_obj
)

# Optional: Install the Python module to a standard location or site-packages
# pybind11_install(PyEOS) # Needs pybind11_DIR to be set for find_package(pybind11) version
# Manual install:
# install(TARGETS PyEOS DESTINATION your_python_module_install_path)