cmake_minimum_required(VERSION 3.15)
project(EquationOfStateProject LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Path to third-party libraries
set(THIRDPARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party) # Corrected path

# Add src subdirectories
add_subdirectory(src/fortran)
add_subdirectory(src/cpp)

# --- BUILD: Python Interface ---
option(EOS_BUILD_PYTHON_BINDINGS "Build Python bindings for EquationOfStateProject" ON)
if(EOS_BUILD_PYTHON_BINDINGS)
    find_package(PythonInterp REQUIRED)
    add_subdirectory(${THIRDPARTY_DIR}/pybind11)
    add_subdirectory(python)
endif()

# --- BUILD: Examples ---
option(EOS_BUILD_EXAMPLES "Build example programs" ON)
if(EOS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# --- BUILD: Tests ---
option(EOS_BUILD_TESTS "Build the tests for EquationOfStateProject" ON)
if(EOS_BUILD_TESTS)
    enable_testing() # Enables CTest
    add_subdirectory(tests)
    # Make doctest.h available globally for test targets
    include_directories(SYSTEM ${THIRDPARTY_DIR}/doctest) # Corrected path
endif()
