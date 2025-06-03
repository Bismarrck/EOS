cmake_minimum_required(VERSION 3.15)
project(EquationOfStateProject LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    # For MinGW gfortran with MSVC C++, BIND(C) symbols might need __declspec(dllimport/dllexport)
    # or specific handling if CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS isn't sufficient for Fortran.
    # For Intel Fortran with MSVC, CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS usually works well.
endif()

# --- Define Library Output Path ---
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Add Source Subdirectories (they define object libraries) ---
add_subdirectory(src/fortran) # Defines eos_fortran_core_obj
add_subdirectory(src/cpp)     # Defines eos_v1_obj

# --- Define the EOS Library (Shared and Static) ---
# Combine C++ and Fortran object files into the final libraries
set(EOS_LIBRARY_OBJECTS $<TARGET_OBJECTS:eos_v1_obj> $<TARGET_OBJECTS:eos_fortran_core_obj>)

# Shared Library
add_library(EquationOfState SHARED ${EOS_LIBRARY_OBJECTS}) # Main target name
set_target_properties(EquationOfState PROPERTIES PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/src/cpp/EquationOfStateV1.h")
# On Windows, this will create EquationOfState.dll and EquationOfState.lib (import library)
# On Linux, EquationOfState.so

# Static Library
add_library(EquationOfState_static STATIC ${EOS_LIBRARY_OBJECTS})
set_target_properties(EquationOfState_static PROPERTIES OUTPUT_NAME EquationOfState) # Keep output name simple if desired
# On Windows, EquationOfState.lib
# On Linux, libEquationOfState.a


# --- Installation Rules (Phase 6) ---
include(GNUInstallDirs) # Provides CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR etc.
install(TARGETS EquationOfState EquationOfState_static
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # For DLLs on Windows
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/EquationOfState # For EquationOfStateV1.h
)
# For the PUBLIC_HEADER to work well, EquationOfStateV1.h should be listed in
# target_sources for the 'EquationOfState' library, or use install(FILES ...)
install(FILES ${PROJECT_SOURCE_DIR}/src/cpp/EquationOfStateV1.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/EquationOfState
)


# --- Build Tests (Optional) ---
option(BUILD_TESTS "Build the project's tests" ON) # Default to ON
if(BUILD_TESTS)
    enable_testing() # For CTest
    add_subdirectory(tests)
endif()

# --- Build Examples (Optional, for later) ---
# option(BUILD_EXAMPLES "Build example programs" ON)
# if(BUILD_EXAMPLES)
#     add_subdirectory(examples)
# endif()