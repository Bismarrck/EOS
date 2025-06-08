cmake_minimum_required(VERSION 3.15)
project(EquationOfStateProject LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Symbol Visibility ---
# For creating a well-behaved shared library, especially on Windows and for C++ classes.
# CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS is a broad approach.
# For more fine-grained control (good practice for libraries not just exposing a C API):
# include(GenerateExportHeader)
# However, for a library primarily exposing a C-API (eos_c_api.cpp) and used by pybind11
# which binds to the C++ class, CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS might be sufficient to start.
# If you were to distribute the C++ headers for others to link against EquationOfStateV1 directly,
# then GenerateExportHeader is better for C++ class exports.
if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    # For MinGW when building shared libraries that link Fortran:
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--export-all-symbols")
    endif()
endif()
# For GCC/Clang, explicitly set visibility for C++ if not using GenerateExportHeader
# and wanting to hide internal symbols. C API functions are usually fine.
# set(CMAKE_CXX_VISIBILITY_PRESET hidden)
# set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# Path to third-party libraries
set(THIRDPARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party) # Corrected path

# Add src subdirectories
add_subdirectory(src/fortran)
add_subdirectory(src/cpp)

# --- Define the Core Shared Library ---
add_library(EOSCore SHARED
        $<TARGET_OBJECTS:EquationOfStateV1_lib>    # Objects from the C++ static lib
        $<TARGET_OBJECTS:eos_fortran_core_obj> # Objects from the Fortran object lib
)

# Set properties for the shared library
set_target_properties(EOSCore PROPERTIES
        OUTPUT_NAME "EOSCore" # Base name: EOSCore.dll, libEOSCore.so
        # If you have project version variables set:
        # VERSION ${PROJECT_VERSION_STRING} # e.g., "1.2.3"
        # SOVERSION ${PROJECT_VERSION_MAJOR} # e.g., "1"
)

# Link the C++ static library
target_link_libraries(EOSCore PUBLIC EquationOfStateV1_lib)

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

# --- Installation Rules ---
include(GNUInstallDirs)

# Install the primary shared library and its C++ headers
install(TARGETS EOSCore
        EXPORT EOSCoreTargets # For CMake package configuration (see below)
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT DevelopmentLibraries
)

install(FILES src/cpp/EquationOfStateV1.h src/cpp/utils/string_utils.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/EOSCore # Install headers to include/EOSCore
        COMPONENT DevelopmentHeaders
)

# Optional: Install the static C++ library if you want to distribute it as well
install(TARGETS EquationOfStateV1_lib
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT DevelopmentLibraries
)

# --- CMake Package Configuration (for find_package) ---
# This makes your library easily findable by other CMake projects.
# 1. Create EOSCoreConfig.cmake.in and EOSCoreConfigVersion.cmake.in in a 'cmake' subdir
#    (Example content for EOSCoreConfig.cmake.in):
#    @PACKAGE_INIT@
#    include("${CMAKE_CURRENT_LIST_DIR}/EOSCoreTargets.cmake")
#    check_required_components(EOSCore)

# 2. Configure and install these files
# include(CMakePackageConfigHelpers) # Provides helper macros
# configure_package_config_file(
#     ${CMAKE_CURRENT_SOURCE_DIR}/cmake/EOSCoreConfig.cmake.in
#     ${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfig.cmake
#     INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/EOSCore # Standard location
#     PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR
# )
# write_basic_package_version_file(
#     ${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfigVersion.cmake
#     VERSION "1.0.0" # Replace with your project version
#     COMPATIBILITY AnyNewerVersion
# )
# install(FILES
#     ${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfig.cmake
#     ${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfigVersion.cmake
#     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/EOSCore
#     COMPONENT DevelopmentFiles
# )
# install(EXPORT EOSCoreTargets # Generated by install(TARGETS EOSCore EXPORT ...)
#     FILE EOSCoreTargets.cmake
#     NAMESPACE EOS:: # Optional namespace for imported targets
#     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/EOSCore
#     COMPONENT DevelopmentFiles
# )
