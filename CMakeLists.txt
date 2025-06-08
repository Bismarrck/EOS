cmake_minimum_required(VERSION 3.15)
project(EquationOfStateProject VERSION 1.0.0 LANGUAGES CXX Fortran) # Added VERSION

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Symbol Visibility & Paths ---
if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
set(THIRDPARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party)

# --- Include Package Helpers ---
include(GNUInstallDirs) # For CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR
include(CMakePackageConfigHelpers) # For configure_package_config_file, write_basic_package_version_file

# --- Source Subdirectories (Define library components) ---
# These add_subdirectory calls ensure EquationOfStateV1_lib (STATIC) and eos_fortran_core_obj (OBJECT)
# are defined. src/cpp/CMakeLists.txt should contain:
#   install(TARGETS EquationOfStateV1_lib EXPORT EOSCoreTargets ARCHIVE ... )
add_subdirectory(src/fortran)
add_subdirectory(src/cpp)

# --- Define the Core Shared Library ---
add_library(EOSCore SHARED
        $<TARGET_OBJECTS:EquationOfStateV1_lib>
        $<TARGET_OBJECTS:eos_fortran_core_obj>
)
set_target_properties(EOSCore PROPERTIES OUTPUT_NAME "EOSCore")
target_link_libraries(EOSCore PUBLIC EquationOfStateV1_lib) # For usage requirements

# This adds EOSCore to the EOSCoreTargets export set AND defines its install rule
install(TARGETS EOSCore
        EXPORT EOSCoreTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT DevelopmentLibraries
)

# 1. Generate EOSCoreTargets.cmake in the build directory for build-tree usage.
#    This uses the targets that have been associated with the 'EOSCoreTargets' export set.
export(EXPORT EOSCoreTargets # Use the EXPORT SET name here
        FILE "${CMAKE_CURRENT_BINARY_DIR}/EOSCoreTargets.cmake"
        NAMESPACE EOS::
)

# --- Configure CMake Package Files for both build-tree and install-tree ---
configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/EOSCoreConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfig.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/EOSCore"
        PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR
)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

# --- Prepare for find_package() in sub-projects (examples, tests) within this build ---
list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})


# --- BUILD Options (Python, Examples, Tests) ---
# These will be processed *after* CMAKE_PREFIX_PATH is updated.
option(EOS_BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
if(EOS_BUILD_PYTHON_BINDINGS)
    find_package(PythonInterp REQUIRED) # Ensure Python interpreter is found before pybind11 subdir
    include_directories(SYSTEM ${THIRDPARTY_DIR}/pybind11/include) # Ensure pybind11 headers are findable
    add_subdirectory(${THIRDPARTY_DIR}/pybind11) # If pybind11 is built as a target
    add_subdirectory(python)
endif()

option(EOS_BUILD_TESTS "Build the tests" ON)
if(EOS_BUILD_TESTS)
    enable_testing()
    include_directories(SYSTEM ${THIRDPARTY_DIR}/doctest)
    add_subdirectory(tests)
endif()

option(EOS_BUILD_EXAMPLES "Build example programs" ON)
if(EOS_BUILD_EXAMPLES)
    add_subdirectory(examples) # examples/cpp_example/CMakeLists.txt will call find_package(EOSCore)
endif()

# --- Installation Rules ---)

# Install C++ headers
install(FILES src/cpp/EquationOfStateV1.h src/cpp/utils/string_utils.h
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/EOSCore"
        COMPONENT DevelopmentHeaders
)

# Install the configured EOSCoreConfig.cmake and EOSCoreConfigVersion.cmake files
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/EOSCoreConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/EOSCore"
        COMPONENT DevelopmentFiles
)

# Install the export set's targets file for the install tree
# This creates <prefix>/lib/cmake/EOSCore/EOSCoreTargets.cmake with install-tree paths
install(EXPORT EOSCoreTargets
        FILE EOSCoreTargets.cmake # Filename within the DESTINATION directory
        NAMESPACE EOS::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/EOSCore"
        COMPONENT DevelopmentFiles
)