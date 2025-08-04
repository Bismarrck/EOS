# tools/CMakeLists.txt
project(EOSTools CXX)
cmake_policy(SET CMP0074 NEW)

# Check if the main library target EOSCore exists from the parent build
if(NOT TARGET EOSCore)
    message(FATAL_ERROR "The 'eos-tool' executable requires the 'EOSCore' library target, which was not found. Ensure the main project is configured first.")
endif()

# Define the command-line tool executable
add_executable(eos-tool eos-tool.cpp)

# Link the tool against our main EOSCore library.
# This provides all the necessary functions from EquationOfStateV1 and the Fortran routines.
target_link_libraries(eos-tool PRIVATE EOSCore)

# HDF5
target_link_libraries(eos-tool PRIVATE HDF5::HDF5)
target_include_directories(eos-tool PUBLIC ${HDF5_INCLUDE_DIRS})

# --- Installation Rule ---
# This command specifies where to install the 'eos-tool' executable.
# It uses GNUInstallDirs to get standard directory names like CMAKE_INSTALL_BINDIR.
if(NOT CMAKE_INSTALL_BINDIR)
    include(GNUInstallDirs)
endif()

install(TARGETS eos-tool
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # Installs to bin/ on Linux/macOS, or root of prefix on Windows
        COMPONENT Applications # Group this installation component
)

# Optional: Add the tool to CTest for a basic "runs without crashing" test
# This is NOT a functional test, just a sanity check.
if(EOS_BUILD_TESTS) # Check if parent scope enabled tests
    add_test(
            NAME EOSToolRuns
            COMMAND eos-tool --help # Running with --help is a simple way to check if it links and starts
    )
    # The test will pass if the command returns an exit code of 0.
    # CLI11 returns 0 when --help is successfully processed.
    set_tests_properties(EOSToolRuns PROPERTIES LABELS "tool")
endif()