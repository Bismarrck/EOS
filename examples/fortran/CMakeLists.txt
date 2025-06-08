cmake_minimum_required(VERSION 3.15)
project(FortranEOSDemo LANGUAGES Fortran)

if(NOT TARGET EOSCore OR NOT TARGET eos_fortran_core_obj OR NOT TARGET EquationOfStateV1_lib)
    message(FATAL_ERROR "FortranEOSDemo: Required library targets not found.")
endif()

# This example compiles the interface module itself.
add_executable(FortranEOSDemo main.f90)

# Link against the build-tree targets
target_link_libraries(FortranEOSDemo PRIVATE
        EOSCore                  # For the C-API symbols
        eos_fortran_core_obj     # For any other Fortran modules used by the interface or main
        # and to propagate .mod include path for module_eos_fortran_interface
)
# Ensure main.f90 can find module_eos_fortran_interface.mod (generated in this target's build dir)
# and module_eos_fortran_interface.f90 can find .mod files from eos_fortran_core_obj
target_include_directories(FortranEOSDemo PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# The INTERFACE include from eos_fortran_core_obj should provide paths to its own .mod files.