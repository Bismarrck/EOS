cmake_minimum_required(VERSION 3.15)
project(FortranEOSDemo LANGUAGES Fortran)

if(NOT TARGET EOSCore OR NOT TARGET eos_fortran_core_obj OR NOT TARGET EquationOfStateV1_lib)
    message(FATAL_ERROR "FortranEOSDemo: Required library targets not found.")
endif()

# This example compiles the interface module itself.
add_executable(FortranEOSDemo main.f90)

if(EOS_LIBRARY_IS_BUILDING_EXAMPLES)
    # Link against the build-tree targets
    target_link_libraries(FortranEOSDemo PRIVATE
            EOSCore                  # For the C-API symbols
            eos_fortran_core_obj     # For any other Fortran modules used by the interface or main
            # and to propagate .mod include path for module_eos_fortran_interface
    )
    # Ensure main.f90 can find module_eos_fortran_interface.mod (generated in this target's build dir)
    # and module_eos_fortran_interface.f90 can find .mod files from eos_fortran_core_obj
    target_include_directories(FortranEOSDemo PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    # The INTERFACE include from eos_fortran_core_obj should provide paths to its own .mod files
else()
    # --- Building as a standalone example (e.g., CI testing installed library) ---
    message(STATUS "FortranEOSDemo: Building against installed EOSCore (using find_package).")

    # CMAKE_PREFIX_PATH should be set by the CI script (or user) to point to the install_dir
    find_package(EOSCore 1.0.0 REQUIRED) # Finds installed library

    # Link against the imported target (e.g., EOS::EOSCore)
    target_link_libraries(FortranEOSDemo PRIVATE EOS::EOSCore)
endif()
