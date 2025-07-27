cmake_minimum_required(VERSION 3.15)
project(FortranEOSDemo LANGUAGES Fortran C CXX)

# Source file for the main example program
set(MAIN_F90_SRC main.f90)

# Sources for the FortranEOSDemo executable
list(APPEND DEMO_SOURCES ${MAIN_F90_SRC})

if(EOS_LIBRARY_IS_BUILDING_EXAMPLES)
    if(NOT TARGET EOSCore OR NOT TARGET eos_fortran_core_obj OR NOT TARGET EquationOfStateV1_lib)
        message(FATAL_ERROR "FortranEOSDemo: Required library targets not found.")
    endif()

    add_executable(FortranEOSDemo ${MAIN_F90_SRC})

    # Link against the build-tree targets
    target_link_libraries(FortranEOSDemo PRIVATE
            EOSCore                  # For the C-API symbols
            eos_fortran_core_obj     # For any other Fortran modules used by the interface or main
            # and to propagate .mod include path for module_eos_fortran_interface
    )
    target_include_directories(FortranEOSDemo PRIVATE ${HDF5_INCLUDE_DIRS})

else()
    # --- Building as a standalone example (e.g., CI testing installed library) ---
    message(STATUS "FortranEOSDemo: Building against installed EOSCore (using find_package).")

    # CMAKE_PREFIX_PATH should be set by the CI script (or user) to the install_dir
    find_package(EOSCore 1.0.0 REQUIRED CONFIG)

    if(NOT DEFINED EOSCore_FORTRAN_MODULE_DIRECTORY)
        # Get the path to the installed Fortran interface module source file
        if(NOT DEFINED EOSCore_FORTRAN_INTERFACE_MODULE_FILE OR NOT EXISTS ${EOSCore_FORTRAN_INTERFACE_MODULE_FILE})
            message(FATAL_ERROR "Installed module_eos_fortran_interface.f90 not found. "
                    "EOSCore_FORTRAN_INTERFACE_MODULE_FILE was not set or path is invalid: "
                    "'${EOSCore_FORTRAN_INTERFACE_MODULE_FILE}'. "
                    "Ensure EOSCoreConfig.cmake defines it and the file was installed.")
        endif()

        # Add the installed Fortran interface module source to be compiled with this example
        list(APPEND DEMO_SOURCES ${EOSCore_FORTRAN_INTERFACE_MODULE_FILE})
    endif ()

    add_executable(FortranEOSDemo ${DEMO_SOURCES})

    # When module_eos_fortran_interface.f90 (from installed source) is compiled as part of this target,
    # its .mod file is generated in CMAKE_CURRENT_BINARY_DIR (e.g., build_example_fortran/CMakeFiles/...).
    # main.f90 needs to find it there for its USE statement.
    target_include_directories(FortranEOSDemo PRIVATE
            ${CMAKE_CURRENT_BINARY_DIR}
            ${HDF5_INCLUDE_DIRS})

    # Link against the imported target EOS::EOSCore for the C-API symbols
    target_link_libraries(FortranEOSDemo PRIVATE EOS::EOSCore)
endif()


# CTest registration (this registers the example for testing when built as part of main project)
if(EOS_BUILD_TESTS AND EOS_LIBRARY_IS_BUILDING_EXAMPLES)
    add_test(NAME EosFortranExampleAsTest COMMAND FortranEOSDemo)
    # Set working directory if main.f90 expects data files relative to its executable location
    # For example, if main.f90 uses data_dir = "../../../eos_data_dir"
    # and FortranEOSDemo is in build/examples/fortran_example/
    set_tests_properties(EosFortranExampleAsTest PROPERTIES
            LABELS "fortran_example"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endif()
