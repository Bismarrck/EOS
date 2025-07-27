# src/fortran/CMakeLists.txt
add_library(eos_fortran_core_obj OBJECT
        module_analytic_eos.f90
        module_global_controls.f90
        module_complicated_eos.f90
        module_tfd_eos.f90
        fortran_c_wrappers.f90
        module_eos_fortran_interface.f90
        module_poly_eos.f90
)

# Set the directory where .mod files will be generated
# This helps standardize the location within the build tree.
# We'll use a directory relative to the project's include path.
set(FORTRAN_MODULE_PATH ${CMAKE_INSTALL_INCLUDEDIR}/EOSCore/fortran_modules)
set_target_properties(eos_fortran_core_obj PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FORTRAN_MODULE_PATH}
        POSITION_INDEPENDENT_CODE ON
)

# For build-tree usage ensure this directory is part of the INTERFACE includes.
target_include_directories(eos_fortran_core_obj INTERFACE
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${FORTRAN_MODULE_PATH}>
        # For installed usage, the consumer will need to add the installed .mod directory.
        # We will provide this path in the EOSCoreConfig.cmake file.
        $<INSTALL_INTERFACE:${FORTRAN_MODULE_PATH}>
)

# Conditional installation of the generated .mod files
if(EOS_INSTALL_FORTRAN_MODULES)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/${FORTRAN_MODULE_PATH}/
            DESTINATION ${FORTRAN_MODULE_PATH}
            COMPONENT FortranModules
    )
endif()