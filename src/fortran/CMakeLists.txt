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

set_target_properties(eos_fortran_core_obj PROPERTIES
        POSITION_INDEPENDENT_CODE ON
)

# To make .mod files easily findable by other Fortran targets
# that might compile against this module's outputs directly
# (like the Fortran tests/examples).
target_include_directories(eos_fortran_core_obj INTERFACE
        ${CMAKE_CURRENT_BINARY_DIR} # This makes .mod files findable
)

# Add Fortran compiler flags if needed
# target_compile_options(eos_fortran_core_obj PRIVATE -g -Wall) # for gfortran