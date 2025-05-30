# src/fortran/CMakeLists.txt
add_library(eos_fortran_core_obj OBJECT
        module_analytic_eos.f90
        module_global_controls.f90
        module_complicated_eos.f90
        module_tfd_eos.f90
        fortran_c_wrappers.f90
)

# Add Fortran compiler flags if needed
# target_compile_options(eos_fortran_core_obj PRIVATE -g -Wall) # for gfortran