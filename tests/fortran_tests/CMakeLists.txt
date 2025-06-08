if(NOT TARGET EquationOfStateV1_lib OR NOT TARGET eos_fortran_core_obj)
    message(FATAL_ERROR "Fortran tests require EquationOfStateV1_lib and eos_fortran_core_obj targets.")
endif()

add_executable(FortranNumericalTest test_eos_numerical_f.f90)

target_link_libraries(FortranNumericalTest PRIVATE
        EquationOfStateV1_lib
        eos_fortran_core_obj
)

add_test(NAME EosFortranNumericalTest COMMAND FortranNumericalTest)
