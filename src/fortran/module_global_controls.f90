! src/fortran/module_global_controls.f90
MODULE module_global_controls
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

    LOGICAL :: complicated_eos_use_old_cold_term = .FALSE.
    LOGICAL :: use_tfd_data_ver1 = .TRUE.  ! Default value

    SAVE complicated_eos_use_old_cold_term, use_tfd_data_ver1

    PUBLIC :: complicated_eos_use_old_cold_term, use_tfd_data_ver1
    ! Add any other global variables here

CONTAINS
    ! No subroutines needed here directly, setters/getters will be in wrappers module
END MODULE module_global_controls