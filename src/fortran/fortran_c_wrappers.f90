! src/fortran/fortran_c_wrappers.f90
MODULE eos_c_wrappers
    USE, INTRINSIC :: ISO_C_BINDING
    USE module_analytic_eos       ! To access air_eos_2000, carbon_eos_2001
    USE module_global_controls    ! To access control variables
    IMPLICIT NONE

CONTAINS

    ! --- Wrappers for Analytic EOS functions ---

    SUBROUTINE c_air_eos_2000(rho, T, P, E, istat_out) &
            BIND(C, name='c_air_eos_2000')
        !DEC$ ATTRIBUTES DLLEXPORT :: c_air_eos_2000
        REAL(C_DOUBLE), VALUE, INTENT(IN) :: rho, T
        REAL(C_DOUBLE), INTENT(OUT) :: P, E
        INTEGER(C_INT), INTENT(OUT) :: istat_out

        CALL air_eos_2000(rho, T, P, E)
        istat_out = 0 ! Success, as original function has no istat
    END SUBROUTINE c_air_eos_2000

    SUBROUTINE c_carbon_eos_2001(rho, T, P, E, istat_out) &
            BIND(C, name='c_carbon_eos_2001')
        !DEC$ ATTRIBUTES DLLEXPORT :: c_carbon_eos_2001
        REAL(C_DOUBLE), VALUE, INTENT(IN) :: rho, T
        REAL(C_DOUBLE), INTENT(OUT) :: P, E
        INTEGER(C_INT), INTENT(OUT) :: istat_out ! This will be the istat from carbon_eos_2001

        CALL carbon_eos_2001(rho, T, P, E, istat_out)
    END SUBROUTINE c_carbon_eos_2001

    ! --- Wrappers for Global Control Variables (Setters) ---

    SUBROUTINE c_set_complicated_eos_use_old_cold_term(value, istat_out) &
            BIND(C, name='c_set_complicated_eos_use_old_cold_term')
        !DEC$ ATTRIBUTES DLLEXPORT :: c_set_complicated_eos_use_old_cold_term
        LOGICAL(C_BOOL), VALUE, INTENT(IN) :: value
        INTEGER(C_INT), INTENT(OUT) :: istat_out

        complicated_eos_use_old_cold_term = value
        istat_out = 0 ! Success
        PRINT *, "Fortran: complicated_eos_use_old_cold_term set to: ", complicated_eos_use_old_cold_term
    END SUBROUTINE c_set_complicated_eos_use_old_cold_term

    SUBROUTINE c_set_use_tfd_data_ver1(value, istat_out) &
            BIND(C, name='c_set_use_tfd_data_ver1')
        !DEC$ ATTRIBUTES DLLEXPORT :: c_set_use_tfd_data_ver1
        LOGICAL(C_BOOL), VALUE, INTENT(IN) :: value
        INTEGER(C_INT), INTENT(OUT) :: istat_out

        use_tfd_data_ver1 = value
        istat_out = 0 ! Success
        PRINT *, "Fortran: use_tfd_data_ver1 set to: ", use_tfd_data_ver1
    END SUBROUTINE c_set_use_tfd_data_ver1

    ! --- (Optional) Wrappers for Global Control Variables (Getters) ---
    ! Example:
    ! FUNCTION c_get_use_tfd_data_ver1(istat_out) RESULT(res) &
    !     BIND(C, name='c_get_use_tfd_data_ver1')
    !     !DEC$ ATTRIBUTES DLLEXPORT :: c_get_use_tfd_data_ver1
    !     INTEGER(C_INT), INTENT(OUT) :: istat_out
    !     LOGICAL(C_BOOL) :: res
    !
    !     res = use_tfd_data_ver1
    !     istat_out = 0 ! Success
    ! END FUNCTION c_get_use_tfd_data_ver1

END MODULE eos_c_wrappers