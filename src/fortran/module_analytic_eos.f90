! src/fortran/module_analytic_eos.f90
MODULE module_analytic_eos
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

CONTAINS

    SUBROUTINE air_eos_2000(rho, T, P, E)
        REAL(C_DOUBLE), INTENT(IN) :: rho, T
        REAL(C_DOUBLE), INTENT(OUT) :: P, E

        ! Placeholder implementation
        PRINT *, "Fortran (air_eos_2000): rho=", rho, " T=", T
        P = rho * T * 1.0_C_DOUBLE  ! Example calculation
        E = rho * T * 1.5_C_DOUBLE  ! Example calculation
    END SUBROUTINE air_eos_2000

    SUBROUTINE carbon_eos_2001(rho, T, P, E, istat)
        REAL(C_DOUBLE), INTENT(IN) :: rho, T
        REAL(C_DOUBLE), INTENT(OUT) :: P, E
        INTEGER(C_INT), INTENT(OUT) :: istat

        ! Placeholder implementation
        PRINT *, "Fortran (carbon_eos_2001): rho=", rho, " T=", T
        IF (rho < 0.0_C_DOUBLE .OR. T < 0.0_C_DOUBLE) THEN
            P = 0.0_C_DOUBLE
            E = 0.0_C_DOUBLE
            istat = 101 ! Invalid input error code
            PRINT *, "Fortran (carbon_eos_2001): Error - invalid input."
            RETURN
        END IF

        P = rho * T * 2.0_C_DOUBLE  ! Example calculation
        E = rho * T * 2.5_C_DOUBLE  ! Example calculation
        istat = 0 ! Success
    END SUBROUTINE carbon_eos_2001

END MODULE module_analytic_eos