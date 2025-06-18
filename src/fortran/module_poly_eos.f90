MODULE module_poly_eos

    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

CONTAINS

    SUBROUTINE poly_eos_compute_kernel(poly_coeffs, num_coeffs, rho, T, P, E, istat)

        REAL(C_DOUBLE), INTENT(IN) :: poly_coeffs(num_coeffs)
        INTEGER(C_INT), INTENT(IN) :: num_coeffs
        REAL(C_DOUBLE), INTENT(IN) :: rho, T
        REAL(C_DOUBLE), INTENT(OUT) :: P, E
        INTEGER(C_INT), INTENT(OUT) :: istat

        P = rho * poly_coeffs(1)
        E = rho * poly_coeffs(2)
        istat = 0

    END SUBROUTINE

END MODULE module_poly_eos