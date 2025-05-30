! src/fortran/module_complicated_eos.f90
MODULE module_complicated_eos
    USE, INTRINSIC :: ISO_C_BINDING
    USE module_tfd_eos, ONLY: TFD_eos_compute ! To call the TFD routine
    IMPLICIT NONE

CONTAINS

    SUBROUTINE Complicated_EOS(params_in, num_params, rho, T, &
            matrix_A_in, N1_A, N2_A, &
            matrix_B_in, N1_B, N2_B, &
            P, E, dPdT, dEdT, dPdrho, &
            istat)
        ! Inputs
        INTEGER(C_INT), INTENT(IN) :: num_params
        REAL(C_DOUBLE), DIMENSION(num_params), INTENT(IN) :: params_in
        REAL(C_DOUBLE), INTENT(IN) :: rho, T
        INTEGER(C_INT), INTENT(IN) :: N1_A, N2_A, N1_B, N2_B
        REAL(C_DOUBLE), DIMENSION(N1_A, N2_A), INTENT(IN) :: matrix_A_in
        REAL(C_DOUBLE), DIMENSION(N1_B, N2_B), INTENT(IN) :: matrix_B_in
        ! Outputs
        REAL(C_DOUBLE), INTENT(OUT) :: P, E, dPdT, dEdT, dPdrho
        INTEGER(C_INT), INTENT(OUT) :: istat

        ! Local variables for TFD results
        REAL(C_DOUBLE) :: tfd_result_x, tfd_result_y
        INTEGER(C_INT) :: tfd_istat

        PRINT *, "Fortran (Complicated_EOS): rho=", rho, " T=", T
        IF (num_params > 0) THEN
            PRINT *, "Fortran (Complicated_EOS): Received ", num_params, " parameters. First param=", params_in(1)
        ELSE
            PRINT *, "Fortran (Complicated_EOS): Received 0 parameters."
            istat = 301 ! Error: No parameters provided
            P = 0.0_C_DOUBLE; E = 0.0_C_DOUBLE; dPdT = 0.0_C_DOUBLE
            dEdT = 0.0_C_DOUBLE; dPdrho = 0.0_C_DOUBLE
            RETURN
        END IF
        IF (num_params < 2) THEN ! Assuming we need at least 2 params for this stub
            PRINT *, "Fortran (Complicated_EOS): Error - not enough parameters (need at least 2)."
            istat = 302 ! Error: Not enough parameters
            P = 0.0_C_DOUBLE; E = 0.0_C_DOUBLE; dPdT = 0.0_C_DOUBLE
            dEdT = 0.0_C_DOUBLE; dPdrho = 0.0_C_DOUBLE
            RETURN
        END IF


        PRINT *, "Fortran (Complicated_EOS): Calling TFD_eos_compute..."
        CALL TFD_eos_compute(rho, T, matrix_A_in, matrix_B_in, N1_A, N2_A, N1_B, N2_B, &
                tfd_result_x, tfd_result_y, tfd_istat)

        IF (tfd_istat /= 0) THEN
            PRINT *, "Fortran (Complicated_EOS): Error from TFD_eos_compute. Istat=", tfd_istat
            istat = 310 + tfd_istat ! Differentiate error source
            P = 0.0_C_DOUBLE; E = 0.0_C_DOUBLE; dPdT = 0.0_C_DOUBLE
            dEdT = 0.0_C_DOUBLE; dPdrho = 0.0_C_DOUBLE
            RETURN
        END IF
        PRINT *, "Fortran (Complicated_EOS): TFD results: X=", tfd_result_x, " Y=", tfd_result_y

        ! Placeholder calculations:
        P = params_in(1) * rho + tfd_result_x
        E = params_in(2) * T   + tfd_result_y
        dPdT   = params_in(1) * 0.1_C_DOUBLE
        dEdT   = params_in(2) * 1.1_C_DOUBLE
        dPdrho = params_in(1)

        istat = 0 ! Success
        PRINT *, "Fortran (Complicated_EOS): Successfully computed P, E, etc."

    END SUBROUTINE Complicated_EOS

END MODULE module_complicated_eos