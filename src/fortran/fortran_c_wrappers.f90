! src/fortran/fortran_c_wrappers.f90
MODULE eos_c_wrappers
    USE, INTRINSIC :: ISO_C_BINDING
    USE module_analytic_eos       ! To access air_eos_2000, carbon_eos_2001
    USE module_global_controls    ! To access control variables
    USE module_tfd_eos            ! To access TFD
    USE module_complicated_eos, ONLY: Complicated_EOS
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

    ! --- Wrapper for TFD EOS compute ---
    SUBROUTINE c_tfd_eos_compute(rho_c, T_c, &
            matrix_a_flat_c, matrix_b_flat_c, &
            N1_A_c, N2_A_c, N1_B_c, N2_B_c, &
            tfd_result_x_c, tfd_result_y_c, &
            istat_out) &
            BIND(C, name='c_tfd_eos_compute')
        !DEC$ ATTRIBUTES DLLEXPORT :: c_tfd_eos_compute

        REAL(C_DOUBLE), VALUE, INTENT(IN) :: rho_c, T_c
        INTEGER(C_INT), VALUE, INTENT(IN) :: N1_A_c, N2_A_c, N1_B_c, N2_B_c

        ! C passes flat arrays. Fortran compute routine expects 2D arrays.
        ! We receive pointers to the start of the flat data.
        REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: matrix_a_flat_c
        REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: matrix_b_flat_c

        REAL(C_DOUBLE), INTENT(OUT) :: tfd_result_x_c, tfd_result_y_c
        INTEGER(C_INT), INTENT(OUT) :: istat_out

        ! --- Reshape flat C arrays to 2D Fortran arrays for the call ---
        ! This requires Fortran 2003+ pointer reshaping or careful indexing if using older Fortran.
        ! With modern Fortran and ISO_C_BINDING, passing the flat array and dimensions to
        ! TFD_eos_compute which declares it as DIMENSION(N1_A, N2_A) is the standard way.
        ! The actual "reshaping" happens by how Fortran interprets the memory based on
        ! the dimensions passed to the TFD_eos_compute subroutine.

        ! Pointers for reshaping (Fortran 2003+)
        ! REAL(C_DOUBLE), DIMENSION(:,:), POINTER :: matrix_A_ptr, matrix_B_ptr
        ! INTEGER, PARAMETER :: A_LBOUND(2) = [1,1]
        ! INTEGER :: A_UBOUND(2)
        ! INTEGER :: B_UBOUND(2)
        !
        ! A_UBOUND = [N1_A_c, N2_A_c]
        ! B_UBOUND = [N1_B_c, N2_B_c]
        !
        ! CALL C_F_POINTER(C_LOC(matrix_a_flat_c(1)), matrix_A_ptr, SHAPE=A_UBOUND) ! C_LOC of first element
        ! CALL C_F_POINTER(C_LOC(matrix_b_flat_c(1)), matrix_B_ptr, SHAPE=B_UBOUND)

        ! Check if dimensions are valid before calling
        IF (N1_A_c <=0 .OR. N2_A_c <=0 .OR. N1_B_c <=0 .OR. N2_B_c <=0) THEN
            istat_out = 202 ! Error code: invalid dimensions passed to C wrapper
            tfd_result_x_c = -1.0_C_DOUBLE
            tfd_result_y_c = -1.0_C_DOUBLE
            RETURN
        END IF

        ! Call the actual computation routine directly.
        ! TFD_eos_compute is designed to take array descriptors if it were an F77 style routine,
        ! or if it's F90+ and expects assumed-shape arrays, the compiler handles it if the interface is correct.
        ! Here, TFD_eos_compute declares matrix_A(N1_A_c, N2_A_c), so passing matrix_a_flat_c(1) (the start)
        ! and the dimensions should work, as Fortran will interpret the contiguous memory based on those dims.
        CALL TFD_eos_compute(rho_c, T_c, &
                matrix_a_flat_c(1:N1_A_c*N2_A_c), & ! Pass the relevant slice of the flat array
                matrix_b_flat_c(1:N1_B_c*N2_B_c), &
                N1_A_c, N2_A_c, N1_B_c, N2_B_c, &
                tfd_result_x_c, tfd_result_y_c, istat_out)
        ! Note: The above way of passing matrix_a_flat_c(1:N1_A_c*N2_A_c) creates a temporary array copy
        ! for an INTENT(IN) argument if the actual argument is not contiguous or simply-contiguous.
        ! A more direct way, if TFD_eos_compute expects a contiguous block starting at matrix_a_flat_c(1),
        ! is to just pass matrix_a_flat_c (or matrix_a_flat_c(1) if the callee expects a scalar start).
        ! The `DIMENSION(N1_A, N2_A)` in TFD_eos_compute should handle interpretation of the contiguous block.
        ! Let's simplify the call assuming TFD_eos_compute handles it:
        ! CALL TFD_eos_compute(rho_c, T_c, matrix_a_flat_c, matrix_b_flat_c, &
        !                      N1_A_c, N2_A_c, N1_B_c, N2_B_c, &
        !                      tfd_result_x_c, tfd_result_y_c, istat_out)
        ! The above might not work directly if TFD_eos_compute expects explicit shape.
        ! The issue is that matrix_a_flat_c is DIMENSION(*).
        !
        ! The most robust way for the BIND(C) wrapper to pass to a modern Fortran routine
        ! that expects assumed-shape arrays (like matrix_A(N1_A, N2_A)) is if that routine
        ! is also BIND(C) or if we use an intermediate procedure.
        !
        ! However, if TFD_eos_compute is defined as:
        ! SUBROUTINE TFD_eos_compute(rho, T, matrix_A_actual(N1_A, N2_A), ...)
        ! Then passing the start of the flat array should be fine, as Fortran knows how to
        ! interpret the contiguous memory using column-major order given the dimensions.

    END SUBROUTINE c_tfd_eos_compute

    SUBROUTINE c_complicated_eos(params_c, num_params_c, rho_c, T_c, &
            matrix_a_flat_c, N1_A_c, N2_A_c, &
            matrix_b_flat_c, N1_B_c, N2_B_c, &
            P_c, E_c, dPdT_c, dEdT_c, dPdrho_c, &
            istat_out) &
            BIND(C, name='c_complicated_eos')
        !DEC$ ATTRIBUTES DLLEXPORT :: c_complicated_eos

        ! Inputs
        INTEGER(C_INT), VALUE, INTENT(IN) :: num_params_c
        REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: params_c ! Flat array of params
        REAL(C_DOUBLE), VALUE, INTENT(IN) :: rho_c, T_c

        INTEGER(C_INT), VALUE, INTENT(IN) :: N1_A_c, N2_A_c, N1_B_c, N2_B_c
        REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: matrix_a_flat_c
        REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: matrix_b_flat_c

        ! Outputs
        REAL(C_DOUBLE), INTENT(OUT) :: P_c, E_c, dPdT_c, dEdT_c, dPdrho_c
        INTEGER(C_INT), INTENT(OUT) :: istat_out

        ! Validate inputs before calling the main routine
        IF (num_params_c <= 0 .OR. N1_A_c <= 0 .OR. N2_A_c <= 0 .OR. &
                N1_B_c <= 0 .OR. N2_B_c <= 0 ) THEN
            istat_out = 399 ! Invalid input to C wrapper
            P_c = 0.0_C_DOUBLE; E_c = 0.0_C_DOUBLE; dPdT_c = 0.0_C_DOUBLE
            dEdT_c = 0.0_C_DOUBLE; dPdrho_c = 0.0_C_DOUBLE
            RETURN
        END IF

        CALL Complicated_EOS(params_c(1:num_params_c), num_params_c, rho_c, T_c, &
                matrix_a_flat_c, N1_A_c, N2_A_c, & ! Pass flat C arrays directly
                matrix_b_flat_c, N1_B_c, N2_B_c, &
                P_c, E_c, dPdT_c, dEdT_c, dPdrho_c, &
                istat_out)
        ! The Complicated_EOS subroutine declares its array arguments (params_in, matrix_A_in, etc.)
        ! with explicit dimensions. Fortran will correctly map the contiguous data from
        ! params_c, matrix_a_flat_c, etc.

    END SUBROUTINE c_complicated_eos

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