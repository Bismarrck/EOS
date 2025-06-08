PROGRAM test_eos_numerical_f
    USE module_eos_fortran_interface
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE
    IMPLICIT NONE
    REAL(C_DOUBLE) :: P, E, dPdT, dEdT, dPdrho
    CHARACTER(LEN=200) :: data_dir
    INTEGER :: istat
    REAL(C_DOUBLE) :: tol = 1.0E-7_C_DOUBLE

    data_dir = "../../../eos_data_dir" ! Adjust

    WRITE(*,*) "Fortran Numerical Validation"
    WRITE(*,*) "----------------------------"

    ! Test Case 1: Analytic Air (matches a row in CSV)
    CALL eos_free() ! Ensure clean state
    istat = eos_initialize(TRIM(data_dir), [INTEGER :: 1])
    IF (istat /= 0) THEN
        WRITE(*,*) "FAIL: Init for TC1"; STOP;
    END IF

    istat = eos_compute(1, 1.2_C_DOUBLE, 300.0_C_DOUBLE, P, E, dPdT, dEdT, dPdrho)
    IF (istat /= 0) THEN
        WRITE(*,*) "FAIL: Compute TC1"; STOP;
    END IF
    WRITE(*,*) "TC1 (Air): P=", P, " E=", E ! (Print all outputs)
    ! Compare with known ref_P=360.0, ref_E=540.0
    IF (ABS(P - 360.0_C_DOUBLE) > tol .OR. ABS(E - 540.0_C_DOUBLE) > tol) THEN
        WRITE(*,*) "FAIL: TC1 values mismatch"
    ELSE
        WRITE(*,*) "PASS: TC1"
    END IF

    ! Test Case 2: Complicated with TFDv1 (matches a row in CSV)
    CALL eos_free()
    istat = eos_set_control_variable("use_tfd_data_ver1", .TRUE.)
    istat = eos_initialize(TRIM(data_dir), [INTEGER :: 10000])
    IF (istat /= 0) THEN
        WRITE(*,*) "FAIL: Init for TC2"; STOP;
    END IF

    istat = eos_compute(10000, 1.8_C_DOUBLE, 1200.0_C_DOUBLE, P, E, dPdT, dEdT, dPdrho)
    IF (istat /= 0) THEN
        WRITE(*,*) "FAIL: Compute TC2"; STOP;
    END IF
    WRITE(*,*) "TC2 (Comp TFDv1): P=", P, " E=", E ! (Print all outputs)
    ! Compare with known ref_P=19.98, ref_E=36120.0
    IF (ABS(P - 19.98_C_DOUBLE) > tol .OR. ABS(E - 36120.0_C_DOUBLE) > tol) THEN
        WRITE(*,*) "FAIL: TC2 values mismatch"
    ELSE
        WRITE(*,*) "PASS: TC2"
    END IF

    CALL eos_free()
END PROGRAM test_eos_numerical_f
