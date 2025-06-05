! examples/fortran_example/main.f90
PROGRAM fortran_eos_example
    USE module_eos_fortran_interface
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE
    IMPLICIT NONE

    CHARACTER(LEN=200) :: data_dir
    INTEGER, DIMENSION(1) :: eos_id_single_check(1) ! For checking a single ID
    INTEGER, DIMENSION(3) :: eos_ids_to_load_arr(3) ! For loading multiple IDs
    INTEGER :: istat, i
    INTEGER :: required_size, actual_packed_size

    ! Fortran allocatable buffer for pack/unpack
    CHARACTER(KIND=C_CHAR, LEN=1), ALLOCATABLE, DIMENSION(:) :: packed_data_f

    REAL(C_DOUBLE) :: rho, temp, P, E, dPdT, dEdT, dPdrho

    WRITE(*,*) "Fortran Example: Starting EOS Test"
    WRITE(*,*) "=================================="

    ! --- Define Data Directory (adjust path as needed relative to executable) ---
    data_dir = "../../../eos_data_dir" ! Example: if exe is in build/examples/fortran_example/
    WRITE(*,*) "Using data_dir: ", TRIM(data_dir)

    ! --- Initialize EOS Instance ---
    ! Correct Fortran array constructor or assignment
    eos_ids_to_load_arr(1) = 1
    eos_ids_to_load_arr(2) = 2
    eos_ids_to_load_arr(3) = 10000

    istat = eos_initialize(TRIM(data_dir), eos_ids_to_load_arr)
    WRITE(*,*) "eos_initialize status: ", istat
    IF (istat /= 0) THEN
        WRITE(*,*) "Initialization failed. Exiting."
        CALL eos_free()
        STOP 1
    END IF
    WRITE(*,*) "EOS Initialized."

    ! --- Check Data Directory ---
    eos_id_single_check(1) = 10000 ! Correct Fortran array assignment
    istat = eos_check_data_dir(data_dir, eos_id_single_check)
    WRITE(*,*) "eos_check_data_dir status (for ID 10000): ", istat
    IF (istat /= 0) THEN
        WRITE(*,*) "Check data dir failed. Exiting."
        STOP 1
    END IF

    ! --- Set a Control Variable ---
    istat = eos_set_control_variable("complicated_eos_use_old_cold_term", .TRUE.)
    WRITE(*,*) "eos_set_control_variable (old_cold_term) status: ", istat
    istat = eos_set_control_variable("use_tfd_data_ver1", .TRUE.)
    WRITE(*,*) "eos_set_control_variable (tfd_ver1) status: ", istat


    ! --- Perform Computations ---
    WRITE(*,*) ""
    WRITE(*,*) "Performing Computations:"
    ! Analytic Air (ID 1)
    rho = 1.2_C_DOUBLE; temp = 300.0_C_DOUBLE
    istat = eos_compute(1, rho, temp, P, E, dPdT, dEdT, dPdrho)
    WRITE(*,'(A,I0,A,I0)') "  EOS ID 1 compute status: ", istat
    IF (istat == 0) THEN
        WRITE(*,'(A, F10.5, A, F12.5, A, F10.5, A, F10.5, A, F10.5)') &
                "    P=", P, " E=", E, " dPdT=", dPdT, " dEdT=", dEdT, " dPdrho=", dPdrho
    END IF

    ! Complicated EOS (ID 10000)
    rho = 1.8_C_DOUBLE; temp = 1200.0_C_DOUBLE
    istat = eos_compute(10000, rho, temp, P, E, dPdT, dEdT, dPdrho)
    WRITE(*,'(A,I0,A,I0)') "  EOS ID 10000 compute status: ", istat
    IF (istat == 0) THEN
        WRITE(*,'(A, F12.5, A, F14.5, A, F10.5, A, F10.5, A, F10.5)') &
                "    P=", P, " E=", E, " dPdT=", dPdT, " dEdT=", dEdT, " dPdrho=", dPdrho
    END IF

    ! --- Pack Data ---
    WRITE(*,*) ""
    WRITE(*,*) "Packing Data:"
    istat = eos_pack_get_size(required_size)
    WRITE(*,*) "eos_pack_get_size status: ", istat, ", Required size: ", required_size
    IF (istat == 0 .AND. required_size > 0) THEN
        ALLOCATE(CHARACTER(KIND=C_CHAR, LEN=1) :: packed_data_f(required_size))
        istat = eos_pack_to_buffer(packed_data_f, required_size, actual_packed_size)
        WRITE(*,*) "eos_pack_to_buffer status: ", istat, ", Actual packed size: ", actual_packed_size
        IF (actual_packed_size /= required_size .AND. istat == 0) THEN
            WRITE(*,*) "Warning: actual packed size differs from required size."
        END IF

        ! --- Unpack Data (into a new conceptual instance, simulated by freeing and re-init from buffer) ---
        IF (istat == 0 .AND. actual_packed_size > 0) THEN
            WRITE(*,*) ""
            WRITE(*,*) "Unpacking Data:"
            CALL eos_free()
            istat = eos_unpack(packed_data_f(1:actual_packed_size))
            WRITE(*,*) "eos_unpack status: ", istat

            IF (istat == 0) THEN
                WRITE(*,*) "Unpack successful. Re-computing ID 10000 to verify:"
                istat = eos_compute(10000, rho, temp, P, E, dPdT, dEdT, dPdrho)
                WRITE(*,'(A,I0,A,I0)') "  EOS ID 10000 (post-unpack) compute status: ", istat
                IF (istat == 0) THEN
                    WRITE(*,'(A, F12.5, A, F14.5, A, F10.5, A, F10.5, A, F10.5)') &
                            "    P=", P, " E=", E, " dPdT=", dPdT, " dEdT=", dEdT, " dPdrho=", dPdrho
                END IF
            END IF
        END IF
        IF (ALLOCATED(packed_data_f)) DEALLOCATE(packed_data_f)
    END IF


    ! --- Free EOS Instance ---
    CALL eos_free()
    WRITE(*,*) ""
    WRITE(*,*) "EOS Instance Freed."
    WRITE(*,*) "=================================="
    WRITE(*,*) "Fortran Example: Test Complete."

END PROGRAM fortran_eos_example