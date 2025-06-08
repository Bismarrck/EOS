! examples/fortran_example/main.f90
PROGRAM fortran_eos_example
    USE module_eos_fortran_interface
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE
    IMPLICIT NONE

    CHARACTER(LEN=256) :: data_dir_arg ! Increased length for path
    CHARACTER(LEN=256) :: default_data_dir
    INTEGER :: num_args, arg_len
    LOGICAL :: data_dir_provided

    INTEGER, DIMENSION(1) :: eos_id_single_check(1)
    INTEGER, DIMENSION(3) :: eos_ids_to_load_arr(3)
    INTEGER :: istat
    INTEGER :: required_size, actual_packed_size
    CHARACTER(KIND=C_CHAR, LEN=1), ALLOCATABLE, DIMENSION(:) :: packed_data_f
    REAL(C_DOUBLE) :: rho, temp, P, E, dPdT, dEdT, dPdrho

    WRITE(*,*) "Fortran Example: Starting EOS Test"
    WRITE(*,*) "=================================="

    ! --- Get Data Directory from Command Line ---
    data_dir_provided = .FALSE.
    num_args = IARGC() ! Get number of command line arguments (gfortran extension)

    IF (num_args >= 1) THEN
        CALL GETARG(1, data_dir_arg) ! Get the first argument (gfortran extension)
        arg_len = len(data_dir_arg)
        IF (arg_len > 0 .AND. arg_len <= LEN(data_dir_arg)) THEN
            data_dir_arg = data_dir_arg(1:arg_len) ! Use actual length
            data_dir_provided = .TRUE.
            WRITE(*,*) "Data directory provided via command line: ", TRIM(data_dir_arg)
        ELSE
            WRITE(*,*) "Warning: Could not retrieve command line argument for data_dir or it was empty."
        END IF
    END IF

    IF (.NOT. data_dir_provided) THEN
        ! Default path if no command-line argument is given
        ! This default needs to be sensible for typical local execution.
        ! For CI, we will always pass the argument.
        default_data_dir = "../../../eos_data_dir" ! Example if exe is in build/examples/fortran_example/
        data_dir_arg = default_data_dir
        WRITE(*,*) "No data directory provided via command line. Using default: ", TRIM(data_dir_arg)
    END IF

    ! --- Initialize EOS Instance ---
    eos_ids_to_load_arr(1) = 1
    eos_ids_to_load_arr(2) = 2
    eos_ids_to_load_arr(3) = 10000
    istat = eos_initialize(TRIM(data_dir_arg), eos_ids_to_load_arr) ! Use TRIM
    WRITE(*,*) "eos_initialize status: ", istat
    IF (istat /= 0) THEN
        WRITE(*,*) "Initialization failed. Exiting."
        CALL eos_free()
        STOP 1
    END IF
    WRITE(*,*) "EOS Initialized using data_dir: ", TRIM(data_dir_arg)

    ! ... (rest of the Fortran example program: set control, compute, pack/unpack, free) ...
    ! Ensure TRIM(data_dir_arg) is used wherever the path is needed if it was used for logging.

    ! --- Set a Control Variable ---
    istat = eos_set_control_variable("complicated_eos_use_old_cold_term", .TRUE.)
    WRITE(*,*) "eos_set_control_variable (old_cold_term) status: ", istat
    istat = eos_set_control_variable("use_tfd_data_ver1", .TRUE.)
    WRITE(*,*) "eos_set_control_variable (tfd_ver1) status: ", istat

    ! --- Perform Computations ---
    WRITE(*,*) ""
    WRITE(*,*) "Performing Computations:"
    rho = 1.2_C_DOUBLE; temp = 300.0_C_DOUBLE
    istat = eos_compute(1, rho, temp, P, E, dPdT, dEdT, dPdrho)
    WRITE(*,'(A,I0,A,I0)') "  EOS ID 1 compute status: ", istat
    IF (istat == 0) THEN
        WRITE(*,'(A, F10.5, A, F12.5, A, F10.5, A, F10.5, A, F10.5)') &
                "    P=", P, " E=", E, " dPdT=", dPdT, " dEdT=", dEdT, " dPdrho=", dPdrho
    END IF

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

    CALL eos_free()
    WRITE(*,*) ""
    WRITE(*,*) "EOS Instance Freed."
    WRITE(*,*) "=================================="
    WRITE(*,*) "Fortran Example: Test Complete."

END PROGRAM fortran_eos_example