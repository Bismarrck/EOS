! src/fortran/module_eos_fortran_interface.f90
MODULE module_eos_fortran_interface
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: eos_initialize, eos_free, eos_set_control_variable, &
            eos_compute, eos_pack_get_size, eos_pack_to_buffer, eos_unpack, &
            eos_check_data_dir

    ! --- C Interface Definitions --- (No changes here, they are correct)
    INTERFACE
        FUNCTION c_eos_initialize_instance(data_dir_c, data_dir_len, eos_ids_c, num_eos_ids) &
                BIND(C, name='eos_c_initialize_instance') RESULT(istat)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*), INTENT(IN) :: data_dir_c
            INTEGER(C_INT), VALUE, INTENT(IN) :: data_dir_len
            INTEGER(C_INT), DIMENSION(*), INTENT(IN) :: eos_ids_c
            INTEGER(C_INT), VALUE, INTENT(IN) :: num_eos_ids
            INTEGER(C_INT) :: istat
        END FUNCTION c_eos_initialize_instance

        SUBROUTINE c_eos_free_instance() BIND(C, name='eos_c_free_instance')
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
        END SUBROUTINE c_eos_free_instance

        FUNCTION c_eos_set_bool_control(control_name_c, name_len, value) &
                BIND(C, name='eos_c_set_bool_control') RESULT(istat)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*), INTENT(IN) :: control_name_c
            INTEGER(C_INT), VALUE, INTENT(IN) :: name_len
            LOGICAL(C_BOOL), VALUE, INTENT(IN) :: value
            INTEGER(C_INT) :: istat
        END FUNCTION c_eos_set_bool_control

        FUNCTION c_eos_compute(eos_id, rho, T, P, E, dPdT, dEdT, dPdrho) &
                BIND(C, name='eos_c_compute') RESULT(istat)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), VALUE, INTENT(IN) :: eos_id
            REAL(C_DOUBLE), VALUE, INTENT(IN) :: rho, T
            REAL(C_DOUBLE), INTENT(OUT) :: P, E, dPdT, dEdT, dPdrho
            INTEGER(C_INT) :: istat
        END FUNCTION c_eos_compute

        FUNCTION c_eos_pack_data_get_size(required_buffer_size) &
                BIND(C, name='eos_c_pack_data_get_size') RESULT(istat)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(OUT) :: required_buffer_size
            INTEGER(C_INT) :: istat
        END FUNCTION c_eos_pack_data_get_size

        FUNCTION c_eos_pack_data_to_buffer(buffer_c, buffer_size_inout) &
                BIND(C, name='eos_c_pack_data_to_buffer') RESULT(istat)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*), INTENT(INOUT) :: buffer_c
            INTEGER(C_INT), INTENT(INOUT) :: buffer_size_inout
            INTEGER(C_INT) :: istat
        END FUNCTION c_eos_pack_data_to_buffer

        FUNCTION c_eos_unpack_data(buffer_c, buffer_size) &
                BIND(C, name='eos_c_unpack_data') RESULT(istat)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*), INTENT(IN) :: buffer_c
            INTEGER(C_INT), VALUE, INTENT(IN) :: buffer_size
            INTEGER(C_INT) :: istat
        END FUNCTION c_eos_unpack_data

        FUNCTION c_eos_check_eos_data_dir(data_dir_c, data_dir_len, eos_ids_c, num_ids_to_check) &
                BIND(C, name='eos_c_check_eos_data_dir') RESULT(istat)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*), INTENT(IN) :: data_dir_c
            INTEGER(C_INT), VALUE, INTENT(IN) :: data_dir_len
            INTEGER(C_INT), DIMENSION(*), INTENT(IN) :: eos_ids_c
            INTEGER(C_INT), VALUE, INTENT(IN) :: num_ids_to_check
            INTEGER(C_INT) :: istat
        END FUNCTION c_eos_check_eos_data_dir
    END INTERFACE

CONTAINS

    ! --- Fortran Wrapper Subroutines/Functions ---

    FUNCTION eos_initialize(data_dir, eos_ids) RESULT(istat)
        CHARACTER(LEN=*), INTENT(IN) :: data_dir
        INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: eos_ids
        INTEGER(C_INT) :: istat
        CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(MAX(1, LEN_TRIM(data_dir))) :: data_dir_c_arr
        INTEGER(C_INT), ALLOCATABLE, TARGET, DIMENSION(:) :: eos_ids_c_arr_alloc
        INTEGER(C_INT), POINTER, DIMENSION(:) :: eos_ids_c_arr_ptr => NULL() ! Use pointer for optional
        INTEGER(C_INT) :: num_eos_ids_val, local_size_eos_ids

        CALL string_to_c_char_array(data_dir, data_dir_c_arr)

        num_eos_ids_val = 0
        IF (PRESENT(eos_ids)) THEN
            local_size_eos_ids = SIZE(eos_ids)
            IF (local_size_eos_ids > 0) THEN
                ALLOCATE(eos_ids_c_arr_alloc(local_size_eos_ids))
                eos_ids_c_arr_alloc = INT(eos_ids(1:local_size_eos_ids), KIND=C_INT)
                eos_ids_c_arr_ptr => eos_ids_c_arr_alloc
                num_eos_ids_val = INT(local_size_eos_ids, KIND=C_INT)
            END IF
        END IF

        IF (ASSOCIATED(eos_ids_c_arr_ptr)) THEN
            istat = c_eos_initialize_instance(data_dir_c_arr, INT(LEN_TRIM(data_dir), C_INT), &
                    eos_ids_c_arr_ptr, num_eos_ids_val)
        ELSE
            ! Pass a conforming null pointer or a dummy if C expects non-null for zero size
            ! For C, if num_eos_ids_val is 0, eos_ids_c can be NULL.
            ! The C API checks for num_eos_ids > 0 before dereferencing eos_ids_arr.
            ! To be absolutely safe for older C standards or if C might still try to access,
            ! we could pass a dummy, but modern C handles NULL fine if size is 0.
            ! Let's pass a non-associated pointer which will be NULL on C side.
            ! Or, more explicitly:
            IF (num_eos_ids_val == 0) THEN
                ! Create a dummy zero-sized array to satisfy non-null pointer expectations if any.
                ! However, the c_eos_initialize_instance will check num_eos_ids.
                ! A non-associated pointer is better.
                ! The C API is robust to eos_ids_arr being NULL if num_eos_ids is 0.
                istat = c_eos_initialize_instance(data_dir_c_arr, INT(LEN_TRIM(data_dir), C_INT), &
                        eos_ids_c_arr_ptr, num_eos_ids_val) ! Pass non-associated pointer
            END IF
        END IF
        IF (ALLOCATED(eos_ids_c_arr_alloc)) DEALLOCATE(eos_ids_c_arr_alloc)

    END FUNCTION eos_initialize

    SUBROUTINE eos_free()
        CALL c_eos_free_instance()
    END SUBROUTINE eos_free

    FUNCTION eos_set_control_variable(control_name, value) RESULT(istat)
        CHARACTER(LEN=*), INTENT(IN) :: control_name
        LOGICAL, INTENT(IN) :: value
        INTEGER(C_INT) :: istat
        CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(MAX(1, LEN_TRIM(control_name))) :: control_name_c_arr

        CALL string_to_c_char_array(control_name, control_name_c_arr)
        istat = c_eos_set_bool_control(control_name_c_arr, INT(LEN_TRIM(control_name), C_INT), &
                LOGICAL(value, C_BOOL))
    END FUNCTION eos_set_control_variable

    FUNCTION eos_compute(eos_id, rho, T, P, E, dPdT, dEdT, dPdrho) RESULT(istat)
        INTEGER, INTENT(IN) :: eos_id
        REAL(C_DOUBLE), INTENT(IN) :: rho, T
        REAL(C_DOUBLE), INTENT(OUT) :: P, E, dPdT, dEdT, dPdrho
        INTEGER(C_INT) :: istat

        istat = c_eos_compute(INT(eos_id, C_INT), rho, T, P, E, dPdT, dEdT, dPdrho)
    END FUNCTION eos_compute

    FUNCTION eos_pack_get_size(required_buffer_size) RESULT(istat)
        INTEGER, INTENT(OUT) :: required_buffer_size
        INTEGER(C_INT) :: istat
        INTEGER(C_INT) :: required_buffer_size_c

        istat = c_eos_pack_data_get_size(required_buffer_size_c)
        required_buffer_size = INT(required_buffer_size_c)
    END FUNCTION eos_pack_get_size

    FUNCTION eos_pack_to_buffer(f_buffer, capacity_bytes, actual_bytes_written) RESULT(istat)
        CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(:), INTENT(INOUT) :: f_buffer
        INTEGER, INTENT(IN) :: capacity_bytes
        INTEGER, INTENT(OUT) :: actual_bytes_written
        INTEGER(C_INT) :: istat
        INTEGER(C_INT) :: buffer_size_c_inout

        IF (SIZE(f_buffer) < capacity_bytes) THEN
            WRITE(*,*) "Fortran (eos_pack_to_buffer) Error: f_buffer size (", SIZE(f_buffer), &
                    ") is less than capacity_bytes (", capacity_bytes, ")."
            istat = -101
            actual_bytes_written = 0
            RETURN
        END IF

        buffer_size_c_inout = INT(capacity_bytes, C_INT)
        istat = c_eos_pack_data_to_buffer(f_buffer(1:capacity_bytes), buffer_size_c_inout)
        actual_bytes_written = INT(buffer_size_c_inout)
    END FUNCTION eos_pack_to_buffer

    FUNCTION eos_unpack(f_buffer) RESULT(istat)
        CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(:), INTENT(IN) :: f_buffer
        INTEGER(C_INT) :: istat

        istat = c_eos_unpack_data(f_buffer, INT(SIZE(f_buffer), C_INT))
    END FUNCTION eos_unpack

    FUNCTION eos_check_data_dir(data_dir, eos_ids_to_check) RESULT(istat)
        CHARACTER(LEN=*), INTENT(IN) :: data_dir
        INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: eos_ids_to_check
        INTEGER(C_INT) :: istat
        CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(MAX(1, LEN_TRIM(data_dir))) :: data_dir_c_arr
        INTEGER(C_INT), ALLOCATABLE, TARGET, DIMENSION(:) :: eos_ids_c_arr_alloc
        INTEGER(C_INT), POINTER, DIMENSION(:) :: eos_ids_c_arr_ptr => NULL()
        INTEGER(C_INT) :: num_ids_val, local_size_eos_ids

        CALL string_to_c_char_array(data_dir, data_dir_c_arr)

        num_ids_val = 0
        IF (PRESENT(eos_ids_to_check)) THEN
            local_size_eos_ids = SIZE(eos_ids_to_check)
            IF (local_size_eos_ids > 0) THEN
                ALLOCATE(eos_ids_c_arr_alloc(local_size_eos_ids))
                eos_ids_c_arr_alloc = INT(eos_ids_to_check(1:local_size_eos_ids), KIND=C_INT)
                eos_ids_c_arr_ptr => eos_ids_c_arr_alloc
                num_ids_val = INT(local_size_eos_ids, KIND=C_INT)
            END IF
        END IF

        IF (ASSOCIATED(eos_ids_c_arr_ptr)) THEN
            istat = c_eos_check_eos_data_dir(data_dir_c_arr, INT(LEN_TRIM(data_dir), C_INT), &
                    eos_ids_c_arr_ptr, num_ids_val)
        ELSE
            istat = c_eos_check_eos_data_dir(data_dir_c_arr, INT(LEN_TRIM(data_dir), C_INT), &
                    eos_ids_c_arr_ptr, num_ids_val) ! Pass non-associated pointer
        END IF
        IF (ALLOCATED(eos_ids_c_arr_alloc)) DEALLOCATE(eos_ids_c_arr_alloc)

    END FUNCTION eos_check_data_dir

    ! Helper subroutine to convert Fortran string to C char array
    SUBROUTINE string_to_c_char_array(f_string, c_char_array)
        CHARACTER(LEN=*), INTENT(IN) :: f_string
        ! c_char_array is assumed-shape here, its size is taken from the actual argument.
        CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(:), INTENT(OUT) :: c_char_array
        INTEGER :: i, L_f_string, L_c_array

        L_f_string = LEN_TRIM(f_string)
        L_c_array = SIZE(c_char_array) ! Now this is valid as c_char_array is assumed-shape

        IF (L_c_array < L_f_string) THEN
            WRITE(*,*) "Error (string_to_c_char_array): Output C_CHAR array (size ", L_c_array, &
                    ") too small for Fortran string (trimmed length ", L_f_string, ")."
            ! Handle error: could stop, or truncate, or return status
            ! For now, let's fill what we can. A real error system would be better.
            DO i = 1, L_c_array
                c_char_array(i) = f_string(i:i)
            END DO
            RETURN
        END IF

        DO i = 1, L_f_string
            c_char_array(i) = f_string(i:i)
        END DO
        ! Fill rest with spaces or nulls if C API relies on it beyond explicit length.
        ! Our C API uses explicit length, so this is not strictly necessary.
        IF (L_f_string < L_c_array) THEN
            c_char_array(L_f_string+1 : L_c_array) = ' ' ! Or C_NULL_CHAR if appropriate
        END IF
    END SUBROUTINE string_to_c_char_array

END MODULE module_eos_fortran_interface