MODULE csv_xs_7g
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPLICIT NONE

    INTEGER, PARAMETER :: n_groups_7g = 7
    INTEGER, PARAMETER :: n_materials_7g = 4

    REAL(DP) :: xs_sigtr_7g(n_groups_7g, n_materials_7g)
    REAL(DP) :: xs_siga_7g(n_groups_7g, n_materials_7g)
    REAL(DP) :: xs_sigf_7g(n_groups_7g, n_materials_7g)
    REAL(DP) :: xs_nusigf_7g(n_groups_7g, n_materials_7g)
    REAL(DP) :: xs_chi_7g(n_groups_7g, n_materials_7g)
    REAL(DP) :: xs_sigs_7g(n_groups_7g, n_groups_7g, n_materials_7g)

CONTAINS

    SUBROUTINE load_7g_xs_from_csv(macro_path, scatter_path)
        CHARACTER(LEN=*), INTENT(IN) :: macro_path, scatter_path
        INTEGER :: m
        REAL(DP) :: s

        xs_sigtr_7g = 0.0_DP
        xs_siga_7g = 0.0_DP
        xs_sigf_7g = 0.0_DP
        xs_nusigf_7g = 0.0_DP
        xs_chi_7g = 0.0_DP
        xs_sigs_7g = 0.0_DP

        CALL read_macro_csv(macro_path)
        CALL read_scatter_csv(scatter_path)

        DO m = 1, n_materials_7g
            s = SUM(xs_chi_7g(:, m))
            IF (s > 0.0_DP) xs_chi_7g(:, m) = xs_chi_7g(:, m) / s
        END DO
    END SUBROUTINE load_7g_xs_from_csv

    SUBROUTINE read_macro_csv(path)
        CHARACTER(LEN=*), INTENT(IN) :: path
        INTEGER :: u, ios, g
        CHARACTER(LEN=512) :: line
        CHARACTER(LEN=64) :: fields(8), xs_type
        REAL(DP) :: uo2, mox, h2o, cr

        OPEN(NEWUNIT=u, FILE=TRIM(path), STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, "ERROR: cannot open ", TRIM(path)
            STOP 1
        END IF

        READ(u, '(A)', IOSTAT=ios) line ! header
        IF (ios /= 0) THEN
            PRINT *, "ERROR: empty file ", TRIM(path)
            STOP 1
        END IF

        DO
            READ(u, '(A)', IOSTAT=ios) line
            IF (ios /= 0) EXIT
            IF (LEN_TRIM(line) == 0) CYCLE

            CALL split_csv_simple(line, fields)
            xs_type = uppercase(TRIM(fields(1)))
            g = parse_group_token(TRIM(fields(2)))
            IF (g < 1 .OR. g > n_groups_7g) CYCLE

            READ(fields(3), *, IOSTAT=ios) uo2
            IF (ios /= 0) CYCLE
            READ(fields(4), *, IOSTAT=ios) mox
            IF (ios /= 0) CYCLE
            READ(fields(5), *, IOSTAT=ios) h2o
            IF (ios /= 0) CYCLE
            READ(fields(6), *, IOSTAT=ios) cr
            IF (ios /= 0) CYCLE

            SELECT CASE (xs_type)
            CASE ('SIGTR')
                xs_sigtr_7g(g, 1) = uo2
                xs_sigtr_7g(g, 2) = mox
                xs_sigtr_7g(g, 3) = h2o
                xs_sigtr_7g(g, 4) = cr
            CASE ('SIGA')
                xs_siga_7g(g, 1) = uo2
                xs_siga_7g(g, 2) = mox
                xs_siga_7g(g, 3) = h2o
                xs_siga_7g(g, 4) = cr
            CASE ('SIGF')
                xs_sigf_7g(g, 1) = uo2
                xs_sigf_7g(g, 2) = mox
                xs_sigf_7g(g, 3) = h2o
                xs_sigf_7g(g, 4) = cr
            CASE ('NUSIGF')
                xs_nusigf_7g(g, 1) = uo2
                xs_nusigf_7g(g, 2) = mox
                xs_nusigf_7g(g, 3) = h2o
                xs_nusigf_7g(g, 4) = cr
            CASE ('CHI')
                xs_chi_7g(g, 1) = uo2
                xs_chi_7g(g, 2) = mox
                xs_chi_7g(g, 3) = h2o
                xs_chi_7g(g, 4) = cr
            END SELECT
        END DO

        CLOSE(u)
    END SUBROUTINE read_macro_csv

    SUBROUTINE read_scatter_csv(path)
        CHARACTER(LEN=*), INTENT(IN) :: path
        INTEGER :: u, ios, g_from, g_to, mat_idx
        REAL(DP) :: val
        CHARACTER(LEN=512) :: line
        CHARACTER(LEN=64) :: fields(8)

        OPEN(NEWUNIT=u, FILE=TRIM(path), STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, "ERROR: cannot open ", TRIM(path)
            STOP 1
        END IF

        READ(u, '(A)', IOSTAT=ios) line ! header
        IF (ios /= 0) THEN
            PRINT *, "ERROR: empty file ", TRIM(path)
            STOP 1
        END IF

        DO
            READ(u, '(A)', IOSTAT=ios) line
            IF (ios /= 0) EXIT
            IF (LEN_TRIM(line) == 0) CYCLE

            CALL split_csv_simple(line, fields)
            mat_idx = material_index(TRIM(fields(1)))
            IF (mat_idx < 1 .OR. mat_idx > n_materials_7g) CYCLE

            READ(fields(2), *, IOSTAT=ios) g_from
            IF (ios /= 0) CYCLE
            READ(fields(3), *, IOSTAT=ios) g_to
            IF (ios /= 0) CYCLE
            READ(fields(4), *, IOSTAT=ios) val
            IF (ios /= 0) CYCLE

            IF (g_from >= 1 .AND. g_from <= n_groups_7g .AND. g_to >= 1 .AND. g_to <= n_groups_7g) THEN
                xs_sigs_7g(g_from, g_to, mat_idx) = val
            END IF
        END DO

        CLOSE(u)
    END SUBROUTINE read_scatter_csv

    INTEGER FUNCTION material_index(name)
        CHARACTER(LEN=*), INTENT(IN) :: name
        CHARACTER(LEN=32) :: up
        up = uppercase(TRIM(name))
        SELECT CASE (TRIM(up))
        CASE ('UO2')
            material_index = 1
        CASE ('MOX')
            material_index = 2
        CASE ('H2O')
            material_index = 3
        CASE ('CR')
            material_index = 4
        CASE DEFAULT
            material_index = -1
        END SELECT
    END FUNCTION material_index

    INTEGER FUNCTION parse_group_token(tok)
        CHARACTER(LEN=*), INTENT(IN) :: tok
        CHARACTER(LEN=32) :: t
        INTEGER :: ios
        parse_group_token = -1
        t = uppercase(TRIM(tok))
        IF (LEN_TRIM(t) >= 2 .AND. t(1:1) == 'G') THEN
            READ(t(2:LEN_TRIM(t)), *, IOSTAT=ios) parse_group_token
            IF (ios /= 0) parse_group_token = -1
        END IF
    END FUNCTION parse_group_token

    SUBROUTINE split_csv_simple(line, fields)
        CHARACTER(LEN=*), INTENT(IN) :: line
        CHARACTER(LEN=*), INTENT(OUT) :: fields(:)
        INTEGER :: i, start_pos, end_pos, n
        CHARACTER(LEN=LEN(line)) :: tmp

        DO i = 1, SIZE(fields)
            fields(i) = ''
        END DO

        tmp = TRIM(line)
        n = 1
        start_pos = 1

        DO i = 1, LEN_TRIM(tmp) + 1
            IF (i > LEN_TRIM(tmp) .OR. tmp(i:i) == ',') THEN
                end_pos = i - 1
                IF (n <= SIZE(fields)) THEN
                    IF (end_pos >= start_pos) THEN
                        fields(n) = ADJUSTL(tmp(start_pos:end_pos))
                    ELSE
                        fields(n) = ''
                    END IF
                END IF
                n = n + 1
                start_pos = i + 1
            END IF
        END DO
    END SUBROUTINE split_csv_simple

    FUNCTION uppercase(s) RESULT(out)
        CHARACTER(LEN=*), INTENT(IN) :: s
        CHARACTER(LEN=LEN(s)) :: out
        INTEGER :: i, c

        out = s
        DO i = 1, LEN(s)
            c = ICHAR(out(i:i))
            IF (c >= ICHAR('a') .AND. c <= ICHAR('z')) out(i:i) = ACHAR(c - 32)
        END DO
    END FUNCTION uppercase

END MODULE csv_xs_7g
