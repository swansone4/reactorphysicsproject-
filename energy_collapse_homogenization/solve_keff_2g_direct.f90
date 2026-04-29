PROGRAM solve_keff_2g_direct
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    USE fd_solvers_1d, ONLY: build_direct_inverse, apply_direct_inverse
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: in_xs = "energy_collapse_homogenization/outputs/two_group_homogenized_xs.csv"
    CHARACTER(LEN=*), PARAMETER :: out_csv = "energy_collapse_homogenization/outputs/two_group_keff_direct.csv"
    INTEGER, PARAMETER :: max_tok = 64
    CHARACTER(LEN=1024) :: line
    CHARACTER(LEN=128) :: tok(max_tok)
    INTEGER :: uin, uout, ios, ntok
    CHARACTER(LEN=64) :: case_id
    REAL(DP) :: k2g

    OPEN(NEWUNIT=uin, FILE=in_xs, STATUS="OLD", ACTION="READ", IOSTAT=ios)
    IF (ios /= 0) STOP 1
    OPEN(NEWUNIT=uout, FILE=out_csv, STATUS="REPLACE", ACTION="WRITE", IOSTAT=ios)
    IF (ios /= 0) STOP 1
    WRITE(uout,'(A)') "case_id,k_eff_2g_direct"

    READ(uin,'(A)',IOSTAT=ios) line ! header
    DO
        READ(uin,'(A)',IOSTAT=ios) line
        IF (ios /= 0) EXIT
        CALL split_csv(line, tok, ntok)
        IF (ntok < 14) CYCLE
        case_id = TRIM(tok(1))
        k2g = solve_case_2g(TRIM(case_id), tok)
        WRITE(uout,'(A,",",ES24.16)') TRIM(case_id), k2g
    END DO

    CLOSE(uin)
    CLOSE(uout)

CONTAINS
    SUBROUTINE split_csv(s, out, nout)
        CHARACTER(LEN=*), INTENT(IN) :: s
        CHARACTER(LEN=*), INTENT(OUT) :: out(:)
        INTEGER, INTENT(OUT) :: nout
        INTEGER :: i, p0, p1, n
        n = 0
        p0 = 1
        DO i = 1, LEN_TRIM(s) + 1
            IF (i > LEN_TRIM(s) .OR. s(i:i) == ",") THEN
                p1 = i - 1
                n = n + 1
                IF (n <= SIZE(out)) out(n) = ADJUSTL(s(p0:p1))
                p0 = i + 1
            END IF
        END DO
        nout = n
    END SUBROUTINE split_csv

    REAL(DP) FUNCTION solve_case_2g(cid, toks) RESULT(k)
        CHARACTER(LEN=*), INTENT(IN) :: cid
        CHARACTER(LEN=*), INTENT(IN) :: toks(:)
        CHARACTER(LEN=256) :: fpath
        INTEGER :: n, i, u, ios, nt
        REAL(DP), ALLOCATABLE :: dx(:), phi1(:), phi2(:), rhs1(:), rhs2(:), inv1(:,:), inv2(:,:)
        REAL(DP), ALLOCATABLE :: al1(:), ad1(:), au1(:), al2(:), ad2(:), au2(:)
        REAL(DP) :: k_old, fiss_old, fiss_new, sigtr1, sigtr2, siga1, siga2, nus1, nus2, s12, s21
        LOGICAL :: ok1, ok2
        CHARACTER(LEN=1024) :: line
        CHARACTER(LEN=128) :: tt(64)

        READ(toks(5), *) sigtr1
        READ(toks(6), *) sigtr2
        READ(toks(7), *) siga1
        READ(toks(8), *) siga2
        READ(toks(9), *) nus1
        READ(toks(10), *) nus2
        READ(toks(12), *) s12
        READ(toks(13), *) s21

        fpath = "energy_collapse_homogenization/outputs/" // TRIM(cid) // "/flux_cells.csv"
        OPEN(NEWUNIT=u, FILE=TRIM(fpath), STATUS="OLD", ACTION="READ", IOSTAT=ios)
        IF (ios /= 0) THEN
            k = -1.0_DP
            RETURN
        END IF
        READ(u,'(A)',IOSTAT=ios) line
        n = 0
        DO
            READ(u,'(A)',IOSTAT=ios) line
            IF (ios /= 0) EXIT
            n = n + 1
        END DO
        REWIND(u)
        READ(u,'(A)',IOSTAT=ios) line
        ALLOCATE(dx(n), phi1(n), phi2(n), rhs1(n), rhs2(n), inv1(n,n), inv2(n,n))
        ALLOCATE(al1(n), ad1(n), au1(n), al2(n), ad2(n), au2(n))

        i = 0
        DO
            READ(u,'(A)',IOSTAT=ios) line
            IF (ios /= 0) EXIT
            i = i + 1
            CALL split_csv(line, tt, nt)
            READ(tt(4),*) dx(i)
        END DO
        CLOSE(u)

        CALL assemble_tridiag(n, dx, 1.0_DP/(3.0_DP*MAX(sigtr1,1.0E-12_DP)), MAX(siga1+s12,1.0E-12_DP), al1, ad1, au1)
        CALL assemble_tridiag(n, dx, 1.0_DP/(3.0_DP*MAX(sigtr2,1.0E-12_DP)), MAX(siga2+s21,1.0E-12_DP), al2, ad2, au2)
        CALL build_direct_inverse(n, al1, ad1, au1, inv1, ok1)
        CALL build_direct_inverse(n, al2, ad2, au2, inv2, ok2)
        IF (.NOT. ok1 .OR. .NOT. ok2) THEN
            k = -1.0_DP
            RETURN
        END IF

        phi1 = 1.0_DP
        phi2 = 1.0_DP
        k = 1.0_DP
        DO i = 1, 500
            fiss_old = SUM(dx * (nus1 * phi1 + nus2 * phi2))
            rhs1 = dx * (s21 * phi2 + (1.0_DP / MAX(k,1.0E-20_DP)) * (nus1 * phi1 + nus2 * phi2))
            rhs2 = dx * (s12 * phi1)
            CALL apply_direct_inverse(n, inv1, rhs1, phi1)
            CALL apply_direct_inverse(n, inv2, rhs2, phi2)
            fiss_new = SUM(dx * (nus1 * phi1 + nus2 * phi2))
            k_old = k
            IF (fiss_old > 1.0E-30_DP .AND. fiss_new > 1.0E-30_DP) k = k * (fiss_new / fiss_old)
            IF (fiss_new > 1.0E-30_DP) THEN
                phi1 = phi1 / fiss_new
                phi2 = phi2 / fiss_new
            END IF
            IF (ABS(k - k_old) / MAX(ABS(k),1.0E-20_DP) < 1.0E-8_DP) EXIT
        END DO
    END FUNCTION solve_case_2g

    SUBROUTINE assemble_tridiag(n, dx, dcell, rrem, al, ad, au)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: dx(n), dcell, rrem
        REAL(DP), INTENT(OUT) :: al(n), ad(n), au(n)
        REAL(DP), ALLOCATABLE :: dD(:), L(:)
        INTEGER :: i
        ALLOCATE(dD(n))
        DO i = 1, n
            dD(i) = dcell / MAX(dx(i), 1.0E-20_DP)
        END DO
        IF (n > 1) ALLOCATE(L(n-1))
        IF (n > 1) THEN
            DO i = 1, n - 1
                L(i) = 2.0_DP * dD(i) * dD(i+1) / MAX(dD(i)+dD(i+1), 1.0E-20_DP)
            END DO
        END IF
        al = 0.0_DP; ad = 0.0_DP; au = 0.0_DP
        IF (n == 1) THEN
            ad(1) = dx(1) * rrem
            RETURN
        END IF
        ad(1) = L(1) + dx(1) * rrem
        au(1) = -L(1)
        DO i = 2, n - 1
            al(i) = -L(i-1)
            au(i) = -L(i)
            ad(i) = L(i-1) + L(i) + dx(i) * rrem
        END DO
        al(n) = -L(n-1)
        ad(n) = L(n-1) + dx(n) * rrem
    END SUBROUTINE assemble_tridiag
END PROGRAM solve_keff_2g_direct

