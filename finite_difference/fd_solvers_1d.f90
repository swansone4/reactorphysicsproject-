! Tridiagonal 1D solvers (Lecture 19): direct matrix inversion + Jacobi / SR(Gauss-Seidel) / SOR.
MODULE fd_solvers_1d
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPLICIT NONE
    INTEGER, PARAMETER, PUBLIC :: solve_direct = 0, solve_jacobi = 1, solve_gs = 2, solve_sor = 3
    INTERFACE
        SUBROUTINE dgetrf(m, n, a, lda, ipiv, info)
            INTEGER, INTENT(IN) :: m, n, lda
            INTEGER, INTENT(OUT) :: ipiv(*), info
            DOUBLE PRECISION, INTENT(INOUT) :: a(lda,*)
        END SUBROUTINE dgetrf
        SUBROUTINE dgetri(n, a, lda, ipiv, work, lwork, info)
            INTEGER, INTENT(IN) :: n, lda, lwork
            INTEGER, INTENT(IN) :: ipiv(*)
            INTEGER, INTENT(OUT) :: info
            DOUBLE PRECISION, INTENT(INOUT) :: a(lda,*), work(*)
        END SUBROUTINE dgetri
    END INTERFACE
CONTAINS

    SUBROUTINE build_direct_inverse(n, al, ad, au, invA, ok)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: al(n), ad(n), au(n)
        REAL(DP), INTENT(OUT) :: invA(n,n)
        LOGICAL, INTENT(OUT) :: ok
        REAL(DP), ALLOCATABLE :: A(:,:), work(:)
        INTEGER, ALLOCATABLE :: ipiv(:)
        INTEGER :: i, info, lwork

        ok = .FALSE.
        IF (n <= 0) RETURN
        ALLOCATE(A(n,n), ipiv(n))
        A = 0.0_DP
        DO i = 1, n
            A(i,i) = ad(i)
            IF (i > 1) A(i,i-1) = al(i)
            IF (i < n) A(i,i+1) = au(i)
        END DO

        CALL dgetrf(n, n, A, n, ipiv, info)
        IF (info /= 0) THEN
            DEALLOCATE(A, ipiv)
            RETURN
        END IF
        lwork = MAX(1, 64*n)
        ALLOCATE(work(lwork))
        CALL dgetri(n, A, n, ipiv, work, lwork, info)
        DEALLOCATE(work)
        IF (info /= 0) THEN
            DEALLOCATE(A, ipiv)
            RETURN
        END IF

        invA = A
        ok = .TRUE.
        DEALLOCATE(A, ipiv)
    END SUBROUTINE build_direct_inverse

    SUBROUTINE apply_direct_inverse(n, invA, b, x)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: invA(n,n), b(n)
        REAL(DP), INTENT(OUT) :: x(n)
        x = MATMUL(invA, b)
    END SUBROUTINE apply_direct_inverse

    SUBROUTINE solve_tridiag_work(n, al0, ad0, au0, rh, x, ok)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: al0(n), ad0(n), au0(n), rh(n)
        REAL(DP), INTENT(OUT) :: x(n)
        LOGICAL, INTENT(OUT) :: ok
        REAL(DP), ALLOCATABLE :: invA(:,:)
        ok = n > 0
        IF (.NOT. ok) RETURN
        ALLOCATE(invA(n,n))
        CALL build_direct_inverse(n, al0, ad0, au0, invA, ok)
        IF (ok) CALL apply_direct_inverse(n, invA, rh, x)
        DEALLOCATE(invA)
    END SUBROUTINE solve_tridiag_work

    ! max norm of r = b - A x, tridiagonal with al(1)=0, au(n)=0
    PURE REAL(DP) FUNCTION tridiad_resid_inf(n, aal, aad, aau, bv, xv) RESULT(r)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: aal(n), aad(n), aau(n), bv(n), xv(n)
        INTEGER :: j1
        REAL(DP) :: t
        r = 0.0_DP
        DO j1 = 1, n
            t = aad(j1) * xv(j1)
            IF (j1 > 1) t = t + aal(j1) * xv(j1 - 1)
            IF (j1 < n) t = t + aau(j1) * xv(j1 + 1)
            t = t - bv(j1)
            r = MAX(r, ABS(t))
        END DO
    END FUNCTION tridiad_resid_inf

    SUBROUTINE iterative_solve( &
        method, n, al, ad, au, bvec, x, max_iter, tol, omega, it_out, res_out)
        INTEGER, INTENT(IN) :: method, n, max_iter
        REAL(DP), INTENT(IN) :: al(n), ad(n), au(n), bvec(n), tol, omega
        REAL(DP), INTENT(INOUT) :: x(n)
        INTEGER, INTENT(OUT) :: it_out
        REAL(DP), INTENT(OUT) :: res_out
        REAL(DP) :: xn(n), t, rres, rch
        INTEGER :: it, j

        it_out = 0
        res_out = 0.0_DP
        IF (n <= 0) RETURN

        DO it = 1, max_iter
            IF (method == solve_jacobi) THEN
                DO j = 1, n
                    t = bvec(j)
                    IF (j > 1) t = t - al(j) * x(j - 1)
                    IF (j < n) t = t - au(j) * x(j + 1)
                    xn(j) = t / ad(j)
                END DO
                x = xn
            ELSE IF (method == solve_gs) THEN
                rch = 0.0_DP
                DO j = 1, n
                    t = bvec(j)
                    IF (j > 1) t = t - al(j) * x(j - 1)
                    IF (j < n) t = t - au(j) * x(j + 1)
                    t = t / ad(j) - x(j)
                    rch = MAX(rch, ABS(t))
                    x(j) = t + x(j)
                END DO
            ELSE IF (method == solve_sor) THEN
                rch = 0.0_DP
                DO j = 1, n
                    t = bvec(j)
                    IF (j > 1) t = t - al(j) * x(j - 1)
                    IF (j < n) t = t - au(j) * x(j + 1)
                    t = t / ad(j) - x(j)
                    rch = MAX(rch, ABS(omega * t))
                    x(j) = x(j) + omega * t
                END DO
            END IF
            it_out = it
            rres = tridiad_resid_inf(n, al, ad, au, bvec, x)
            res_out = rres
            IF (rres < tol) RETURN
        END DO
    END SUBROUTINE iterative_solve

END MODULE fd_solvers_1d
