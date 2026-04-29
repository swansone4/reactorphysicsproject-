! 7-group 1D multigroup diffusion eigenvalue: finite difference from D&H-style lectures
!   (3-point with harmonic d_{i,i+1}, albedo/reflector boundaries, power iteration for k).
!   CLI:  diffusion_solver_7g  [ direct | jacobi | gs | sor ]  [ omega for sor; default 1.5 ].
!   Direct = matrix inversion. Iterative methods: Jacobi / SR(Gauss-Seidel) / SOR.
PROGRAM diffusion_solver_7g
    ! 64-bit floats; standard library name for the kind
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    ! Reads input_file.txt, mesh, materials, 2G stuff we ignore, and 7G boundary flags
    USE read_input_file
    ! Fills x_left, x_right, cell volumes from the repeating lattice in the input
    USE mesh_generation
    ! Loads macroscopic 7G cross sections from the CSVs in seven_groups/
    USE csv_xs_7g
    ! Tridiagonal: direct matrix inversion, Jacobi, SR(Gauss-Seidel), SOR
    USE fd_solvers_1d
    IMPLICIT NONE
    ! Fixed number of energy groups from the XS module (must match the CSVs)
    INTEGER, PARAMETER :: ng = n_groups_7g

    ! set_num: which ConfigSets row to use; n_cells: total fine cells; h, it_k, it_in: loop indices; mth: which solver
    INTEGER :: set_num, n_cells, h, it_k, it_in, mth
    ! Per group: 0=albedo, 1=reflector on left and right; shape (ng) from the input
    INTEGER :: left_bc(ng), right_bc(ng)
    ! Return code from optional Python plotting call
    INTEGER :: plot_x
    ! SOR weight if mth = sor; tolerances: outer k, inner multigroup, inner linear tridiagonal
    REAL(DP) :: omega, tol_k, tol_in, tol_solve
    ! Counters: max power steps, max inner (scatter) iters, max Jacobi/SOR sweeps
    INTEGER :: k_max, in_max, solve_max, it_step

    ! Cell geometry: left/right edges, width, center; all length n_cells
    REAL(DP), ALLOCATABLE :: x_left(:), x_right(:), dx(:), x_center(:)
    ! Per cell, per group: total/sigma_a/nu*sigma_f/chi; sigs g_in -> g_out; d_cell = diff coef; R_rem = removal
    REAL(DP), ALLOCATABLE :: sigt(:,:), siga(:,:), nusigf(:,:), chi(:,:)
    ! Scatter matrix sigs(i, g_from, g_to) in code storage used in build_rhs
    REAL(DP), ALLOCATABLE :: sigs(:,:,:)
    ! Diffusion coefficient 1/(3*Sigma_tr) and diagonal removal (absorption + net outscatter) per D&H multigroup
    REAL(DP), ALLOCATABLE :: d_cell(:,:)
    REAL(DP), ALLOCATABLE :: R_rem(:,:)
    ! phi: flux; phi_oldk: buffer for the inner iteration convergence test; rhs_g: right-hand for one group
    REAL(DP), ALLOCATABLE :: phi(:,:), phi_oldk(:,:), rhs_g(:)
    REAL(DP), ALLOCATABLE :: phi_prev_g(:)
    ! Tridiagonal al, ad, au for one group
    REAL(DP), ALLOCATABLE :: al(:), ad(:), au(:), w_ad(:)
    REAL(DP), ALLOCATABLE :: inv_direct(:,:,:)
    LOGICAL :: inv_ready(ng)
    ! k_eff, k_old: k before/after this outer step; fiss_old, fiss_new: integrated fission; various diagnostics
    REAL(DP) :: k_eff, k_old, fiss_old, fiss_new, k_fiss_over_rem, rem_int, fiss_at_end, k_candidate
    ! t1: last linear residual (iterative); tridiag_worst: max t1 seen, printed once; err_in, err_k, f0
    REAL(DP) :: t1, tridiag_worst, err_in, err_k, f0
    ! Current energy group in inner loops; last-residual over groups after convergence
    INTEGER :: ghi
    ! After convergence: one pass over all groups, max of |b - A*phi|_inf, should be ~tol or roundoff for direct
    REAL(DP) :: max_res_tdi
    REAL(DP) :: power_target_w, energy_per_fission_j, fission_rate_target
    REAL(DP) :: fission_rate_model, flux_scale_factor
    LOGICAL :: ok, conv_k, power_norm_applied
    ! All CSV/figures land here (path relative to repo root)
    CHARACTER(LEN=*), PARAMETER :: out_dir = "finite_difference/outputs"
    CHARACTER(LEN=*), PARAMETER :: macro_xs = "seven_groups/macroscopic_cross_sections.csv"
    CHARACTER(LEN=*), PARAMETER :: scat_xs = "seven_groups/scattering_matrix.csv"
    ! Optional plots after the run; ignore failure if no Python
    CHARACTER(LEN=*), PARAMETER :: plot_cmd = "python finite_difference/plot_fd_flux_7g.py"

    ! mth: 0=direct matrix inversion, 1=Jacobi, 2=SR(Gauss-Seidel), 3=SOR; omega for SOR
    CALL get_cli(mth, omega)
    ! Outer (k) iterations. Power update can be slow; 1e-3 rel k often converges in O(1e3) iters, not 1e-5
    k_max = 20000
    ! Inner loop: repeated group sweeps until group fluxes stop moving (fission+scatter)
    in_max = 200
    ! Per group linear solve: Jacobi/SOR can need many sweeps; keep this high so we are not the limiting factor
    solve_max = 2000000
    ! Outer: |dk/k| stop (relative). 1e-3 = 0.1% in k. Tighter (1e-4) needs large k_max and more CPU.
    tol_k = 1.0E-3_DP
    ! Inner (multigroup) iteration: max change in any cell, any group, between sweeps
    tol_in = 1.0E-5_DP
    ! Stopping for Jacobi/GS/SOR on |b-Ax|_inf. 1e-8 is far too strict (millions of sweeps/group, looks like "never");
    ! 1e-4 is the usual "good enough" lecture-style residual and is fast. Thomas is not affected.
    tol_solve = 1.0E-4_DP
    ! k starts at 1.0; actual eigenvalue is whatever the power iteration converges to
    k_eff = 1.0_DP

    ! So CSV writes do not fail on missing directory
    CALL ensure_out_dir()
    ! Fills global module variables: config, mat IDs, etc.
    CALL read_input('input_file.txt')
    ! If your input uses more than one geometry set, this wires them up
    CALL init_configurations()
    ! Populates the xs_* arrays in csv_xs_7g
    CALL load_7g_xs_from_csv(macro_xs, scat_xs)

    IF (mth == solve_direct) THEN
        tol_k = 1.0E-2_DP
    END IF

    ! Input is 0-based "config" index; we need 1-based "set" index for the mesh routines
    set_num = config + 1
    ! Guard so we do not read past the configured number of sets
    IF (set_num < 1 .OR. set_num > configs) THEN
        PRINT *, "ERROR: config out of range. Config =", config, " Configs =", configs
        STOP 1
    END IF
    ! Pull reflect/vacuum flags for each of the 7 groups (from input / module)
    CALL build_7g_bc(left_bc, right_bc)

    ! Total number of fine mesh cells in this 1D layout (sum over segments, fuel vs not)
    n_cells = count_mesh_cells_7g(set_num)
    ! Nothing to solve on an empty domain
    IF (n_cells < 1) THEN
        PRINT *, "ERROR: zero mesh cells"
        STOP 1
    END IF

    ! One block so array shapes stay consistent
    ALLOCATE(x_left(n_cells), x_right(n_cells), dx(n_cells), x_center(n_cells))
    ALLOCATE(sigt(n_cells,ng), siga(n_cells,ng), nusigf(n_cells,ng), chi(n_cells,ng), sigs(n_cells,ng,ng))
    ALLOCATE(d_cell(n_cells,ng), R_rem(n_cells,ng), phi(n_cells,ng), phi_oldk(n_cells,ng), rhs_g(n_cells))
    ALLOCATE(phi_prev_g(n_cells))
    ALLOCATE(al(n_cells), ad(n_cells), au(n_cells), w_ad(n_cells))
    ALLOCATE(inv_direct(n_cells,n_cells,ng))

    ! Walk the lattice: for each input segment, add mpwr or mpfr cells, map MatID to CSV column, set XS on each cell
    CALL build_mesh_xs_7g(set_num, n_cells, x_left, x_right, dx, x_center, sigt, siga, nusigf, chi, sigs, R_rem, d_cell)

    ! Start from a positive guess so the first fission is nonzero; flat is enough for this k-eigen problem
    phi = 1.0_DP
    k_eff = 1.0_DP
    inv_ready = .FALSE.
    ! Total fission rate ~ integral nu*Sigma_f * phi; scale phi so fission = 1 (D&H-style normalization, removes arbitrary scale)
    f0 = fission_production(n_cells, nusigf, phi, dx)
    ! Avoid div by zero if all nu*sigf = 0 (degenerate case)
    IF (f0 > 1.0E-30_DP) THEN
        phi = phi / f0
    END IF

    ! Stays false if we never satisfy err_k < tol_k
    conv_k = .FALSE.
    ! ------------- outer loop: k-eigenvalue (power) iteration -----------------
    DO it_k = 1, k_max
        tridiag_worst = 0.0_DP
        ! Fission "source" as seen at this k before the inner (scatter) solve; used in k update ratio
        fiss_old = fission_production(n_cells, nusigf, phi, dx)
        k_old = k_eff
        ! We use phi_oldk to measure inner convergence (flux change sweep-to-sweep)
        phi_oldk = phi

        ! ---------- inner loop: fixed k, converge multigroup balance (fission+down/up scatter) ----------
        ! Start the inner as "infinite" error so we at least do one pass
        err_in = huge(1.0_DP)
        DO it_in = 1, in_max
            ! This will be max |phi-phi_old| this sweep
            err_in = 0.0_DP
            ! One full multigroup pass: g = 1..7
            DO ghi = 1, ng
                ! RHS: inscatter from other groups + (chi/k)*fission, times dx (finite-volume source)
                CALL build_rhs(ghi, n_cells, dx, nusigf, chi, sigs, phi, k_eff, rhs_g)
                ! A matrix: 3-point diffusion + R_rem*dx on the diagonal, Lecture 18 boundaries
                CALL assemble_tridiag( &
                    n_cells, d_cell(:, ghi), dx, R_rem(:, ghi), left_bc(ghi), right_bc(ghi), al, ad, au)
                ! Same (al,ad,au) and rhs_g for all methods; only the linear solver changes
                IF (mth == solve_direct) THEN
                    phi_prev_g = phi(:, ghi)
                    IF (.NOT. inv_ready(ghi)) THEN
                        CALL build_direct_inverse(n_cells, al, ad, au, inv_direct(:,:,ghi), ok)
                        IF (.NOT. ok) THEN
                            PRINT *, "ERROR: direct matrix inversion failed for group", ghi
                            STOP 1
                        END IF
                        inv_ready(ghi) = .TRUE.
                    END IF
                    CALL apply_direct_inverse(n_cells, inv_direct(:,:,ghi), rhs_g, phi(:, ghi))
                    phi(:, ghi) = 0.35_DP * phi(:, ghi) + 0.65_DP * phi_prev_g
                ELSE
                    it_step = 0
                    ! Exit residual in t1; no per-group prints here (1e-8+spam was unusable; see tridiag_worst at end)
                    IF (mth == solve_sor) THEN
                        CALL iterative_solve( &
                            solve_sor, n_cells, al, ad, au, rhs_g, phi(:, ghi), solve_max, tol_solve, &
                            omega, it_step, t1)
                    ELSE IF (mth == solve_gs) THEN
                        CALL iterative_solve( &
                            solve_gs, n_cells, al, ad, au, rhs_g, phi(:, ghi), solve_max, tol_solve, 1.0_DP, it_step, t1)
                    ELSE
                        CALL iterative_solve( &
                            solve_jacobi, n_cells, al, ad, au, rhs_g, phi(:, ghi), solve_max, tol_solve, 1.0_DP, it_step, t1)
                    END IF
                    tridiag_worst = MAX(tridiag_worst, t1)
                END IF
            END DO
            ! After a full pass, how much did the flux change vs start of this inner? (L∞ across cells and groups)
            DO h = 1, ng
                err_in = MAX(err_in, MAXVAL(ABS(phi(:, h) - phi_oldk(:, h))))
            END DO
            ! For the next pass, the old is now the current
            phi_oldk = phi
            ! If scatter+within-group is converged, exit inner; else repeat with updated phi
            IF (err_in < tol_in) EXIT
        END DO
        ! ---------- end inner; phi is now consistent for this k in each group (within tolerances) -----------

        ! New fission with updated phi; ratio updates k
        fiss_new = fission_production(n_cells, nusigf, phi, dx)
        IF (fiss_old > 1.0E-30_DP .AND. fiss_new > 1.0E-30_DP) THEN
            k_candidate = k_old * (fiss_new / fiss_old)
            IF (mth == solve_direct) THEN
                ! Direct inversion is aggressive on early outer iterations; damp k updates.
                k_candidate = k_old * MIN(MAX(k_candidate / k_old, 0.80_DP), 1.25_DP)
                k_eff = 0.65_DP * k_old + 0.35_DP * k_candidate
            ELSE
                k_eff = k_candidate
            END IF
            IF (k_eff <= 0.0_DP .OR. k_eff /= k_eff) k_eff = k_old
        END IF
        ! Renormalize phi each outer so the problem stays O(1) in source size (fission = 1 again)
        IF (fiss_new > 1.0E-30_DP .AND. fiss_new == fiss_new) THEN
            phi = phi / fiss_new
        END IF
        ! Relative k change; small means we can stop
        err_k = ABS((k_eff - k_old) / MAX(ABS(k_eff), 1.0E-20_DP))
        IF (err_k < tol_k) THEN
            conv_k = .TRUE.
            EXIT
        END IF
    END DO

    ! One line for iterative: worst |b-Ax| in the last k step (if any) — no per-group spam
    IF (mth /= solve_direct .AND. tridiag_worst > 0.0_DP) THEN
        WRITE (*, *) "Largest tridiagonal res (this last k step) =", tridiag_worst, "  (tol =", tol_solve, ")"
    END IF

    ! After the multigroup (inner) and k (outer) stop, b(phi) and phi are a fixed point only to tol_in.
    ! So max_g |b_g(phi)-A_g*phi| can still be O(tol_in) when you rebuild b from the final phi (not 1e-8).
    CALL final_tridiag_resid_max(max_res_tdi)

    ! Convert converged direct-solve shape to absolute magnitude using requested thermal power.
    power_target_w = 3565.0E6_DP
    energy_per_fission_j = 200.0E6_DP * 1.602176634E-19_DP
    CALL apply_power_normalization( &
        n_cells, dx, nusigf, phi, mth, power_target_w, energy_per_fission_j, &
        fission_rate_target, fission_rate_model, flux_scale_factor, power_norm_applied)

    ! P_fis and "removal" integral: used only for the text diagnostic, not the same as k (diffusion leakage missing in denom)
    fiss_at_end = fission_production(n_cells, nusigf, phi, dx)
    rem_int = removal_leakage_volume(n_cells, R_rem, phi, dx)
    IF (rem_int > 1.0E-30_DP) THEN
        k_fiss_over_rem = fiss_at_end / rem_int
    ELSE
        k_fiss_over_rem = 0.0_DP
    END IF

    ! One row per x with all groups; for plotting
    CALL write_flux_csv(n_cells, x_center, phi, out_dir)
    ! Denser per-cell + summary.txt, including the final tridiagonal residual
    CALL write_flux_cells_csv( &
        n_cells, x_left, x_right, dx, phi, sigt, R_rem, out_dir, k_eff, mth, omega, it_k, conv_k, &
        k_fiss_over_rem, fiss_at_end, rem_int, max_res_tdi, &
        power_target_w, energy_per_fission_j, fission_rate_target, fission_rate_model, flux_scale_factor, power_norm_applied)
    CALL execute_command_line(plot_cmd, WAIT=.TRUE., EXITSTAT=plot_x)
    IF (plot_x /= 0) WRITE (*, *) "Warning: plot command exit status =", plot_x

    WRITE (*, *) "Check |b(phi)-A*phi|_inf (final phi,k, all groups) =", max_res_tdi, "  (coupled; expect <<1 if tol_in tight)"
    PRINT *, "Finite-difference 7G diffusion done."
    PRINT *, "k_eff (power iter) =", k_eff, "  outer_its =", it_k, "  k_conv =", conv_k
    PRINT *, "k_diag P_fis/R_int =", k_fiss_over_rem, "  (R*phi*dx integrated; diffusion not in denom)"

CONTAINS
    ! Parse the optional command line: which inner tridiagonal method and, for SOR, omega
    SUBROUTINE get_cli(mout, omega)
        INTEGER, INTENT(OUT) :: mout
        REAL(DP), INTENT(OUT) :: omega
        INTEGER :: n
        CHARACTER(LEN=32) :: a1, a2
        ! how many args after the executable name
        n = command_argument_count()
        ! default: direct matrix inversion
        mout = solve_direct
        omega = 1.5_DP
        ! No args: stay with default, print a hint
        IF (n < 1) THEN
            PRINT *, "No CLI args: using direct matrix inversion solver. Optional: [direct|jacobi|gs|sor] [omega]"
            RETURN
        END IF
        ! First char only (case-insensitive in spirit) picks the method; loose matching for student typing
        CALL get_command_argument(1, a1)
        IF (a1(1:1) == "j" .OR. a1(1:1) == "J") mout = solve_jacobi
        IF (a1(1:1) == "g" .OR. a1(1:1) == "G") mout = solve_gs
        IF (a1(1:1) == "s" .OR. a1(1:1) == "S") mout = solve_sor
        ! second arg: SOR weight
        IF (n >= 2) THEN
            CALL get_command_argument(2, a2)
            READ (a2, *, ERR=1) omega
        END IF
1       CONTINUE
    END SUBROUTINE get_cli

    ! Create outputs/ if it does not exist, using a tiny one-liner Python
    SUBROUTINE ensure_out_dir()
        CALL execute_command_line( &
            "python -c ""from pathlib import Path; Path('finite_difference/outputs').mkdir(parents=True, exist_ok=True)""", &
            WAIT=.TRUE., EXITSTAT=plot_x)
    END SUBROUTINE ensure_out_dir

    SUBROUTINE apply_power_normalization( &
        n, dxc, nus, ph, method, p_target_w, e_fis_j, f_target, f_model, scale, did_apply)
        INTEGER, INTENT(IN) :: n, method
        REAL(DP), INTENT(IN) :: dxc(n), nus(n,ng), p_target_w, e_fis_j
        REAL(DP), INTENT(INOUT) :: ph(n,ng)
        REAL(DP), INTENT(OUT) :: f_target, f_model, scale
        LOGICAL, INTENT(OUT) :: did_apply

        did_apply = .FALSE.
        f_target = 0.0_DP
        f_model = 0.0_DP
        scale = 1.0_DP
        IF (method /= solve_direct) RETURN
        IF (e_fis_j <= 0.0_DP .OR. p_target_w <= 0.0_DP) RETURN

        f_target = p_target_w / e_fis_j
        f_model = fission_production(n, nus, ph, dxc)
        IF (f_model <= 1.0E-30_DP) RETURN
        scale = f_target / f_model
        ph = ph * scale
        did_apply = .TRUE.
    END SUBROUTINE apply_power_normalization

    ! R_rem * phi * dx volume integral (1D "per unit area"); not the Laplace leakage, just removal term
    PURE FUNCTION removal_leakage_volume(n, rrm, ph, dxc) RESULT(r)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: rrm(n, ng), ph(n, ng), dxc(n)
        REAL(DP) :: r
        INTEGER :: i, g1
        r = 0.0_DP
        DO i = 1, n
            DO g1 = 1, ng
                r = r + rrm(i, g1) * ph(i, g1) * dxc(i)
            END DO
        END DO
    END FUNCTION removal_leakage_volume

    ! Fission "production" = sum over cells of dx * sum_h nu*Sigma_f_h*phi_h (all groups contribute)
    PURE FUNCTION fission_production(n, nus, ph, dxc) RESULT(f)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: nus(n,ng), ph(n,ng), dxc(n)
        REAL(DP) :: f, tsum
        INTEGER :: i, h
        f = 0.0_DP
        DO i = 1, n
            tsum = 0.0_DP
            ! Sum fission in all groups at this cell
            DO h = 1, ng
                tsum = tsum + nus(i, h) * ph(i, h)
            END DO
            f = f + tsum * dxc(i)
        END DO
    END FUNCTION fission_production

    ! Copy 7G boundary type arrays from the global module into left_bc/right_bc as integers
    SUBROUTINE build_7g_bc(lb, rb)
        INTEGER, INTENT(OUT) :: lb(ng), rb(ng)
        INTEGER :: g1
        ! Should have been set when input was read
        IF (.NOT. ALLOCATED(bound_l_7g) .OR. .NOT. ALLOCATED(bound_r_7g)) STOP 1
        DO g1 = 1, ng
            lb(g1) = NINT(bound_l_7g(g1))
            rb(g1) = NINT(bound_r_7g(g1))
        END DO
    END SUBROUTINE build_7g_bc

    ! Count fine cells: each input region adds mpwr cells (MatID=2) or mpfr (otherwise)
    INTEGER FUNCTION count_mesh_cells_7g(seti) RESULT(nc)
        INTEGER, INTENT(IN) :: seti
        INTEGER :: p, m
        nc = 0
        DO p = 1, set_len(seti)
            m = MatID(seti, p)
            IF (m == 2) THEN
                nc = nc + mpwr
            ELSE
                nc = nc + mpfr
            END IF
        END DO
    END FUNCTION count_mesh_cells_7g

    ! Fill cell edges, width, D, and removal, using MatID+1 to pick CSV material column: 1=UO2,2=MOX,3=H2O,4=CR
    SUBROUTINE build_mesh_xs_7g( &
        seti, n, xl, xr, dxc, xc, st, sa, nus, khi, sgs, rrm, dcel)
        INTEGER, INTENT(IN) :: seti, n
        REAL(DP), INTENT(OUT) :: xl(n), xr(n), dxc(n), xc(n)
        REAL(DP), INTENT(OUT) :: st(n,ng), sa(n,ng), nus(n,ng), khi(n,ng), sgs(n,ng,ng), rrm(n,ng), dcel(n,ng)
        ! pos: segment in the user lattice; cell_i: running global index; m: CSV column id (mid+1)
        INTEGER :: pos, mesh_i, cell_i, mid, m, g1, g2
        ! runx: left edge of current cell; wgeom, lasm, scx: geometry/segment weights; t_out: outscatter sum; pitch/diam from input
        REAL(DP) :: runx, wgeom, wsum, lasm, scx, pitch_dp, diam_dp, chs, t_out
        ! n_pin: count of pin-diameter style segments; npl: how many subcells this segment
        INTEGER :: n_pin, npl

        pitch_dp = REAL(rod_pitch, DP)
        diam_dp = REAL(rod_diameter, DP)
        n_pin = 0
        ! How many 0,1,3 (pin) type segments, for the lattice "length" scaling
        DO pos = 1, set_len(seti)
            mid = MatID(seti, pos)
            IF (mid == 0 .OR. mid == 1 .OR. mid == 3) n_pin = n_pin + 1
        END DO
        wsum = 0.0_DP
        ! Total "weight" of all segw pieces (so we can normalize lengths)
        DO pos = 1, set_len(seti)
            wsum = wsum + segw(MatID(seti, pos), pitch_dp, diam_dp)
        END DO
        ! lasm: physical length of one repeating "assembly" in cm for this problem (how the notes pin pitch works)
        IF (n_pin > 0) THEN
            lasm = REAL(n_pin, DP) * pitch_dp
        ELSE
            lasm = MAX(wsum, pitch_dp)
        END IF
        ! Scale: map sum of segw to lasm, so the cell widths in cm are consistent
        IF (wsum > 1.0E-15_DP) THEN
            scx = lasm / wsum
        ELSE
            scx = 1.0_DP
        END IF

        cell_i = 0
        runx = 0.0_DP
        ! For each user segment, subdivide into npl fine cells, assign that material to each subcell
        DO pos = 1, set_len(seti)
            mid = MatID(seti, pos)
            ! Fuel strip vs not uses different subcell count from input
            IF (mid == 2) THEN
                npl = mpwr
            ELSE
                npl = mpfr
            END IF
            ! width of the segment in cm, after scaling; split equally among the npl subcells
            wgeom = segw(mid, pitch_dp, diam_dp) * scx
            m = mid + 1
            DO mesh_i = 1, npl
                cell_i = cell_i + 1
                dxc(cell_i) = wgeom / REAL(npl, DP)
                xl(cell_i) = runx
                xr(cell_i) = runx + dxc(cell_i)
                ! Cell center: midpoint, where we will tabulate the flux
                xc(cell_i) = 0.5_DP * (xl(cell_i) + xr(cell_i))
                runx = xr(cell_i)
                ! chs will normalize chi to sum 1
                chs = 0.0_DP
                ! Pull 7G constants from the CSV arrays for this material
                DO g1 = 1, ng
                    st(cell_i, g1) = MAX(xs_sigtr_7g(g1, m), 1.0E-15_DP)
                    sa(cell_i, g1) = MAX(xs_siga_7g(g1, m), 0.0_DP)
                    nus(cell_i, g1) = MAX(xs_nusigf_7g(g1, m), 0.0_DP)
                    khi(cell_i, g1) = MAX(xs_chi_7g(g1, m), 0.0_DP)
                    chs = chs + khi(cell_i, g1)
                    DO g2 = 1, ng
                        sgs(cell_i, g1, g2) = MAX(xs_sigs_7g(g1, g2, m), 0.0_DP)
                    END DO
                END DO
                ! Fission spectrum must be a true pdf over groups
                IF (chs > 0.0_DP) khi(cell_i, :) = khi(cell_i, :) / chs
                ! rrm: diagonal "removal" in the multigroup sense (see comment below)
                DO g1 = 1, ng
                    t_out = 0.0_DP
                    DO g2 = 1, ng
                        IF (g2 /= g1) t_out = t_out + sgs(cell_i, g1, g2)
                    END DO
                    ! (Sigma_t - Sigma_s,gg) = Siga + sum of outscatter; small floor keeps D and A well-behaved
                    rrm(cell_i, g1) = MAX(sa(cell_i, g1) + t_out, 1.0E-8_DP * st(cell_i, g1))
                    ! P1 diffusion: D = 1/(3 Sigma_tr)
                    dcel(cell_i, g1) = 1.0_DP / (3.0_DP * st(cell_i, g1))
                END DO
            END DO
        END DO
    END SUBROUTINE build_mesh_xs_7g

    ! For the assembly pitch layout: "gap" (fuel) vs "pin" width, used in the 1D mapping to cm
    PURE REAL(DP) FUNCTION segw(aid, pitch, diade)
        INTEGER, INTENT(IN) :: aid
        REAL(DP), INTENT(IN) :: pitch, diade
        IF (aid == 2) THEN
            segw = MAX(pitch - diade, 0.0_DP)
        ELSE
            segw = diade
        END IF
    END FUNCTION segw

    ! Finite-volume right-hand: dx * [ inscatter + chi_g * (1/k) * fission at all h ]
    SUBROUTINE build_rhs( &
        g1, n, dxc, nus, khi, sgs, ph, kcur, rvec)
        INTEGER, INTENT(IN) :: g1, n
        REAL(DP), INTENT(IN) :: dxc(n), nus(n,ng), khi(n,ng), sgs(n,ng,ng), ph(n,ng)
        ! Current k in the fission source (outer iteration)
        REAL(DP), INTENT(IN) :: kcur
        REAL(DP), INTENT(OUT) :: rvec(n)
        INTEGER :: i, h1
        ! s_in: inscatter; s_f: total fission; kinv: 1/k
        REAL(DP) :: s_in, s_f, kinv
        ! 1/k from the eigenvalue form
        kinv = 1.0_DP / MAX(kcur, 1.0E-20_DP)
        DO i = 1, n
            s_in = 0.0_DP
            ! h1 /= g1: inscatter from other groups into g1
            DO h1 = 1, ng
                IF (h1 == g1) CYCLE
                s_in = s_in + sgs(i, h1, g1) * ph(i, h1)
            END DO
            s_f = 0.0_DP
            ! Total fission in all groups at (i) times chi in g1 / k
            DO h1 = 1, ng
                s_f = s_f + nus(i, h1) * ph(i, h1)
            END DO
            s_f = khi(i, g1) * kinv * s_f
            ! Volume weight dx: consistent with the Laplacian discretization
            rvec(i) = dxc(i) * (s_in + s_f)
        END DO
    END SUBROUTINE build_rhs

    ! Harmonic face coupling for different D, different dx: lecture d_{i,i+1} formula
    PURE REAL(DP) FUNCTION d_couple(d1, d2) RESULT(c)
        REAL(DP), INTENT(IN) :: d1, d2
        c = 2.0_DP * d1 * d2 / MAX(d1 + d2, 1.0E-20_DP)
    END FUNCTION d_couple

    ! Lecture 18-3: beta( alpha ), with alpha=0 (vac) or 1 (ref); dedge = D/dx at the boundary face
    PURE REAL(DP) FUNCTION albedo_beta(ibound, dedge) RESULT(bta)
        INTEGER, INTENT(IN) :: ibound
        REAL(DP), INTENT(IN) :: dedge
        REAL(DP) :: afa, t
        ! Degenerate mesh
        IF (dedge < 1.0E-20_DP) THEN
            bta = 0.0_DP
            RETURN
        END IF
        afa = 0.0_DP
        IF (ibound == 1) afa = 1.0_DP
        ! (1- alpha)/(1+alpha) is 1 for vac, 0 for reflect; then 1/4D term as in the slide
        t = (1.0_DP - afa) / MAX(1.0_DP + afa, 1.0E-20_DP) * 0.25_DP / dedge
        bta = 1.0_DP / (1.0_DP + t)
    END FUNCTION albedo_beta

    ! Assemble a_l, a_d, a_u for the global unknowns phi_1..phi_n for one group. Units are consistent with RHS* ~ dx*source
    SUBROUTINE assemble_tridiag( &
        n, dcel, dxc, R_rem, lbc, rbc, aall, aadd, aauu)
        INTEGER, INTENT(IN) :: n, lbc, rbc
        REAL(DP), INTENT(IN) :: dcel(n), dxc(n), R_rem(n)
        REAL(DP), DIMENSION(n), INTENT(OUT) :: aall, aadd, aauu
        INTEGER :: k1
        ! dD: Dg/dx in each cell; L_if: face coupling between k1 and k1+1; bL,bR: albedo betas; remk: dx*R
        REAL(DP) :: dD(n), L_if(n-1), bL, bR, remk, tlo, thi, L12, Lnm

        ! d_i = D_i / Delta x_i
        DO k1 = 1, n
            dD(k1) = dcel(k1) / MAX(dxc(k1), 1.0E-20_DP)
        END DO
        ! Interior face links d_{i,i+1} = 2 d_i d_{i+1} / (d_i + d_{i+1})
        DO k1 = 1, n - 1
            L_if(k1) = d_couple(dD(k1), dD(k1 + 1))
        END DO
        bL = albedo_beta(lbc, dD(1))
        bR = albedo_beta(rbc, dD(n))

        aall = 0.0_DP
        aauu = 0.0_DP
        aadd = 0.0_DP
        aall(1) = 0.0_DP
        ! Single cell case: two boundary terms + removal, no off-diagonal neighbor
        IF (n == 1) THEN
            aauu(1) = 0.0_DP
            aadd(1) = 2.0_DP * dD(1) * (1.0_DP - bL) + 2.0_DP * dD(1) * (1.0_DP - bR) + dxc(1) * R_rem(1)
            RETURN
        END IF
        L12 = L_if(1)
        Lnm = L_if(n - 1)
        ! k = 1 row: albedo (1-beta) term, link to 2, plus removal*dx
        remk = dxc(1) * R_rem(1)
        aadd(1) = 2.0_DP * dD(1) * (1.0_DP - bL) + L12 + remk
        aauu(1) = -L12
        ! k = n row: no upper neighbor; link from n-1, right boundary, removal
        remk = dxc(n) * R_rem(n)
        aall(n) = -Lnm
        aauu(n) = 0.0_DP
        aadd(n) = Lnm + 2.0_DP * dD(n) * (1.0_DP - bR) + remk
        ! Interior: standard 3-point + dx*R on the diagonal
        DO k1 = 2, n - 1
            tlo = L_if(k1 - 1)
            thi = L_if(k1)
            remk = dxc(k1) * R_rem(k1)
            aall(k1) = -tlo
            aauu(k1) = -thi
            aadd(k1) = tlo + thi + remk
        END DO
    END SUBROUTINE assemble_tridiag

    ! Rebuild each group's A and b from final phi, k, and take |b - A*phi|_inf — sanity check for any solver
    SUBROUTINE final_tridiag_resid_max(rout)
        REAL(DP), INTENT(OUT) :: rout
        INTEGER :: g1
        rout = 0.0_DP
        DO g1 = 1, ng
            CALL build_rhs(g1, n_cells, dx, nusigf, chi, sigs, phi, k_eff, rhs_g)
            CALL assemble_tridiag( &
                n_cells, d_cell(:, g1), dx, R_rem(:, g1), left_bc(g1), right_bc(g1), al, ad, au)
            rout = MAX(rout, tridiad_resid_inf(n_cells, al, ad, au, rhs_g, phi(:, g1)))
        END DO
    END SUBROUTINE final_tridiag_resid_max

    ! CSV: one row per cell center, columns are group fluxes, for the supplied Python plotter
    SUBROUTINE write_flux_csv(n, xc, fmc, odir)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: xc(n), fmc(n, ng)
        CHARACTER(LEN=*), INTENT(IN) :: odir
        INTEGER :: u, ios, j
        CHARACTER(LEN=256) :: fn
        fn = TRIM(odir) // '/flux_fundamental.csv'
        OPEN(NEWUNIT=u, FILE=TRIM(fn), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, "ERROR: could not open ", TRIM(fn)
            RETURN
        END IF
        WRITE (u, '(A)') "cell,x_cm,g1,g2,g3,g4,g5,g6,g7"
        DO j = 1, n
            WRITE (u, '(*(G0, :, ","))') j, xc(j), fmc(j,1), fmc(j,2), fmc(j,3), fmc(j,4), fmc(j,5), fmc(j,6), fmc(j,7)
        END DO
        CLOSE(u)
    END SUBROUTINE write_flux_csv

    ! Cell-by-cell flux+XS, plus human-readable run summary; tridiag_res: max |b-A*phi|_inf for final (k,phi)
    SUBROUTINE write_flux_cells_csv( &
        n, xl, xr, dxc, fmc, st, srem, odir, keff, mth, om, oit, didconv, kfr, fend, rint, tridiag_res, &
        p_target_w, e_fis_j, f_target, f_model, f_scale, p_applied)
        INTEGER, INTENT(IN) :: n, mth, oit
        REAL(DP), INTENT(IN) :: xl(n), xr(n), dxc(n), fmc(n,ng), st(n,ng), srem(n,ng)
        REAL(DP), INTENT(IN) :: keff, om, kfr, fend, rint, tridiag_res
        REAL(DP), INTENT(IN) :: p_target_w, e_fis_j, f_target, f_model, f_scale
        CHARACTER(LEN=*), INTENT(IN) :: odir
        LOGICAL, INTENT(IN) :: didconv, p_applied
        INTEGER :: u, ios, j
        CHARACTER(LEN=256) :: fsum
        fsum = TRIM(odir) // "/flux_cells.csv"
        OPEN(NEWUNIT=u, FILE=TRIM(fsum), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
        IF (ios == 0) THEN
            WRITE (u, '(A)') "cell,x_left_cm,x_right_cm,dx_cm,flux_g1,flux_g2,flux_g3,flux_g4,flux_g5,flux_g6,flux_g7," // &
                "sigTR_g1,sigTR_g2,sigTR_g3,sigTR_g4,sigTR_g5,sigTR_g6,sigTR_g7," // &
                "Rrem_g1,Rrem_g2,Rrem_g3,Rrem_g4,Rrem_g5,Rrem_g6,Rrem_g7"
            DO j = 1, n
                WRITE (u, '(*(G0,:,","))') j, xl(j), xr(j), dxc(j), &
                    fmc(j,1),fmc(j,2),fmc(j,3),fmc(j,4),fmc(j,5),fmc(j,6),fmc(j,7), &
                    st(j,1),st(j,2),st(j,3),st(j,4),st(j,5),st(j,6),st(j,7), &
                    srem(j,1),srem(j,2),srem(j,3),srem(j,4),srem(j,5),srem(j,6),srem(j,7)
            END DO
            CLOSE(u)
        END IF
        fsum = TRIM(odir) // "/summary.txt"
        OPEN(NEWUNIT=u, FILE=TRIM(fsum), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
        IF (ios == 0) THEN
            WRITE (u, '(A,I0)') "config_set_0based = ", config
            WRITE (u, '(A,ES20.10)') "k_eff_fd = ", keff
            WRITE (u, '(A,L1)') "k_converged = ", didconv
            WRITE (u, '(A,I0)') "outer_iterations = ", oit
            WRITE (u, '(A,I0)') "tridiag_solver = ", mth
            WRITE (u, '(A,ES20.10)') "sor_omega = ", om
            WRITE (u, '(A)') "method: 0=direct matrix inversion 1=Jacobi 2=SR(Gauss-Seidel) 3=SOR"
            WRITE (u, '(A,I0)') "total_mesh_cells = ", n
            WRITE (u, '(A)') "fd_formulation: 3-point d_{i,i+1} + removal + Lecture18 albedo/reflect on boundaries; power iter k"
            ! One-group A*phi vs b(phi) with final k,phi (fixed-point slop ~ tol_in, not tol_solve)
            WRITE (u, '(A,ES20.10)') "max_inf_norm_b_minus_Aphi_final_phi = ", tridiag_res
            WRITE (u, '(A,ES20.10)') "tol_solve = ", tol_solve
            WRITE (u, '(A,ES20.10)') "tol_k = ", tol_k
            WRITE (u, '(A,ES20.10)') "tol_in = ", tol_in
            WRITE (u, '(A,ES20.10)') "k_diag_fission_over_removal_Rphi = ", kfr
            WRITE (u, '(A,ES20.10)') "integral_P_fis = ", fend
            WRITE (u, '(A,ES20.10)') "integral_R_phi_dx = ", rint
            WRITE (u, '(A,L1)') "power_normalization_applied = ", p_applied
            WRITE (u, '(A,ES20.10)') "power_target_w = ", p_target_w
            WRITE (u, '(A,ES20.10)') "energy_per_fission_j = ", e_fis_j
            WRITE (u, '(A,ES20.10)') "fission_rate_target_1_per_s = ", f_target
            WRITE (u, '(A,ES20.10)') "fission_rate_model_pre_scale = ", f_model
            WRITE (u, '(A,ES20.10)') "flux_scale_factor = ", f_scale
            WRITE (u, '(A)') "power_norm_assumption: 1D chain treated as representative single-assembly proxy."
            WRITE (u, '(A)') "k_diag note: denominator is sum R_g*phi_g*dx only; edge diffusion leakage is separate."
            CLOSE (u)
        END IF
    END SUBROUTINE write_flux_cells_csv
END PROGRAM diffusion_solver_7g
