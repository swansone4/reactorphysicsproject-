PROGRAM monte_carlo_solver_7g
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    USE read_input_file
    USE mesh_generation
    USE csv_xs_7g
    IMPLICIT NONE

    INTEGER, PARAMETER :: max_steps_per_history = 4000
    INTEGER, PARAMETER :: ng = n_groups_7g

    INTEGER :: set_number, n_pos, n_cells
    INTEGER :: generation_index, history_index, pos_index, cell_index
    INTEGER :: n_points, birth_cell, birth_group
    INTEGER :: n_active_generations, summary_unit
    INTEGER :: step_count
    INTEGER :: left_bc(ng), right_bc(ng)
    INTEGER :: ig
    INTEGER :: plot_exitstat

    REAL(DP), ALLOCATABLE :: x_left(:), x_right(:), dx(:), x_center(:)
    REAL(DP), ALLOCATABLE :: sigma_total(:,:), sigma_absorption(:,:), sigma_fission(:,:), nu_bar(:,:)
    REAL(DP), ALLOCATABLE :: chi_birth(:,:), sigma_scatter(:,:,:), sigma_scatter_total(:,:)
    REAL(DP), ALLOCATABLE :: source_pdf(:,:), source_next(:,:), track_length(:,:), k_history(:)
    REAL(DP), ALLOCATABLE :: flux_sum(:,:), flux_fundamental(:,:)
    REAL(DP), ALLOCATABLE :: flux_int(:)

    REAL(DP) :: x, mu, free_path, remaining_path, distance_to_edge, distance_move
    REAL(DP) :: rand_u, p_capture, p_fission, p_scatter, p_sum
    REAL(DP) :: sigma_capture_only, nu_emit
    REAL(DP) :: fission_neutrons_produced, source_sum
    REAL(DP) :: k_mean, k_std, k_last

    LOGICAL :: alive
    INTEGER :: neutron_group

    CHARACTER(LEN=*), PARAMETER :: output_dir = "seven_groups/outputs/monte_carlo"
    CHARACTER(LEN=*), PARAMETER :: macro_xs_csv = "seven_groups/macroscopic_cross_sections.csv"
    CHARACTER(LEN=*), PARAMETER :: scatter_xs_csv = "seven_groups/scattering_matrix.csv"
    CHARACTER(LEN=*), PARAMETER :: plot_cmd = "python seven_groups/plot_monte_carlo_flux_7g.py"

    CALL ensure_output_dir_exists()
    CALL RANDOM_SEED()

    CALL read_input('input_file.txt')
    CALL init_configurations()
    CALL load_7g_xs_from_csv(macro_xs_csv, scatter_xs_csv)

    CALL build_7g_boundaries(left_bc, right_bc)

    set_number = config + 1
    IF (set_number < 1 .OR. set_number > configs) THEN
        PRINT *, "ERROR: Config index out of range. Config=", config, " Configs=", configs
        STOP 1
    END IF

    n_pos = set_len(set_number)
    n_cells = 0
    DO pos_index = 1, n_pos
        IF (MatID(set_number, pos_index) == 2) THEN
            n_points = mpwr
        ELSE
            n_points = mpfr
        END IF
        n_cells = n_cells + n_points
    END DO

    ALLOCATE(x_left(n_cells), x_right(n_cells), dx(n_cells), x_center(n_cells))
    ALLOCATE(sigma_total(n_cells, ng), sigma_absorption(n_cells, ng), sigma_fission(n_cells, ng), nu_bar(n_cells, ng))
    ALLOCATE(chi_birth(n_cells, ng), sigma_scatter(n_cells, ng, ng), sigma_scatter_total(n_cells, ng))
    ALLOCATE(source_pdf(n_cells, ng), source_next(n_cells, ng), track_length(n_cells, ng))
    ALLOCATE(k_history(generations), flux_sum(n_cells, ng), flux_fundamental(n_cells, ng), flux_int(ng))

    CALL build_geometry_and_xs_7g( &
        set_number, n_cells, x_left, x_right, dx, x_center, &
        sigma_total, sigma_absorption, sigma_fission, nu_bar, chi_birth, sigma_scatter, sigma_scatter_total)

    source_pdf = 0.0_DP
    DO cell_index = 1, n_cells
        IF (MAXVAL(sigma_fission(cell_index, :)) > 0.0_DP) THEN
            source_pdf(cell_index, :) = chi_birth(cell_index, :)
        END IF
    END DO
    IF (SUM(source_pdf) <= 0.0_DP) source_pdf(:, 1) = 1.0_DP
    source_pdf = source_pdf / SUM(source_pdf)

    flux_sum = 0.0_DP

    DO generation_index = 1, generations
        track_length = 0.0_DP
        source_next = 0.0_DP
        fission_neutrons_produced = 0.0_DP

        DO history_index = 1, histories
            CALL sample_birth_state(source_pdf, birth_cell, birth_group)
            CALL RANDOM_NUMBER(rand_u)
            x = x_left(birth_cell) + rand_u * dx(birth_cell)
            cell_index = birth_cell
            neutron_group = birth_group
            alive = .TRUE.

            step_count = 0
            DO WHILE (alive .AND. step_count < max_steps_per_history)
                step_count = step_count + 1
                CALL RANDOM_NUMBER(rand_u)
                mu = 2.0_DP * rand_u - 1.0_DP
                IF (ABS(mu) < 1.0E-12_DP) mu = 1.0E-12_DP

                CALL RANDOM_NUMBER(rand_u)
                free_path = -LOG(MAX(rand_u, 1.0E-15_DP)) / MAX(sigma_total(cell_index, neutron_group), 1.0E-15_DP)
                remaining_path = free_path

                DO WHILE (alive .AND. remaining_path > 1.0E-12_DP)
                    IF (mu > 0.0_DP) THEN
                        distance_to_edge = (x_right(cell_index) - x) / mu
                    ELSE
                        distance_to_edge = (x - x_left(cell_index)) / ABS(mu)
                    END IF
                    distance_to_edge = MAX(distance_to_edge, 0.0_DP)

                    IF (distance_to_edge >= remaining_path) THEN
                        distance_move = remaining_path
                        track_length(cell_index, neutron_group) = track_length(cell_index, neutron_group) + distance_move
                        x = x + mu * distance_move
                        remaining_path = 0.0_DP
                    ELSE
                        distance_move = distance_to_edge
                        track_length(cell_index, neutron_group) = track_length(cell_index, neutron_group) + distance_move
                        x = x + mu * distance_move
                        remaining_path = remaining_path - distance_move

                        IF (mu > 0.0_DP) THEN
                            IF (cell_index == n_cells) THEN
                                IF (right_bc(neutron_group) == 1) THEN
                                    mu = -mu
                                    x = x_right(n_cells) - 1.0E-12_DP
                                ELSE
                                    alive = .FALSE.
                                END IF
                            ELSE
                                cell_index = cell_index + 1
                                x = MIN(x, x_right(cell_index) - 1.0E-12_DP)
                            END IF
                        ELSE
                            IF (cell_index == 1) THEN
                                IF (left_bc(neutron_group) == 1) THEN
                                    mu = -mu
                                    x = x_left(1) + 1.0E-12_DP
                                ELSE
                                    alive = .FALSE.
                                END IF
                            ELSE
                                cell_index = cell_index - 1
                                x = MAX(x, x_left(cell_index) + 1.0E-12_DP)
                            END IF
                        END IF
                    END IF
                END DO

                IF (.NOT. alive) EXIT

                sigma_capture_only = MAX(sigma_absorption(cell_index, neutron_group) - &
                    sigma_fission(cell_index, neutron_group), 0.0_DP)
                p_capture = sigma_capture_only / sigma_total(cell_index, neutron_group)
                p_fission = sigma_fission(cell_index, neutron_group) / sigma_total(cell_index, neutron_group)
                p_scatter = sigma_scatter_total(cell_index, neutron_group) / sigma_total(cell_index, neutron_group)

                p_capture = MAX(0.0_DP, p_capture)
                p_fission = MAX(0.0_DP, p_fission)
                p_scatter = MAX(0.0_DP, p_scatter)
                p_sum = p_capture + p_fission + p_scatter
                IF (p_sum > 1.0_DP) THEN
                    p_capture = p_capture / p_sum
                    p_fission = p_fission / p_sum
                    p_scatter = p_scatter / p_sum
                END IF

                CALL RANDOM_NUMBER(rand_u)
                IF (rand_u < p_capture) THEN
                    alive = .FALSE.
                ELSE IF (rand_u < p_capture + p_fission) THEN
                    nu_emit = nu_bar(cell_index, neutron_group)
                    source_next(cell_index, :) = source_next(cell_index, :) + nu_emit * chi_birth(cell_index, :)
                    fission_neutrons_produced = fission_neutrons_produced + nu_emit
                    alive = .FALSE.
                ELSE IF (rand_u < p_capture + p_fission + p_scatter) THEN
                    CALL sample_scatter_group(sigma_scatter(cell_index, neutron_group, :), neutron_group)
                ELSE
                    ! leftover weight -> stay in same group
                END IF
            END DO
        END DO

        k_history(generation_index) = fission_neutrons_produced / REAL(histories, DP)
        IF (generation_index > skip) flux_sum = flux_sum + track_length

        source_sum = SUM(source_next)
        IF (source_sum > 0.0_DP) THEN
            source_pdf = source_next / source_sum
        ELSE
            source_pdf = 1.0_DP / REAL(n_cells * ng, DP)
        END IF
    END DO

    n_active_generations = MAX(generations - skip, 1)
    k_mean = SUM(k_history(skip+1:generations)) / REAL(n_active_generations, DP)
    k_std = SQRT(MAX(0.0_DP, SUM((k_history(skip+1:generations) - k_mean)**2) / REAL(n_active_generations, DP)))
    k_last = k_history(generations)

    DO cell_index = 1, n_cells
        DO ig = 1, ng
            flux_fundamental(cell_index, ig) = flux_sum(cell_index, ig) / &
                (REAL(histories * n_active_generations, DP) * dx(cell_index))
        END DO
    END DO

    DO ig = 1, ng
        flux_int(ig) = SUM(flux_fundamental(:, ig) * dx(:))
    END DO

    CALL write_flux_fundamental_csv_7g(n_cells, x_center, flux_fundamental)
    CALL write_flux_cells_csv_7g(n_cells, x_left, x_right, dx, flux_fundamental, &
        sigma_total, sigma_absorption, sigma_fission, nu_bar, sigma_scatter_total)
    CALL write_summary_7g(summary_unit, n_pos, n_cells, left_bc, right_bc, &
        k_last, k_mean, k_std, flux_int)

    PRINT *, "7-group Monte Carlo solver complete."
    PRINT *, "Set (0-based):", config
    PRINT *, "Cells:", n_cells, " Generations:", generations, " Histories:", histories
    PRINT *, "k_eff (last) =", k_last
    PRINT *, "k_eff (mean after skip) =", k_mean
    PRINT *, "k_eff std =", k_std
    PRINT *, "Wrote seven_groups/outputs/monte_carlo/summary.txt and flux CSV files."

    CALL execute_command_line(plot_cmd, WAIT=.TRUE., EXITSTAT=plot_exitstat)
    IF (plot_exitstat /= 0) THEN
        PRINT *, "Warning: 7g plot command failed (EXITSTAT=", plot_exitstat, ")."
    END IF

CONTAINS

    SUBROUTINE ensure_output_dir_exists()
        INTEGER :: stat_mkdir
        CALL execute_command_line( &
            "python -c ""from pathlib import Path; Path('seven_groups/outputs/monte_carlo').mkdir(parents=True, exist_ok=True)""", &
            WAIT=.TRUE., EXITSTAT=stat_mkdir)
        IF (stat_mkdir /= 0) THEN
            PRINT *, "Warning: could not ensure output directory exists."
        END IF
    END SUBROUTINE ensure_output_dir_exists

    SUBROUTINE build_7g_boundaries(left, right)
        INTEGER, INTENT(OUT) :: left(ng), right(ng)
        INTEGER :: g

        IF (.NOT. ALLOCATED(bound_l_7g) .OR. .NOT. ALLOCATED(bound_r_7g)) THEN
            PRINT *, "ERROR: BoundL_7G / BoundR_7G missing in input_file.txt parsing."
            STOP 1
        END IF
        IF (SIZE(bound_l_7g) < ng .OR. SIZE(bound_r_7g) < ng) THEN
            PRINT *, "ERROR: BoundL_7G / BoundR_7G must have 7 entries."
            STOP 1
        END IF

        DO g = 1, ng
            left(g) = NINT(bound_l_7g(g))
            right(g) = NINT(bound_r_7g(g))
        END DO

        DO g = 1, ng
            IF ((left(g) /= 0 .AND. left(g) /= 1) .OR. (right(g) /= 0 .AND. right(g) /= 1)) THEN
                PRINT *, "ERROR: boundaries must be 0 or 1. Group=", g, " left=", left(g), " right=", right(g)
                STOP 1
            END IF
        END DO
    END SUBROUTINE build_7g_boundaries

    SUBROUTINE build_geometry_and_xs_7g( &
        set_num, total_cells, x_left, x_right, dx, x_center, &
        sigma_t, sigma_a, sigma_f, nu, chi, sigma_s, sigma_s_total)
        INTEGER, INTENT(IN) :: set_num, total_cells
        REAL(DP), INTENT(OUT) :: x_left(total_cells), x_right(total_cells), dx(total_cells), x_center(total_cells)
        REAL(DP), INTENT(OUT) :: sigma_t(total_cells, ng), sigma_a(total_cells, ng), sigma_f(total_cells, ng), nu(total_cells, ng)
        REAL(DP), INTENT(OUT) :: chi(total_cells, ng), sigma_s(total_cells, ng, ng), sigma_s_total(total_cells, ng)

        INTEGER :: pos_i, mesh_i, n_points_local, cell_i, mid_i, m, g, h, n_pin_cells
        REAL(DP) :: running_x, d_length, pitch_dp, diam_dp, w_geom, w_sum, l_asm, scale_x, chi_sum

        pitch_dp = REAL(rod_pitch, DP)
        diam_dp = REAL(rod_diameter, DP)

        n_pin_cells = 0
        DO pos_i = 1, set_len(set_num)
            mid_i = MatID(set_num, pos_i)
            IF (mid_i == 0 .OR. mid_i == 1 .OR. mid_i == 3) n_pin_cells = n_pin_cells + 1
        END DO

        w_sum = 0.0_DP
        DO pos_i = 1, set_len(set_num)
            w_sum = w_sum + segment_width_x_cm(MatID(set_num, pos_i), pitch_dp, diam_dp)
        END DO

        IF (n_pin_cells > 0) THEN
            l_asm = REAL(n_pin_cells, DP) * pitch_dp
        ELSE
            l_asm = MAX(w_sum, pitch_dp)
        END IF

        IF (w_sum > 1.0E-15_DP) THEN
            scale_x = l_asm / w_sum
        ELSE
            scale_x = 1.0_DP
        END IF

        cell_i = 0
        running_x = 0.0_DP
        DO pos_i = 1, set_len(set_num)
            IF (MatID(set_num, pos_i) == 2) THEN
                n_points_local = mpwr
            ELSE
                n_points_local = mpfr
            END IF
            w_geom = segment_width_x_cm(MatID(set_num, pos_i), pitch_dp, diam_dp) * scale_x
            d_length = w_geom / REAL(n_points_local, DP)
            mid_i = MatID(set_num, pos_i)
            m = mid_i + 1
            IF (m < 1 .OR. m > n_materials_7g) THEN
                PRINT *, "ERROR: MatID out of range at pos=", pos_i, " mid=", mid_i
                STOP 1
            END IF

            DO mesh_i = 1, n_points_local
                cell_i = cell_i + 1
                x_left(cell_i) = running_x
                x_right(cell_i) = running_x + d_length
                dx(cell_i) = d_length
                x_center(cell_i) = 0.5_DP * (x_left(cell_i) + x_right(cell_i))
                running_x = x_right(cell_i)

                chi_sum = 0.0_DP
                DO g = 1, ng
                    sigma_t(cell_i, g) = MAX(xs_sigtr_7g(g, m), 1.0E-15_DP)
                    sigma_a(cell_i, g) = MAX(xs_siga_7g(g, m), 0.0_DP)
                    sigma_f(cell_i, g) = MAX(xs_sigf_7g(g, m), 0.0_DP)
                    IF (sigma_f(cell_i, g) > 1.0E-15_DP) THEN
                        nu(cell_i, g) = xs_nusigf_7g(g, m) / sigma_f(cell_i, g)
                    ELSE
                        nu(cell_i, g) = 0.0_DP
                    END IF
                    chi(cell_i, g) = MAX(xs_chi_7g(g, m), 0.0_DP)
                    chi_sum = chi_sum + chi(cell_i, g)

                    sigma_s_total(cell_i, g) = 0.0_DP
                    DO h = 1, ng
                        sigma_s(cell_i, g, h) = MAX(xs_sigs_7g(g, h, m), 0.0_DP)
                        sigma_s_total(cell_i, g) = sigma_s_total(cell_i, g) + sigma_s(cell_i, g, h)
                    END DO
                END DO
                IF (chi_sum > 0.0_DP) chi(cell_i, :) = chi(cell_i, :) / chi_sum
            END DO
        END DO
    END SUBROUTINE build_geometry_and_xs_7g

    PURE FUNCTION segment_width_x_cm(mid, pitch, diam) RESULT(w)
        INTEGER, INTENT(IN) :: mid
        REAL(DP), INTENT(IN) :: pitch, diam
        REAL(DP) :: w
        IF (mid == 2) THEN
            w = MAX(pitch - diam, 0.0_DP)
        ELSE
            w = diam
        END IF
    END FUNCTION segment_width_x_cm

    SUBROUTINE sample_scatter_group(scatter_row, sampled_group)
        REAL(DP), INTENT(IN) :: scatter_row(:)
        INTEGER, INTENT(INOUT) :: sampled_group
        REAL(DP) :: r, accum, row_sum
        INTEGER :: h

        row_sum = SUM(scatter_row)
        IF (row_sum <= 0.0_DP) RETURN
        CALL RANDOM_NUMBER(r)
        r = r * row_sum
        accum = 0.0_DP
        DO h = 1, SIZE(scatter_row)
            accum = accum + scatter_row(h)
            IF (r <= accum) THEN
                sampled_group = h
                RETURN
            END IF
        END DO
        sampled_group = SIZE(scatter_row)
    END SUBROUTINE sample_scatter_group

    SUBROUTINE sample_birth_state(pdf, sampled_cell, sampled_group)
        REAL(DP), INTENT(IN) :: pdf(:, :)
        INTEGER, INTENT(OUT) :: sampled_cell, sampled_group
        REAL(DP) :: r
        INTEGER :: i, g

        CALL RANDOM_NUMBER(r)
        DO i = 1, SIZE(pdf, 1)
            DO g = 1, SIZE(pdf, 2)
                r = r - pdf(i, g)
                IF (r <= 0.0_DP) THEN
                    sampled_cell = i
                    sampled_group = g
                    RETURN
                END IF
            END DO
        END DO
        sampled_cell = SIZE(pdf, 1)
        sampled_group = 1
    END SUBROUTINE sample_birth_state

    SUBROUTINE write_flux_fundamental_csv_7g(n, xc, flux_mc)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: xc(n), flux_mc(n, ng)
        INTEGER :: u, i, ios

        OPEN(NEWUNIT=u, FILE=output_dir // '/flux_fundamental.csv', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, "ERROR: could not open flux_fundamental.csv"
            STOP 1
        END IF

        WRITE(u, '(A)') "cell,x_cm,g1,g2,g3,g4,g5,g6,g7"
        DO i = 1, n
            WRITE(u, '(*(g0,:,","))') i, xc(i), flux_mc(i,1), flux_mc(i,2), flux_mc(i,3), flux_mc(i,4), &
                flux_mc(i,5), flux_mc(i,6), flux_mc(i,7)
        END DO
        CLOSE(u)
    END SUBROUTINE write_flux_fundamental_csv_7g

    SUBROUTINE write_flux_cells_csv_7g(n, xl, xr, dxc, flux_mc, st, sa, sf, nu, sst)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: xl(n), xr(n), dxc(n), flux_mc(n, ng), st(n, ng), sa(n, ng), sf(n, ng), nu(n, ng), sst(n, ng)
        INTEGER :: u, i, ios

        OPEN(NEWUNIT=u, FILE=output_dir // '/flux_cells.csv', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, "ERROR: could not open flux_cells.csv"
            STOP 1
        END IF

        WRITE(u, '(A)') "cell,x_left_cm,x_right_cm,dx_cm,flux_g1,flux_g2,flux_g3,flux_g4,flux_g5,flux_g6,flux_g7," // &
            "sigTR_g1,sigTR_g2,sigTR_g3,sigTR_g4,sigTR_g5,sigTR_g6,sigTR_g7," // &
            "sigA_g1,sigA_g2,sigA_g3,sigA_g4,sigA_g5,sigA_g6,sigA_g7," // &
            "sigF_g1,sigF_g2,sigF_g3,sigF_g4,sigF_g5,sigF_g6,sigF_g7," // &
            "nu_g1,nu_g2,nu_g3,nu_g4,nu_g5,nu_g6,nu_g7," // &
            "sigS_tot_g1,sigS_tot_g2,sigS_tot_g3,sigS_tot_g4,sigS_tot_g5,sigS_tot_g6,sigS_tot_g7"

        DO i = 1, n
            WRITE(u, '(*(g0,:,","))') i, xl(i), xr(i), dxc(i), &
                flux_mc(i,1), flux_mc(i,2), flux_mc(i,3), flux_mc(i,4), flux_mc(i,5), flux_mc(i,6), flux_mc(i,7), &
                st(i,1), st(i,2), st(i,3), st(i,4), st(i,5), st(i,6), st(i,7), &
                sa(i,1), sa(i,2), sa(i,3), sa(i,4), sa(i,5), sa(i,6), sa(i,7), &
                sf(i,1), sf(i,2), sf(i,3), sf(i,4), sf(i,5), sf(i,6), sf(i,7), &
                nu(i,1), nu(i,2), nu(i,3), nu(i,4), nu(i,5), nu(i,6), nu(i,7), &
                sst(i,1), sst(i,2), sst(i,3), sst(i,4), sst(i,5), sst(i,6), sst(i,7)
        END DO
        CLOSE(u)
    END SUBROUTINE write_flux_cells_csv_7g

    SUBROUTINE write_summary_7g(u, n_pos, n_cells, lbc, rbc, k_last, k_mean, k_std, flux_int)
        INTEGER, INTENT(INOUT) :: u
        INTEGER, INTENT(IN) :: n_pos, n_cells, lbc(ng), rbc(ng)
        REAL(DP), INTENT(IN) :: k_last, k_mean, k_std, flux_int(ng)
        INTEGER :: g

        OPEN(NEWUNIT=u, FILE=output_dir // '/summary.txt', STATUS='REPLACE', ACTION='WRITE')
        WRITE(u, '(A,I0)') "config_set_0based = ", config
        WRITE(u, '(A,I0)') "positions_in_set = ", n_pos
        WRITE(u, '(A,I0)') "total_mesh_cells = ", n_cells
        WRITE(u, '(A,I0)') "generations = ", generations
        WRITE(u, '(A,I0)') "histories = ", histories
        WRITE(u, '(A,I0)') "skip = ", skip
        WRITE(u, '(A,7(I0,1X))') "bound_left_g1_to_g7 = ", lbc
        WRITE(u, '(A,7(I0,1X))') "bound_right_g1_to_g7 = ", rbc
        WRITE(u, '(A)') "sigma_collision = SigTR from CSV; scatter from full 7x7 matrix CSV"
        WRITE(u, '(A,F12.6)') "k_eff_mc_last = ", k_last
        WRITE(u, '(A,F12.6)') "k_eff_mc_mean = ", k_mean
        WRITE(u, '(A,F12.6)') "k_eff_mc_std = ", k_std
        DO g = 1, ng
            WRITE(u, '(A,I0,A,F12.6)') "flux_integral_g", g, " = ", flux_int(g)
        END DO
        CLOSE(u)
    END SUBROUTINE write_summary_7g

END PROGRAM monte_carlo_solver_7g
