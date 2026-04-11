PROGRAM monte_carlo_solver
    USE read_input_file
    USE mesh_generation
    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER, PARAMETER :: max_steps_per_history = 4000
    INTEGER, PARAMETER :: number_of_energy_groups = 2

    INTEGER :: set_number, n_pos, n_cells
    INTEGER :: generation_index, history_index, pos_index, cell_index
    INTEGER :: n_points, birth_cell
    INTEGER :: n_active_generations, summary_unit
    INTEGER :: step_count

    REAL(dp), ALLOCATABLE :: x_left(:), x_right(:), dx(:), x_center(:)
    REAL(dp), ALLOCATABLE :: sigma_total(:,:), sigma_in_scatter(:,:), sigma_absorption(:,:), sigma_fission(:,:), nu_bar(:,:)
    REAL(dp), ALLOCATABLE :: sigma_down_scatter(:)
    REAL(dp), ALLOCATABLE :: source_pdf(:), source_next(:), track_length(:), k_history(:)

    REAL(dp) :: x, mu, free_path, remaining_path, distance_to_edge, distance_move
    REAL(dp) :: rand_u, p_capture, p_fission, p_in_scatter, p_down_scatter
    REAL(dp) :: sigma_capture_only
    REAL(dp) :: fission_neutrons_produced, source_sum
    REAL(dp) :: k_mean, k_std, k_last
    REAL(dp) :: ref_k_table_A, ref_eta_nu_sigf_over_sig_capture

    LOGICAL :: alive
    INTEGER :: neutron_group

    CHARACTER(LEN=*), PARAMETER :: output_dir = "outputs/monte_carlo"
    CHARACTER(LEN=*), PARAMETER :: mkdir_cmd = "mkdir -p outputs/monte_carlo"

    ! Make sure output folder exists before OPEN.
    CALL execute_command_line(mkdir_cmd, WAIT=.TRUE.)
    CALL RANDOM_SEED()

    CALL read_input('input_file.txt')
    CALL init_configurations()

    set_number = config + 1
    IF (set_number < 1 .OR. set_number > configs) THEN
        PRINT *, "ERROR: Config index out of range. Config=", config, " Configs=", configs
        STOP 1
    END IF

    CALL map_XS_to_mesh(set_number)

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
    ALLOCATE(sigma_total(n_cells, number_of_energy_groups))
    ALLOCATE(sigma_in_scatter(n_cells, number_of_energy_groups))
    ALLOCATE(sigma_absorption(n_cells, number_of_energy_groups))
    ALLOCATE(sigma_fission(n_cells, number_of_energy_groups))
    ALLOCATE(nu_bar(n_cells, number_of_energy_groups))
    ALLOCATE(sigma_down_scatter(n_cells))
    ALLOCATE(source_pdf(n_cells), source_next(n_cells), track_length(n_cells))
    ALLOCATE(k_history(generations))

    CALL build_geometry_and_xs( &
        set_number, n_cells, x_left, x_right, dx, x_center, &
        sigma_total, sigma_in_scatter, sigma_down_scatter, sigma_absorption, sigma_fission, nu_bar)

    ! Quick check vs handout Table 1: tabulated SigA already counts fission, so k ~ nu*SigF/SigA on fuel g2.
    ref_k_table_A = -1.0_dp
    ref_eta_nu_sigf_over_sig_capture = -1.0_dp
    DO cell_index = 1, n_cells
        IF (sigma_fission(cell_index, 2) > 1.0E-12_dp) THEN
            ref_k_table_A = nu_bar(cell_index, 2) * sigma_fission(cell_index, 2) / &
                MAX(sigma_absorption(cell_index, 2), 1.0E-15_dp)
            ref_eta_nu_sigf_over_sig_capture = nu_bar(cell_index, 2) * sigma_fission(cell_index, 2) / &
                MAX(sigma_absorption(cell_index, 2) - sigma_fission(cell_index, 2), 1.0E-15_dp)
            EXIT
        END IF
    END DO

    source_pdf = 0.0_dp
    DO cell_index = 1, n_cells
        IF (MAXVAL(sigma_fission(cell_index,:)) > 0.0_dp) source_pdf(cell_index) = 1.0_dp
    END DO
    IF (SUM(source_pdf) <= 0.0_dp) source_pdf = 1.0_dp
    source_pdf = source_pdf / SUM(source_pdf)

    ! Superhistories: each generation uses fission sites from the last one as births.
    DO generation_index = 1, generations
        track_length = 0.0_dp
        source_next = 0.0_dp
        fission_neutrons_produced = 0.0_dp

        DO history_index = 1, histories
            CALL sample_birth_cell(source_pdf, birth_cell)
            CALL RANDOM_NUMBER(rand_u)
            x = x_left(birth_cell) + rand_u * dx(birth_cell)
            cell_index = birth_cell
            neutron_group = 1
            alive = .TRUE.

            step_count = 0
            DO WHILE (alive .AND. step_count < max_steps_per_history)
                step_count = step_count + 1
                CALL RANDOM_NUMBER(rand_u)
                mu = 2.0_dp * rand_u - 1.0_dp
                IF (ABS(mu) < 1.0E-12_dp) mu = 1.0E-12_dp

                CALL RANDOM_NUMBER(rand_u)
                ! Exponential flight; Sigma_t = SigTR from input.
                free_path = -LOG(MAX(rand_u, 1.0E-15_dp)) / MAX(sigma_total(cell_index, neutron_group), 1.0E-15_dp)
                remaining_path = free_path

                DO WHILE (alive .AND. remaining_path > 1.0E-12_dp)
                    IF (mu > 0.0_dp) THEN
                        distance_to_edge = (x_right(cell_index) - x) / mu
                    ELSE
                        distance_to_edge = (x - x_left(cell_index)) / ABS(mu)
                    END IF
                    distance_to_edge = MAX(distance_to_edge, 0.0_dp)

                    IF (distance_to_edge >= remaining_path) THEN
                        distance_move = remaining_path
                        track_length(cell_index) = track_length(cell_index) + distance_move
                        x = x + mu * distance_move
                        remaining_path = 0.0_dp
                    ELSE
                        distance_move = distance_to_edge
                        track_length(cell_index) = track_length(cell_index) + distance_move
                        x = x + mu * distance_move
                        remaining_path = remaining_path - distance_move

                        ! Reflecting ends (symmetry / zero net current at assembly edge).
                        IF (mu > 0.0_dp) THEN
                            IF (cell_index == n_cells) THEN
                                mu = -mu
                                x = x_right(n_cells) - 1.0E-12_dp
                            ELSE
                                cell_index = cell_index + 1
                                x = MIN(x, x_right(cell_index) - 1.0E-12_dp)
                            END IF
                        ELSE
                            IF (cell_index == 1) THEN
                                mu = -mu
                                x = x_left(1) + 1.0E-12_dp
                            ELSE
                                cell_index = cell_index - 1
                                x = MAX(x, x_left(cell_index) + 1.0E-12_dp)
                            END IF
                        END IF
                    END IF
                END DO

                IF (.NOT. alive) EXIT

                ! Analog collision: split SigTR into capture (SigA minus fission), fission, scatter.
                sigma_capture_only = MAX(sigma_absorption(cell_index, neutron_group) - &
                    sigma_fission(cell_index, neutron_group), 0.0_dp)
                p_capture = sigma_capture_only / sigma_total(cell_index, neutron_group)
                p_fission = sigma_fission(cell_index, neutron_group) / sigma_total(cell_index, neutron_group)
                p_in_scatter = sigma_in_scatter(cell_index, neutron_group) / sigma_total(cell_index, neutron_group)
                IF (neutron_group == 1) THEN
                    p_down_scatter = sigma_down_scatter(cell_index) / sigma_total(cell_index, neutron_group)
                ELSE
                    p_down_scatter = 0.0_dp
                END IF
                p_capture = MAX(0.0_dp, p_capture)
                p_fission = MAX(0.0_dp, p_fission)
                p_in_scatter = MAX(0.0_dp, p_in_scatter)
                p_down_scatter = MAX(0.0_dp, p_down_scatter)

                CALL RANDOM_NUMBER(rand_u)
                IF (rand_u < p_capture) THEN
                    alive = .FALSE.
                ELSE IF (rand_u < p_capture + p_fission) THEN
                    source_next(cell_index) = source_next(cell_index) + nu_bar(cell_index, neutron_group)
                    fission_neutrons_produced = fission_neutrons_produced + nu_bar(cell_index, neutron_group)
                    alive = .FALSE.
                ELSE IF (rand_u < p_capture + p_fission + p_in_scatter) THEN
                    ! in-group scatter: continue with same group
                ELSE IF (rand_u < p_capture + p_fission + p_in_scatter + p_down_scatter) THEN
                    IF (neutron_group == 1) neutron_group = 2
                ELSE
                    ! any remainder probability is treated as scatter continuation
                END IF
            END DO
        END DO

        k_history(generation_index) = fission_neutrons_produced / REAL(histories, dp)

        source_sum = SUM(source_next)
        IF (source_sum > 0.0_dp) THEN
            source_pdf = source_next / source_sum
        ELSE
            source_pdf = 1.0_dp / REAL(n_cells, dp)
        END IF
    END DO

    n_active_generations = MAX(generations - skip, 1)
    k_mean = SUM(k_history(skip+1:generations)) / REAL(n_active_generations, dp)
    k_std = SQRT(MAX(0.0_dp, SUM((k_history(skip+1:generations) - k_mean)**2) / REAL(n_active_generations, dp)))
    k_last = k_history(generations)

    OPEN(NEWUNIT=summary_unit, FILE=output_dir // '/summary.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE(summary_unit,'(A,I0)') "config_set_0based = ", config
    WRITE(summary_unit,'(A,I0)') "positions_in_set = ", n_pos
    WRITE(summary_unit,'(A,I0)') "total_mesh_cells = ", n_cells
    WRITE(summary_unit,'(A,I0)') "generations = ", generations
    WRITE(summary_unit,'(A,I0)') "histories = ", histories
    WRITE(summary_unit,'(A,I0)') "skip = ", skip
    WRITE(summary_unit,'(A)') "boundary_outer = reflective_J0"
    WRITE(summary_unit,'(A)') "sigma_collision = SigTR; capture_prob uses (SigA-SigF)+ per Table 1"
    IF (ref_k_table_A >= 0.0_dp) THEN
        WRITE(summary_unit,'(A,F12.6)') "ref_k_nu_sigf_over_sigA_table = ", ref_k_table_A
        WRITE(summary_unit,'(A,F12.6)') "ref_eta_nu_sigf_over_sig_capture = ", ref_eta_nu_sigf_over_sig_capture
    END IF
    WRITE(summary_unit,'(A,F12.6)') "k_eff_mc_last = ", k_last
    WRITE(summary_unit,'(A,F12.6)') "k_eff_mc_mean = ", k_mean
    WRITE(summary_unit,'(A,F12.6)') "k_eff_mc_std = ", k_std
    CLOSE(summary_unit)

    PRINT *, "Monte Carlo solver complete."
    PRINT *, "Set (0-based):", config
    PRINT *, "Cells:", n_cells, " Generations:", generations, " Histories:", histories
    IF (ref_k_table_A >= 0.0_dp) THEN
        PRINT *, "Ref k (nu*SigF/SigA from table) =", ref_k_table_A
        PRINT *, "Ref eta (nu*SigF/(SigA-SigF)) =", ref_eta_nu_sigf_over_sig_capture
    END IF
    PRINT *, "k_eff (last) =", k_last
    PRINT *, "k_eff (mean after skip) =", k_mean
    PRINT *, "k_eff std =", k_std
    PRINT *, "Wrote outputs/monte_carlo/summary.txt"

CONTAINS

    SUBROUTINE build_geometry_and_xs( &
        set_num, total_cells, x_left, x_right, dx, x_center, &
        sigma_t, sigma_is, sigma_ds, sigma_a, sigma_f, nu)
        INTEGER, INTENT(IN) :: set_num, total_cells
        REAL(dp), INTENT(OUT) :: x_left(total_cells), x_right(total_cells), dx(total_cells), x_center(total_cells)
        REAL(dp), INTENT(OUT) :: sigma_t(total_cells, number_of_energy_groups)
        REAL(dp), INTENT(OUT) :: sigma_is(total_cells, number_of_energy_groups)
        REAL(dp), INTENT(OUT) :: sigma_ds(total_cells)
        REAL(dp), INTENT(OUT) :: sigma_a(total_cells, number_of_energy_groups)
        REAL(dp), INTENT(OUT) :: sigma_f(total_cells, number_of_energy_groups)
        REAL(dp), INTENT(OUT) :: nu(total_cells, number_of_energy_groups)

        INTEGER :: pos_i, mesh_i, n_points_local, cell_i
        REAL(dp) :: running_x, d_length

        cell_i = 0
        running_x = 0.0_dp
        DO pos_i = 1, set_len(set_num)
            IF (MatID(set_num, pos_i) == 2) THEN
                n_points_local = mpwr
            ELSE
                n_points_local = mpfr
            END IF
            d_length = rod_diameter / REAL(n_points_local, dp)

            DO mesh_i = 1, n_points_local
                cell_i = cell_i + 1
                x_left(cell_i) = running_x
                x_right(cell_i) = running_x + d_length
                dx(cell_i) = d_length
                x_center(cell_i) = 0.5_dp * (x_left(cell_i) + x_right(cell_i))
                running_x = x_right(cell_i)

                sigma_is(cell_i,1) = REAL(mesh_sigIS(pos_i, mesh_i, 1), dp)
                sigma_is(cell_i,2) = REAL(mesh_sigIS(pos_i, mesh_i, 2), dp)
                sigma_ds(cell_i) = REAL(mesh_sigDS(pos_i, mesh_i, 1), dp)
                sigma_a(cell_i,1) = REAL(mesh_sigA(pos_i, mesh_i, 1), dp)
                sigma_a(cell_i,2) = REAL(mesh_sigA(pos_i, mesh_i, 2), dp)
                sigma_f(cell_i,1) = REAL(mesh_sigF(pos_i, mesh_i, 1), dp)
                sigma_f(cell_i,2) = REAL(mesh_sigF(pos_i, mesh_i, 2), dp)
                nu(cell_i,1) = REAL(mesh_nuT(pos_i, mesh_i, 1), dp)
                nu(cell_i,2) = REAL(mesh_nuT(pos_i, mesh_i, 2), dp)

                ! Sigma_t for transport is SigTR (first column in the XS table).
                sigma_t(cell_i,1) = REAL(mesh_sigTR(pos_i, mesh_i, 1), dp)
                sigma_t(cell_i,2) = REAL(mesh_sigTR(pos_i, mesh_i, 2), dp)
                sigma_t(cell_i,1) = MAX(sigma_t(cell_i,1), 1.0E-15_dp)
                sigma_t(cell_i,2) = MAX(sigma_t(cell_i,2), 1.0E-15_dp)
            END DO
        END DO
    END SUBROUTINE build_geometry_and_xs

    ! Inverse-CDF sample from discrete source_pdf over cells.
    SUBROUTINE sample_birth_cell(source_pdf, sampled_cell)
        REAL(dp), INTENT(IN) :: source_pdf(:)
        INTEGER, INTENT(OUT) :: sampled_cell
        REAL(dp) :: r
        INTEGER :: i

        CALL RANDOM_NUMBER(r)
        DO i = 1, SIZE(source_pdf)
            r = r - source_pdf(i)
            IF (r <= 0.0_dp) THEN
                sampled_cell = i
                RETURN
            END IF
        END DO
        sampled_cell = SIZE(source_pdf)
    END SUBROUTINE sample_birth_cell

END PROGRAM monte_carlo_solver
