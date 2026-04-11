PROGRAM monte_carlo_bare_slab_benchmark
    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER, PARAMETER :: number_of_cells = 80
    INTEGER, PARAMETER :: number_of_generations = 800
    INTEGER, PARAMETER :: number_of_histories = 50000
    INTEGER, PARAMETER :: skip_generations = 40
    INTEGER, PARAMETER :: max_steps_per_history = 4000

    REAL(dp), PARAMETER :: slab_width_cm = 8.0_dp
    REAL(dp), PARAMETER :: sigma_total = 0.88525_dp
    REAL(dp), PARAMETER :: sigma_absorption = 0.11746_dp
    REAL(dp), PARAMETER :: sigma_fission = 0.08307_dp
    REAL(dp), PARAMETER :: nu_bar = 2.4695_dp
    REAL(dp), PARAMETER :: sigma_scatter = sigma_total - sigma_absorption

    REAL(dp), PARAMETER :: dx = slab_width_cm / REAL(number_of_cells, dp)

    INTEGER :: cell_index, generation_index, history_index
    INTEGER :: neutron_cell
    REAL(dp) :: cell_left(number_of_cells), cell_right(number_of_cells), cell_center(number_of_cells)
    REAL(dp) :: source_pdf(number_of_cells), source_next(number_of_cells)
    REAL(dp) :: track_length(number_of_cells), flux(number_of_cells), flux_fundamental(number_of_cells)
    REAL(dp) :: source_fundamental(number_of_cells)
    REAL(dp) :: k_history(number_of_generations), k_mc, k_mean, k_std, k_sq_sum, k_sum
    REAL(dp) :: x, mu, rand_u, free_path, remaining_path, distance_to_edge, distance_move
    REAL(dp) :: p_cap, p_fis
    REAL(dp) :: fission_neutrons_produced
    REAL(dp) :: source_sum
    LOGICAL :: alive
    INTEGER :: step_counter

    ! Diffusion reference: k_inf, buckling B^2=(pi/a)^2, k_eff = k_inf/(1+B^2 L^2), etc.
    REAL(dp) :: k_infinite
    REAL(dp) :: diffusion_coefficient
    REAL(dp) :: migration_area
    REAL(dp) :: geometric_buckling_sq
    REAL(dp) :: extrapolation_distance_cm
    REAL(dp) :: effective_slab_width_cm
    REAL(dp) :: k_analytic
    REAL(dp) :: shape_analytic(number_of_cells)
    INTEGER :: summary_unit

    CHARACTER(LEN=*), PARAMETER :: output_dir = "outputs"
    CHARACTER(LEN=*), PARAMETER :: mkdir_cmd = "mkdir -p outputs"

    CALL execute_command_line(mkdir_cmd, WAIT=.TRUE.)
    CALL RANDOM_SEED()

    p_cap = (sigma_absorption - sigma_fission) / sigma_total
    p_fis = sigma_fission / sigma_total

    DO cell_index = 1, number_of_cells
        cell_left(cell_index) = (REAL(cell_index - 1, dp)) * dx
        cell_right(cell_index) = (REAL(cell_index, dp)) * dx
        cell_center(cell_index) = 0.5_dp * (cell_left(cell_index) + cell_right(cell_index))
    END DO

    ! Uniform initial fission source guess.
    source_pdf = 1.0_dp / REAL(number_of_cells, dp)
    flux_fundamental = 0.0_dp
    source_fundamental = 0.0_dp
    k_sum = 0.0_dp
    k_sq_sum = 0.0_dp

    DO generation_index = 1, number_of_generations
        track_length = 0.0_dp
        source_next = 0.0_dp
        fission_neutrons_produced = 0.0_dp

        DO history_index = 1, number_of_histories
            CALL sample_birth_cell(source_pdf, neutron_cell)
            CALL RANDOM_NUMBER(rand_u)
            x = cell_left(neutron_cell) + rand_u * dx
            alive = .TRUE.
            step_counter = 0

            DO WHILE (alive .AND. step_counter < max_steps_per_history)
                step_counter = step_counter + 1

                CALL RANDOM_NUMBER(rand_u)
                mu = 2.0_dp * rand_u - 1.0_dp
                IF (ABS(mu) < 1.0E-12_dp) mu = SIGN(1.0E-12_dp, mu + 1.0E-13_dp)

                CALL RANDOM_NUMBER(rand_u)
                free_path = -LOG(MAX(rand_u, 1.0E-15_dp)) / sigma_total
                remaining_path = free_path

                DO WHILE (alive .AND. remaining_path > 1.0E-14_dp)
                    IF (mu > 0.0_dp) THEN
                        distance_to_edge = (cell_right(neutron_cell) - x) / mu
                    ELSE
                        distance_to_edge = (x - cell_left(neutron_cell)) / ABS(mu)
                    END IF
                    distance_to_edge = MAX(distance_to_edge, 0.0_dp)

                    IF (distance_to_edge >= remaining_path) THEN
                        ! Collision in current cell.
                        distance_move = remaining_path
                        track_length(neutron_cell) = track_length(neutron_cell) + distance_move
                        x = x + mu * distance_move
                        remaining_path = 0.0_dp
                    ELSE
                        ! Cross into adjacent cell (or leak at domain edge).
                        distance_move = distance_to_edge
                        track_length(neutron_cell) = track_length(neutron_cell) + distance_move
                        x = x + mu * distance_move
                        remaining_path = remaining_path - distance_move

                        IF (mu > 0.0_dp) THEN
                            IF (neutron_cell == number_of_cells) THEN
                                ! Right vacuum boundary -> leakage/termination.
                                alive = .FALSE.
                            ELSE
                                neutron_cell = neutron_cell + 1
                                x = MIN(x, cell_right(neutron_cell) - 1.0E-14_dp)
                            END IF
                        ELSE
                            IF (neutron_cell == 1) THEN
                                ! Left vacuum boundary -> leakage/termination.
                                alive = .FALSE.
                            ELSE
                                neutron_cell = neutron_cell - 1
                                x = MAX(x, cell_left(neutron_cell) + 1.0E-14_dp)
                            END IF
                        END IF
                    END IF
                END DO

                IF (.NOT. alive) EXIT

                ! Collision type CDF: capture, fission, scatter.
                CALL RANDOM_NUMBER(rand_u)
                IF (rand_u < p_cap) THEN
                    alive = .FALSE.
                ELSE IF (rand_u < p_cap + p_fis) THEN
                    source_next(neutron_cell) = source_next(neutron_cell) + nu_bar
                    fission_neutrons_produced = fission_neutrons_produced + nu_bar
                    alive = .FALSE.
                ELSE
                    ! Scatter: continue history with new direction/path.
                END IF
            END DO
        END DO

        ! Track-length flux estimator (shape metric for comparison).
        DO cell_index = 1, number_of_cells
            flux(cell_index) = track_length(cell_index) / (REAL(number_of_histories, dp) * dx)
        END DO

        ! Monte Carlo generation keff estimate.
        k_mc = fission_neutrons_produced / REAL(number_of_histories, dp)
        k_history(generation_index) = k_mc

        source_sum = SUM(source_next)
        IF (source_sum > 0.0_dp) THEN
            source_pdf = source_next / source_sum
        ELSE
            source_pdf = 1.0_dp / REAL(number_of_cells, dp)
        END IF

        IF (generation_index > skip_generations) THEN
            flux_fundamental = flux_fundamental + flux
            source_fundamental = source_fundamental + source_pdf
            k_sum = k_sum + k_mc
            k_sq_sum = k_sq_sum + k_mc * k_mc
        END IF
    END DO

    flux_fundamental = flux_fundamental / REAL(MAX(number_of_generations - skip_generations, 1), dp)
    source_fundamental = source_fundamental / REAL(MAX(number_of_generations - skip_generations, 1), dp)
    k_mean = k_sum / REAL(MAX(number_of_generations - skip_generations, 1), dp)
    k_std = SQRT(MAX(0.0_dp, k_sq_sum / REAL(MAX(number_of_generations - skip_generations, 1), dp) - k_mean * k_mean))

    ! Analytic diffusion buckling estimate.
    k_infinite = (nu_bar * sigma_fission) / sigma_absorption
    diffusion_coefficient = 1.0_dp / (3.0_dp * sigma_total)
    migration_area = diffusion_coefficient / sigma_absorption
    ! Vacuum diffusion correction using extrapolated boundaries:
    !   d_extrap ≈ 2.131*D,   a_eff = a + 2*d_extrap,   B^2 = (pi/a_eff)^2
    extrapolation_distance_cm = 2.131_dp * diffusion_coefficient
    effective_slab_width_cm = slab_width_cm + 2.0_dp * extrapolation_distance_cm
    geometric_buckling_sq = (ACOS(-1.0_dp) / effective_slab_width_cm)**2
    k_analytic = k_infinite / (1.0_dp + migration_area * geometric_buckling_sq)

    ! Analytic fundamental mode shape with extrapolated vacuum boundaries:
    ! phi(x) ~ cos(B*(x-a/2)), with B = pi/a_eff
    DO cell_index = 1, number_of_cells
        shape_analytic(cell_index) = COS((ACOS(-1.0_dp) / effective_slab_width_cm) * &
                                         (cell_center(cell_index) - 0.5_dp * slab_width_cm))
        shape_analytic(cell_index) = MAX(shape_analytic(cell_index), 0.0_dp)
    END DO
    CALL write_k_history(output_dir // "/k_history.csv", k_history)
    CALL write_flux_and_shape(output_dir // "/flux_vs_analytic_shape.csv", cell_center, flux_fundamental, shape_analytic)
    CALL write_source(output_dir // "/source_fundamental.csv", cell_center, source_fundamental)

    OPEN(NEWUNIT=summary_unit, FILE=output_dir // "/summary.txt", STATUS="REPLACE", ACTION="WRITE")
    WRITE(summary_unit,'(A,F12.6)') "slab_width_cm = ", slab_width_cm
    WRITE(summary_unit,'(A,F12.6)') "sigma_total = ", sigma_total
    WRITE(summary_unit,'(A,F12.6)') "sigma_absorption = ", sigma_absorption
    WRITE(summary_unit,'(A,F12.6)') "sigma_fission = ", sigma_fission
    WRITE(summary_unit,'(A,F12.6)') "nu_bar = ", nu_bar
    WRITE(summary_unit,'(A,F12.6)') "k_inf = ", k_infinite
    WRITE(summary_unit,'(A,F12.6)') "L2 = ", migration_area
    WRITE(summary_unit,'(A,F12.6)') "D = ", diffusion_coefficient
    WRITE(summary_unit,'(A,F12.6)') "extrapolation_distance_cm = ", extrapolation_distance_cm
    WRITE(summary_unit,'(A,F12.6)') "effective_slab_width_cm = ", effective_slab_width_cm
    WRITE(summary_unit,'(A,F12.6)') "B2 = ", geometric_buckling_sq
    WRITE(summary_unit,'(A,F12.6)') "k_eff_analytic_buckling = ", k_analytic
    WRITE(summary_unit,'(A,F12.6)') "k_eff_mc_mean = ", k_mean
    WRITE(summary_unit,'(A,F12.6)') "k_eff_mc_std = ", k_std
    WRITE(summary_unit,'(A,F12.6)') "k_eff_mc_last = ", k_history(number_of_generations)
    WRITE(summary_unit,'(A,F12.6)') "relative_error_percent = ", 100.0_dp * (k_mean - k_analytic) / k_analytic
    CLOSE(summary_unit)

    PRINT *, "Bare slab benchmark complete."
    PRINT *, "Wrote bare_slab_benchmark/outputs/summary.txt"
    PRINT *, "Wrote bare_slab_benchmark/outputs/k_history.csv"
    PRINT *, "Wrote bare_slab_benchmark/outputs/flux_vs_analytic_shape.csv"
    PRINT *, "Wrote bare_slab_benchmark/outputs/source_fundamental.csv"

CONTAINS

    SUBROUTINE sample_birth_cell(source_pdf, sampled_cell)
        REAL(dp), INTENT(IN) :: source_pdf(number_of_cells)
        INTEGER, INTENT(OUT) :: sampled_cell
        REAL(dp) :: r
        INTEGER :: i

        CALL RANDOM_NUMBER(r)
        DO i = 1, number_of_cells
            r = r - source_pdf(i)
            IF (r <= 0.0_dp) THEN
                sampled_cell = i
                RETURN
            END IF
        END DO
        sampled_cell = number_of_cells
    END SUBROUTINE sample_birth_cell

    SUBROUTINE write_k_history(path, k_hist)
        CHARACTER(LEN=*), INTENT(IN) :: path
        REAL(dp), INTENT(IN) :: k_hist(number_of_generations)
        INTEGER :: u, i

        OPEN(NEWUNIT=u, FILE=path, STATUS="REPLACE", ACTION="WRITE")
        WRITE(u,'(A)') "generation,k_eff"
        DO i = 1, number_of_generations
            WRITE(u,'(I0,A,F12.6)') i, ",", k_hist(i)
        END DO
        CLOSE(u)
    END SUBROUTINE write_k_history

    SUBROUTINE write_flux_and_shape(path, x_cm, flux_mc, flux_shape_analytic)
        CHARACTER(LEN=*), INTENT(IN) :: path
        REAL(dp), INTENT(IN) :: x_cm(number_of_cells), flux_mc(number_of_cells), flux_shape_analytic(number_of_cells)
        INTEGER :: u, i

        OPEN(NEWUNIT=u, FILE=path, STATUS="REPLACE", ACTION="WRITE")
        WRITE(u,'(A)') "cell,x_cm,flux_mc,flux_shape_analytic"
        DO i = 1, number_of_cells
            WRITE(u,'(I0,A,F12.6,A,F12.6,A,F12.6)') i, ",", x_cm(i), ",", flux_mc(i), ",", flux_shape_analytic(i)
        END DO
        CLOSE(u)
    END SUBROUTINE write_flux_and_shape

    SUBROUTINE write_source(path, x_cm, source_pdf)
        CHARACTER(LEN=*), INTENT(IN) :: path
        REAL(dp), INTENT(IN) :: x_cm(number_of_cells), source_pdf(number_of_cells)
        INTEGER :: u, i

        OPEN(NEWUNIT=u, FILE=path, STATUS="REPLACE", ACTION="WRITE")
        WRITE(u,'(A)') "cell,x_cm,source_pdf"
        DO i = 1, number_of_cells
            WRITE(u,'(I0,A,F12.6,A,F12.6)') i, ",", x_cm(i), ",", source_pdf(i)
        END DO
        CLOSE(u)
    END SUBROUTINE write_source

END PROGRAM monte_carlo_bare_slab_benchmark
