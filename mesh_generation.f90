MODULE mesh_generation 
    USE read_input_file
    IMPLICIT NONE 

    !===================================================================================
    ! 
    ! Here, we read the parameters from the input file reader (read_input_file.f90) 
    ! and then use the parameters to generate the mesh. 
    !
    ! We must model (for Part I of the project) the fuel assembly which is visually represented in 
    ! \images\fuel_assembly_config.png
    ! 
    ! Our mesh is defined at the bottom of the input file in the ConfigSets section.
    !
    ! To create a new mesh, simply create a new "Set = N" where N is the next available set number
    ! and then define the material IDs 
    !
    !===================================================================================
     
     ! Here we declare some important variables 
     LOGICAL :: plot_visuals = .TRUE.

     INTEGER :: num_sets ! Number of configuration sets
     INTEGER, ALLOCATABLE :: MatID(:,:) ! Material IDs for positions [set, pos]
     INTEGER, ALLOCATABLE :: set_len(:) ! number of positions for each set
     REAL, ALLOCATABLE :: mesh(:,:) ! Mesh points for positions [pos, mesh_index]

     !===================================================================================
     !
     ! Cross section data mapped to each position in a configuration set.
     ! These arrays are populated by map_XS_to_meshes() using:
     !   - MatID(set_number, pos_index)
     !   - test_case (from input file)
     !   - xs_*(:,:,:) arrays built in read_input_file.f90
     !
     ! Dimensions: (position, energy_group)
     !
     !===================================================================================
    REAL, ALLOCATABLE :: pos_sigTR(:,:), pos_sigIS(:,:), pos_sigDS(:,:)
    REAL, ALLOCATABLE :: pos_sigA(:,:),  pos_sigF(:,:),  pos_nuT(:,:), pos_chiT(:,:)
    REAL, ALLOCATABLE :: mesh_sigTR(:,:,:), mesh_sigIS(:,:,:), mesh_sigDS(:,:,:)
    REAL, ALLOCATABLE :: mesh_sigA(:,:,:),  mesh_sigF(:,:,:),  mesh_nuT(:,:,:), mesh_chiT(:,:,:)

CONTAINS

    !===========================================================
    ! Initialize configuration sets
    ! Each rod is assigned a material ID (0 = UO2, 1 = MOX, 2 = H2O, 3 = CR)
    !-----------------------------------------------------------
    SUBROUTINE init_configurations()
        INTEGER :: set_index, max_len

        num_sets = configs

        IF (.NOT. ALLOCATED(config_set_len)) THEN
            PRINT *, "ERROR: ConfigSets were not loaded. Check input_file.txt formatting."
            STOP 1
        END IF

        IF (SIZE(config_set_len) /= num_sets) THEN
            PRINT *, "ERROR: ConfigSets size mismatch vs Configs."
            STOP 1
        END IF

        IF (MAXVAL(config_set_len) <= 0) THEN
            PRINT *, "ERROR: ConfigSets present but empty."
            STOP 1
        END IF

        max_len = MAXVAL(config_set_len)

        IF (ALLOCATED(MatID)) DEALLOCATE(MatID)
        IF (ALLOCATED(set_len)) DEALLOCATE(set_len)
        ALLOCATE(MatID(num_sets, max_len))
        ALLOCATE(set_len(num_sets))

        MatID = -1
        set_len = config_set_len

        IF (ALLOCATED(config_set_matid)) THEN
            MatID(:, :) = config_set_matid(:, :)
        ELSE
            PRINT *, "ERROR: ConfigSets lengths exist but MatID data missing."
            STOP 1
        END IF

        DO set_index = 1, num_sets
            IF (set_len(set_index) <= 0) THEN
                PRINT *, "ERROR: Config set ", set_index, " has no MatID entries."
                STOP 1
            END IF
        END DO
    END SUBROUTINE init_configurations

    !===========================================================
    ! Map cross sections to each mesh position for a given set.
    !
    ! This builds a simple "material -> XS column" mapping so later solvers
    ! (finite difference / Monte Carlo) can query cross sections per position.
    !-----------------------------------------------------------
    SUBROUTINE map_XS_to_mesh(set_number)
        INTEGER, INTENT(IN) :: set_number
        ! pos_index = index of the position in the configuration set
        ! mesh_index = index of mesh point in the rod
        ! g = energy group index 
        ! n_pos = number of positions in the configuration set
        ! n_points = number of mesh points in this rod
        ! mid = material ID (0 = UO2, 1 = MOX, 2 = H2O, 3 = CR)
        INTEGER :: pos_index, mesh_index, g, n_pos, n_points, mid, case_idx

        n_pos = set_len(set_number)

        ! Input TestCase = 0..cases-1 maps to xs_* first dimension 1..cases (Fortran index = test_case + 1).
        case_idx = test_case + 1
        IF (case_idx < 1 .OR. case_idx > cases) THEN
            PRINT *, "ERROR: test_case out of range for loaded XS cases. test_case=", test_case, " cases=", cases
            STOP 1
        END IF

        IF (ALLOCATED(pos_sigTR)) DEALLOCATE(pos_sigTR)
        IF (ALLOCATED(pos_sigIS)) DEALLOCATE(pos_sigIS)
        IF (ALLOCATED(pos_sigDS)) DEALLOCATE(pos_sigDS)
        IF (ALLOCATED(pos_sigA))  DEALLOCATE(pos_sigA)
        IF (ALLOCATED(pos_sigF))  DEALLOCATE(pos_sigF)
        IF (ALLOCATED(pos_nuT))   DEALLOCATE(pos_nuT)
        IF (ALLOCATED(pos_chiT))  DEALLOCATE(pos_chiT)
        IF (ALLOCATED(mesh_sigTR)) DEALLOCATE(mesh_sigTR)
        IF (ALLOCATED(mesh_sigIS)) DEALLOCATE(mesh_sigIS)
        IF (ALLOCATED(mesh_sigDS)) DEALLOCATE(mesh_sigDS)
        IF (ALLOCATED(mesh_sigA))  DEALLOCATE(mesh_sigA)
        IF (ALLOCATED(mesh_sigF))  DEALLOCATE(mesh_sigF)
        IF (ALLOCATED(mesh_nuT))   DEALLOCATE(mesh_nuT)
        IF (ALLOCATED(mesh_chiT))  DEALLOCATE(mesh_chiT)

        ALLOCATE(pos_sigTR(n_pos, energy_groups))
        ALLOCATE(pos_sigIS(n_pos, energy_groups))
        ALLOCATE(pos_sigDS(n_pos, energy_groups))
        ALLOCATE(pos_sigA (n_pos, energy_groups))
        ALLOCATE(pos_sigF (n_pos, energy_groups))
        ALLOCATE(pos_nuT  (n_pos, energy_groups))
        ALLOCATE(pos_chiT (n_pos, energy_groups))
        ALLOCATE(mesh_sigTR(n_pos, MAX(mpfr, mpwr), energy_groups))
        ALLOCATE(mesh_sigIS(n_pos, MAX(mpfr, mpwr), energy_groups))
        ALLOCATE(mesh_sigDS(n_pos, MAX(mpfr, mpwr), energy_groups))
        ALLOCATE(mesh_sigA (n_pos, MAX(mpfr, mpwr), energy_groups))
        ALLOCATE(mesh_sigF (n_pos, MAX(mpfr, mpwr), energy_groups))
        ALLOCATE(mesh_nuT  (n_pos, MAX(mpfr, mpwr), energy_groups))
        ALLOCATE(mesh_chiT (n_pos, MAX(mpfr, mpwr), energy_groups))

        mesh_sigTR = 0.0
        mesh_sigIS = 0.0
        mesh_sigDS = 0.0
        mesh_sigA  = 0.0
        mesh_sigF  = 0.0
        mesh_nuT   = 0.0
        mesh_chiT  = 0.0

        DO pos_index = 1, n_pos
            mid = MatID(set_number, pos_index)
            IF (mid < 0 .OR. mid >= mattypes) THEN
                PRINT *, "ERROR: MatID out of range at pos=", pos_index, " mid=", mid, " MatTypes=", mattypes
                STOP 1
            END IF

            IF (mid == 0 .OR. mid == 1 .OR. mid == 3) THEN
                n_points = mpfr
            ELSE
                n_points = mpwr
            END IF

            DO g = 1, energy_groups
                pos_sigTR(pos_index, g) = xs_sigTR(case_idx, g, mid + 1)
                pos_sigIS(pos_index, g) = xs_sigIS(case_idx, g, mid + 1)
                pos_sigDS(pos_index, g) = xs_sigDS(case_idx, g, mid + 1)
                pos_sigA (pos_index, g) = xs_sigA (case_idx, g, mid + 1)
                pos_sigF (pos_index, g) = xs_sigF (case_idx, g, mid + 1)
                pos_nuT  (pos_index, g) = xs_nuT (case_idx, g, mid + 1)
                pos_chiT (pos_index, g) = xs_chiT(case_idx, g, mid + 1)

                DO mesh_index = 1, n_points
                    mesh_sigTR(pos_index, mesh_index, g) = pos_sigTR(pos_index, g)
                    mesh_sigIS(pos_index, mesh_index, g) = pos_sigIS(pos_index, g)
                    mesh_sigDS(pos_index, mesh_index, g) = pos_sigDS(pos_index, g)
                    mesh_sigA (pos_index, mesh_index, g) = pos_sigA (pos_index, g)
                    mesh_sigF (pos_index, mesh_index, g) = pos_sigF (pos_index, g)
                    mesh_nuT  (pos_index, mesh_index, g) = pos_nuT  (pos_index, g)
                    mesh_chiT (pos_index, mesh_index, g) = pos_chiT (pos_index, g)
                END DO
            END DO
        END DO
    END SUBROUTINE map_XS_to_mesh

    ! Backward-compatible wrapper for existing call sites.
    SUBROUTINE map_XS_to_meshes(set_number)
        INTEGER, INTENT(IN) :: set_number
        CALL map_XS_to_mesh(set_number)
    END SUBROUTINE map_XS_to_meshes

    !===========================================================
    ! Writes one CSV for a single XS library case: xs_*(case_idx, g, matcol) with
    ! matcol = 1..mattypes (same indexing as map_XS_to_meshes uses via matid+1).
    ! Enable with VerifyXS = 1 in the input file. Case selection:
    !   - Optional argument selected_xs_case: if >= 0, use that case index (0=A,1=B,...);
    !     if present and -1, use input TestCase.
    !   - If omitted: use VerifyXSCase from input when >= 0, else TestCase.
    !-----------------------------------------------------------
    SUBROUTINE verify_XS_mapping(selected_xs_case)
        INTEGER, INTENT(IN), OPTIONAL :: selected_xs_case
        INTEGER :: u, ios, mcol, g, case_idx, case_input
        CHARACTER(LEN=64) :: fname
        CHARACTER(LEN=16) :: mname

        IF (.NOT. verify_xs_request) RETURN
        IF (.NOT. ALLOCATED(xs_sigTR)) THEN
            PRINT *, "verify_XS_mapping: XS arrays not loaded; skipping."
            RETURN
        END IF

        IF (PRESENT(selected_xs_case)) THEN
            IF (selected_xs_case < 0) THEN
                case_input = test_case
            ELSE
                case_input = selected_xs_case
            END IF
        ELSE IF (verify_xs_case_choice < 0) THEN
            case_input = test_case
        ELSE
            case_input = verify_xs_case_choice
        END IF

        case_idx = case_input + 1
        IF (case_idx < 1 .OR. case_idx > cases) THEN
            PRINT *, "verify_XS_mapping: case out of range. case_input=", case_input, " cases=", cases
            RETURN
        END IF

        CALL execute_command_line("mkdir -p outputs/xs_mapping", WAIT=.TRUE.)
        WRITE(fname, '(A,I0,A)') "outputs/xs_mapping/xs_mapping_verify_case_", case_input, ".csv"
        OPEN(NEWUNIT=u, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, "verify_XS_mapping: could not open ", TRIM(fname)
            RETURN
        END IF

        WRITE(u, '(A)') "source,xs_case_index_0based,xs_fortran_case_index,material_id," // &
            "material_name,energy_group,sigTR,sigIS,sigDS,sigA,sigF,nuT,chiT"

        DO mcol = 1, mattypes
            SELECT CASE (mcol - 1)
            CASE (0); mname = "UO2"
            CASE (1); mname = "MOX"
            CASE (2); mname = "H2O"
            CASE (3); mname = "CR"
            CASE DEFAULT
                WRITE(mname, '(A,I0)') "M", mcol - 1
            END SELECT
            DO g = 1, energy_groups
                WRITE(u, '(A,I0,A,I0,A,I0,A,A,A,I0,7(A,F0.8))') &
                    "xs_library,", case_input, ",", case_idx, ",", mcol - 1, ",", &
                    TRIM(ADJUSTL(mname)), ",", g, &
                    ",", xs_sigTR(case_idx, g, mcol), &
                    ",", xs_sigIS(case_idx, g, mcol), &
                    ",", xs_sigDS(case_idx, g, mcol), &
                    ",", xs_sigA (case_idx, g, mcol), &
                    ",", xs_sigF (case_idx, g, mcol), &
                    ",", xs_nuT (case_idx, g, mcol), &
                    ",", xs_chiT(case_idx, g, mcol)
            END DO
        END DO

        CLOSE(u)
        PRINT *, "verify_XS_mapping: wrote ", TRIM(fname), " (XS case ", case_input, " only)"
    END SUBROUTINE verify_XS_mapping

    !===========================================================
    ! Dump a lightweight CSV for Python visualization.
    !-----------------------------------------------------------
    SUBROUTINE write_mesh_dump(set_number)
        INTEGER, INTENT(IN) :: set_number
        INTEGER :: pos_index, n_points, n_pos, u, ios
        INTEGER :: row, col, cols_in_row, pos_in_set
        INTEGER :: dump_nrows, dump_ncols
        CHARACTER(LEN=128) :: filename

        n_pos = set_len(set_number)
        CALL execute_command_line("mkdir -p outputs/mesh_dumps", WAIT=.TRUE.)
        WRITE(filename,'(A,I0,A)') "outputs/mesh_dumps/mesh_dump_set_", set_number-1, ".csv"

        OPEN(NEWUNIT=u, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, "ERROR: Could not open dump file:", TRIM(filename)
            STOP 1
        END IF

        WRITE(u,'(A)', IOSTAT=ios) &
            "pos,row,col,nrows,ncols,rod_pitch,rod_diameter,matid,n_points,r_first,r_last"

        row = 1
        col = 0
        cols_in_row = 0
        pos_in_set = 0
        dump_nrows = 1
        dump_ncols = n_pos
        IF (ALLOCATED(config_set_nrows)) dump_nrows = config_set_nrows(set_number)
        IF (ALLOCATED(config_set_row_len)) dump_ncols = MAXVAL(config_set_row_len(set_number, 1:dump_nrows))

        DO pos_index = 1, n_pos
            IF (MatID(set_number, pos_index) == 0 .OR. MatID(set_number, pos_index) == 1 .OR. &
                MatID(set_number, pos_index) == 3) THEN
                n_points = mpfr
            ELSE
                n_points = mpwr
            END IF

            ! Use the actual ConfigSets row structure (each MatID line = a row)
            pos_in_set = pos_in_set + 1
            IF (ALLOCATED(config_set_row_len) .AND. ALLOCATED(config_set_nrows)) THEN
                IF (row > config_set_nrows(set_number)) THEN
                    PRINT *, "ERROR: ConfigSets row structure mismatch for set", set_number
                    STOP 1
                END IF
                IF (cols_in_row == 0) cols_in_row = config_set_row_len(set_number, row)
                col = col + 1
                IF (col > cols_in_row) THEN
                    row = row + 1
                    col = 1
                    IF (row > config_set_nrows(set_number)) THEN
                        PRINT *, "ERROR: ConfigSets row structure mismatch for set", set_number
                        STOP 1
                    END IF
                    cols_in_row = config_set_row_len(set_number, row)
                END IF
            ELSE
                row = 1
                col = pos_index
            END IF

            WRITE(u,'(I0,A,I0,A,I0,A,I0,A,I0,A,F0.6,A,F0.6,A,I0,A,I0,A,F0.6,A,F0.6)') &
                pos_index, ",", row, ",", col, ",", dump_nrows, ",", dump_ncols, ",", &
                rod_pitch, ",", rod_diameter, ",", &
                MatID(set_number, pos_index), ",", n_points, ",", &
                mesh(pos_index,1), ",", mesh(pos_index,n_points)
        END DO

        CLOSE(u)
    END SUBROUTINE write_mesh_dump

    !===========================================================
    ! Call Python to generate visualization PNGs.
    !-----------------------------------------------------------
    SUBROUTINE run_mesh_visualization()
        CALL execute_command_line("python scripts/mesh_visualization.py", WAIT=.TRUE.)
    END SUBROUTINE run_mesh_visualization

    !===========================================================
    ! Generate mesh points for all rods in a configuration set, these are the center points of each mesh segment. 
    ! The number of mesh points is determined by the material type (MPFR for fuel, MPWR for water) as the project requirements specify. 
    !-----------------------------------------------------------
    SUBROUTINE generate_mesh(set_number)
        INTEGER, INTENT(IN) :: set_number  ! Configuration set to use
        INTEGER :: pos_index, mesh_index, n_points, n_pos
        REAL :: d_length                    ! Length of each mesh segment

        IF (set_number < 1 .OR. set_number > num_sets) THEN
            PRINT *, "ERROR: Requested set out of range:", set_number
            STOP 1
        END IF

        n_pos = set_len(set_number)

        IF (ALLOCATED(mesh)) DEALLOCATE(mesh)
        ALLOCATE(mesh(n_pos, MAX(mpfr, mpwr)))

        ! Loop over each position
        DO pos_index = 1, n_pos

            ! Determine number of mesh points for this rod
            IF (MatID(set_number, pos_index) == 0 .OR. MatID(set_number, pos_index) == 1 .OR. &
                MatID(set_number, pos_index) == 3) THEN
                ! Fuel rod uses MPFR mesh points, we determine the number of points here by dividing the diameter by MPFR
                n_points = mpfr
                d_length = rod_diameter / REAL(mpfr)
            ELSE
                ! Water rod uses MPWR mesh points, we determine the number of points here by dividing the diameter by MPWR
                n_points = mpwr
                d_length = rod_diameter / REAL(mpwr)
            END IF

            ! Generate mesh points (center of each segment)
            DO mesh_index = 1, n_points
                mesh(pos_index, mesh_index) = d_length * (mesh_index - 0.5)
            END DO

        END DO
    END SUBROUTINE generate_mesh

    !===========================================================
    ! Print the mesh for a configuration set
    !-----------------------------------------------------------
    SUBROUTINE print_mesh(set_number)
        INTEGER, INTENT(IN) :: set_number
        INTEGER :: pos_index, n_points, n_pos
        INTEGER :: n_pin_cells, mid_i
        REAL :: total_config_length

        n_pin_cells = 0
        DO pos_index = 1, set_len(set_number)
            mid_i = MatID(set_number, pos_index)
            IF (mid_i == 0 .OR. mid_i == 1 .OR. mid_i == 3) n_pin_cells = n_pin_cells + 1
        END DO
        total_config_length = REAL(n_pin_cells) * rod_pitch

        PRINT *, "=== Mesh for Configuration Set ", set_number - 1 , " ==="
        PRINT *, "The total length of this configuration (x - cm) is:", total_config_length
        PRINT *, "There are a total of ", set_len(set_number), " positions in this set."
        PRINT *, "(n) corresponds to the number of mesh points for this position."
        PRINT *, "(r) corresponds to the radius of each mesh point for each axial position."
        PRINT *, "(M) corresponds to the material ID. 0 = UO2, 1 = MOX, 2 = H2O, 3 = CR"

        n_pos = set_len(set_number)

        ! Minimal per-position summary (material, points, first/last radius)
        DO pos_index = 1, n_pos
            IF (MatID(set_number, pos_index) == 0 .OR. MatID(set_number, pos_index) == 1 .OR. &
                MatID(set_number, pos_index) == 3) THEN
                n_points = mpfr
            ELSE
                n_points = mpwr
            END IF

            WRITE(*,'(A,I4,A,I2,A,I2,A,F8.4,A,F8.4,A)') &
                "Pos ", pos_index, " | M=", MatID(set_number, pos_index), &
                " | n=", n_points, " | r=[", mesh(pos_index,1), ", ", mesh(pos_index,n_points), "]"
        END DO

        PRINT *, "==============================================="
    END SUBROUTINE print_mesh

END MODULE mesh_generation