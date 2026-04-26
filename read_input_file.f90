MODULE read_input_file
    IMPLICIT NONE

    !===================================================================================
    ! Basic input file variables and parameters
    !===================================================================================

    ! solution variable 
    INTEGER :: solution ! method for solving the problem (0 = Finite difference, to be created), (1 = Monte Carlo) 
    ! test case variable 
    INTEGER :: test_case 
    ! config variable 
    INTEGER :: config 
    ! Analk variable
    INTEGER :: analk
    ! cases variable 
    INTEGER :: cases
    ! configs variable 
    INTEGER :: configs
    ! mattypes variable 
    INTEGER :: mattypes
    ! energy groups variable (default 2 if omitted from input)
    INTEGER :: energy_groups = 2
    ! solver variable 
    INTEGER :: solver
    ! generations variable 
    INTEGER :: generations
    ! histories variable 
    INTEGER :: histories 
    ! skip variable 
    INTEGER :: skip 
    ! If .TRUE., write xs_mapping_verify_case_<n>.csv (see mesh_generation::verify_XS_mapping).
    LOGICAL :: verify_xs_request = .FALSE.
    ! XS case to dump (0=A,1=B,...); if < 0, use current TestCase from input (unless verify_XS_mapping is called with an explicit case).
    INTEGER :: verify_xs_case_choice = -1
    ! number of assemblies variable 
    INTEGER :: num_assemblies
    ! number of rods variable 
    INTEGER :: num_rods
    ! rod diameter value 
    REAL :: rod_diameter
    ! rod pitch value 
    REAL :: rod_pitch
    ! number of Meshes Per Fuel Rod (MPFR) variable 
    INTEGER :: mpfr
    ! number of Meshes Per Water Rod (MPWR). MPWR should be an even number 
    INTEGER :: mpwr
    ! left / right boundary (0 = vacuum, 1 = reflective), length = EnergyGroups
    REAL, ALLOCATABLE :: bound_l(:)
    REAL, ALLOCATABLE :: bound_r(:)
    ! Explicit seven-group boundaries for seven_groups MC (optional in input_file.txt).
    REAL, ALLOCATABLE :: bound_l_7g(:)
    REAL, ALLOCATABLE :: bound_r_7g(:)

    !===================================================================================
    ! Configuration sets (from ConfigSets block)
    !===================================================================================
    INTEGER, ALLOCATABLE :: config_set_len(:)        ! number of MatID entries per set (flattened)
    INTEGER, ALLOCATABLE :: config_set_matid(:,:)    ! [set_index(1..configs), 1..max_len], padded with -1
    INTEGER, ALLOCATABLE :: config_set_nrows(:)      ! number of MatID lines (rows) per set
    INTEGER, ALLOCATABLE :: config_set_row_len(:,:)  ! row lengths per set [set, row]

    !===================================================================================
    ! Now, the cross-section (XS) data arrays 
    !===================================================================================

    ! Defining the mattype variable for indexing the cross section data arrays 
    INTEGER :: mattype

    !===================================================================================
    ! Cross section data for: 
    ! Energy group = 1, Test case = A (case = 0) 
    ! Declaring the length of all of these to be 4 since thats how it is in the input file, can be changed however 
    !===================================================================================

    REAL :: sigTR_1_A(4) ! Transport cross section (sigTR) for energy group 1 for test case A
    REAL :: sigIS_1_A(4) ! inelastic scatterring cross section (sigIS) for energy group 1
    REAL :: sigDS_1_A(4) ! downscattering cross section (sigDS) for energy group 1
    REAL :: sigA_1_A(4)  ! absorption cross section (sigA) for energy group 1
    REAL :: sigF_1_A(4)  ! fission cross section (sigF) for energy group 1
    REAL :: nuT_1_A(4)   ! average number of neutrons emitted per fission (nuT) for energy group 1
    REAL :: chiT_1_A(4)  ! total neutron fission specture (chiT) for energy group 1

    !===================================================================================
    ! Cross section data for: 
    ! Energy group = 2, Test case = A (case = 0)
    !===================================================================================

    REAL :: sigTR_2_A(4) ! transport cross section (sigTR) for energy group 2 for test case A
    REAL :: sigIS_2_A(4) ! inelastic scatterring cross section (sigIS) for energy group 2
    REAL :: sigDS_2_A(4) ! downscattering cross section (sigDS) for energy group 2
    REAL :: sigA_2_A(4) ! absorption cross section (sigA) for energy group 2
    REAL :: sigF_2_A(4) ! fission cross section (sigF) for energy group 2
    REAL :: nuT_2_A(4)  ! average number of neutrons emitted per fission (nuT) for energy group 2
    REAL :: chiT_2_A(4) ! total neutron fission specture (chiT) for energy group 2

    !===================================================================================
    ! Cross section data for: 
    ! Energy group = 1, Test case = B (case = 1)
    !===================================================================================

    REAL :: sigTR_1_B(4) ! transport cross section (sigTR) for energy group 1 for test case B
    REAL :: sigIS_1_B(4) ! inelastic scatterring cross section (sigIS) for energy group 1
    REAL :: sigDS_1_B(4) ! downscattering cross section (sigDS) for energy group 1
    REAL :: sigA_1_B(4) ! absorption cross section (sigA) for energy group 1
    REAL :: sigF_1_B(4) ! fission cross section (sigF) for energy group 1
    REAL :: nuT_1_B(4) ! average number of neutrons emitted per fission (nuT) for energy group 1
    REAL :: chiT_1_B(4) ! total neutron fission specture (chiT) for energy group 1

    !===================================================================================
    ! Cross section data for: 
    ! Energy group = 2, Test case = B (case = 1)
    !===================================================================================

    REAL :: sigTR_2_B(4) ! transport cross section (sigTR) for energy group 2 for test case B
    REAL :: sigIS_2_B(4) ! inelastic scatterring cross section (sigIS) for energy group 2
    REAL :: sigDS_2_B(4) ! downscattering cross section (sigDS) for energy group 2
    REAL :: sigA_2_B(4) ! absorption cross section (sigA) for energy group 2
    REAL :: sigF_2_B(4) ! fission cross section (sigF) for energy group 2
    REAL :: nuT_2_B(4) ! average number of neutrons emitted per fission (nuT) for energy group 2
    REAL :: chiT_2_B(4) ! total neutron fission specture (chiT) for energy group 2

    !===================================================================================
    ! Cross section data for: 
    ! Energy group = 1, Test case = C (case = 2)
    !===================================================================================

    REAL :: sigTR_1_C(4) ! transport cross section (sigTR) for energy group 1 for test case C
    REAL :: sigIS_1_C(4) ! inelastic scatterring cross section (sigIS) for energy group 1
    REAL :: sigDS_1_C(4) ! downscattering cross section (sigDS) for energy group 1 
    REAL :: sigA_1_C(4) ! absorption cross section (sigA) for energy group 1
    REAL :: sigF_1_C(4) ! fission cross section (sigF) for energy group 1
    REAL :: nuT_1_C(4) ! average number of neutrons emitted per fission (nuT) for energy group 1
    REAL :: chiT_1_C(4) ! total neutron fission specture (chiT) for energy group 1

    !===================================================================================
    ! Cross section data for: 
    ! Energy group = 2, Test case = C (case = 2)
    !===================================================================================

    REAL :: sigTR_2_C(4) ! transport cross section (sigTR) for energy group 2 for test case C
    REAL :: sigIS_2_C(4) ! inelastic scatterring cross section (sigIS) for energy group 2
    REAL :: sigDS_2_C(4) ! downscattering cross section (sigDS) for energy group 2
    REAL :: sigA_2_C(4) ! absorption cross section (sigA) for energy group 2
    REAL :: sigF_2_C(4) ! fission cross section (sigF) for energy group 2
    REAL :: nuT_2_C(4) ! average number of neutrons emitted per fission (nuT) for energy group 2
    REAL :: chiT_2_C(4) ! total neutron fission specture (chiT) for energy group 2

    !===================================================================================
    ! Cross section data mapped by (case, energy group, material column)
    !
    ! After reading the XSData block, we map each column to a specific material.
    ! Column 1 -> matID = 0 (UO2), column 2 -> matID = 1 (MOX), etc.
    !
    ! Access pattern (example):
    !   xs_sigTR(test_case+1, 1, matID+1)  ! SigTR for selected test_case, group 1, material matID
    !===================================================================================
    REAL, ALLOCATABLE :: xs_sigTR(:,:,:)
    REAL, ALLOCATABLE :: xs_sigIS(:,:,:)
    REAL, ALLOCATABLE :: xs_sigDS(:,:,:)
    REAL, ALLOCATABLE :: xs_sigA(:,:,:)
    REAL, ALLOCATABLE :: xs_sigF(:,:,:)
    REAL, ALLOCATABLE :: xs_nuT(:,:,:)
    REAL, ALLOCATABLE :: xs_chiT(:,:,:)

CONTAINS

    !===================================================================================
    ! Helper: parse integers from the RHS of a line like "MatID = 2 0 2 ..."
    ! Returns up to SIZE(values) integers in values(1:n_values).
    !===================================================================================
    SUBROUTINE parse_int_values(rhs, values, n_values)
        CHARACTER(LEN=*), INTENT(IN) :: rhs
        INTEGER, INTENT(OUT) :: values(:)
        INTEGER, INTENT(OUT) :: n_values

        INTEGER :: i, n, token_start, token_end, ios
        CHARACTER(LEN=64) :: token
        LOGICAL :: in_token

        values = 0
        n_values = 0
        n = LEN_TRIM(rhs)
        in_token = .FALSE.
        token_start = 0

        DO i = 1, n + 1
            IF (i <= n) THEN
                IF (rhs(i:i) /= ' ' .AND. rhs(i:i) /= CHAR(9)) THEN
                    IF (.NOT. in_token) THEN
                        in_token = .TRUE.
                        token_start = i
                    END IF
                ELSE
                    IF (in_token) THEN
                        token_end = i - 1
                        token = rhs(token_start:token_end)
                        READ(token, *, IOSTAT=ios) values(n_values + 1)
                        IF (ios == 0) n_values = n_values + 1
                        IF (n_values >= SIZE(values)) RETURN
                        in_token = .FALSE.
                    END IF
                END IF
            ELSE
                IF (in_token) THEN
                    token_end = n
                    token = rhs(token_start:token_end)
                    READ(token, *, IOSTAT=ios) values(n_values + 1)
                    IF (ios == 0) n_values = n_values + 1
                END IF
            END IF
        END DO
    END SUBROUTINE parse_int_values

    !===================================================================================
    ! Helper: find '=' and return RHS (adjusted/trimmed)
    !===================================================================================
    SUBROUTINE get_rhs_after_equals(line, rhs)
        CHARACTER(LEN=*), INTENT(IN) :: line
        CHARACTER(LEN=*), INTENT(OUT) :: rhs
        INTEGER :: p

        rhs = ''
        p = INDEX(line, '=')
        IF (p > 0 .AND. p < LEN(line)) THEN
            rhs = ADJUSTL(line(p+1:))
        END IF
    END SUBROUTINE get_rhs_after_equals

    !===================================================================================
    ! Subroutine: read_input
    ! Now, I've heard of using Namelist which reads the input file in the native ForTRAN format
    ! But, I'm more used to doing a line by line read.
    !===================================================================================
    SUBROUTINE read_input(filename)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=256) :: line
        CHARACTER(LEN=256) :: adj
        CHARACTER(LEN=32)  :: keyword
        CHARACTER(LEN=1)   :: equals_sign
        INTEGER :: ios
        INTEGER :: ios2, set_id, set_idx, max_len, nvals
        INTEGER, ALLOCATABLE :: counts(:)
        INTEGER, ALLOCATABLE :: fill_pos(:)
        CHARACTER(LEN=256) :: rhs
        INTEGER :: tmp_vals(2048)
        LOGICAL :: found_configsets
        INTEGER :: required_size
        INTEGER, ALLOCATABLE :: new_counts(:)
        INTEGER, ALLOCATABLE :: row_counts(:)
        INTEGER, ALLOCATABLE :: row_lens(:,:)
        INTEGER, ALLOCATABLE :: new_row_counts(:)
        INTEGER, ALLOCATABLE :: new_row_lens(:,:)
        INTEGER :: row_idx, max_rows, max_row_len
        INTEGER :: verify_xs_file_flag
        LOGICAL :: got_bound_l_7g, got_bound_r_7g

        OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='READ')

        ! Read Solution
        DO
            READ(10, '(A)') line
            IF (INDEX(line,"Solution") == 1) THEN
                READ(line, *) keyword, equals_sign, solution
                EXIT
            END IF
        END DO

        ! Read TestCase
        REWIND(10)
        DO
            READ(10, '(A)') line
            IF (INDEX(line,"TestCase") == 1) THEN
                READ(line, *) keyword, equals_sign, test_case
                EXIT
            END IF
        END DO

        ! Read other scalar inputs (one pass; Configs before Config). EnergyGroups must be
        ! read before boundaries so BoundL/BoundR can be allocated to length EnergyGroups.
        REWIND(10)
        DO
            READ(10, '(A)', IOSTAT=ios) line
            IF (ios /= 0) EXIT
            adj = ADJUSTL(line)
            IF (LEN_TRIM(adj) == 0) CYCLE
            IF (adj(1:1) == '#') CYCLE
            IF (INDEX(adj, "Configs") == 1) THEN
                READ(line, *) keyword, equals_sign, configs
            ELSE IF (INDEX(adj, "Config") == 1 .AND. INDEX(adj, "ConfigSets") == 0) THEN
                READ(line, *) keyword, equals_sign, config
            ELSE IF (INDEX(adj, "EnergyGroups") == 1) THEN
                READ(line, *) keyword, equals_sign, energy_groups
            ELSE IF (INDEX(adj, "MatTypes") == 1) THEN
                READ(line, *) keyword, equals_sign, mattypes
            ELSE IF (INDEX(adj, "Cases") == 1) THEN
                READ(line, *) keyword, equals_sign, cases
            ELSE IF (INDEX(adj, "Analk") == 1) THEN
                READ(line, *) keyword, equals_sign, analk
            ELSE IF (INDEX(adj, "solver") == 1) THEN
                READ(line, *) keyword, equals_sign, solver
            ELSE IF (INDEX(adj, "Generations") == 1) THEN
                READ(line, *) keyword, equals_sign, generations
            ELSE IF (INDEX(adj, "Histories") == 1) THEN
                READ(line, *) keyword, equals_sign, histories
            ELSE IF (INDEX(adj, "VerifyXSCase") == 1) THEN
                READ(line, *) keyword, equals_sign, verify_xs_case_choice
            ELSE IF (INDEX(adj, "VerifyXS") == 1) THEN
                READ(line, *) keyword, equals_sign, verify_xs_file_flag
                verify_xs_request = (verify_xs_file_flag /= 0)
            ELSE IF (INDEX(adj, "Skip") == 1) THEN
                READ(line, *) keyword, equals_sign, skip
            ELSE IF (INDEX(adj, "NumAss") == 1) THEN
                READ(line, *) keyword, equals_sign, num_assemblies
            ELSE IF (INDEX(adj, "NumRods") == 1) THEN
                READ(line, *) keyword, equals_sign, num_rods
            ELSE IF (INDEX(adj, "RodPitch") == 1) THEN
                READ(line, *) keyword, equals_sign, rod_pitch
            ELSE IF (INDEX(adj, "RodDia") == 1) THEN
                READ(line, *) keyword, equals_sign, rod_diameter
            ELSE IF (INDEX(adj, "MPFR") == 1) THEN
                READ(line, *) keyword, equals_sign, mpfr
            ELSE IF (INDEX(adj, "MPWR") == 1) THEN
                READ(line, *) keyword, equals_sign, mpwr
            END IF
        END DO

        IF (energy_groups <= 0) energy_groups = 2

        IF (ALLOCATED(bound_l)) DEALLOCATE(bound_l)
        IF (ALLOCATED(bound_r)) DEALLOCATE(bound_r)
        IF (ALLOCATED(bound_l_7g)) DEALLOCATE(bound_l_7g)
        IF (ALLOCATED(bound_r_7g)) DEALLOCATE(bound_r_7g)
        ALLOCATE(bound_l(energy_groups))
        ALLOCATE(bound_r(energy_groups))
        ALLOCATE(bound_l_7g(7))
        ALLOCATE(bound_r_7g(7))
        bound_l = 1.0
        bound_r = 1.0
        bound_l_7g = 1.0
        bound_r_7g = 1.0
        got_bound_l_7g = .FALSE.
        got_bound_r_7g = .FALSE.
        REWIND(10)
        DO
            READ(10, '(A)', IOSTAT=ios) line
            IF (ios /= 0) EXIT
            adj = ADJUSTL(line)
            IF (LEN_TRIM(adj) == 0) CYCLE
            IF (adj(1:1) == '#') CYCLE
            IF (INDEX(adj, "BoundL_7G") == 1) THEN
                READ(line, *, IOSTAT=ios) keyword, equals_sign, bound_l_7g
                IF (ios /= 0) THEN
                    PRINT *, "ERROR: BoundL_7G must provide 7 numeric values."
                    STOP 1
                END IF
                got_bound_l_7g = .TRUE.
            ELSE IF (INDEX(adj, "BoundR_7G") == 1) THEN
                READ(line, *, IOSTAT=ios) keyword, equals_sign, bound_r_7g
                IF (ios /= 0) THEN
                    PRINT *, "ERROR: BoundR_7G must provide 7 numeric values."
                    STOP 1
                END IF
                got_bound_r_7g = .TRUE.
            ELSE IF (INDEX(adj, "BoundL") == 1) THEN
                READ(line, *) keyword, equals_sign, bound_l
            ELSE IF (INDEX(adj, "BoundR") == 1) THEN
                READ(line, *) keyword, equals_sign, bound_r
            END IF
        END DO

        ! If explicit 7G boundaries were omitted, derive defaults from BoundL/BoundR.
        IF (.NOT. got_bound_l_7g) THEN
            IF (SIZE(bound_l) >= 1) bound_l_7g(1) = bound_l(1)
            IF (SIZE(bound_l) >= 2) THEN
                bound_l_7g(2:7) = bound_l(2)
            ELSE
                bound_l_7g(2:7) = bound_l_7g(1)
            END IF
        END IF
        IF (.NOT. got_bound_r_7g) THEN
            IF (SIZE(bound_r) >= 1) bound_r_7g(1) = bound_r(1)
            IF (SIZE(bound_r) >= 2) THEN
                bound_r_7g(2:7) = bound_r(2)
            ELSE
                bound_r_7g(2:7) = bound_r_7g(1)
            END IF
        END IF

        IF (energy_groups == 7) THEN
            PRINT *, "EnergyGroups = 7: use the seven-group MC in seven_groups/ (see seven_groups/README.txt)."
            STOP 1
        ELSE IF (energy_groups /= 2) THEN
            PRINT *, "ERROR: EnergyGroups must be 2 for this project (main + root monte_carlo). Got:", energy_groups
            STOP 1
        END IF

        ! Read XS Data
        REWIND(10)
        DO
            READ(10, '(A)') line
            IF (INDEX(line,"XSData") > 0) THEN
                DO
                    READ(10, '(A)') line
                    IF (INDEX(line,"END") > 0) EXIT
                    IF (line(1:1) == '#') CYCLE

                    IF (INDEX(line,"case = 0") > 0) THEN
                        DO
                            READ(10, '(A)') line
                            IF (INDEX(line,"SigTR") > 0) EXIT
                        END DO
                        READ(line, *) keyword, equals_sign, sigTR_1_A, sigTR_2_A
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigIS_1_A, sigIS_2_A
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigDS_1_A, sigDS_2_A
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigA_1_A, sigA_2_A
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigF_1_A, sigF_2_A
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, nuT_1_A, nuT_2_A
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, chiT_1_A, chiT_2_A
                    END IF

                    IF (INDEX(line,"case = 1") > 0) THEN
                        DO
                            READ(10, '(A)') line
                            IF (INDEX(line,"SigTR") > 0) EXIT
                        END DO
                        READ(line, *) keyword, equals_sign, sigTR_1_B, sigTR_2_B
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigIS_1_B, sigIS_2_B
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigDS_1_B, sigDS_2_B
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigA_1_B, sigA_2_B
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigF_1_B, sigF_2_B
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, nuT_1_B, nuT_2_B
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, chiT_1_B, chiT_2_B
                    END IF

                    IF (INDEX(line,"case = 2") > 0) THEN
                        DO
                            READ(10, '(A)') line
                            IF (INDEX(line,"SigTR") > 0) EXIT
                        END DO
                        READ(line, *) keyword, equals_sign, sigTR_1_C, sigTR_2_C
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigIS_1_C, sigIS_2_C
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigDS_1_C, sigDS_2_C
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigA_1_C, sigA_2_C
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, sigF_1_C, sigF_2_C
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, nuT_1_C, nuT_2_C
                        READ(10, '(A)') line
                        READ(line, *) keyword, equals_sign, chiT_1_C, chiT_2_C
                    END IF
                END DO
                EXIT
            END IF
        END DO

        ! Assign each column of the cross-section arrays to a specific material 
        ! The structure of each Case is such that: 
        ! ---- Case A ----
        ! SigTR Group 1 =  0.248960003      0.247979999      0.238230005      0.100000001    
        ! SigIS Group 1 =  0.226980001      0.224820003      0.219040006       0.00000000    
        ! SigDS Group 1 =   1.10200001E-02   8.86000041E-03   2.77600009E-02   0.00000000    
        ! SigA  Group 1 =   1.16299996E-02   1.49499997E-02   1.59799997E-02  0.100000001    
        ! SigF  Group 1 =   3.94000020E-03   4.00000019E-03   0.00000000       0.00000000    
        ! nuT   Group 1 =   2.46950006       2.84100008       0.00000000       0.00000000    
        ! chiT  Group 1 =   1.00000000       1.00000000       0.00000000       0.00000000    
        ! SigTR Group 2 =  0.885249972       1.12071002       1.22518003       1.10000002    
        ! [...] 
        ! Where column one is the cross-sectional data for matID = 0 or UO2, column two is for matID = 1 or MOX, etc.

        !===================================================================================
        ! Build a (case, energy group, material column) mapping so we can directly index by
        ! matID from the geometry. This retains the "column = material" structure from the
        ! input file, but stores it in a way the solver can use.
        !===================================================================================

        !===================================================================================
        !
        ! This structure will allow us to identify the cross-sections for each energy group for any given case per materal
        ! To find, for example the transport cross section for case 1, energy group 1, and material 2 (MOX):
        ! xs_sigTR(1, 1, 2)
        ! 
        ! So, the first paramter is the case, the second is the energy group, and the third is the material
        !
        !===================================================================================

        IF (mattypes <= 0) THEN
            PRINT *, "ERROR: MatTypes must be > 0 to map XSData columns."
            STOP 1
        END IF
        IF (cases <= 0) THEN
            PRINT *, "ERROR: Cases must be > 0 to map XSData cases."
            STOP 1
        END IF

        IF (ALLOCATED(xs_sigTR)) DEALLOCATE(xs_sigTR)
        IF (ALLOCATED(xs_sigIS)) DEALLOCATE(xs_sigIS)
        IF (ALLOCATED(xs_sigDS)) DEALLOCATE(xs_sigDS)
        IF (ALLOCATED(xs_sigA))  DEALLOCATE(xs_sigA)
        IF (ALLOCATED(xs_sigF))  DEALLOCATE(xs_sigF)
        IF (ALLOCATED(xs_nuT))   DEALLOCATE(xs_nuT)
        IF (ALLOCATED(xs_chiT))  DEALLOCATE(xs_chiT)

        ALLOCATE(xs_sigTR(cases, energy_groups, mattypes))
        ALLOCATE(xs_sigIS(cases, energy_groups, mattypes))
        ALLOCATE(xs_sigDS(cases, energy_groups, mattypes))
        ALLOCATE(xs_sigA(cases,  energy_groups, mattypes))
        ALLOCATE(xs_sigF(cases,  energy_groups, mattypes))
        ALLOCATE(xs_nuT(cases,   energy_groups, mattypes))
        ALLOCATE(xs_chiT(cases,  energy_groups, mattypes))

        xs_sigTR = 0.0
        xs_sigIS = 0.0
        xs_sigDS = 0.0
        xs_sigA  = 0.0
        xs_sigF  = 0.0
        xs_nuT   = 0.0
        xs_chiT  = 0.0

        ! Case A (case = 0 -> index 1)
        IF (cases >= 1) THEN
            xs_sigTR(1,1,:) = sigTR_1_A(1:mattypes) ! transport cross section array for case A, EG 1, 
            xs_sigIS(1,1,:) = sigIS_1_A(1:mattypes)
            xs_sigDS(1,1,:) = sigDS_1_A(1:mattypes)
            xs_sigA (1,1,:) = sigA_1_A(1:mattypes)
            xs_sigF (1,1,:) = sigF_1_A(1:mattypes)
            xs_nuT  (1,1,:) = nuT_1_A(1:mattypes)
            xs_chiT (1,1,:) = chiT_1_A(1:mattypes)

            xs_sigTR(1,2,:) = sigTR_2_A(1:mattypes)
            xs_sigIS(1,2,:) = sigIS_2_A(1:mattypes)
            xs_sigDS(1,2,:) = sigDS_2_A(1:mattypes)
            xs_sigA (1,2,:) = sigA_2_A(1:mattypes)
            xs_sigF (1,2,:) = sigF_2_A(1:mattypes)
            xs_nuT  (1,2,:) = nuT_2_A(1:mattypes)
            xs_chiT (1,2,:) = chiT_2_A(1:mattypes)
        END IF

        ! Case B (case = 1 -> index 2)
        IF (cases >= 2) THEN
            xs_sigTR(2,1,:) = sigTR_1_B(1:mattypes)
            xs_sigIS(2,1,:) = sigIS_1_B(1:mattypes)
            xs_sigDS(2,1,:) = sigDS_1_B(1:mattypes)
            xs_sigA (2,1,:) = sigA_1_B(1:mattypes)
            xs_sigF (2,1,:) = sigF_1_B(1:mattypes)
            xs_nuT  (2,1,:) = nuT_1_B(1:mattypes)
            xs_chiT (2,1,:) = chiT_1_B(1:mattypes)

            xs_sigTR(2,2,:) = sigTR_2_B(1:mattypes)
            xs_sigIS(2,2,:) = sigIS_2_B(1:mattypes)
            xs_sigDS(2,2,:) = sigDS_2_B(1:mattypes)
            xs_sigA (2,2,:) = sigA_2_B(1:mattypes)
            xs_sigF (2,2,:) = sigF_2_B(1:mattypes)
            xs_nuT  (2,2,:) = nuT_2_B(1:mattypes)
            xs_chiT (2,2,:) = chiT_2_B(1:mattypes)
        END IF

        ! Case C (case = 2 -> index 3)
        IF (cases >= 3) THEN
            xs_sigTR(3,1,:) = sigTR_1_C(1:mattypes)
            xs_sigIS(3,1,:) = sigIS_1_C(1:mattypes)
            xs_sigDS(3,1,:) = sigDS_1_C(1:mattypes)
            xs_sigA (3,1,:) = sigA_1_C(1:mattypes)
            xs_sigF (3,1,:) = sigF_1_C(1:mattypes)
            xs_nuT  (3,1,:) = nuT_1_C(1:mattypes)
            xs_chiT (3,1,:) = chiT_1_C(1:mattypes)

            xs_sigTR(3,2,:) = sigTR_2_C(1:mattypes)
            xs_sigIS(3,2,:) = sigIS_2_C(1:mattypes)
            xs_sigDS(3,2,:) = sigDS_2_C(1:mattypes)
            xs_sigA (3,2,:) = sigA_2_C(1:mattypes)
            xs_sigF (3,2,:) = sigF_2_C(1:mattypes)
            xs_nuT  (3,2,:) = nuT_2_C(1:mattypes)
            xs_chiT (3,2,:) = chiT_2_C(1:mattypes)
        END IF

        ! Example (UO2, matID=0) for verification:
        ! xs_sigTR(test_case+1, 1, 0+1) is SigTR for group 1, selected case, UO2 column

        ! Read ConfigSets block (Set = N, MatID = ... lines, END terminator)
        ! Modular behavior: determine number of sets from the block itself.
        found_configsets = .FALSE.
        IF (ALLOCATED(counts)) DEALLOCATE(counts)
        ALLOCATE(counts(0))
        counts = 0
        IF (ALLOCATED(row_counts)) DEALLOCATE(row_counts)
        ALLOCATE(row_counts(0))
        row_counts = 0
        IF (ALLOCATED(row_lens)) DEALLOCATE(row_lens)
        ALLOCATE(row_lens(0,0))

        ! Pass 1: determine how many MatID ints per set (flattened across MatID lines)
        REWIND(10)
        DO
            READ(10, '(A)', IOSTAT=ios2) line
            IF (ios2 /= 0) EXIT
            adj = ADJUSTL(line)
            IF (LEN_TRIM(adj) == 0) CYCLE
            IF (adj(1:1) == '#') CYCLE
            IF (INDEX(adj, "ConfigSets") == 1) THEN
                found_configsets = .TRUE.
                EXIT
            END IF
        END DO

        set_idx = 0
        IF (found_configsets) THEN
            DO
                READ(10, '(A)', IOSTAT=ios2) line
                IF (ios2 /= 0) EXIT
                adj = ADJUSTL(line)
                IF (LEN_TRIM(adj) == 0) CYCLE
                IF (adj(1:1) == '#') CYCLE
                IF (INDEX(adj, "END") > 0) EXIT

                IF (INDEX(adj, "Set") == 1) THEN
                    ! Set IDs in input start at 0; map to 1-based indexing internally
                    READ(line, *, IOSTAT=ios2) keyword, equals_sign, set_id
                    IF (ios2 == 0 .AND. set_id >= 0) THEN
                        set_idx = set_id + 1
                        required_size = set_idx
                        IF (SIZE(counts) < required_size) THEN
                            IF (ALLOCATED(new_counts)) DEALLOCATE(new_counts)
                            ALLOCATE(new_counts(required_size))
                            new_counts = 0
                            IF (SIZE(counts) > 0) new_counts(1:SIZE(counts)) = counts
                            CALL MOVE_ALLOC(new_counts, counts)

                            IF (ALLOCATED(new_row_counts)) DEALLOCATE(new_row_counts)
                            ALLOCATE(new_row_counts(required_size))
                            new_row_counts = 0
                            IF (SIZE(row_counts) > 0) new_row_counts(1:SIZE(row_counts)) = row_counts
                            CALL MOVE_ALLOC(new_row_counts, row_counts)
                        END IF
                    ELSE
                        set_idx = 0
                    END IF
                ELSE IF (INDEX(adj, "MatID") == 1 .AND. set_idx >= 1 .AND. set_idx <= SIZE(counts)) THEN
                    CALL get_rhs_after_equals(line, rhs)
                    CALL parse_int_values(rhs, tmp_vals, nvals)
                    counts(set_idx) = counts(set_idx) + nvals

                    ! Track row breaks (one MatID line = one row)
                    row_counts(set_idx) = row_counts(set_idx) + 1
                    row_idx = row_counts(set_idx)
                    IF (row_idx == 1) THEN
                        ! First row seen for this set; lazily size row_lens to [sets, some rows]
                        ! We'll finalize sizes after pass 1 by copying into config_set_row_len.
                    END IF

                    ! Grow row_lens matrix if needed (rows dimension).
                    IF (SIZE(row_lens,1) < SIZE(counts) .OR. SIZE(row_lens,2) < row_idx) THEN
                        max_row_len = SIZE(row_lens,1)
                        max_rows = SIZE(row_lens,2)
                        IF (ALLOCATED(new_row_lens)) DEALLOCATE(new_row_lens)
                        ALLOCATE(new_row_lens(SIZE(counts), MAX(max_rows, row_idx)))
                        new_row_lens = 0
                        IF (max_row_len > 0 .AND. max_rows > 0) THEN
                            new_row_lens(1:max_row_len,1:max_rows) = row_lens
                        END IF
                        CALL MOVE_ALLOC(new_row_lens, row_lens)
                    END IF
                    row_lens(set_idx, row_idx) = nvals
                END IF
            END DO
        END IF

        IF (found_configsets .AND. SIZE(counts) > 0) THEN
            max_len = MAXVAL(counts)
            IF (ALLOCATED(config_set_len)) DEALLOCATE(config_set_len)
            IF (ALLOCATED(config_set_matid)) DEALLOCATE(config_set_matid)
            ALLOCATE(config_set_len(SIZE(counts)))
            config_set_len = counts

            ! Override configs so downstream code automatically sees the new sets
            configs = SIZE(counts)

            ! Store row structure: number of MatID lines and their lengths
            IF (ALLOCATED(config_set_nrows)) DEALLOCATE(config_set_nrows)
            IF (ALLOCATED(config_set_row_len)) DEALLOCATE(config_set_row_len)
            ALLOCATE(config_set_nrows(configs))
            config_set_nrows = row_counts
            max_rows = MAXVAL(config_set_nrows)
            IF (max_rows > 0) THEN
                ALLOCATE(config_set_row_len(configs, max_rows))
                config_set_row_len = 0
                IF (SIZE(row_lens,1) >= configs .AND. SIZE(row_lens,2) >= max_rows) THEN
                    config_set_row_len(:, :) = row_lens(:, 1:max_rows)
                END IF
            END IF

            IF (max_len > 0) THEN
                ALLOCATE(config_set_matid(configs, max_len))
                config_set_matid = -1
            END IF

            ! Pass 2: fill config_set_matid
            IF (max_len > 0) THEN
                IF (ALLOCATED(fill_pos)) DEALLOCATE(fill_pos)
                ALLOCATE(fill_pos(configs))
                fill_pos = 0

                REWIND(10)
                DO
                    READ(10, '(A)', IOSTAT=ios2) line
                    IF (ios2 /= 0) EXIT
                    adj = ADJUSTL(line)
                    IF (LEN_TRIM(adj) == 0) CYCLE
                    IF (adj(1:1) == '#') CYCLE
                    IF (INDEX(adj, "ConfigSets") == 1) EXIT
                END DO

                set_idx = 0
                DO
                    READ(10, '(A)', IOSTAT=ios2) line
                    IF (ios2 /= 0) EXIT
                    adj = ADJUSTL(line)
                    IF (LEN_TRIM(adj) == 0) CYCLE
                    IF (adj(1:1) == '#') CYCLE
                    IF (INDEX(adj, "END") > 0) EXIT

                    IF (INDEX(adj, "Set") == 1) THEN
                        READ(line, *, IOSTAT=ios2) keyword, equals_sign, set_id
                        IF (ios2 == 0 .AND. set_id >= 0) THEN
                            set_idx = set_id + 1
                        ELSE
                            set_idx = 0
                        END IF
                    ELSE IF (INDEX(adj, "MatID") == 1 .AND. set_idx >= 1 .AND. set_idx <= configs) THEN
                        CALL get_rhs_after_equals(line, rhs)
                        CALL parse_int_values(rhs, tmp_vals, nvals)
                        IF (nvals > 0) THEN
                            IF (fill_pos(set_idx) + nvals <= max_len) THEN
                                config_set_matid(set_idx, fill_pos(set_idx)+1:fill_pos(set_idx)+nvals) = tmp_vals(1:nvals)
                                fill_pos(set_idx) = fill_pos(set_idx) + nvals
                            END IF
                        END IF
                    END IF
                END DO
            END IF
        END IF

        IF (ALLOCATED(fill_pos)) DEALLOCATE(fill_pos)
        IF (ALLOCATED(counts)) DEALLOCATE(counts)
        IF (ALLOCATED(new_counts)) DEALLOCATE(new_counts)
        IF (ALLOCATED(row_counts)) DEALLOCATE(row_counts)
        IF (ALLOCATED(row_lens)) DEALLOCATE(row_lens)
        IF (ALLOCATED(new_row_counts)) DEALLOCATE(new_row_counts)
        IF (ALLOCATED(new_row_lens)) DEALLOCATE(new_row_lens)

        CLOSE(10)
    END SUBROUTINE read_input

    !===================================================================================
    ! Subroutine: print_input
    !===================================================================================
    SUBROUTINE print_input()
        PRINT *, "================ INPUT DATA ================"

        PRINT *, "Solution =", solution
        PRINT *, "Test Case =", test_case
        PRINT *, "Config =", config
        PRINT *, "Analk =", analk
        PRINT *, "Cases =", cases
        PRINT *, "Configs =", configs
        PRINT *, "MatTypes =", mattypes
        PRINT *, "Energy Groups =", energy_groups
        PRINT *, "Solver =", solver
        PRINT *, "Generations =", generations
        PRINT *, "Histories =", histories
        PRINT *, "Skip =", skip
        PRINT *, "VerifyXS (0=no CSV, 1=write xs_mapping_verify_case_*.csv) =", verify_xs_request
        PRINT *, "VerifyXSCase (XS A=0,B=1,...; negative uses TestCase) =", verify_xs_case_choice
        PRINT *, "Number of Assemblies =", num_assemblies
        PRINT *, "Number of Rods =", num_rods
        PRINT *, "Rod Diameter =", rod_diameter
        PRINT *, "Rod Pitch =", rod_pitch
        PRINT *, "MPFR =", mpfr
        PRINT *, "MPWR =", mpwr

        PRINT *, "Bound Left =", bound_l
        PRINT *, "Bound Right =", bound_r
        PRINT *, "Bound Left 7G =", bound_l_7g
        PRINT *, "Bound Right 7G =", bound_r_7g

        PRINT *, "================ XS DATA ================"

        PRINT *, "---- Case A ----"
        PRINT *, "SigTR Group 1 =", sigTR_1_A
        PRINT *, "SigIS Group 1 =", sigIS_1_A
        PRINT *, "SigDS Group 1 =", sigDS_1_A
        PRINT *, "SigA  Group 1 =", sigA_1_A
        PRINT *, "SigF  Group 1 =", sigF_1_A
        PRINT *, "nuT   Group 1 =", nuT_1_A
        PRINT *, "chiT  Group 1 =", chiT_1_A

        PRINT *, "SigTR Group 2 =", sigTR_2_A
        PRINT *, "SigIS Group 2 =", sigIS_2_A
        PRINT *, "SigDS Group 2 =", sigDS_2_A
        PRINT *, "SigA  Group 2 =", sigA_2_A
        PRINT *, "SigF  Group 2 =", sigF_2_A
        PRINT *, "nuT   Group 2 =", nuT_2_A
        PRINT *, "chiT  Group 2 =", chiT_2_A

        PRINT *, "---- Case B ----"
        PRINT *, "SigTR Group 1 =", sigTR_1_B
        PRINT *, "SigIS Group 1 =", sigIS_1_B
        PRINT *, "SigDS Group 1 =", sigDS_1_B
        PRINT *, "SigA  Group 1 =", sigA_1_B
        PRINT *, "SigF  Group 1 =", sigF_1_B
        PRINT *, "nuT   Group 1 =", nuT_1_B
        PRINT *, "chiT  Group 1 =", chiT_1_B

        PRINT *, "SigTR Group 2 =", sigTR_2_B
        PRINT *, "SigIS Group 2 =", sigIS_2_B
        PRINT *, "SigDS Group 2 =", sigDS_2_B
        PRINT *, "SigA  Group 2 =", sigA_2_B
        PRINT *, "SigF  Group 2 =", sigF_2_B
        PRINT *, "nuT   Group 2 =", nuT_2_B
        PRINT *, "chiT  Group 2 =", chiT_2_B

        PRINT *, "---- Case C ----"
        PRINT *, "SigTR Group 1 =", sigTR_1_C
        PRINT *, "SigIS Group 1 =", sigIS_1_C
        PRINT *, "SigDS Group 1 =", sigDS_1_C
        PRINT *, "SigA  Group 1 =", sigA_1_C
        PRINT *, "SigF  Group 1 =", sigF_1_C
        PRINT *, "nuT   Group 1 =", nuT_1_C
        PRINT *, "chiT  Group 1 =", chiT_1_C

        PRINT *, "SigTR Group 2 =", sigTR_2_C
        PRINT *, "SigIS Group 2 =", sigIS_2_C
        PRINT *, "SigDS Group 2 =", sigDS_2_C
        PRINT *, "SigA  Group 2 =", sigA_2_C
        PRINT *, "SigF  Group 2 =", sigF_2_C
        PRINT *, "nuT   Group 2 =", nuT_2_C
        PRINT *, "chiT  Group 2 =", chiT_2_C

        PRINT *, "=========================================="
    END SUBROUTINE print_input

END MODULE read_input_file