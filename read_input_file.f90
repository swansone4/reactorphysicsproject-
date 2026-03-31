MODULE read_input_file
    IMPLICIT NONE

    !===================================================================================
    ! Basic input file variables and parameters
    !===================================================================================

    ! solution variable 
    INTEGER :: solution ! method for solving the problem (0 = Finite difference), (1 = Monte Carlo) 
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
    ! energy groups variable
    INTEGER :: energy_groups
    ! solver variable 
    INTEGER :: solver
    ! generations variable 
    INTEGER :: generations
    ! histories variable 
    INTEGER :: histories 
    ! skip variable 
    INTEGER :: skip 
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
    ! left boundary condition value 
    REAL :: bound_l(2)
    ! right boundary condition value    
    REAL :: bound_r(2)

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

CONTAINS

    !===================================================================================
    ! Subroutine: read_input
    !===================================================================================
    SUBROUTINE read_input(filename)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=256) :: line
        CHARACTER(LEN=256) :: adj
        CHARACTER(LEN=32)  :: keyword
        CHARACTER(LEN=1)   :: equals_sign
        INTEGER :: ios

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

        ! Read Boundaries
        REWIND(10)
        DO
            READ(10, '(A)') line
            IF (INDEX(line,"BoundL") == 1) THEN
                READ(line, *) keyword, equals_sign, bound_l
            END IF
            IF (INDEX(line,"BoundR") == 1) THEN
                READ(line, *) keyword, equals_sign, bound_r
                EXIT
            END IF
        END DO

        ! Read other scalar inputs (one pass; Configs before Config)
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
        PRINT *, "Number of Assemblies =", num_assemblies
        PRINT *, "Number of Rods =", num_rods
        PRINT *, "Rod Diameter =", rod_diameter
        PRINT *, "Rod Pitch =", rod_pitch
        PRINT *, "MPFR =", mpfr
        PRINT *, "MPWR =", mpwr

        PRINT *, "Bound Left =", bound_l
        PRINT *, "Bound Right =", bound_r

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