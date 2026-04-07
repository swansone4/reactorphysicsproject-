PROGRAM main
! this is the main program that runs the entire code. It reads the input file and prints the data.
    USE read_input_file
    USE mesh_generation
    IMPLICIT NONE

    ! boolean to ask if the user wants to print the input data
    LOGICAL :: print_input_data
    INTEGER :: set_num
    print_input_data = .TRUE.

    ! Set plot_visuals in mesh_generation.f90 if desired.
    ! LOGICAL :: plot_visuals = .TRUE.

    !===================================================================================
    ! Read input file; optional XS verify CSV (input VerifyXS); optional echo of inputs
    !===================================================================================
    CALL read_input('input_file.txt')
    IF (verify_xs_request) THEN
        CALL verify_XS_mapping()
    ELSE
        PRINT *, "XS verify CSV not written: set VerifyXS = 1 in input_file.txt to create xs_mapping_verify_case_*.csv."
    END IF
    IF (print_input_data) THEN
        CALL print_input()
    ELSE
        PRINT *, 'Input echo skipped (print_input_data = .FALSE. in main.f90).'
    END IF
    
    ! Now, we will initialize the configuration sets and generate the mesh points 
    CALL init_configurations()

    ! loop over all of the configuration sets, which can be found in the input file at the bottom 
    DO set_num = 1, configs 
        ! Map XS columns to materials for this configuration set
        CALL map_XS_to_mesh(set_num)

        ! generate the mesh for this iteration of the loop
        CALL generate_mesh(set_num)

        ! print the configuration set to the terminal to make sure everything looks good 
        CALL print_mesh(set_num)

        IF (plot_visuals) THEN
            CALL write_mesh_dump(set_num)
        END IF
    END DO 

    IF (plot_visuals) THEN
        CALL run_mesh_visualization()
    END IF
END PROGRAM main