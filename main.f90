PROGRAM main
    USE read_input_file
    USE mesh_generation
    IMPLICIT NONE

    LOGICAL :: print_input_data
    INTEGER :: set_num

    print_input_data = .TRUE.

    ! Turn off mesh dumps / Python plots by setting plot_visuals = .FALSE. in mesh_generation if you want.
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
    
    CALL init_configurations()

    DO set_num = 1, configs
        CALL map_XS_to_mesh(set_num)
        CALL generate_mesh(set_num)
        CALL print_mesh(set_num)

        IF (plot_visuals) THEN
            CALL write_mesh_dump(set_num)
        END IF
    END DO 

    IF (plot_visuals) THEN
        CALL run_mesh_visualization()
    END IF
END PROGRAM main