PROGRAM main
! this is the main program that runs the entire code. It reads the input file and prints the data.
    USE read_input_file
    IMPLICIT NONE

    ! boolean to ask if the user wants to print the input data
    LOGICAL :: print_input_data
    print_input_data = .FALSE.

    !===================================================================================
    ! Read input file and print data
    !===================================================================================
    IF (print_input_data) THEN
        CALL read_input('input_file.txt')
        CALL print_input()
    ELSE ! if the user has print_input_data to false, let them know 
        PRINT *, 'Input data will not be printed. Set print_input_data to true in main.f90 to print the input data.'
    END IF
    
END PROGRAM main