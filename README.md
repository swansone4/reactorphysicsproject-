# FDDS input reader (Fortran)

Small program that reads `input_file.txt` (main settings, boundaries, and `XSData` cross sections) and prints the values loaded into module `read_input_file`.

## Requirements

- [GNU Fortran](https://gcc.gnu.org/fortran/) (`gfortran`) on your PATH.

## Build

From this directory:

```bash
gfortran -o main_program read_input_file.f90 main.f90
```

Compile the module source before `main.f90` so the `read_input_file` module is available to the main program.

## Run

```bash
./main_program
```

The program opens `input_file.txt` in the **current working directory**. Run the command from the project folder (or pass a path only if you change the code to accept a different filename; the default in `main.f90` is `input_file.txt`).

## What it does

1. **`main.f90`** sets `print_input_data` and calls `read_input('input_file.txt')` then `print_input()`.
2. **`read_input_file.f90`** declares the input variables and implements:
   - Line-based reads for `Solution`, `TestCase`, `BoundL`, `BoundR`, and the scalar block (`Config`, `Configs`, `Cases`, etc.—see the source for the exact names).
   - Parsing of the `XSData` … `END` block: for each `case = 0/1/2`, lines are advanced until a `SigTR` data line (skipping comment lines that sit between the case label and the data).
   - **`Config`** is distinguished from **`ConfigSets`** so the geometry block title does not get parsed as `Config`.

Cross-section arrays are fixed length `(4)` per energy group/material column, matching the sample input layout.
