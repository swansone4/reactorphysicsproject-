Seven-group Monte Carlo (minimal guide)
======================================

What is here:

- `csv_xs_7g.f90`: reads the two CSV libraries into 7-group XS arrays.
- `monte_carlo_solver_7g.f90`: runs the 1D Monte Carlo (same style as main MC, but 7 groups).
- `plot_monte_carlo_flux_7g.py`: makes the 7-group plots from solver CSV outputs.
- `macroscopic_cross_sections.csv` and `scattering_matrix.csv`: 7-group XS inputs.

Input source:

- Geometry and run controls still come from `input_file.txt` (`ConfigSets`, `Config`, `MPFR`, `MPWR`, `RodPitch`, `RodDia`, `Generations`, `Histories`, `Skip`, `BoundL_7G`, `BoundR_7G`).

Build + run from repo root:

    gfortran -O2 -o seven_groups/monte_carlo_solver_7g.exe read_input_file.f90 mesh_generation.f90 seven_groups/csv_xs_7g.f90 seven_groups/monte_carlo_solver_7g.f90
    seven_groups/monte_carlo_solver_7g.exe

Outputs:

- `seven_groups/outputs/monte_carlo/summary.txt`
- `seven_groups/outputs/monte_carlo/flux_cells.csv`
- `seven_groups/outputs/monte_carlo/flux_fundamental.csv`
- `seven_groups/outputs/monte_carlo/mc_flux_g1_to_g7.png` (all-group flux)
- `seven_groups/outputs/monte_carlo/mc_cell_center_flux_g1_to_g7.png` (cell-center view)
- `seven_groups/outputs/monte_carlo/mc_cell_edge_current_g1_to_g7.png` (cell-edge current from Fick's law with a simple finite difference between cell centers)
