# ENU 6106 — reactor physics project

## What you need

- `gfortran` (GNU Fortran)
- Python 3 with `matplotlib` if you run the plotting scripts

## Main program (mesh + XS check)

Build and run from the repo root (needs `input_file.txt` in the working directory):

```bash
gfortran -o main_program read_input_file.f90 mesh_generation.f90 main.f90
./main_program
```

- Set `VerifyXS = 1` in `input_file.txt` to write `outputs/xs_mapping/xs_mapping_verify_case_*.csv`.
- With `plot_visuals` true in `mesh_generation.f90`, mesh CSVs go to `outputs/mesh_dumps/` and `scripts/mesh_visualization.py` is called for PNGs under `outputs/plots/mesh/`.

Optional XS-over-mesh figure (after you have mesh dumps):

```bash
python3 scripts/XS_mesh_mapped_visualization.py
```

## Assembly Monte Carlo

Uses the same input reader and `Config` / `TestCase` as above.

```bash
gfortran -o monte_carlo_solver read_input_file.f90 mesh_generation.f90 monte_carlo_solver.f90
./monte_carlo_solver
```

Summary: `outputs/monte_carlo/summary.txt`.

## Bare slab benchmark (homogeneous slab, vacuum BCs)

```bash
cd bare_slab_benchmark
gfortran -O2 -o monte_carlo_bare_slab_benchmark monte_carlo_bare_slab_benchmark.f90
./monte_carlo_bare_slab_benchmark
python3 plot_bare_slab_benchmark.py
```

Outputs appear under `bare_slab_benchmark/outputs/` (created when you run).

## Layout

| `read_input_file.f90` | Parses `input_file.txt` |
| `mesh_generation.f90` | Config sets, XS on mesh, mesh dump, verify CSV |
| `main.f90` | Driver for mesh / XS echo |
| `monte_carlo_solver.f90` | 1D MC for the assembly geometry |
| `bare_slab_benchmark/` | Standalone slab vs diffusion reference |
| `scripts/` | Python visualization helpers |
| `images/fuel_assembly_config.png` | Schematic for the assembly layout |
