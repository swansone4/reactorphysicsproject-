# ENU 6106 — reactor physics project

Course project for a **1D fuel assembly chain**: mesh + XS mapping, **2-group MC**, **7-group MC** (CSVs in `seven_groups/`), **7-group diffusion** (FD solver in `finite_difference/`), a small **energy collapse / homogenization** script, and a **bare slab** MC benchmark. Anything that uses the lattice geometry reads **`input_file.txt`** from whatever directory you run the exe from (almost always the **repo root**).

Start with the folder table, then the numbered build/run sections.

---

## What you need installed

| Tool | Used for |
|------|----------|
| **GNU Fortran (`gfortran`)** | All Fortran drivers and solvers |
| **LAPACK + BLAS** (linked as `-llapack -lblas`) | `finite_difference/diffusion_solver_7g.f90` when built as the direct / LAPACK driver (see `finite_difference/README.txt`) |
| **Python 3** with **matplotlib** (and usually **numpy**) | Plotting and orchestration scripts under `scripts/`, `seven_groups/`, `finite_difference/`, `energy_collapse_homogenization/`, `bare_slab_benchmark/` |

**Windows:** If you use **MSYS2 MinGW**, install a `gfortran` toolchain and LAPACK/BLAS for that environment, then run the same `gfortran ... -llapack -lblas` commands from a shell where those libraries are on the link path. Executable names in this doc use `.\` or backslashes where that matches PowerShell; from Git Bash you can use `./` like Linux.

---

## How the repository is organized (navigation map)

**Which folder answers which question:**

| Location | Role |
|----------|------|
| **`input_file.txt`** | Single configuration file: geometry, materials, boundary flags, optional verification toggles. **Must be in the working directory** when you run any Fortran program that calls `read_input_file`. |
| **`read_input_file.f90`** | Parses `input_file.txt` into Fortran modules used everywhere. |
| **`mesh_generation.f90`** | Builds the 1D mesh, maps materials, optional mesh dumps and XS verification CSVs. |
| **`main.f90`** | Small driver: run this when you only want mesh / XS echo (no transport solve in this executable). |
| **`monte_carlo_solver.f90`** | **Two-group** Monte Carlo for the assembly chain (not the seven-group CSV library). |
| **`seven_groups/`** | **Seven-group** library: `macroscopic_cross_sections.csv`, `scattering_matrix.csv`, `csv_xs_7g.f90`, `monte_carlo_solver_7g.f90`, plotting scripts, and **`seven_groups/outputs/monte_carlo/`** for MC outputs. |
| **`finite_difference/`** | **Seven-group diffusion** (FD): `diffusion_solver_7g.f90`, `fd_solvers_1d.f90`, `plot_fd_flux_7g.py`, helper scripts, and **`finite_difference/outputs/`** for FD CSVs/plots. See **`finite_difference/README.txt`** for theory and build lines. |
| **`energy_collapse_homogenization/`** | Python pipeline that snapshots FD cases, collapses groups G1–4 and G5–7, homogenizes XS, calls Fortran **`solve_keff_2g_direct.f90`** for a direct 2-group \(k\), and writes plots/CSVs under **`energy_collapse_homogenization/outputs/`**. See **`energy_collapse_homogenization/README.txt`**. |
| **`bare_slab_benchmark/`** | Standalone homogeneous slab MC vs reference; self-contained under that folder. |
| **`scripts/`** | Shared Python helpers (e.g. mesh visualization, XS-over-mesh figures). |
| **`outputs/`** (repo root) | Default location for **main / mesh** dumps and related plots when driven from `mesh_generation.f90` (e.g. `outputs/mesh_dumps/`, `outputs/plots/mesh/`). |
| **`images/`** | Static figures such as the fuel assembly schematic. |
| **`additional_plots/`** | Extra figures produced outside the main output trees (e.g. MC \(k_\text{eff}\) vs generation); treat as **documentation / results**, not as inputs to the solvers. |

**Important distinction:** **`outputs/monte_carlo/`** (under the repo root) is associated with the **two-group** `monte_carlo_solver` workflow. **Seven-group** Monte Carlo writes under **`seven_groups/outputs/monte_carlo/`**. Do not mix the two when comparing summaries or flux files.

---

## The central configuration file: `input_file.txt`

- All assembly-style Fortran programs expect **`input_file.txt` in the current directory** (almost always the repo root).
- **`VerifyXS = 1`** triggers extra XS mapping verification CSVs under `outputs/xs_mapping/` (see comments in the input file and `mesh_generation.f90`).
- **`Config` / `ConfigSets` / material IDs** drive which fuel/water layout is built; the seven-group CSVs and FD solver use the same geometry as the mesh code. For detailed boundary and \(k\) interpretation notes for FD vs MC, read **`finite_difference/README.txt`** (especially the paragraph on **`k_diag_fission_over_removal_Rphi`** vs true \(k_\text{eff}\)).

---

## 1) Main program — mesh generation and XS check

**Purpose:** Build the mesh, optionally dump it, run optional XS verification — **no** full MC transport in this executable.

From the **repository root**:

```bash
gfortran -o main_program read_input_file.f90 mesh_generation.f90 main.f90
./main_program
```

On Windows (PowerShell), after build:

```powershell
.\main_program.exe
```

- With **`plot_visuals`** enabled in `mesh_generation.f90`, mesh CSVs go to **`outputs/mesh_dumps/`**, and **`scripts/mesh_visualization.py`** writes PNGs under **`outputs/plots/mesh/`**.
- XS-over-mesh figure (needs mesh dumps already on disk):

```bash
python scripts/XS_mesh_mapped_visualization.py
```

---

## 2) Two-group assembly Monte Carlo

**Purpose:** 1D MC for the assembly geometry using the **two-group** data path wired through `read_input_file` / `mesh_generation` (not the `seven_groups` CSV files).

From the **repository root**:

```bash
gfortran -o monte_carlo_solver read_input_file.f90 mesh_generation.f90 monte_carlo_solver.f90
./monte_carlo_solver
```

Summary and related outputs are written under **`outputs/monte_carlo/`** (for example **`summary.txt`**).

---

## 3) Seven-group assembly Monte Carlo

**Purpose:** Same geometry and `input_file.txt` controls, but **macroscopic cross sections and scattering** come from the CSVs in **`seven_groups/`**.

**Build** (from repo root):

```bash
gfortran -O2 -o seven_groups/monte_carlo_solver_7g.exe read_input_file.f90 mesh_generation.f90 seven_groups/csv_xs_7g.f90 seven_groups/monte_carlo_solver_7g.f90
```

**Run** (from repo root, so `input_file.txt` is found):

```bash
seven_groups/monte_carlo_solver_7g.exe
```

**Inputs:**

- `input_file.txt` (repo root)
- `seven_groups/macroscopic_cross_sections.csv`
- `seven_groups/scattering_matrix.csv`

**Typical outputs** (under **`seven_groups/outputs/monte_carlo/`**):

- `summary.txt`
- `flux_cells.csv`, `flux_fundamental.csv`
- `mc_flux_g1_to_g7.png`, `mc_cell_center_flux_g1_to_g7.png`, `mc_cell_edge_current_g1_to_g7.png`

**Plotting / sweeps:**

- `seven_groups/plot_monte_carlo_flux_7g.py` — flux plots from the MC CSVs under `seven_groups/outputs/monte_carlo/`; run from the repo root.
- `seven_groups/plot_Monte_Carlo_criticality_and_time.py` — from the repo root, sweeps **Generations** and **Histories** in `input_file.txt`, runs **`seven_groups/monte_carlo_solver_7g.exe`** each time, then writes **`seven_groups/outputs/monte_carlo/criticality_and_time_sweep.csv`** and **`criticality_and_time_vs_histories_times_generations.png`**. Edit the sweep lists inside the script if you want different points.

---

## 4) Seven-group finite-difference diffusion

**Purpose:** Lecture-style **multigroup diffusion** on the same 1D mesh and **seven-group** CSV library as in `seven_groups/`, with reflective/vacuum boundaries as documented in code and notes.

**Authoritative instructions** (build flags, solver modes `direct` / `jacobi` / `gs` / `sor`, output file names, and **how to interpret \(k\)** vs MC) live in:

**`finite_difference/README.txt`**

Typical **build** (from repo root, requires LAPACK/BLAS):

```bash
gfortran -O2 -o finite_difference/diffusion_solver_7g.exe read_input_file.f90 mesh_generation.f90 seven_groups/csv_xs_7g.f90 finite_difference/fd_solvers_1d.f90 finite_difference/diffusion_solver_7g.f90 -llapack -lblas
```

**Run** (from repo root):

```bash
finite_difference/diffusion_solver_7g.exe
finite_difference/diffusion_solver_7g.exe gs
finite_difference/diffusion_solver_7g.exe sor 1.3
```

**Outputs:** under **`finite_difference/outputs/`** (see `finite_difference/README.txt` for the exact list). The Fortran driver can shell out to **`finite_difference/plot_fd_flux_7g.py`** after a run when Python is on `PATH`.

**Which `.exe` the pipelines use:** `energy_collapse_homogenization/run_energy_collapse_homogenization.py` (and `finite_difference/energy_group_analysis/run_energy_group_analysis.py`) pick the first file that exists, in this order: **`diffusion_solver_7g_lapack.exe`**, **`diffusion_solver_7g_lapack_outputs.exe`**, **`diffusion_solver_7g.exe`**, all under **`finite_difference/`**.

---

## 5) Energy collapse and 2-group homogenization

**Purpose:** Automate **three material layouts** (MOX strip, UO2 strip, half-and-half), run the FD snapshot step, **collapse** G1–4 and G5–7, **homogenize** cross sections over the whole configuration, and solve a **2-group \(k_\text{eff}\)** in Fortran.

Read **`energy_collapse_homogenization/README.txt`** for definitions and the exact output tree.

**Run** (from repo root):

```bash
python energy_collapse_homogenization/run_energy_collapse_homogenization.py
```

The script orchestrates case edits to `input_file.txt`, calls the FD executable under **`finite_difference/`**, and compiles/runs **`energy_collapse_homogenization/solve_keff_2g_direct.f90`** as needed. Aggregated and per-case artifacts land under **`energy_collapse_homogenization/outputs/`**.

---

## 6) Bare slab benchmark

**Purpose:** Small standalone benchmark (homogeneous slab, vacuum BCs) with its own Fortran source and Python plot script — **does not** use `input_file.txt`.

```bash
cd bare_slab_benchmark
gfortran -O2 -o monte_carlo_bare_slab_benchmark monte_carlo_bare_slab_benchmark.f90
./monte_carlo_bare_slab_benchmark
python plot_bare_slab_benchmark.py
```

Outputs appear under **`bare_slab_benchmark/outputs/`** when you run the programs.

---

## Quick reference: where outputs go

| Program / script | Main output directory |
|------------------|------------------------|
| `main_program` | `outputs/mesh_dumps/`, `outputs/plots/mesh/`, optional `outputs/xs_mapping/` |
| `monte_carlo_solver` (2-group) | `outputs/monte_carlo/` |
| `seven_groups/monte_carlo_solver_7g.exe` | `seven_groups/outputs/monte_carlo/` |
| `finite_difference/diffusion_solver_7g.exe` | `finite_difference/outputs/` |
| `run_energy_collapse_homogenization.py` | `energy_collapse_homogenization/outputs/` |
| `bare_slab_benchmark` | `bare_slab_benchmark/outputs/` |

Rerunning a solver **overwrites** its usual CSV/PNG outputs. For git, either commit snapshots you care about or leave generated folders untracked.

---

## Version control and build artifacts

- **`.gitignore`** ignores Fortran intermediates (`*.o`, `*.mod`, …), **`*.exe`**, and Python **`__pycache__/`** / **`*.pyc`**. It does **not** ignore the **`finite_difference/`** source tree.
- After a fresh clone you **rebuild** the Fortran executables with the commands in this file and in **`finite_difference/README.txt`**.

