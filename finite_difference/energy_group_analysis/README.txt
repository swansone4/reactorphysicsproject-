Energy group analysis (direct matrix inversion)
================================================

This folder holds labeled seven-group flux plots for three MatID layouts
(same rod count and mesh as `input_file.txt`: MPFR / MPWR unchanged).

Plots
-----
1. `01_mox_alternating_2-1-2-1.png` — H2O–MOX alternating (MatID 2–1–2–1).
2. `02_uo2_alternating_2-0-2-0.png` — H2O–UO2 alternating (MatID 2–0–2–0), uranium-oxide fuel.
3. `03_half_uo2_half_mox.png` — First half of lattice 2–0 pattern, second half 2–1 pattern;
   total MatID length unchanged.

Regenerate
----------
From the repository root (with `finite_difference/diffusion_solver_7g_lapack.exe` built):

  python finite_difference/energy_group_analysis/run_energy_group_analysis.py

Solver writes CSV under `finite_difference/flux_magnitude_analysis_power/`; the script
copies the seven-group profile into PNGs here with explicit titles.
