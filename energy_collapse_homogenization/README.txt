Cross-Section Homogenization and 2-Group Energy Collapse
========================================================

This folder contains a reproducible workflow for collapsing 7-group direct
finite-difference results into 2 groups and homogenizing cross sections for
three configurations:

1) MOX 2-1-2-1
2) UO2 2-0-2-0
3) Half MOX (left) / Half UO2 (right)

Collapse definition
-------------------
- Group 1 (fast): G1-G4
- Group 2 (thermal): G5-G7

Homogenization scope
--------------------
Whole 1D configuration for each case (volume-and-flux weighted).

Main script
-----------
Run from repo root:

  python energy_collapse_homogenization/run_energy_collapse_homogenization.py

Outputs
-------
Under energy_collapse_homogenization/outputs:
- case_summaries_7g.csv
- two_group_homogenized_xs.csv
- two_group_keff_direct.csv
- 01_mox/, 02_uo2/, 03_half/ (snapshots + per-case 2-group plot + JSON)
- plots/overlay_group1_fast.png
- plots/overlay_group2_thermal.png

Notes
-----
- The script temporarily edits only the MatID line in input_file.txt and restores it.
- It uses the available finite_difference solver executable under finite_difference/.
- The 2-group eigenvalue solve is done in Fortran by:
    energy_collapse_homogenization/solve_keff_2g_direct.f90
  Python is used for orchestration (case setup/IO) and plotting only.

