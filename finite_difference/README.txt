Finite-difference 7-group diffusion (Lecture 17–19)
==================================================

Formulation
-----------
- 3-point finite-difference: coupling d_{i,i+1} = 2 d_i d_{i+1}/(d_i + d_{i+1}) with
  d_i = D_i / \Delta x_i, D = 1/(3 \Sigma_tr).
- Removal per group: R = \Sigma_a + \sum_{g' \ne g} \Sigma_{s,g \to g'} (inscatter to g
  and fission stay on the RHS; power iteration for k, Lecture 19).
- Boundaries: reflective / vacuum (albedo) from Lecture 18-3, same as the boxed notes: first and
  last row use 2 d (1 - \beta) on the end faces plus interior coupling d_{i-1,i}, d_{i,i+1};
  d = D / \Delta x, and
  \beta = 1 / ( 1 + (1/(4 d)) (1 - \alpha)/(1 + \alpha) )  with  \alpha = 0  (vacuum)  or  1  (reflect)  from
  BoundL_7G / BoundR_7G.  For \alpha = 1, (1 - \alpha)/(1 + \alpha) = 0, so  \beta = 1
  and 2 d (1 - \beta) = 0 on that face, matching a_{11} / a_{NN}  in the lecture.

Data / mesh
---------
- Same 7g CSV libraries as `seven_groups`: `macroscopic_cross_sections.csv`,
  `scattering_matrix.csv`.
- Same geometry: `input_file.txt` (Config, MPFR, MPWR, ConfigSets, RodPitch, RodDia,
  BoundL_7G, BoundR_7G, ...).

Tridiagonal solve (Lecture 19)
---------------------------
Jacobi / G--S / SOR stop on the tridiagonal max residual ||b - A x||_inf (see fd_solvers_1d);
the direct method is matrix inversion (\phi = A^{-1}Q), as in Lecture 19. Large meshes may need
many SOR sweeps; power iteration k stopping uses a relative 1.0E-3 in diffusion_solver_7g.f90.

Command-line first argument:
  (omit)  or  direct   — direct matrix inversion
  jacobi
  gauss   /  gs
  sor
Optional second: SOR \omega (default 1.5) when using sor.

Build (from repository root)
----------------------------
gfortran -O2 -o finite_difference/diffusion_solver_7g.exe \
  read_input_file.f90 mesh_generation.f90 seven_groups/csv_xs_7g.f90 \
  finite_difference/fd_solvers_1d.f90 finite_difference/diffusion_solver_7g.f90 \
  -llapack -lblas

Run (from repository root; uses input_file.txt in cwd)
------------------------------------------------------
finite_difference\\diffusion_solver_7g.exe
finite_difference\\diffusion_solver_7g.exe gs
finite_difference\\diffusion_solver_7g.exe sor 1.3

The executable writes CSV under `finite_difference/outputs/`, then runs
`python finite_difference/plot_fd_flux_7g.py` if available.

Outputs
-------
- `finite_difference/outputs/flux_fundamental.csv`
- `finite_difference/outputs/flux_cells.csv`
- `finite_difference/outputs/summary.txt`
- `finite_difference/outputs/fd_*.png` (after plot script)

K_eff vs Monte Carlo / 2G analytic
----------------------------------
- k in `summary.txt` as `k_eff_fd`  comes from  **power iteration** on the multigroup
  diffusion  eigenvalue (Lecture~19), using the **seven_group**  CSV macro + scatter data.
- The line  `k_diag_fission_over_removal_Rphi`  is  NOT  the same as k: it is the ratio
  of integrated  \nu\Sigma_f \phi  to integrated  R\phi  with  R  =  removal; it **does
  not**  include  diffusion  leakage  through  edges  in  the  denominator, so do not
  compare  it  to  k  or  to  MC  directly.
- 7G  finite-difference  (P1, 1D  assembly  chain)  is  not  the  same  physics
  as  7G  Monte  Carlo  transport;  k  can  differ  from  k_MC  by  order  \sim0.1--0.2
  in  a  strong  water/fuel  lattice, which  is  not  by  itself  evidence  of  a  coding
  error:  the  D\&H  stencil  and  albedo  above  match  the  lecture.  A  self-contained
  check  of  the  *diffusion*  solver  is  a  simple 1D/2G  case  with  a  **closed-form**
  diffusion  k,  not 7G  transport  k.
- 2-group  analytic  k  (e.g.  with  `Analk`  and  the 2G  `XSData`  in  the  main
  input)  use  a  **different  library**  from  the 7G  `seven_groups`  CSVs.  You  should
  **not**  expect  k_eff_fd  (~0.9--1.2  in  7G  diffusion)  to  line  up  with  2G
  1.26/1.12;  7G  vs 2G  and  diffusion  vs  transport  all  move  k.
-  **ConfigSets**  MatID  in  the  main  `input_file.txt`  (UO2=0,  MOX=1,  H2O=2,  CR=3):
  **H2O/fuel  “UO2”  strip  **: 2,0,2,0,…;  **H2O/fuel  “MOX”  strip  **: 2,1,2,1,…  (same
  `NumRods`  count).  Use  `Config  =  0`  or  `1`  when  two  `Set`  lines  are  present.
