"""
Mesh sensitivity (UO2): hold MPWR fixed (from input_file.txt), vary MPFR only.
Each run also sets Generations, Histories, and Skip to fixed values (see module constants).
Writes: outputs/monte_carlo/mesh_sensitivity_study_UO2.png
"""
import re
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 69 positions: alternating water (2) / UO2 (0).
_UO2_MATID = " ".join(str(2 if i % 2 == 0 else 0) for i in range(69))

# Vary only meshes per fuel rod (MPFR). MPWR comes from your input file (unchanged across runs).
_MPFR_VALUES = list(range(3, 16))  # 3 … 15

# Monte Carlo run controls: fixed for every sweep point (matches documented “ideal” MC setup).
_MC_GENERATIONS = 100
_MC_HISTORIES = 1000
_MC_SKIP = 5


def _parse_mpwr(text: str) -> int:
    m = re.search(r"(?m)^\s*MPWR\s*=\s*(\d+)", text)
    if not m:
        raise ValueError("Could not parse MPWR from input_file.txt")
    return int(m.group(1))


def _read_summary(path: Path) -> dict[str, float]:
    vals: dict[str, float] = {}
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        k = k.strip()
        try:
            vals[k] = float(v.strip())
        except ValueError:
            continue
    return vals


def _patch_input(text: str, mpfr: int, mpwr: int) -> str:
    text = re.sub(r"(?m)^\s*MPFR\s*=\s*[^\n#]*", f"MPFR = {mpfr}", text, count=1)
    text = re.sub(r"(?m)^\s*MPWR\s*=\s*[^\n#]*", f"MPWR = {mpwr}", text, count=1)
    text = re.sub(
        r"(?m)^\s*Generations\s*=\s*[^\n#]*",
        f"Generations = {_MC_GENERATIONS}",
        text,
        count=1,
    )
    text = re.sub(
        r"(?m)^\s*Histories\s*=\s*[^\n#]*",
        f"Histories = {_MC_HISTORIES}",
        text,
        count=1,
    )
    text = re.sub(r"(?m)^\s*Skip\s*=\s*[^\n#]*", f"Skip = {_MC_SKIP}", text, count=1)
    text = re.sub(
        r"(?m)^\s*MatID\s*=\s*[^\n]*",
        f"    MatID = {_UO2_MATID}",
        text,
        count=1,
    )
    return text


def main() -> int:
    repo_root = Path(".").resolve()
    input_path = repo_root / "input_file.txt"
    out_dir = repo_root / "outputs" / "monte_carlo"
    out_dir.mkdir(parents=True, exist_ok=True)

    exe = repo_root / "monte_carlo_solver.exe"
    if not exe.is_file():
        exe = repo_root / "monte_carlo_solver"
    if not exe.is_file():
        raise FileNotFoundError("Build monte_carlo_solver first.")

    original = input_path.read_text(encoding="utf-8", errors="replace")
    mpwr_fixed = _parse_mpwr(original)

    n_cells_list: list[float] = []
    k_list: list[float] = []

    try:
        n_total = len(_MPFR_VALUES)
        print(
            f"MPWR fixed at {mpwr_fixed} (from input); sweeping MPFR only; "
            f"MC: Generations={_MC_GENERATIONS}, Histories={_MC_HISTORIES}, Skip={_MC_SKIP}.",
            flush=True,
        )
        for idx, mpfr in enumerate(_MPFR_VALUES, start=1):
            print(f"Mesh sweep {idx}/{n_total}: MPFR={mpfr} MPWR={mpwr_fixed}", flush=True)
            input_path.write_text(_patch_input(original, mpfr, mpwr_fixed), encoding="utf-8")
            subprocess.run([str(exe)], check=True, cwd=str(repo_root))
            summary = _read_summary(out_dir / "summary.txt")
            n_cells_list.append(summary["total_mesh_cells"])
            k_list.append(summary["k_eff_mc_mean"])
    finally:
        input_path.write_text(original, encoding="utf-8")

    mpfr_x = np.array(_MPFR_VALUES, dtype=float)
    n_cells_arr = np.array(n_cells_list, dtype=float)

    # Matplotlib passes scalars or 1D arrays; np.interp matches that and must not be wrapped in float().
    def _mpfr_to_n_cells(mpfr):
        return np.interp(np.asarray(mpfr, dtype=float), mpfr_x, n_cells_arr)

    def _n_cells_to_mpfr(n_cells):
        return np.interp(np.asarray(n_cells, dtype=float), n_cells_arr, mpfr_x)

    mc_legend = (
        f"Generations = {_MC_GENERATIONS}\n"
        f"Histories = {_MC_HISTORIES}\n"
        f"Skip = {_MC_SKIP}"
    )

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.plot(
        mpfr_x,
        k_list,
        marker="o",
        lw=1.6,
        color="#1f55a4",
        markersize=6,
        label=mc_legend,
    )
    ax.legend(loc="best", fontsize=9, framealpha=0.95)
    ax.set_xlabel("Meshes per fuel rod (MPFR)")
    ax.set_ylabel(r"$k_\mathrm{eff}$ (mean, post-skip)")
    ax.set_title(f"Mesh sensitivity (UO2): MPWR = {mpwr_fixed} fixed, MPFR varied")
    ax.grid(alpha=0.25)

    secax = ax.secondary_xaxis("top", functions=(_mpfr_to_n_cells, _n_cells_to_mpfr))
    secax.set_xlabel("Total mesh cells")

    fig.tight_layout()

    out_png = out_dir / "mesh_sensitivity_study_UO2.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
