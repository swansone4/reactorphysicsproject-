from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


GROUP_FAST = [1, 2, 3, 4]
GROUP_THERMAL = [5, 6, 7]
MAT_NAME = {0: "UO2", 1: "MOX", 2: "H2O", 3: "CR"}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_input(input_text: str) -> tuple[int, int, list[int]]:
    mpfr = int(re.search(r"(?m)^\s*MPFR\s*=\s*(\d+)\s*$", input_text).group(1))
    mpwr = int(re.search(r"(?m)^\s*MPWR\s*=\s*(\d+)\s*$", input_text).group(1))
    m = re.search(r"(?m)^\s*MatID\s*=\s*(.+)$", input_text)
    matid = [int(x) for x in m.group(1).split()]
    return mpfr, mpwr, matid


def replace_matid(input_text: str, new_matid: list[int]) -> str:
    return re.sub(
        r"(?m)^(\s*MatID\s*=\s*).*$",
        lambda m: m.group(1) + " ".join(str(x) for x in new_matid),
        input_text,
    )


def matid_case_mox_2121(n: int) -> list[int]:
    return [2 if i % 2 == 0 else 1 for i in range(n)]


def matid_case_uo2_2020(n: int) -> list[int]:
    return [2 if i % 2 == 0 else 0 for i in range(n)]


def matid_case_half_mox_half_uo2(n: int) -> list[int]:
    n1 = n // 2
    n2 = n - n1
    left = [2 if i % 2 == 0 else 1 for i in range(n1)]   # MOX-left
    right = [2 if i % 2 == 0 else 0 for i in range(n2)]  # UO2-right
    return left + right


def choose_solver_exe(root: Path) -> Path:
    for name in (
        "diffusion_solver_7g_lapack.exe",
        "diffusion_solver_7g_lapack_outputs.exe",
        "diffusion_solver_7g.exe",
    ):
        p = root / "finite_difference" / name
        if p.exists():
            return p
    raise FileNotFoundError("No finite_difference diffusion solver executable found.")


def build_2g_fortran_solver(root: Path, env: dict) -> Path:
    exe = root / "energy_collapse_homogenization" / "solve_keff_2g_direct.exe"
    src = root / "energy_collapse_homogenization" / "solve_keff_2g_direct.f90"
    solver_mod = root / "finite_difference" / "fd_solvers_1d.f90"
    if exe.exists() and exe.stat().st_mtime >= src.stat().st_mtime:
        return exe
    cmd = [
        "gfortran",
        "-O2",
        "-o",
        str(exe),
        str(solver_mod),
        str(src),
        "-LC:\\msys64\\mingw64\\lib",
        "-llapack",
        "-lblas",
    ]
    subprocess.run(cmd, cwd=str(root), env=env, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return exe


def read_macro_xs(path: Path) -> dict[str, np.ndarray]:
    mats = ["UO2", "MOX", "H2O", "CR"]
    idx = {m: i for i, m in enumerate(mats)}
    sigtr = np.zeros((7, 4))
    siga = np.zeros((7, 4))
    nusigf = np.zeros((7, 4))
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            g = int(row["Group"].replace("G", "")) - 1
            rtype = row["Type"].strip().upper()
            vals = np.array([float(row[m]) for m in mats], dtype=float)
            if rtype == "SIGTR":
                sigtr[g, :] = vals
            elif rtype == "SIGA":
                siga[g, :] = vals
            elif rtype == "NUSIGF":
                nusigf[g, :] = vals
    return {"sigtr": sigtr, "siga": siga, "nusigf": nusigf}


def read_scatter_xs(path: Path) -> np.ndarray:
    mats = ["UO2", "MOX", "H2O", "CR"]
    idx = {m: i for i, m in enumerate(mats)}
    s = np.zeros((4, 7, 7), dtype=float)
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            m = idx[row["Material"].strip()]
            g_from = int(row["G_from"]) - 1
            g_to = int(row["G_to"]) - 1
            s[m, g_from, g_to] = float(row["Value"])
    return s


def build_cell_materials(matid: list[int], mpfr: int, mpwr: int) -> list[int]:
    mats: list[int] = []
    for m in matid:
        n = mpwr if m == 2 else mpfr
        mats.extend([m] * n)
    return mats


def load_flux_cells(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    x = np.array([0.5 * (float(r["x_left_cm"]) + float(r["x_right_cm"])) for r in rows], dtype=float)
    dx = np.array([float(r["dx_cm"]) for r in rows], dtype=float)
    phi = np.array([[float(r[f"flux_g{i}"]) for i in range(1, 8)] for r in rows], dtype=float)
    return x, dx, phi


def parse_summary(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def collapse_and_homogenize(
    phi: np.ndarray,
    dx: np.ndarray,
    cell_mats: list[int],
    macro: dict[str, np.ndarray],
    scatter: np.ndarray,
) -> dict[str, float]:
    cell_mats_arr = np.array(cell_mats, dtype=int)
    sigtr = np.zeros_like(phi)
    siga = np.zeros_like(phi)
    nusigf = np.zeros_like(phi)
    for i, m in enumerate(cell_mats_arr):
        sigtr[i, :] = macro["sigtr"][:, m]
        siga[i, :] = macro["siga"][:, m]
        nusigf[i, :] = macro["nusigf"][:, m]

    def homogenized(group: list[int], xs: np.ndarray) -> float:
        gi = np.array([g - 1 for g in group], dtype=int)
        num = np.sum(dx[:, None] * xs[:, gi] * phi[:, gi])
        den = np.sum(dx[:, None] * phi[:, gi])
        return float(num / max(den, 1.0e-30))

    def homog_scatter(g_from: list[int], g_to: list[int]) -> float:
        gf = [g - 1 for g in g_from]
        gt = [g - 1 for g in g_to]
        num = 0.0
        den = 0.0
        for i, m in enumerate(cell_mats_arr):
            for g in gf:
                w = dx[i] * phi[i, g]
                den += w
                num += w * float(np.sum(scatter[m, g, gt]))
        return num / max(den, 1.0e-30)

    phi_fast_cell = np.sum(phi[:, [g - 1 for g in GROUP_FAST]], axis=1)
    phi_therm_cell = np.sum(phi[:, [g - 1 for g in GROUP_THERMAL]], axis=1)

    return {
        "Sigma_tr_G1": homogenized(GROUP_FAST, sigtr),
        "Sigma_tr_G2": homogenized(GROUP_THERMAL, sigtr),
        "Sigma_a_G1": homogenized(GROUP_FAST, siga),
        "Sigma_a_G2": homogenized(GROUP_THERMAL, siga),
        "nuSigma_f_G1": homogenized(GROUP_FAST, nusigf),
        "nuSigma_f_G2": homogenized(GROUP_THERMAL, nusigf),
        "Sigma_s_G1_to_G1": homog_scatter(GROUP_FAST, GROUP_FAST),
        "Sigma_s_G1_to_G2": homog_scatter(GROUP_FAST, GROUP_THERMAL),
        "Sigma_s_G2_to_G1": homog_scatter(GROUP_THERMAL, GROUP_FAST),
        "Sigma_s_G2_to_G2": homog_scatter(GROUP_THERMAL, GROUP_THERMAL),
        "phi_fast_avg": float(np.sum(dx * phi_fast_cell) / max(np.sum(dx), 1.0e-30)),
        "phi_thermal_avg": float(np.sum(dx * phi_therm_cell) / max(np.sum(dx), 1.0e-30)),
    }


def load_2g_keff(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            out[row["case_id"]] = float(row["k_eff_2g_direct"])
    return out


def plot_two_group_flux(x: np.ndarray, phi: np.ndarray, out_png: Path, title: str) -> None:
    fast = np.sum(phi[:, [g - 1 for g in GROUP_FAST]], axis=1)
    therm = np.sum(phi[:, [g - 1 for g in GROUP_THERMAL]], axis=1)
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.plot(x, fast, lw=1.8, label="Collapsed Group 1 (G1-G4)")
    ax.plot(x, therm, lw=1.8, label="Collapsed Group 2 (G5-G7)")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel(r"Collapsed scalar flux, $\phi_G$ (n cm$^{-2}$ s$^{-1}$)")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_overlay(cases: list[dict], out_dir: Path) -> None:
    fig1, ax1 = plt.subplots(figsize=(9.2, 5.0))
    fig2, ax2 = plt.subplots(figsize=(9.2, 5.0))
    for c in cases:
        x = c["x"]
        phi = c["phi"]
        fast = np.sum(phi[:, [g - 1 for g in GROUP_FAST]], axis=1)
        therm = np.sum(phi[:, [g - 1 for g in GROUP_THERMAL]], axis=1)
        ax1.plot(x, fast, lw=1.6, label=c["case_name"])
        ax2.plot(x, therm, lw=1.6, label=c["case_name"])
    ax1.set_xlabel("x (cm)")
    ax1.set_ylabel(r"$\phi_{G1}$ (n cm$^{-2}$ s$^{-1}$)")
    ax1.set_title("Collapsed Group 1 (G1-G4) flux comparison")
    ax1.grid(alpha=0.25)
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig(out_dir / "overlay_group1_fast.png", dpi=180)
    plt.close(fig1)

    ax2.set_xlabel("x (cm)")
    ax2.set_ylabel(r"$\phi_{G2}$ (n cm$^{-2}$ s$^{-1}$)")
    ax2.set_title("Collapsed Group 2 (G5-G7) flux comparison")
    ax2.grid(alpha=0.25)
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig(out_dir / "overlay_group2_thermal.png", dpi=180)
    plt.close(fig2)


def main() -> int:
    root = repo_root()
    work_dir = root / "energy_collapse_homogenization"
    out_dir = work_dir / "outputs"
    plots_dir = out_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    input_path = root / "input_file.txt"
    input_original = input_path.read_text(encoding="utf-8")
    mpfr, mpwr, matid_orig = parse_input(input_original)
    n_ids = len(matid_orig)

    macro = read_macro_xs(root / "seven_groups" / "macroscopic_cross_sections.csv")
    scatter = read_scatter_xs(root / "seven_groups" / "scattering_matrix.csv")
    solver_exe = choose_solver_exe(root)

    env = dict(**subprocess.os.environ)
    env["PATH"] = r"C:\msys64\mingw64\bin;" + env.get("PATH", "")

    cases = [
        ("01_mox", "MOX 2-1-2-1", matid_case_mox_2121(n_ids)),
        ("02_uo2", "UO2 2-0-2-0", matid_case_uo2_2020(n_ids)),
        ("03_half", "Half MOX left / Half UO2 right", matid_case_half_mox_half_uo2(n_ids)),
    ]

    case_summary_rows: list[dict] = []
    xs_rows: list[dict] = []
    overlay_data: list[dict] = []

    try:
        for case_id, case_name, matid in cases:
            input_path.write_text(replace_matid(input_original, matid), encoding="utf-8")
            subprocess.run(
                [str(solver_exe)],
                cwd=str(root),
                env=env,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )

            src_dir = root / "finite_difference" / "outputs"
            case_dir = out_dir / case_id
            case_dir.mkdir(parents=True, exist_ok=True)
            for fn in ("flux_fundamental.csv", "flux_cells.csv", "summary.txt"):
                shutil.copy2(src_dir / fn, case_dir / fn)

            x, dx, phi = load_flux_cells(case_dir / "flux_cells.csv")
            summary = parse_summary(case_dir / "summary.txt")
            k_eff = float(summary["k_eff_fd"])
            total_cells = int(summary["total_mesh_cells"])
            model_length_cm = float(np.sum(dx))

            homog = collapse_and_homogenize(phi, dx, build_cell_materials(matid, mpfr, mpwr), macro, scatter)
            homog_row = {"case_id": case_id, "case_name": case_name, "k_eff_7g": k_eff, "k_eff_2g_direct": 0.0}
            homog_row.update(homog)
            xs_rows.append(homog_row)

            case_summary_rows.append(
                {
                    "case_id": case_id,
                    "case_name": case_name,
                    "mpfr": mpfr,
                    "mpwr": mpwr,
                    "matid_length": len(matid),
                    "total_cells": total_cells,
                    "model_length_cm": model_length_cm,
                    "k_eff_7g": k_eff,
                    "k_eff_2g_direct": 0.0,
                    "matid_string": " ".join(str(v) for v in matid),
                }
            )

            plot_two_group_flux(
                x,
                phi,
                case_dir / "collapsed_two_group_flux.png",
                f"{case_name} — collapsed 2-group flux (G1-4 / G5-7)",
            )
            overlay_data.append({"case_name": case_name, "x": x, "phi": phi})

            with (case_dir / "homogenized_2g_xs.json").open("w", encoding="utf-8") as f:
                json.dump(homog_row, f, indent=2)

        plot_overlay(overlay_data, plots_dir)
    finally:
        input_path.write_text(input_original, encoding="utf-8")

    with (out_dir / "case_summaries_7g.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "case_id",
                "case_name",
                "mpfr",
                "mpwr",
                "matid_length",
                "total_cells",
                "model_length_cm",
                "k_eff_7g",
                "k_eff_2g_direct",
                "matid_string",
            ],
        )
        w.writeheader()
        w.writerows(case_summary_rows)

    with (out_dir / "two_group_homogenized_xs.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "case_id",
                "case_name",
                "k_eff_7g",
                "k_eff_2g_direct",
                "Sigma_tr_G1",
                "Sigma_tr_G2",
                "Sigma_a_G1",
                "Sigma_a_G2",
                "nuSigma_f_G1",
                "nuSigma_f_G2",
                "Sigma_s_G1_to_G1",
                "Sigma_s_G1_to_G2",
                "Sigma_s_G2_to_G1",
                "Sigma_s_G2_to_G2",
                "phi_fast_avg",
                "phi_thermal_avg",
            ],
        )
        w.writeheader()
        w.writerows(xs_rows)

    solver_2g_exe = build_2g_fortran_solver(root, env)
    subprocess.run([str(solver_2g_exe)], cwd=str(root), env=env, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    k2g_map = load_2g_keff(out_dir / "two_group_keff_direct.csv")

    for row in case_summary_rows:
        row["k_eff_2g_direct"] = k2g_map.get(row["case_id"], float("nan"))
    for row in xs_rows:
        row["k_eff_2g_direct"] = k2g_map.get(row["case_id"], float("nan"))
    for case_id, _, _ in cases:
        case_json = out_dir / case_id / "homogenized_2g_xs.json"
        data = json.loads(case_json.read_text(encoding="utf-8"))
        data["k_eff_2g_direct"] = k2g_map.get(case_id, float("nan"))
        case_json.write_text(json.dumps(data, indent=2), encoding="utf-8")

    with (out_dir / "case_summaries_7g.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "case_id",
                "case_name",
                "mpfr",
                "mpwr",
                "matid_length",
                "total_cells",
                "model_length_cm",
                "k_eff_7g",
                "k_eff_2g_direct",
                "matid_string",
            ],
        )
        w.writeheader()
        w.writerows(case_summary_rows)

    with (out_dir / "two_group_homogenized_xs.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "case_id",
                "case_name",
                "k_eff_7g",
                "k_eff_2g_direct",
                "Sigma_tr_G1",
                "Sigma_tr_G2",
                "Sigma_a_G1",
                "Sigma_a_G2",
                "nuSigma_f_G1",
                "nuSigma_f_G2",
                "Sigma_s_G1_to_G1",
                "Sigma_s_G1_to_G2",
                "Sigma_s_G2_to_G1",
                "Sigma_s_G2_to_G2",
                "phi_fast_avg",
                "phi_thermal_avg",
            ],
        )
        w.writeheader()
        w.writerows(xs_rows)

    print(f"Wrote {out_dir / 'case_summaries_7g.csv'}")
    print(f"Wrote {out_dir / 'two_group_homogenized_xs.csv'}")
    print(f"Wrote {out_dir / 'two_group_keff_direct.csv'}")
    print(f"Wrote plots under {plots_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

