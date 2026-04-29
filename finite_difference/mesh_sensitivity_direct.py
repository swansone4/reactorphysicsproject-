import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


THERMAL_GROUP = "G7"


def parse_scalar(input_text: str, key: str) -> float:
    m = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s*$", input_text)
    if not m:
        raise ValueError(f"Could not parse '{key}' from input_file.txt")
    return float(m.group(1))


def parse_config0_counts(input_text: str) -> tuple[int, int]:
    pat = re.compile(r"(?ms)^\s*Set\s*=\s*0\s*$.*?^\s*MatID\s*=\s*(.*?)\s*$")
    m = pat.search(input_text)
    if not m:
        raise ValueError("Could not find Set=0 MatID line in input_file.txt")
    ids = [int(x) for x in re.findall(r"-?\d+", m.group(1))]
    n_water = sum(1 for v in ids if v == 2)
    return len(ids) - n_water, n_water


def read_thermal_sigtr(xs_csv: Path) -> tuple[float, float]:
    sigtr_uo2 = None
    sigtr_h2o = None
    with xs_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            if row["Type"] == "SigTR" and row["Group"] == THERMAL_GROUP:
                sigtr_uo2 = float(row["UO2"])
                sigtr_h2o = float(row["H2O"])
                break
    if sigtr_uo2 is None or sigtr_h2o is None:
        raise ValueError(f"Missing SigTR,{THERMAL_GROUP} in {xs_csv}")
    return sigtr_uo2, sigtr_h2o


def total_meshes(n_fuel_like: int, n_water: int, mpfr: int, mpwr: int) -> int:
    return n_fuel_like * mpfr + n_water * mpwr


def plot_mesh_over_mfp(rows: list[dict], out_png: Path) -> None:
    labels = [f"{int(r['mpfr'])}/{int(r['mpwr'])}" for r in rows]
    y = [r["mesh_over_mfp_control"] for r in rows]
    x = list(range(len(rows)))
    meshes = [r["cells"] for r in rows]

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.plot(x, y, marker="o", linewidth=1.5)
    ax.axhline(1.0, linestyle="--", linewidth=1.0, color="tab:red", label=r"Target = 1")
    ax.set_xlabel("MPFR / MPWR")
    ax.set_ylabel(r"Mesh size / thermal mean free path")
    ax.set_title("Direct Matrix Inversion Thermal Neutron Mesh Criteria")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=8)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(x)
    ax_top.set_xticklabels([str(m) for m in meshes], fontsize=8)
    ax_top.set_xlabel("Total number of meshes")

    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    input_path = repo / "input_file.txt"
    xs_csv = repo / "seven_groups" / "macroscopic_cross_sections.csv"
    out_csv = repo / "finite_difference" / "outputs" / "mesh_sensitivity_thermal_mfp_results.csv"
    out_png = repo / "finite_difference" / "outputs" / "mesh_sensitivity_thermal_mfp_ratio.png"

    mesh_cases = [
        (5, 2), (6, 2), (7, 2), (8, 2), (9, 2), (10, 2),
        (11, 4), (12, 4), (13, 4), (14, 4), (15, 4), (16, 4),
        (17, 6), (18, 6), (19, 6), (20, 6), (21, 8), (22, 8), (23, 8), (24, 8),
    ]

    input_text = input_path.read_text(encoding="utf-8")
    rod_diameter = parse_scalar(input_text, "RodDia")
    n_fuel_like, n_water = parse_config0_counts(input_text)
    sigtr_uo2_g7, sigtr_h2o_g7 = read_thermal_sigtr(xs_csv)

    lambda_th_f = 1.0 / sigtr_uo2_g7
    lambda_th_w = 1.0 / sigtr_h2o_g7

    rows = []
    for mpfr, mpwr in mesh_cases:
        h_f = rod_diameter / mpfr
        h_w = rod_diameter / mpwr
        ratio_f = h_f / lambda_th_f
        ratio_w = h_w / lambda_th_w
        control = max(ratio_f, ratio_w)
        rows.append(
            {
                "mpfr": mpfr,
                "mpwr": mpwr,
                "mesh_fraction": mpfr / mpwr,
                "cells": total_meshes(n_fuel_like, n_water, mpfr, mpwr),
                "h_f_cm": h_f,
                "h_w_cm": h_w,
                "thermal_group": THERMAL_GROUP,
                "sigtr_thermal_fuel_cm_inv": sigtr_uo2_g7,
                "sigtr_thermal_water_cm_inv": sigtr_h2o_g7,
                "lambda_thermal_fuel_cm": lambda_th_f,
                "lambda_thermal_water_cm": lambda_th_w,
                "mesh_over_mfp_fuel": ratio_f,
                "mesh_over_mfp_water": ratio_w,
                "mesh_over_mfp_control": control,
                "abs_to_one": abs(control - 1.0),
            }
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "mpfr", "mpwr", "mesh_fraction", "cells", "h_f_cm", "h_w_cm", "thermal_group",
                "sigtr_thermal_fuel_cm_inv", "sigtr_thermal_water_cm_inv",
                "lambda_thermal_fuel_cm", "lambda_thermal_water_cm",
                "mesh_over_mfp_fuel", "mesh_over_mfp_water", "mesh_over_mfp_control", "abs_to_one",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    plot_mesh_over_mfp(rows, out_png)
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_png}")
    best = min(rows, key=lambda r: r["abs_to_one"])
    print(
        "Best-by-target (closest to 1): "
        f"MPFR/MPWR={best['mpfr']}/{best['mpwr']}, "
        f"control_ratio={best['mesh_over_mfp_control']:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
