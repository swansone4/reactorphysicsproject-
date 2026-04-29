from __future__ import annotations

# Runs FD three ways (MOX strip / UO2 strip / half-half), flux PNGs land in energy_group_analysis/.
# Needs input_file.txt at repo root, seven_groups CSVs, and a built diffusion_solver_7g*.exe.

import csv
import re
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def parse_mpfr_mpwr(text: str) -> tuple[int, int]:
    mpfr = int(re.search(r"(?m)^\s*MPFR\s*=\s*(\d+)\s*$", text).group(1))
    mpwr = int(re.search(r"(?m)^\s*MPWR\s*=\s*(\d+)\s*$", text).group(1))
    return mpfr, mpwr


def replace_matid(text: str, matid_tokens: list[str]) -> str:
    line = " ".join(matid_tokens)
    return re.sub(r"(?m)^(\s*MatID\s*=\s*).*$", rf"\g<1>{line}", text)


def matid_2_1(n: int) -> list[str]:
    return ["2" if i % 2 == 0 else "1" for i in range(n)]


def matid_2_0(n: int) -> list[str]:
    return ["2" if i % 2 == 0 else "0" for i in range(n)]


def matid_half_uo2_half_mox(n: int) -> list[str]:
    n1 = n // 2
    n2 = n - n1
    first = ["2" if i % 2 == 0 else "1" for i in range(n1)]
    second = ["2" if i % 2 == 0 else "0" for i in range(n2)]
    return first + second


def find_solver_exe(root: Path) -> Path:
    for name in (
        "diffusion_solver_7g_lapack.exe",
        "diffusion_solver_7g_lapack_outputs.exe",
        "diffusion_solver_7g.exe",
    ):
        p = root / "finite_difference" / name
        if p.is_file():
            return p
    raise FileNotFoundError(
        "No diffusion_solver executable found under finite_difference/. "
        "Build with gfortran + LAPACK (see finite_difference/README.txt)."
    )


def latest_flux_csv(root: Path) -> Path:
    candidates = [
        root / "finite_difference" / "outputs" / "flux_fundamental.csv",
        root / "finite_difference" / "flux_magnitude_analysis_power" / "flux_fundamental.csv",
    ]
    existing = [p for p in candidates if p.is_file()]
    if not existing:
        raise FileNotFoundError("No flux_fundamental.csv found in expected solver output folders.")
    return max(existing, key=lambda p: p.stat().st_mtime)


def plot_flux_fundamental(
    csv_path: Path,
    out_png: Path,
    title: str,
    subtitle: str,
) -> None:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"Empty {csv_path}")
    x = np.array([float(r["x_cm"]) for r in rows], dtype=float)
    group_keys = [f"g{i}" for i in range(1, 8)]
    phi = np.array([[float(r[k]) for k in group_keys] for r in rows], dtype=float)
    colors = plt.cm.viridis([i / 6.0 for i in range(7)])

    fig, ax = plt.subplots(figsize=(10.0, 5.8))
    for i in range(7):
        ax.plot(x, phi[:, i], lw=1.7, color=colors[i], label=f"Group {i + 1}")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel(r"Neutron scalar flux, $\phi_g$ (n cm$^{-2}$ s$^{-1}$)")
    ax.set_title(title + "\n" + subtitle)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper center", ncol=4, fontsize=9, frameon=True)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> int:
    root = repo_root()
    input_path = root / "input_file.txt"
    analysis_dir = root / "finite_difference" / "energy_group_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    original = input_path.read_text(encoding="utf-8")
    mpfr, mpwr = parse_mpfr_mpwr(original)
    m = re.search(r"(?m)^\s*MatID\s*=\s*(.+)$", original)
    if not m:
        print("ERROR: MatID line not found", file=sys.stderr)
        return 1
    n_ids = len(m.group(1).split())
    if n_ids < 1:
        print("ERROR: empty MatID", file=sys.stderr)
        return 1

    exe = find_solver_exe(root)
    env = dict(**subprocess.os.environ)
    env["PATH"] = r"C:\msys64\mingw64\bin;" + env.get("PATH", "")

    cases: list[tuple[str, list[str], str, str]] = [
        (
            "01_mox_2-1-2-1.png",
            matid_2_1(n_ids),
            "Direct matrix inversion — seven-group flux",
            f"MOX assembly (MatID 2–1–2–1), MPFR={mpfr}, MPWR={mpwr}",
        ),
        (
            "02_uo2_2-0-2-0.png",
            matid_2_0(n_ids),
            "Direct matrix inversion — seven-group flux",
            f"UO₂ assembly (MatID 2–0–2–0), MPFR={mpfr}, MPWR={mpwr}",
        ),
        (
            "03_half_uo2_half_mox.png",
            matid_half_uo2_half_mox(n_ids),
            "Direct matrix inversion — seven-group flux",
            f"Half MOX / half UO₂ along lattice (first {n_ids // 2} positions 2–1, remainder 2–0), MPFR={mpfr}, MPWR={mpwr}",
        ),
    ]

    try:
        for fname, tokens, t1, t2 in cases:
            new_text = replace_matid(original, tokens)
            input_path.write_text(new_text, encoding="utf-8")
            subprocess.run(
                [str(exe)],
                cwd=str(root),
                env=env,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            csv_flux = latest_flux_csv(root)
            plot_flux_fundamental(
                csv_flux,
                analysis_dir / fname,
                t1,
                t2,
            )
            print(f"Wrote {analysis_dir / fname}")
    finally:
        input_path.write_text(original, encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
