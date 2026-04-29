import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _read_input_mesh_params(repo_root: Path):
    # scrape MPFR, MPWR, MatID from input_file.txt
    text = (repo_root / "input_file.txt").read_text(encoding="utf-8", errors="replace")
    mpfr = mpwr = None
    for line in text.splitlines():
        s = line.split("#", 1)[0].strip()
        if s.upper().startswith("MPFR"):
            mpfr = int(float(s.split("=", 1)[1].strip().split()[0]))
        elif s.upper().startswith("MPWR"):
            mpwr = int(float(s.split("=", 1)[1].strip().split()[0]))
    m = re.search(r"^\s*MatID\s*=\s*(.+)$", text, re.MULTILINE | re.IGNORECASE)
    if not m or mpfr is None or mpwr is None:
        raise ValueError("Could not parse MPFR/MPWR/MatID from input_file.txt")
    mats = [int(x) for x in m.group(1).split()]
    return mpfr, mpwr, mats


def _flux_group_columns(fieldnames):
    # (group number, column name) sorted, from g1 g2 ... in the csv header
    if not fieldnames:
        return []
    out = []
    for name in fieldnames:
        m = re.fullmatch(r"g(\d+)", name.strip(), flags=re.IGNORECASE)
        if m:
            out.append((int(m.group(1)), name.strip()))
    out.sort(key=lambda t: t[0])
    return out


def _cells_per_pos(mid: int, mpfr: int, mpwr: int) -> int:
    return mpwr if mid == 2 else mpfr


def _fuel_segment_averages(rows, mat_ids, mpfr, mpwr, col_a, col_b):
    # one point per fuel rod (mat 0/1): mean flux in that mpfr chunk, x = mean x_cm in chunk
    n = len(rows)
    i0 = 0
    xc, y1, y2 = [], [], []
    for mid in mat_ids:
        nc = _cells_per_pos(mid, mpfr, mpwr)
        if i0 + nc > n:
            raise ValueError(
                f"MatID layout ({len(mat_ids)} positions) does not match flux rows ({n} cells)."
            )
        chunk = rows[i0 : i0 + nc]
        i0 += nc
        if mid not in (0, 1):
            continue
        xcm = [float(r["x_cm"]) for r in chunk]
        ga = [float(r[col_a]) for r in chunk]
        gb = [float(r[col_b]) for r in chunk]
        xc.append(sum(xcm) / len(xcm))
        y1.append(sum(ga) / len(ga))
        y2.append(sum(gb) / len(gb))
    if i0 != n:
        raise ValueError(f"Cell count mismatch: mapped {i0} cells, CSV has {n}.")
    return xc, y1, y2


def _read_flux_fundamental(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} is empty.")

    groups = _flux_group_columns(list(fieldnames) if fieldnames else [])
    if len(groups) < 1:
        raise ValueError(f"{path} has no g1..gN flux columns.")
    g_first = groups[0][1]
    g_last = groups[-1][1]

    x = [float(r["x_cm"]) for r in rows]
    phi_first = [float(r[g_first]) for r in rows]
    phi_last = [float(r[g_last]) for r in rows]
    return rows, x, phi_first, phi_last, len(groups), g_first, g_last


def _fuel_tag(mat_ids: list[int]) -> str:
    fuels = [m for m in mat_ids if m in (0, 1)]
    if fuels and all(m == 0 for m in fuels):
        return "UO2"
    if fuels and all(m == 1 for m in fuels):
        return "MOX"
    return "fuel_mix"


def main() -> int:
    repo_root = Path(".").resolve()
    out_dir = repo_root / "outputs" / "monte_carlo"
    out_dir.mkdir(parents=True, exist_ok=True)

    mc_csv = out_dir / "flux_fundamental.csv"
    if not mc_csv.is_file():
        raise FileNotFoundError(f"Missing {mc_csv}. Run monte_carlo_solver first.")

    mpfr, mpwr, mat_ids = _read_input_mesh_params(repo_root)
    tag = sys.argv[1] if len(sys.argv) > 1 else _fuel_tag(mat_ids)
    rows, x, phi_fast, phi_therm, n_groups, g_first, g_last = _read_flux_fundamental(mc_csv)
    xc, av1, av2 = _fuel_segment_averages(rows, mat_ids, mpfr, mpwr, g_first, g_last)

    fig, ax = plt.subplots(figsize=(9.0, 5.2))

    ax.plot(x, phi_fast, lw=1.9, color="#1f77b4", label=f"Monte Carlo flux ({g_first})")
    ax.plot(x, phi_therm, lw=1.9, color="#ff7f0e", label=f"Monte Carlo flux ({g_last})")
    ax.scatter(
        xc,
        av1,
        s=36,
        color="#1f77b4",
        edgecolors="k",
        linewidths=0.6,
        zorder=5,
        label=f"Mean {g_first} per fuel segment",
    )
    ax.scatter(
        xc,
        av2,
        s=36,
        color="#ff7f0e",
        edgecolors="k",
        linewidths=0.6,
        zorder=5,
        label=f"Mean {g_last} per fuel segment",
    )

    ax.set_xlabel("x (cm)")
    ax.set_ylabel("Scalar flux (Monte Carlo)")
    title = f"Monte Carlo flux ({g_first} & {g_last}) — {tag}"
    if n_groups > 2:
        title = f"Monte Carlo flux ({n_groups} groups: {g_first} & {g_last}) — {tag}"
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()

    out_png = out_dir / f"{tag}-mc_{g_first}_{g_last}.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
