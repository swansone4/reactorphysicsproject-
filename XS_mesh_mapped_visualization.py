import csv
import os
import re
from itertools import groupby
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


XS_NAMES = ["sigTR", "sigIS", "sigDS", "sigA", "sigF", "nuT", "chiT"]
# Legend swatches (solid, readable).
MATERIAL_COLORS = {0: "#1f77b4", 1: "#d62728", 2: "#008b8b", 3: "#5c5c5c"}
# Panel backgrounds: distinct pastels so g1 (blue) / g2 (orange) lines stay visible.
MATERIAL_BG = {0: "#b8d4f0", 1: "#f4c2c2", 2: "#9ee8e8", 3: "#c5c5ce"}
MATERIAL_NAMES = {0: "UO2", 1: "MOX", 2: "H2O", 3: "CR"}
# Strong enough to read; still below line zorder.
MATERIAL_BG_ALPHA = 0.55


def read_csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        dr = csv.DictReader(f)
        if dr.fieldnames:
            dr.fieldnames = [h.strip() for h in dr.fieldnames]
        rows = []
        for row in dr:
            rows.append(
                {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items() if k}
            )
        return rows


def _input_file_float(path: Path, line_start_pattern: str):
    if not path.is_file():
        return None
    rx = re.compile(line_start_pattern, re.I)
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        head = line.split("#", 1)[0].strip()
        if head and rx.match(head):
            try:
                return float(head.split("=", 1)[1].strip())
            except ValueError:
                return None
    return None


def resolve_mesh_float(rows, col: str, env_var: str, input_pattern: str, dump_path: str):
    if rows and rows[0].get(col) not in (None, ""):
        return float(rows[0][col])
    v = os.environ.get(env_var)
    if v:
        print(f"Note: {dump_path}: no {col} in CSV; using {env_var}={v}")
        return float(v)
    v = _input_file_float(Path("input_file.txt"), input_pattern)
    if v is not None:
        print(f"Note: {dump_path}: no {col} in CSV; using input_file.txt ({input_pattern})={v}")
        return v
    return None


def segment_width_cm(mid: int, pitch: float, diameter: float) -> float:
    """Same as mesh_visualization: water = gap (P−D), pin = D."""
    return max(pitch - diameter, 0.0) if mid == 2 else diameter


def _n_pins(row_dicts) -> int:
    return sum(1 for r in row_dicts if int(r["matid"]) in (0, 1, 3))


def latest_xs_verify_file():
    files = sorted(Path(".").glob("xs_mapping_verify_case_*.csv"), key=lambda p: p.stat().st_mtime)
    return files[-1] if files else None


def load_xs_library(xs_path):
    xs_rows = read_csv(xs_path)
    xs_map = {}
    for row in xs_rows:
        mid = int(row["material_id"])
        g = int(row["energy_group"])
        xs_map[(mid, g)] = {name: float(row[name]) for name in XS_NAMES}
    return xs_map


def flatten_x_cm_and_spans(sorted_rows, pitch: float, diameter: float):
    """
    Same 1D x(cm) convention as mesh_visualization.py: per config row, slabs scaled to L = N×P;
    mesh point k is at x0 + (k+0.5)/npt * w. Multiple rows are concatenated along x.
    Returns (x_cm_list, rod_spans) where rod_spans is (x_lo, x_hi, matid) per rod in order.
    """
    x_mesh = []
    rod_spans = []
    x_global = 0.0

    for _row_id, group in groupby(sorted_rows, key=lambda r: int(r["row"])):
        rd = list(group)
        wgeom = [segment_width_cm(int(r["matid"]), pitch, diameter) for r in rd]
        gtot = sum(wgeom)
        n_pin = _n_pins(rd)
        L = n_pin * pitch if n_pin else (gtot if gtot > 0.0 else pitch)
        s = L / gtot if gtot > 1e-12 else 1.0
        x0 = x_global

        for r, wg in zip(rd, wgeom):
            w = wg * s
            npt = int(r["n_points"])
            mid = int(r["matid"])
            x_lo = x0
            for k in range(npt):
                x_mesh.append(x0 + (k + 0.5) / npt * w)
            x0 += w
            rod_spans.append((x_lo, x0, mid))

        x_global = x0

    return x_mesh, rod_spans


def build_meshwise_records(sorted_rows, xs_map, x_cm_list):
    records = []
    i = 0
    for row in sorted_rows:
        pos = int(row["pos"])
        mid = int(row["matid"])
        n_points = int(row["n_points"])

        for mesh_idx in range(1, n_points + 1):
            rec = {"pos": pos, "mesh": mesh_idx, "matid": mid, "x_cm": x_cm_list[i]}
            i += 1
            for g in (1, 2):
                for xs_name in XS_NAMES:
                    rec[f"{xs_name}_g{g}"] = xs_map[(mid, g)][xs_name]
            records.append(rec)
    return records


def add_xs_panel(ax, records, xs_name):
    x_vals = [r["x_cm"] for r in records]
    g1 = [r[f"{xs_name}_g1"] for r in records]
    g2 = [r[f"{xs_name}_g2"] for r in records]

    ax.plot(x_vals, g1, linewidth=2.0, label="g1", zorder=3, color="#0d47a1")
    ax.plot(x_vals, g2, linewidth=2.0, label="g2", zorder=3, color="#e65100")
    ax.set_ylabel(xs_name)
    ax.grid(alpha=0.25, zorder=1)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92)


def draw_rod_boundaries(axes, rod_spans):
    for x_lo, x_hi, _mid in rod_spans:
        for ax in axes:
            ax.axvline(x_hi, color="k", alpha=0.12, linewidth=0.7)


def add_material_background(ax, rod_spans):
    for x_lo, x_hi, mid in rod_spans:
        ax.axvspan(
            x_lo,
            x_hi,
            facecolor=MATERIAL_BG[mid],
            alpha=MATERIAL_BG_ALPHA,
            zorder=0,
            linewidth=0,
        )


def plot_set(dump_path, xs_map, output_dir):
    dump_rows = read_csv(dump_path)
    if not dump_rows:
        print(f"Skip {dump_path.name}: empty dump.")
        return

    need = {"row", "col", "matid", "n_points", "pos"}
    miss = need - set(dump_rows[0])
    if miss:
        print(f"Skip {dump_path.name}: missing columns {sorted(miss)}.")
        return

    dp = str(dump_path)
    pitch = resolve_mesh_float(dump_rows, "rod_pitch", "MESH_ROD_PITCH", r"RodPitch\s*=", dp)
    diameter = resolve_mesh_float(dump_rows, "rod_diameter", "MESH_ROD_DIAMETER", r"RodDia\s*=", dp)
    if pitch is None or diameter is None:
        print(f"Skip {dump_path.name}: need rod_pitch and rod_diameter in CSV or input_file.txt.")
        return

    sorted_rows = sorted(dump_rows, key=lambda r: (int(r["row"]), int(r["col"])))
    x_cm_list, rod_spans = flatten_x_cm_and_spans(sorted_rows, pitch, diameter)
    records = build_meshwise_records(sorted_rows, xs_map, x_cm_list)

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 1.0, 1.0])

    xs_axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[2, 0]),
        fig.add_subplot(gs[2, 1]),
        fig.add_subplot(gs[3, 0]),
    ]
    empty_ax = fig.add_subplot(gs[3, 1])
    empty_ax.axis("off")

    for ax, xs_name in zip(xs_axes, XS_NAMES):
        add_material_background(ax, rod_spans)
        add_xs_panel(ax, records, xs_name)

    draw_rod_boundaries(xs_axes, rod_spans)

    x_max = rod_spans[-1][1] if rod_spans else 1.0
    for ax in xs_axes:
        ax.set_xlim(0.0, x_max)

    legend_handles = [
        Patch(facecolor=MATERIAL_COLORS[m], edgecolor="0.35", linewidth=0.8, label=f"{m}={MATERIAL_NAMES[m]}")
        for m in [0, 1, 2, 3]
    ]
    xs_axes[-1].set_xlabel("x (cm) — 1D problem direction")

    # Reserve top margin: title high, legend clearly below (figure coordinates).
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.82])
    fig.subplots_adjust(top=0.82)
    fig.suptitle(f"{dump_path.stem}: meshwise XS mapping", y=0.98, fontsize=12)
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.915),
        bbox_transform=fig.transFigure,
        ncol=4,
        fontsize=9,
        frameon=True,
        fancybox=False,
        edgecolor="0.45",
        facecolor="0.98",
        framealpha=0.95,
    )
    out_path = output_dir / f"xs_mesh_mapped_{dump_path.stem}.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    xs_path = latest_xs_verify_file()
    if xs_path is None:
        print("No xs_mapping_verify_case_*.csv found. Run main_program with VerifyXS = 1 first.")
        return 1

    xs_map = load_xs_library(xs_path)
    dumps = sorted(Path(".").glob("mesh_dump_set_*.csv"))
    if not dumps:
        print("No mesh_dump_set_*.csv found. Run main_program first.")
        return 1

    output_dir = Path("XS_mesh_mapped_visualizations")
    output_dir.mkdir(exist_ok=True)

    for dump in dumps:
        plot_set(dump, xs_map, output_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
