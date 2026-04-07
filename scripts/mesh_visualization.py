#================================================================================
#1D mesh plotter for mesh_dump_set_*.csv.
#
#Reads Fortran dumps, resolves P (pitch) and D (rod diameter) from CSV or fallbacks,
#draws each MatID as a horizontal slab (gap = P−D, pin = D), scales the chain to span
#[0, N×P] with N = count of pins (MatID 0/1/3), writes outputs/plots/mesh/*.png.
#===============================================================================

import csv
import glob
import os
import re
from pathlib import Path


def read_rows(path: str) -> list:
    """CSV → list of dicts; utf-8-sig strips BOM; strip header keys and string values."""
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        dr = csv.DictReader(f)
        if dr.fieldnames:
            dr.fieldnames = [h.strip() for h in dr.fieldnames]
        return [
            {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items() if k}
            for row in dr
        ]


def _input_file_float(path: Path, line_start_pattern: str):
    """First 'Keyword = float' line in input (ignore # comments); None if missing/bad."""
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


def resolve_mesh_float(rows, col: str, env_var: str, input_pattern: str, input_label: str, dump_path: str):
    """Scalar from CSV column, else env, else input_file.txt; print Note when not from CSV."""
    if rows and rows[0].get(col) not in (None, ""):
        return float(rows[0][col])
    v = os.environ.get(env_var)
    if v:
        print(f"Note: {dump_path}: no {col} in CSV; using {env_var}={v}")
        return float(v)
    v = _input_file_float(Path("input_file.txt"), input_pattern)
    if v is not None:
        print(f"Note: {dump_path}: no {col} in CSV; using input_file.txt {input_label}={v}")
        return v
    return None


def segment_width_cm(mid: int, pitch: float, diameter: float) -> float:
    """Slab width along x: water = gap (P−D), solid pin = D."""
    return max(pitch - diameter, 0.0) if mid == 2 else diameter


def _n_pins(row_dicts) -> int:
    return sum(1 for r in row_dicts if int(r["matid"]) in (0, 1, 3))


def _target_length_cm(row_dicts, pitch: float, diameter: float) -> float:
    """L = N×P for N pins; if no pins, sum segment widths (fallback)."""
    n = _n_pins(row_dicts)
    if n > 0:
        return n * pitch
    g = sum(segment_width_cm(int(r["matid"]), pitch, diameter) for r in row_dicts)
    return g if g > 0 else pitch


def _rows_for_config_row(all_rows, irow: int) -> list:
    return sorted((r for r in all_rows if int(r["row"]) == irow), key=lambda r: int(r["col"]))


def main() -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.ticker import FormatStrFormatter, MultipleLocator

    dumps = sorted(glob.glob("outputs/mesh_dumps/mesh_dump_set_*.csv"))
    if not dumps:
        dumps = sorted(glob.glob("mesh_dump_set_*.csv"))
    if not dumps:
        print("No mesh dump files found (expected outputs/mesh_dumps/mesh_dump_set_*.csv).")
        return 2

    out = Path("outputs/plots/mesh")
    out.mkdir(parents=True, exist_ok=True)

    # MatID → face color; legend text is "id = name"
    col = {0: "#1f77b4", 1: "#d62728", 2: "#17becf", 3: "#7f7f7f"}
    name = {0: "UO2", 1: "MOX", 2: "H2O", 3: "CR"}
    legend_handles = [Patch(facecolor=col[m], edgecolor="none", label=f"{m} = {name[m]}") for m in col]

    need = {"row", "col", "matid", "n_points", "ncols", "pos"}

    for path in dumps:
        rows = read_rows(path)
        if not rows:
            print(f"Skip {path}: empty.")
            continue
        miss = need - set(rows[0])
        if miss:
            print(f"Skip {path}: missing {sorted(miss)}.")
            continue

        # Geometry: must match Fortran (pitch/diameter in dump or same fallbacks as before)
        P = resolve_mesh_float(rows, "rod_pitch", "MESH_ROD_PITCH", r"RodPitch\s*=", "RodPitch", path)
        if P is None:
            print(f"Skip {path}: need rod_pitch (CSV / MESH_ROD_PITCH / input RodPitch).")
            continue
        D = resolve_mesh_float(rows, "rod_diameter", "MESH_ROD_DIAMETER", r"RodDia\s*=", "RodDia", path)
        if D is None:
            print(f"Skip {path}: need rod_diameter (CSV / MESH_ROD_DIAMETER / input RodDia).")
            continue

        n_cfg_rows = max(int(r["row"]) for r in rows)
        x_max = max(_target_length_cm(_rows_for_config_row(rows, k), P, D) for k in range(1, n_cfg_rows + 1))

        fig_w = max(14.0, 0.55 * x_max + 2.6)
        fig_h = max(3.4, 2.3 * n_cfg_rows + 1.5)
        fig, axes = plt.subplots(n_cfg_rows, 1, sharex=True, figsize=(fig_w, fig_h))
        axes = [axes] if n_cfg_rows == 1 else list(axes)

        dy = 0.20  # half-thickness of the horizontal “ribbon” (y is only visual)

        for ir in range(1, n_cfg_rows + 1):
            ax = axes[ir - 1]
            rd = _rows_for_config_row(rows, ir)
            # True slab widths in cm (gap vs pin); sum may exceed L if an extra trailing gap exists in MatID list.
            wgeom = [segment_width_cm(int(r["matid"]), P, D) for r in rd]
            gtot = sum(wgeom)
            n_pin = _n_pins(rd)
            # Target span along x = N×P (Fortran print_mesh); scale compresses/stretches slabs to fit [0, L].
            L = n_pin * P if n_pin else (gtot if gtot > 0 else P)
            s = L / gtot if gtot > 1e-12 else 1.0
            x0 = 0.0

            for r, wg in zip(rd, wgeom):
                pos, mid, npt = int(r["pos"]), int(r["matid"]), int(r["n_points"])
                w = wg * s
                x1, xc = x0 + w, x0 + 0.5 * w
                fc = col.get(mid, "#cccccc")
                ax.add_patch(
                    Rectangle((x0, -dy), w, 2 * dy, facecolor=fc, edgecolor="none", linewidth=0, alpha=0.45, zorder=1)
                )
                if npt >= 1:
                    ax.scatter(
                        [x0 + (k + 0.5) / npt * w for k in range(npt)],
                        [0.0] * npt,
                        s=14,
                        c=fc,
                        edgecolors="k",
                        linewidths=0.25,
                        zorder=3,
                    )
                ax.text(xc, 0.29, str(pos), ha="center", va="bottom", fontsize=9, color="0.20", fontweight="bold")
                x0 = x1

            ax.set_ylim(-0.30, 0.40)
            ax.set_yticks([])
            ax.set_ylabel(f"row {ir}", rotation=0, labelpad=30, va="center", fontsize=10)
            ax.axhline(0.0, color="0.75", linewidth=0.6, zorder=0)

        for ax in axes:
            ax.set_xlim(0.0, x_max)
            ax.xaxis.set_major_locator(MultipleLocator(1.0))
            ax.xaxis.set_minor_locator(MultipleLocator(0.25))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            ax.grid(True, axis="x", which="major", linestyle=":", linewidth=0.5, alpha=0.45)
            ax.grid(True, axis="x", which="minor", linestyle=":", linewidth=0.3, alpha=0.22)
            ax.tick_params(axis="x", labelsize=9)

        axes[-1].set_xlabel("x (cm) — 1D problem direction", fontsize=11, labelpad=4)
        fig.subplots_adjust(top=0.64, bottom=0.14, left=0.08, right=0.98, hspace=0.32)
        fig.suptitle(f"{Path(path).stem} — 1D x-direction mesh view", y=0.97, fontsize=13, fontweight="bold")
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.90),
            ncol=4,
            frameon=True,
            fancybox=False,
            edgecolor="0.5",
            facecolor="0.98",
            fontsize=10,
            handlelength=1.8,
            columnspacing=1.8,
            borderpad=0.6,
        )

        png = out / f"{Path(path).stem}.png"
        fig.savefig(png, dpi=190)
        plt.close(fig)
        print(f"Wrote {png}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
