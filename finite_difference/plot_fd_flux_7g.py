import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    out_dir = repo_root / "finite_difference" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_flux = out_dir / "flux_fundamental.csv"
    csv_cells = out_dir / "flux_cells.csv"

    if not csv_flux.is_file():
        raise FileNotFoundError(f"Missing {csv_flux}. Run finite_difference/diffusion_solver_7g first.")
    if not csv_cells.is_file():
        raise FileNotFoundError(f"Missing {csv_cells}. Run finite_difference/diffusion_solver_7g first.")

    with csv_flux.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"{csv_flux} is empty.")
    with csv_cells.open("r", encoding="utf-8-sig", newline="") as f:
        rows_cells = list(csv.DictReader(f))
    if not rows_cells:
        raise ValueError(f"{csv_cells} is empty.")

    x_center = np.array([float(r["x_cm"]) for r in rows], dtype=float)
    group_keys = [f"g{i}" for i in range(1, 8)]
    for key in group_keys:
        if key not in rows[0]:
            raise ValueError(f"{csv_flux} is missing expected column '{key}'.")

    phi = np.array(
        [[float(r[key]) for key in group_keys] for r in rows],
        dtype=float,
    )
    colors = plt.cm.viridis([i / 6.0 for i in range(7)])

    fig, ax = plt.subplots(figsize=(10.0, 5.8))
    for i in range(7):
        ax.plot(x_center, phi[:, i], lw=1.7, color=colors[i], label=f"Group {i + 1}")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel(r"Neutron scalar flux, $\phi_g$ (n cm$^{-2}$ s$^{-1}$)")
    ax.set_title("Seven-group finite-difference diffusion flux (Groups 1-7)")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper center", ncol=4, fontsize=9, frameon=True)
    fig.tight_layout()
    out_png = out_dir / "fd_flux_g1_to_g7.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    fig2, ax2 = plt.subplots(figsize=(10.0, 5.8))
    markevery = max(1, len(x_center) // 40)
    for i in range(7):
        ax2.plot(x_center, phi[:, i], lw=1.0, color=colors[i], alpha=0.55)
        ax2.plot(
            x_center,
            phi[:, i],
            linestyle="none",
            marker="o",
            markersize=2.8,
            markevery=markevery,
            color=colors[i],
            label=f"Group {i + 1}",
        )
    ax2.set_xlabel("Cell center x (cm)")
    ax2.set_ylabel(r"Cell-center neutron scalar flux, $\phi_g$ (n cm$^{-2}$ s$^{-1}$)")
    ax2.set_title("Seven-group cell-center flux (FD)")
    ax2.grid(alpha=0.25)
    ax2.legend(loc="upper center", ncol=4, fontsize=9, frameon=True)
    fig2.tight_layout()
    out_centers = out_dir / "fd_cell_center_flux_g1_to_g7.png"
    fig2.savefig(out_centers, dpi=200)
    plt.close(fig2)

    sigTR = np.array(
        [[float(r[f"sigTR_g{g}"]) for g in range(1, 8)] for r in rows_cells],
        dtype=float,
    )
    sigTR = np.maximum(sigTR, 1.0e-15)
    dx = np.array([float(r["dx_cm"]) for r in rows_cells], dtype=float)
    dx = np.maximum(dx, 1.0e-15)
    D_cell = 1.0 / (3.0 * sigTR)
    d_cell = D_cell / dx[:, None]
    coupling = 2.0 * d_cell[:-1, :] * d_cell[1:, :] / np.maximum(
        d_cell[:-1, :] + d_cell[1:, :], 1.0e-15
    )
    J_edge = -coupling * (phi[1:, :] - phi[:-1, :])
    x_inner_edges = 0.5 * (x_center[:-1] + x_center[1:])

    fig3, ax3 = plt.subplots(figsize=(10.0, 5.8))
    ax3.axhline(0.0, color="k", linestyle="--", lw=0.8, alpha=0.5)
    for i in range(7):
        ax3.plot(x_inner_edges, J_edge[:, i], lw=1.3, color=colors[i], label=f"Group {i + 1}")
    ax3.set_xlabel("Cell edge x (cm)")
    ax3.set_ylabel(r"Cell-edge neutron current, $J(x)$ (n cm$^{-2}$ s$^{-1}$)")
    ax3.set_title(
        r"FD cell-edge current: "
        r"$J_{i+1/2,g}=-\frac{2d_{i,g}d_{i+1,g}}{d_{i,g}+d_{i+1,g}}\left(\phi_{i+1,g}-\phi_{i,g}\right)$"
    )
    ax3.grid(alpha=0.25)
    ax3.legend(loc="upper center", ncol=4, fontsize=9, frameon=True)
    fig3.tight_layout()
    out_edges = out_dir / "fd_cell_edge_current_g1_to_g7.png"
    fig3.savefig(out_edges, dpi=200)
    plt.close(fig3)

    print(f"Wrote {out_png}")
    print(f"Wrote {out_centers}")
    print(f"Wrote {out_edges}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
