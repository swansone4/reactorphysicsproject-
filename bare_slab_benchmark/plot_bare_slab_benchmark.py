import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ROOT = Path(__file__).resolve().parent
_OUT = _ROOT / "outputs"


def read_csv_rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def parse_summary(path: Path):
    vals = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            vals[k.strip()] = float(v.strip())
    return vals


def main():
    _OUT.mkdir(parents=True, exist_ok=True)

    k_rows = read_csv_rows(_OUT / "k_history.csv")
    flux_rows = read_csv_rows(_OUT / "flux_vs_analytic_shape.csv")
    summary = parse_summary(_OUT / "summary.txt")

    generation = [int(r["generation"]) for r in k_rows]
    k_eff = [float(r["k_eff"]) for r in k_rows]

    x = [float(r["x_cm"]) for r in flux_rows]
    flux_mc = [float(r["flux_mc"]) for r in flux_rows]
    flux_shape_analytic = [float(r["flux_shape_analytic"]) for r in flux_rows]

    k_analytic = summary["k_eff_analytic_buckling"]
    k_mc_mean = summary["k_eff_mc_mean"]

    fig1, ax1 = plt.subplots(figsize=(8.5, 4.8))
    ax1.plot(generation, k_eff, lw=1.6, label="MC k_eff per generation")
    ax1.axhline(k_analytic, color="tab:red", linestyle="--", lw=1.4, label=f"Analytic k_eff={k_analytic:.4f}")
    ax1.axhline(k_mc_mean, color="tab:green", linestyle="-.", lw=1.4, label=f"MC mean k_eff={k_mc_mean:.4f}")
    ax1.set_xlabel("Generation")
    ax1.set_ylabel("k_eff")
    ax1.set_title("Bare Slab Benchmark: k_eff Convergence")
    ax1.grid(alpha=0.25)
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig(_OUT / "benchmark_k_eff.png", dpi=180)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(8.5, 4.8))
    ax2.plot(x, flux_mc, lw=1.8, label="MC scalar flux")
    ax2.plot(x, flux_shape_analytic, lw=1.8, linestyle="--", label="Analytic cosine flux")
    ax2.set_xlabel("x (cm)")
    ax2.set_ylabel("Scalar flux, phi (1/cm^2 per source neutron)")
    ax2.set_title("Bare Slab Benchmark: Flux Shape")
    ax2.grid(alpha=0.25)
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig(_OUT / "benchmark_flux_shape.png", dpi=180)
    plt.close(fig2)

    print(f"Wrote {_OUT / 'benchmark_k_eff.png'}")
    print(f"Wrote {_OUT / 'benchmark_flux_shape.png'}")


if __name__ == "__main__":
    main()
