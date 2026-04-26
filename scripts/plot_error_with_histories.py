import csv
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _replace_scalar(lines: list[str], key: str, value: int) -> list[str]:
    out = []
    replaced = False
    for ln in lines:
        if ln.lstrip().startswith(f"{key}"):
            out.append(f"{key} = {value}\n")
            replaced = True
        else:
            out.append(ln)
    if not replaced:
        out.append(f"{key} = {value}\n")
    return out


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


def main() -> int:
    repo_root = Path(".").resolve()
    input_path = repo_root / "input_file.txt"
    out_dir = repo_root / "outputs" / "monte_carlo"
    out_dir.mkdir(parents=True, exist_ok=True)

    exe = repo_root / "monte_carlo_solver.exe"
    if not exe.is_file():
        exe = repo_root / "monte_carlo_solver"
    if not exe.exists():
        raise FileNotFoundError("Missing monte_carlo_solver executable. Build it first.")

    histories_list = [200, 400, 800, 1600, 3200, 6400]

    original = input_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    rows_out = []
    try:
        for h in histories_list:
            edited = _replace_scalar(original, "Histories", h)
            input_path.write_text("".join(edited), encoding="utf-8")

            # Run solver (writes outputs/monte_carlo/summary.txt)
            subprocess.run([str(exe)], check=True)

            summary = _read_summary(out_dir / "summary.txt")
            gen = int(summary.get("generations", 0))
            skip = int(summary.get("skip", 0))
            n_active = max(gen - skip, 1)
            k_mean = float(summary.get("k_eff_mc_mean", float("nan")))
            k_std = float(summary.get("k_eff_mc_std", float("nan")))
            k_sem = k_std / (n_active**0.5)

            rows_out.append(
                {
                    "histories": h,
                    "inv_sqrt_histories": 1.0 / (h**0.5),
                    "k_eff_mean": k_mean,
                    "k_eff_std": k_std,
                    "k_eff_sem": k_sem,
                }
            )
    finally:
        input_path.write_text("".join(original), encoding="utf-8")

    csv_path = out_dir / "error_vs_histories.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["histories", "inv_sqrt_histories", "k_eff_mean", "k_eff_std", "k_eff_sem"],
        )
        w.writeheader()
        w.writerows(rows_out)

    x = [r["inv_sqrt_histories"] for r in rows_out]
    err = [r["k_eff_sem"] for r in rows_out]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(x, err, marker="o", lw=1.6)
    ax.set_xlabel(r"$1/\sqrt{\mathrm{Histories}}$")
    ax.set_ylabel(r"$k_\mathrm{eff}$ standard error of mean (post-skip)")
    ax.set_title("Monte Carlo k-effective mean error vs histories")
    ax.grid(alpha=0.25)
    fig.tight_layout()

    out_png = out_dir / "error_vs_histories.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Wrote {csv_path}")
    print(f"Wrote {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

