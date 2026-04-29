from __future__ import annotations

import csv
import os
import re
import subprocess
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _replace_scalar(lines: list[str], key: str, value: int) -> list[str]:
    out: list[str] = []
    replaced = False
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=")
    for line in lines:
        if pattern.match(line):
            out.append(f"{key} = {value}\n")
            replaced = True
        else:
            out.append(line)
    if not replaced:
        out.append(f"{key} = {value}\n")
    return out


def _read_keff_mc_mean(summary_path: Path) -> float:
    txt = summary_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"k_eff_mc_mean\s*=\s*([0-9Ee+\-.]+)", txt)
    if not m:
        raise RuntimeError(f"Could not parse k_eff_mc_mean from {summary_path}")
    return float(m.group(1))


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    input_path = repo_root / "input_file.txt"
    out_dir = repo_root / "seven_groups" / "outputs" / "monte_carlo"
    out_dir.mkdir(parents=True, exist_ok=True)

    exe = repo_root / "seven_groups" / "monte_carlo_solver_7g.exe"
    if not exe.is_file():
        raise FileNotFoundError(
            "Missing seven_groups/monte_carlo_solver_7g.exe. Build it first."
        )

    generations_list = [20, 40, 60, 80, 100, 120]
    histories_list = [500, 1000, 1500, 2000, 3000]
    sweep_points: list[tuple[int, int]] = [
        (g, h) for g in generations_list for h in histories_list
    ]

    mpfr = 7
    mpwr = 2

    original = input_path.read_text(encoding="utf-8", errors="replace").splitlines(
        keepends=True
    )

    rows: list[dict[str, float]] = []
    try:
        for generations, histories in sweep_points:
            edited = original
            edited = _replace_scalar(edited, "Generations", generations)
            edited = _replace_scalar(edited, "Histories", histories)
            edited = _replace_scalar(edited, "MPFR", mpfr)
            edited = _replace_scalar(edited, "MPWR", mpwr)
            input_path.write_text("".join(edited), encoding="utf-8")

            t0 = time.perf_counter()
            env = os.environ.copy()
            env["PATH"] = r"C:\msys64\mingw64\bin" + os.pathsep + env.get("PATH", "")
            subprocess.run([str(exe)], cwd=repo_root, check=True, env=env)
            elapsed_s = time.perf_counter() - t0

            summary = out_dir / "summary.txt"
            keff_mean = _read_keff_mc_mean(summary)
            hg = generations * histories
            rows.append(
                {
                    "generations": float(generations),
                    "histories": float(histories),
                    "histories_times_generations": float(hg),
                    "k_eff_mc_mean": keff_mean,
                    "wall_time_s": elapsed_s,
                }
            )
    finally:
        input_path.write_text("".join(original), encoding="utf-8")

    rows.sort(key=lambda r: r["histories_times_generations"])

    csv_out = out_dir / "criticality_and_time_sweep.csv"
    with csv_out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "generations",
                "histories",
                "histories_times_generations",
                "k_eff_mc_mean",
                "wall_time_s",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    x = [r["histories_times_generations"] for r in rows]
    y_keff = [r["k_eff_mc_mean"] for r in rows]
    y_time = [r["wall_time_s"] for r in rows]

    fig, ax_left = plt.subplots(figsize=(10.5, 5.8))
    ax_right = ax_left.twinx()

    line_keff = ax_left.plot(
        x, y_keff, marker="o", lw=1.8, color="#1f77b4", label=r"$k_{\mathrm{eff}}$"
    )
    line_time = ax_right.plot(
        x, y_time, marker="s", lw=1.8, color="#d62728", label="Wall time (s)"
    )

    ax_left.set_xlabel("Histories x Generations")
    ax_left.set_ylabel(r"Multiplication factor, $k_{\mathrm{eff}}$", color="#1f77b4")
    ax_right.set_ylabel("Wall time (s)", color="#d62728")
    ax_left.tick_params(axis="y", labelcolor="#1f77b4")
    ax_right.tick_params(axis="y", labelcolor="#d62728")

    title = "7g MC: k and wall time vs histories×generations (UO2 checker, MPFR=7 MPWR=2)"
    ax_left.set_title(title)
    ax_left.grid(alpha=0.25)

    lines = line_keff + line_time
    labels = [l.get_label() for l in lines]
    ax_left.legend(lines, labels, loc="best")

    fig.tight_layout()
    png_out = out_dir / "criticality_and_time_vs_histories_times_generations.png"
    fig.savefig(png_out, dpi=220)
    plt.close(fig)

    print(f"Wrote {csv_out}")
    print(f"Wrote {png_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
