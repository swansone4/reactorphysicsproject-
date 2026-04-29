from __future__ import annotations

# Quick FD runs: vac right BC, two different half-core MatID strings, then put input_file.txt back how it was.

import os
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input_file.txt"
EXE = ROOT / "finite_difference" / "diffusion_solver_7g_lapack_outputs.exe"
OUT = ROOT / "finite_difference" / "outputs"


def matid_uo2_then_mox() -> str:
    uo2 = []
    for _ in range(17):
        uo2 += [2, 0]
    mox = []
    for _ in range(17):
        mox += [2, 1]
    return " ".join(str(x) for x in uo2 + mox)


def matid_mox_then_uo2() -> str:
    mox = []
    for _ in range(17):
        mox += [2, 1]
    uo2 = []
    for _ in range(17):
        uo2 += [2, 0]
    return " ".join(str(x) for x in mox + uo2)


def replace_line(text: str, pattern: str, new_line: str) -> str:
    if not re.search(pattern, text, flags=re.MULTILINE):
        raise RuntimeError(f"Pattern not found: {pattern!r}")
    return re.sub(pattern, new_line, text, count=1, flags=re.MULTILINE)


def read_keff(summary: Path) -> str:
    for line in summary.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.strip().startswith("k_eff_fd"):
            parts = line.split("=", 1)
            if len(parts) == 2:
                return parts[1].strip()
    return "n/a"


def run_case(name: str, matid: str, original: str) -> str:
    t = original
    t = replace_line(
        t,
        r"^\s*BoundL_7G\s*=.*$",
        "BoundL_7G =    1.0      1.0         1.0         1.0         1.0         1.0         1.0",
    )
    t = replace_line(
        t,
        r"^\s*BoundR_7G\s*=.*$",
        "BoundR_7G =    0.0      0.0         0.0         0.0         0.0         0.0         0.0",
    )
    t = replace_line(t, r"^\s*MatID\s*=.*$", f"    MatID = {matid}")
    INPUT.write_text(t, encoding="utf-8", newline="\n")

    env = os.environ.copy()
    env["PATH"] = r"C:\msys64\mingw64\bin" + os.pathsep + env.get("PATH", "")
    subprocess.run([str(EXE), "0"], cwd=ROOT, env=env, check=True)

    stem = f"fd_vacuum_{name}"
    moves = [
        ("fd_flux_g1_to_g7.png", f"{stem}_flux_g1_to_g7.png"),
        ("fd_cell_center_flux_g1_to_g7.png", f"{stem}_cell_center_flux_g1_to_g7.png"),
        ("fd_cell_edge_current_g1_to_g7.png", f"{stem}_cell_edge_current_g1_to_g7.png"),
    ]
    for src_name, dst_name in moves:
        shutil.move(str(OUT / src_name), str(OUT / dst_name))

    return read_keff(OUT / "summary.txt")


def main() -> int:
    if not EXE.is_file():
        raise FileNotFoundError(f"Missing {EXE}")
    original = INPUT.read_text(encoding="utf-8")

    results: list[tuple[str, str]] = []
    results.append(("UO2left_MOXright", run_case("UO2left_MOXright", matid_uo2_then_mox(), original)))
    # second run: restore BoundR from previous? run_case sets full text from original each time
    results.append(("MOXleft_UO2right", run_case("MOXleft_UO2right", matid_mox_then_uo2(), original)))

    INPUT.write_text(original, encoding="utf-8", newline="\n")
    print("Restored input_file.txt")
    for tag, k in results:
        print(f"{tag} k_eff_fd = {k}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
