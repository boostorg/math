#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "matplotlib>=3.8",
# ]
# ///

import argparse
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TRACE = (
    REPOSITORY_ROOT
    / "build"
    / "vdb"
    / "van_den_bos_unit_square_convergence.csv"
)
DEFAULT_OUTPUT = DEFAULT_TRACE.with_suffix(".png")


def parse_value(s):
    s = s.strip()
    if s.startswith("(") and s.endswith(")") and "," in s:
        a, b = s[1:-1].split(",", 1)
        return complex(float(a), float(b))
    return float(s)

def parse_trace(path):
    cases = {}
    current = None

    for raw in Path(path).read_text().splitlines():
        if raw.startswith("VDB_CASE,"):
            match = re.fullmatch(r"VDB_CASE,name=([^,]+),exact=(.+)", raw)
            if not match:
                raise ValueError(f"cannot parse case record: {raw}")
            current = match.group(1)
            cases[current] = {
                "exact": parse_value(match.group(2)),
                "levels": [],
            }
        elif raw.startswith("VDB_LEVEL,") and current is not None:
            match = re.fullmatch(
                r"VDB_LEVEL,level=(\d+),degree=(\d+),order=(\d+),"
                r"nodes=(\d+),evaluations=(\d+),value=(\([^)]*\)|[^,]+),"
                r"delta=([^,]+),error_estimate=(.+)",
                raw,
            )
            if not match:
                raise ValueError(f"cannot parse level record: {raw}")

            value = parse_value(match.group(6))
            exact = cases[current]["exact"]
            actual_error = abs(value - exact)
            error_text = match.group(8)
            cases[current]["levels"].append({
                "level": int(match.group(1)),
                "degree": int(match.group(2)),
                "order": int(match.group(3)),
                "nodes": int(match.group(4)),
                "evaluations": int(match.group(5)),
                "value": value,
                "actual_error": actual_error,
                "estimated_error": float(error_text)
                    if error_text != "nan"
                    else math.nan,
            })
        elif raw.startswith("VDB_END,"):
            current = None

    return cases

def main():
    ap = argparse.ArgumentParser(
        description="Plot van den Bos unit-square convergence traces.",
    )
    ap.add_argument(
        "trace",
        nargs="?",
        type=Path,
        default=DEFAULT_TRACE,
        help=f"trace CSV (default: {DEFAULT_TRACE})",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output image (default: {DEFAULT_OUTPUT})",
    )
    ap.add_argument(
        "--separate",
        action="store_true",
        help="also write one convergence plot per integrand",
    )
    args = ap.parse_args()

    cases = parse_trace(args.trace)
    if not cases:
        raise ValueError(f"no van den Bos trace records found in {args.trace}")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots()

    for name, case in cases.items():
        xs = [row["evaluations"] for row in case["levels"]]
        ys = [max(row["actual_error"], 1e-18) for row in case["levels"]]
        ax.plot(xs, ys, marker="o", label=name)

    ax.set_yscale("log")
    ax.set_xlabel("Function evaluations")
    ax.set_ylabel("Absolute error")
    ax.set_title("van den Bos unit-square convergence")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    plt.close(fig)
    print(f"wrote {args.out}")

    if args.separate:
        out = args.out
        for name, case in cases.items():
            fig, ax = plt.subplots()
            xs = [row["evaluations"] for row in case["levels"]]
            ys = [max(row["actual_error"], 1e-18) for row in case["levels"]]
            estimates = [
                max(row["estimated_error"], 1e-18)
                if math.isfinite(row["estimated_error"])
                else math.nan
                for row in case["levels"]
            ]
            ax.plot(xs, ys, marker="o", label="actual error")
            ax.plot(xs, estimates, marker="x", label="estimated error")
            ax.set_yscale("log")
            ax.set_xlabel("Function evaluations")
            ax.set_ylabel("Absolute error")
            ax.set_title(name)
            ax.grid(True, which="both", alpha=0.25)
            ax.legend()
            fig.tight_layout()
            fig.savefig(
                out.with_name(f"{out.stem}_{name}{out.suffix}"),
                dpi=180,
            )
            plt.close(fig)

if __name__ == "__main__":
    main()
