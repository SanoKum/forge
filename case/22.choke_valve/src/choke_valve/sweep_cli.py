from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt

from .geometry import load_geometry_input
from .sweep import build_offset_sweep


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize throat metrics over an offset sweep")
    parser.add_argument("geometry", type=Path, help="Path to geometry JSON or CSV")
    parser.add_argument("--offsets", type=float, nargs="+", required=True, help="Offset values to evaluate")
    parser.add_argument("--plot", type=Path, required=True, help="Output plot path")
    parser.add_argument("--csv", type=Path, help="Optional CSV output path for the sweep summary")
    parser.add_argument("--x-min", type=float, help="Optional sampling start x")
    parser.add_argument("--x-max", type=float, help="Optional sampling end x")
    parser.add_argument("--dx", type=float, help="Optional sampling pitch")
    args = parser.parse_args()

    geometry = load_geometry_input(args.geometry, x_min=args.x_min, x_max=args.x_max, dx=args.dx)
    results = build_offset_sweep(geometry, offsets=list(args.offsets))

    if args.csv is not None:
        _write_csv(results, args.csv)
    _plot_sweep(results, args.plot)
    return 0


def _write_csv(results, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["offset", "throat_x", "throat_area", "minimum_gap", "blocked"])
        for result in results:
            writer.writerow([result.offset, result.throat_x, result.throat_area, result.minimum_gap, result.blocked])


def _plot_sweep(results, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    offsets = [result.offset for result in results]
    throat_area = [result.throat_area for result in results]
    minimum_gap = [result.minimum_gap for result in results]
    throat_x = [result.throat_x for result in results]
    blocked = [result.blocked for result in results]

    figure, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(offsets, throat_area, color="#1d4ed8", marker="o", linewidth=2.0)
    axes[0].set_ylabel("Throat Area")
    axes[0].grid(True, alpha=0.25)
    axes[0].set_title("Offset Sweep Summary")

    axes[1].plot(offsets, minimum_gap, color="#9a3412", marker="o", linewidth=2.0)
    axes[1].set_ylabel("Minimum Gap")
    axes[1].grid(True, alpha=0.25)

    axes[2].plot(offsets, throat_x, color="#0f766e", marker="o", linewidth=2.0)
    axes[2].set_ylabel("Throat x")
    axes[2].set_xlabel("Offset")
    axes[2].grid(True, alpha=0.25)

    for axis in axes:
        for offset, is_blocked in zip(offsets, blocked):
            if is_blocked:
                axis.axvline(offset, color="#dc2626", linestyle=":", linewidth=1.0)

    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    raise SystemExit(main())