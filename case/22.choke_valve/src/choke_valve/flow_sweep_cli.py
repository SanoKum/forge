from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt

from .flow import DEFAULT_METHANE
from .flow_sweep import build_flow_offset_sweep
from .geometry import load_geometry_input


REGIME_CODES = {
    "subsonic": 0,
    "choked_supersonic": 1,
    "choked_internal_shock": 2,
    "choked_external_shock": 3,
    "choked_internal_shock_unresolved": 4,
}

TITLE_FONT_SIZE = 22
LABEL_FONT_SIZE = 18
TICK_FONT_SIZE = 14


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize quasi-1D flow results over an offset sweep")
    parser.add_argument("geometry", type=Path, help="Path to geometry JSON or CSV")
    parser.add_argument("--offsets", type=float, nargs="+", required=True, help="Offset values to evaluate")
    parser.add_argument("--p0", type=float, required=True, help="Inlet stagnation pressure")
    parser.add_argument("--t0", type=float, required=True, help="Inlet stagnation temperature")
    parser.add_argument("--pb", type=float, required=True, help="Back pressure")
    parser.add_argument("--plot", type=Path, required=True, help="Output plot path")
    parser.add_argument("--mass-flow-plot", type=Path, help="Optional opening-vs-mass-flow plot path")
    parser.add_argument("--csv", type=Path, help="Optional CSV output path")
    parser.add_argument("--x-min", type=float, help="Optional sampling start x")
    parser.add_argument("--x-max", type=float, help="Optional sampling end x")
    parser.add_argument("--dx", type=float, help="Optional sampling pitch")
    args = parser.parse_args()

    geometry = load_geometry_input(args.geometry, x_min=args.x_min, x_max=args.x_max, dx=args.dx)
    results = build_flow_offset_sweep(
        geometry,
        offsets=list(args.offsets),
        p0_in=args.p0,
        t0_in=args.t0,
        p_back=args.pb,
        gas=DEFAULT_METHANE,
    )

    if args.csv is not None:
        _write_csv(results, args.csv)
    _plot_sweep(results, args.plot)
    if args.mass_flow_plot is not None:
        _plot_mass_flow_by_opening(results, args.mass_flow_plot)
    return 0


def _write_csv(results, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "offset_mm",
                "regime",
                "choked",
                "shock_x_mm",
                "throat_x_mm",
                "throat_area_mm2",
                "mass_flow_kg_s",
                "exit_mach",
                "exit_pressure",
                "critical_back_pressure",
                "ideal_exit_pressure",
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.offset,
                    result.regime,
                    result.choked,
                    result.shock_x,
                    result.throat_x,
                    result.throat_area,
                    result.mass_flow,
                    result.exit_mach,
                    result.exit_pressure,
                    result.critical_back_pressure,
                    result.ideal_exit_pressure,
                ]
            )


def _plot_sweep(results, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    openings = [abs(result.offset) for result in results]
    mass_flow = [result.mass_flow for result in results]
    shock_x = [result.shock_x for result in results]
    exit_mach = [result.exit_mach for result in results]
    regime_codes = [REGIME_CODES[result.regime] for result in results]

    figure, axes = plt.subplots(4, 1, figsize=(13, 16), sharex=True)

    axes[0].plot(openings, mass_flow, color="#1d4ed8", marker="o", linewidth=2.0)
    axes[0].set_ylabel("Mass Flow [kg/s]", fontsize=LABEL_FONT_SIZE)
    axes[0].set_title("Flow Opening Sweep Summary", fontsize=TITLE_FONT_SIZE)
    axes[0].set_ylim(bottom=0.0)
    axes[0].grid(True, alpha=0.25)
    axes[0].tick_params(labelsize=TICK_FONT_SIZE)

    axes[1].plot(openings, exit_mach, color="#0f766e", marker="o", linewidth=2.0)
    axes[1].axhline(1.0, color="#64748b", linewidth=1.0, linestyle="--")
    axes[1].set_ylabel("Exit Mach", fontsize=LABEL_FONT_SIZE)
    axes[1].grid(True, alpha=0.25)
    axes[1].tick_params(labelsize=TICK_FONT_SIZE)

    axes[2].plot(openings, [value if value is not None else float("nan") for value in shock_x], color="#dc2626", marker="o", linewidth=2.0)
    axes[2].set_ylabel("Shock x [mm]", fontsize=LABEL_FONT_SIZE)
    axes[2].grid(True, alpha=0.25)
    axes[2].tick_params(labelsize=TICK_FONT_SIZE)

    axes[3].step(openings, regime_codes, where="mid", color="#7c3aed", linewidth=2.0)
    axes[3].scatter(openings, regime_codes, color="#7c3aed", s=40)
    axes[3].set_ylabel("Regime", fontsize=LABEL_FONT_SIZE)
    axes[3].set_xlabel("Opening [mm]", fontsize=LABEL_FONT_SIZE)
    axes[3].set_yticks(list(REGIME_CODES.values()), list(REGIME_CODES.keys()))
    axes[3].grid(True, alpha=0.25)
    axes[3].tick_params(labelsize=TICK_FONT_SIZE)

    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def _plot_mass_flow_by_opening(results, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    openings = [abs(result.offset) for result in results]
    mass_flow = [result.mass_flow for result in results]

    figure, axis = plt.subplots(figsize=(12, 7))
    axis.plot(openings, mass_flow, color="#1d4ed8", marker="o", linewidth=2.5)
    axis.set_title("Mass Flow vs Opening", fontsize=TITLE_FONT_SIZE)
    axis.set_xlabel("Opening [mm]", fontsize=LABEL_FONT_SIZE)
    axis.set_ylabel("Mass Flow [kg/s]", fontsize=LABEL_FONT_SIZE)
    axis.set_ylim(bottom=0.0)
    axis.grid(True, alpha=0.25)
    axis.tick_params(labelsize=TICK_FONT_SIZE)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    raise SystemExit(main())