from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .flow import solve_isentropic_flow
from .geometry import build_area_profile, load_geometry_input
from .plotting import plot_results


def main() -> int:
    parser = argparse.ArgumentParser(description="Build axisymmetric area distribution from meridional points")
    parser.add_argument("geometry", type=Path, help="Path to geometry JSON or CSV")
    parser.add_argument("--offset", type=float, default=0.0, help="Needle axial shift delta")
    parser.add_argument("--output", type=Path, help="Optional CSV output path")
    parser.add_argument("--plot", type=Path, help="Optional plot output path, e.g. result.png")
    parser.add_argument("--x-min", type=float, help="Optional sampling start x; useful for CSV input")
    parser.add_argument("--x-max", type=float, help="Optional sampling end x; useful for CSV input")
    parser.add_argument("--dx", type=float, help="Optional sampling pitch; useful for CSV input")
    parser.add_argument("--p0", type=float, help="Inlet stagnation pressure")
    parser.add_argument("--t0", type=float, help="Inlet stagnation temperature")
    parser.add_argument("--pb", type=float, help="Back pressure")
    args = parser.parse_args()

    geometry = load_geometry_input(args.geometry, x_min=args.x_min, x_max=args.x_max, dx=args.dx)
    profile = build_area_profile(geometry, offset=args.offset)

    writer = _csv_writer(args.output)
    writer.writerow(["x", "wall_radius", "needle_radius", "gap", "area"])
    for sample in profile.samples:
        writer.writerow([sample.x, sample.wall_radius, sample.needle_radius, sample.gap, sample.area])

    print(f"throat_x={profile.throat_x}", file=sys.stderr)
    print(f"throat_area={profile.throat_area}", file=sys.stderr)
    print(f"minimum_gap={profile.minimum_gap}", file=sys.stderr)
    print(f"blocked={profile.blocked}", file=sys.stderr)

    solution = None
    if any(value is not None for value in (args.p0, args.t0, args.pb)):
        if None in (args.p0, args.t0, args.pb):
            parser.error("--p0, --t0, and --pb must be provided together")
        solution = solve_isentropic_flow(profile, p0_in=args.p0, t0_in=args.t0, p_back=args.pb)
        print(f"regime={solution.regime}", file=sys.stderr)
        print(f"choked={solution.choked}", file=sys.stderr)
        print(f"shock_candidate={solution.shock_candidate}", file=sys.stderr)
        print(f"shock_x={solution.shock_x}", file=sys.stderr)
        print(f"mass_flow_kg_s={solution.mass_flow}", file=sys.stderr)
        print(f"critical_back_pressure={solution.critical_back_pressure}", file=sys.stderr)
        print(f"ideal_exit_pressure={solution.ideal_exit_pressure}", file=sys.stderr)
        print(f"exit_shock_pressure={solution.exit_shock_pressure}", file=sys.stderr)
        print(f"exit_pressure={solution.exit_pressure}", file=sys.stderr)
        print(f"exit_mach={solution.stations[-1].mach}", file=sys.stderr)

    if args.plot is not None:
        plot_results(profile, solution, args.plot)
        print(f"plot_path={args.plot}", file=sys.stderr)
    return 0


def _csv_writer(output_path: Path | None) -> csv.writer:
    if output_path is None:
        return csv.writer(sys.stdout)
    handle = output_path.open("w", encoding="utf-8", newline="")
    return csv.writer(handle)


if __name__ == "__main__":
    raise SystemExit(main())