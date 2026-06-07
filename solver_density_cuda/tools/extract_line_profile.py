#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path

import h5py
import numpy as np


DEFAULT_ABS_TOL = 1.0e-2
DEFAULT_REL_TOL = 1.0e-4


def _nearest_axis_level(values: np.ndarray, target: float) -> float:
    unique = np.unique(np.round(values.astype(np.float64), 10))
    return float(unique[np.argmin(np.abs(unique - target))])


def extract_line(mesh_path: Path, result_path: Path, y_target: float, output_path: Path, z_target: float | None) -> dict:
    with h5py.File(mesh_path, "r") as mesh_h5, h5py.File(result_path, "r") as result_h5:
        centroids = mesh_h5["CELLS/centCoords"][:].reshape(-1, 3)

        x = centroids[:, 0].astype(np.float64)
        y = centroids[:, 1].astype(np.float64)
        z = centroids[:, 2].astype(np.float64)

        sampled_y = _nearest_axis_level(y, y_target)
        sampled_z = _nearest_axis_level(z, z_target if z_target is not None else float(np.median(z)))

        mask = np.isclose(y, sampled_y) & np.isclose(z, sampled_z)
        indices = np.flatnonzero(mask)
        if indices.size == 0:
            raise RuntimeError(f"No cells found on line y={sampled_y}, z={sampled_z}")

        order = np.argsort(x[indices])
        sorted_indices = indices[order]

        ux = result_h5["VALUE/Ux"][:][sorted_indices].astype(np.float64)
        uy = result_h5["VALUE/Uy"][:][sorted_indices].astype(np.float64)
        uz = result_h5["VALUE/Uz"][:][sorted_indices].astype(np.float64)
        ro = result_h5["VALUE/ro"][:][sorted_indices].astype(np.float64)
        pressure = result_h5["VALUE/P"][:][sorted_indices].astype(np.float64)
        velocity_mag = np.sqrt(ux * ux + uy * uy + uz * uz)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "x",
                "y",
                "z",
                "Ux",
                "Uy",
                "Uz",
                "U_mag",
                "ro",
                "P",
            ])
            for idx, ux_i, uy_i, uz_i, umag_i, ro_i, p_i in zip(sorted_indices, ux, uy, uz, velocity_mag, ro, pressure):
                writer.writerow([
                    f"{x[idx]:.10g}",
                    f"{y[idx]:.10g}",
                    f"{z[idx]:.10g}",
                    f"{ux_i:.10g}",
                    f"{uy_i:.10g}",
                    f"{uz_i:.10g}",
                    f"{umag_i:.10g}",
                    f"{ro_i:.10g}",
                    f"{p_i:.10g}",
                ])

    return {
        "points": int(sorted_indices.size),
        "sampled_y": sampled_y,
        "sampled_z": sampled_z,
        "x_min": float(x[sorted_indices].min()),
        "x_max": float(x[sorted_indices].max()),
    }


def compare_profiles(baseline_path: Path, candidate_path: Path, abs_tol: float, rel_tol: float) -> tuple[bool, list[str]]:
    baseline = np.genfromtxt(baseline_path, delimiter=",", names=True)
    candidate = np.genfromtxt(candidate_path, delimiter=",", names=True)

    if baseline.shape != candidate.shape:
        return False, [f"Row count mismatch: baseline={baseline.shape}, candidate={candidate.shape}"]

    messages: list[str] = []
    ok = True
    for field in ["Ux", "Uy", "Uz", "U_mag", "ro", "P"]:
        diff = np.abs(candidate[field] - baseline[field])
        scale = np.maximum(np.abs(baseline[field]), abs_tol)
        rel = diff / scale
        max_abs = float(np.max(diff))
        max_rel = float(np.max(rel))
        messages.append(
            f"{field}: max_abs_diff={max_abs:.6e}, max_nondim_rel_diff={max_rel:.6e}, scale_floor={abs_tol:.6e}"
        )
        if not np.all(rel <= rel_tol):
            ok = False

    return ok, messages


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract a line profile from Forge HDF5 output.")
    parser.add_argument("--mesh", required=True, help="Path to the mesh HDF5 file")
    parser.add_argument("--result", required=True, help="Path to the result HDF5 file")
    parser.add_argument("--y", required=True, type=float, help="Target y coordinate for the line")
    parser.add_argument("--output", required=True, help="CSV file to write")
    parser.add_argument("--z", type=float, default=None, help="Optional target z coordinate")
    parser.add_argument("--compare", help="Optional baseline CSV to compare against")
    parser.add_argument(
        "--abs-tol",
        type=float,
        default=DEFAULT_ABS_TOL,
        help="Floor used to nondimensionalize relative differences near zero baseline values",
    )
    parser.add_argument("--rel-tol", type=float, default=DEFAULT_REL_TOL, help="Maximum nondimensional relative difference for profile comparison")
    args = parser.parse_args()

    summary = extract_line(Path(args.mesh), Path(args.result), args.y, Path(args.output), args.z)
    print(
        "Extracted line profile:",
        f"points={summary['points']}",
        f"sampled_y={summary['sampled_y']}",
        f"sampled_z={summary['sampled_z']}",
        f"x_range=[{summary['x_min']}, {summary['x_max']}]",
    )

    if args.compare:
        ok, messages = compare_profiles(Path(args.compare), Path(args.output), args.abs_tol, args.rel_tol)
        for message in messages:
            print(message)
        if not ok:
            print("Profile comparison: FAILED")
            return 1
        print("Profile comparison: PASSED")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())