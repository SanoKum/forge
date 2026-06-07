#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np


DEFAULT_FIELD_ORDER = [
    "Ps",
    "Ts",
    "ro",
    "Pt",
    "Tt",
    "Ux",
    "Uy",
    "Uz",
    "roUx",
    "roUy",
    "roUz",
    "roe",
    "twall_x",
    "twall_y",
    "twall_z",
    "ypls",
]

XDMF_MIXED_NODE_COUNT = {
    2: 2,
    4: 3,
    5: 4,
}


def _parse_fields(requested: str | None, available: set[str]) -> list[str]:
    if requested is None:
        ordered = [field for field in DEFAULT_FIELD_ORDER if field in available]
        extras = sorted(available - set(ordered))
        return ordered + extras

    fields = [field.strip() for field in requested.split(",") if field.strip()]
    missing = [field for field in fields if field not in available]
    if missing:
        raise ValueError(f"Requested fields are not available: {', '.join(missing)}")
    return fields


def _mixed_element_centroids(coords: np.ndarray, connectivity: np.ndarray) -> np.ndarray:
    centroids: list[np.ndarray] = []
    index = 0
    while index < connectivity.size:
        element_type = int(connectivity[index])
        node_count = XDMF_MIXED_NODE_COUNT.get(element_type)
        if node_count is None:
            raise ValueError(f"Unsupported XDMF mixed element type: {element_type}")
        node_ids = connectivity[index + 1 : index + 1 + node_count].astype(np.int64)
        centroids.append(coords[node_ids].mean(axis=0))
        index += 1 + node_count

    return np.asarray(centroids, dtype=np.float64)


def extract_wall_surface(input_path: Path, output_path: Path, fields: str | None, sort_axis: str | None) -> dict[str, float | int | str]:
    with h5py.File(input_path, "r") as h5_file:
        coords = h5_file["MESH/COORD"][:].reshape(-1, 3).astype(np.float64)
        connectivity = h5_file["MESH/CONNE"][:]
        value_group = h5_file["VALUE"]
        available_fields = set(value_group.keys())
        selected_fields = _parse_fields(fields, available_fields)

        centroids = _mixed_element_centroids(coords, connectivity)
        row_count = centroids.shape[0]

        if selected_fields:
            first_field_size = int(value_group[selected_fields[0]].shape[0])
            if first_field_size != row_count:
                raise ValueError(
                    f"Centroid count ({row_count}) does not match VALUE/{selected_fields[0]} length ({first_field_size})"
                )

        order = np.arange(row_count)
        if sort_axis is not None:
            axis_index = {"x": 0, "y": 1, "z": 2}[sort_axis]
            order = np.argsort(centroids[:, axis_index], kind="stable")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["x", "y", "z", *selected_fields])
            for cell_index in order:
                row = [
                    f"{centroids[cell_index, 0]:.10g}",
                    f"{centroids[cell_index, 1]:.10g}",
                    f"{centroids[cell_index, 2]:.10g}",
                ]
                for field in selected_fields:
                    row.append(f"{float(value_group[field][cell_index]):.10g}")
                writer.writerow(row)

    return {
        "rows": row_count,
        "fields": len(selected_fields),
        "output": str(output_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Export Forge wall-surface HDF5 results to CSV.")
    parser.add_argument("input", help="Path to a res_wall_*.h5 file")
    parser.add_argument("-o", "--output", help="Output CSV path. Defaults to <input>.csv")
    parser.add_argument(
        "--fields",
        help="Comma-separated VALUE fields to export. Defaults to common wall fields plus any remaining datasets.",
    )
    parser.add_argument(
        "--sort-axis",
        choices=["x", "y", "z", "none"],
        default="x",
        help="Sort CSV rows by centroid coordinate along the chosen axis. Use none to keep original order.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path.with_suffix(".csv")
    sort_axis = None if args.sort_axis == "none" else args.sort_axis

    summary = extract_wall_surface(input_path, output_path, args.fields, sort_axis)
    print(
        "Exported wall surface CSV:",
        f"rows={summary['rows']}",
        f"fields={summary['fields']}",
        f"output={summary['output']}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())