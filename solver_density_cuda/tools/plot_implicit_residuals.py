#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


EQUATIONS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"]
PHASE_ORDER = {
    "outer_begin": 0,
    "inner_begin": 1,
    "inner_iter": 2,
    "outer_end": 3,
}


def _load_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise RuntimeError(f"No rows found in {csv_path}")
    return rows


def _filter_rows(rows: list[dict[str, str]], phases: set[str] | None) -> list[dict[str, str]]:
    if phases is None:
        return rows
    filtered = [row for row in rows if row.get("phase", "") in phases]
    if not filtered:
        raise RuntimeError(f"No rows matched phases: {sorted(phases)}")
    return filtered


def _total_step_values(rows: list[dict[str, str]]) -> list[int]:
    phase_offset = {
        "outer_begin": 0,
        "inner_begin": 1,
        "inner_iter": 1,
        "outer_end": 2,
    }

    max_inner_per_step: dict[int, int] = {}
    for row in rows:
        step = int(row.get("step", "0"))
        inner = int(row.get("inner", "-1"))
        if inner >= 0:
            max_inner_per_step[step] = max(max_inner_per_step.get(step, -1), inner)

    default_inner_count = max(max_inner_per_step.values(), default=0) + 1
    x_values: list[int] = []

    for row in rows:
        step = int(row.get("step", "0"))
        inner = int(row.get("inner", "-1"))
        phase = row.get("phase", "")
        inner_count = max_inner_per_step.get(step, default_inner_count - 1) + 1
        base = step * (inner_count + 2)

        if phase == "outer_begin":
            x_values.append(base)
        elif phase == "outer_end":
            x_values.append(base + inner_count + 1)
        else:
            x_values.append(base + phase_offset.get(phase, 1) + max(inner, 0))

    return x_values


def _column_values(rows: list[dict[str, str]], column: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        raw = row.get(column)
        if raw is None or raw == "":
            values.append(float("nan"))
            continue
        values.append(float(raw))
    return values


def _residual_column(rows: list[dict[str, str]], equation: str) -> str:
    rms_column = f"rms_{equation}"
    if rms_column in rows[0]:
        return rms_column
    return f"abs_{equation}"


def _correction_column(rows: list[dict[str, str]], equation: str) -> str | None:
    rms_column = f"rms_dq_{equation}"
    if rms_column in rows[0]:
        return rms_column
    return None


def _relative_values(rows: list[dict[str, str]], equation: str, relative_mode: str) -> list[float]:
    column = f"rel_{relative_mode}_{equation}"
    if column in rows[0]:
        return _column_values(rows, column)

    residual_column = _residual_column(rows, equation)
    residual_values = _column_values(rows, residual_column)
    eps = 1.0e-30

    if relative_mode == "init3":
        references = [
            float(row[residual_column])
            for row in rows
            if row.get("phase") == "outer_begin"
        ]
        if not references:
            raise RuntimeError("No outer_begin rows found to compute init3 relative residual")
        count = min(3, len(references))
        reference = sum(references[:count]) / count
        scale = max(reference, eps)
        return [value / scale for value in residual_values]

    if relative_mode == "inner":
        step_reference: dict[int, float] = {}
        for row in rows:
            if row.get("phase") != "inner_begin":
                continue
            step_reference[int(row["step"])] = float(row[residual_column])

        values: list[float] = []
        for row, value in zip(rows, residual_values):
            step = int(row["step"])
            reference = step_reference.get(step)
            if reference is None:
                reference = value
            values.append(value / max(reference, eps))
        return values

    raise RuntimeError(f"Unsupported relative mode: {relative_mode}")


def _relative_correction_values(rows: list[dict[str, str]], equation: str) -> list[float]:
    correction_column = _correction_column(rows, equation)
    if correction_column is None:
        return [float("nan")] * len(rows)

    correction_values = _column_values(rows, correction_column)
    eps = 1.0e-30
    step_reference: dict[int, float] = {}

    for row, value in zip(rows, correction_values):
        if row.get("phase") != "inner_iter":
            continue
        step = int(row.get("step", "0"))
        if step not in step_reference and value > 0.0:
            step_reference[step] = value

    values: list[float] = []
    for row, value in zip(rows, correction_values):
        phase = row.get("phase", "")
        if phase in {"outer_begin", "inner_begin"}:
            values.append(float("nan"))
            continue

        step = int(row.get("step", "0"))
        reference = step_reference.get(step)
        if reference is None:
            values.append(float("nan"))
            continue

        values.append(value / max(reference, eps))

    return values


def _has_correction_columns(rows: list[dict[str, str]], equations: list[str]) -> bool:
    return any(_correction_column(rows, equation) is not None for equation in equations)


def _equation_has_data(rows: list[dict[str, str]], equation: str) -> bool:
    header = rows[0]
    return f"rms_{equation}" in header or f"abs_{equation}" in header


def _parse_equations(raw: str | None) -> list[str]:
    if not raw:
        return EQUATIONS[:]
    equations = [item.strip() for item in raw.split(",") if item.strip()]
    invalid = [item for item in equations if item not in EQUATIONS]
    if invalid:
        raise RuntimeError(f"Unsupported equations: {', '.join(invalid)}")
    return equations


def _parse_phases(raw: str | None) -> set[str] | None:
    if not raw:
        return None
    phases = {item.strip() for item in raw.split(",") if item.strip()}
    invalid = [item for item in phases if item not in PHASE_ORDER]
    if invalid:
        raise RuntimeError(f"Unsupported phases: {', '.join(sorted(invalid))}")
    return phases


def _sort_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    def key(row: dict[str, str]) -> tuple[int, int, int]:
        step = int(row.get("step", "0"))
        inner = int(row.get("inner", "-1"))
        phase_name = row.get("phase", "")

        if phase_name == "outer_begin":
            return step, 0, 0
        if phase_name == "inner_begin":
            return step, 1, 0
        if phase_name == "inner_iter":
            return step, 2, max(inner, 0)
        if phase_name == "outer_end":
            return step, 3, 0

        return step, 99, max(inner, -1)

    return sorted(rows, key=key)


def plot_history(
    csv_path: Path,
    output_path: Path,
    equations: list[str],
    phases: set[str] | None,
    relative_mode: str,
    title: str | None,
) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to plot residual history") from exc

    rows = _sort_rows(_filter_rows(_load_rows(csv_path), phases))
    x_values = _total_step_values(rows)

    missing = [eq for eq in equations if not _equation_has_data(rows, eq)]
    if missing:
        print(f"Skipping equations absent from CSV: {', '.join(missing)}")
    equations = [eq for eq in equations if _equation_has_data(rows, eq)]
    if not equations:
        raise RuntimeError("No requested equations have residual columns in the CSV")

    show_corrections = _has_correction_columns(rows, equations)
    if show_corrections:
        figure, axes = plt.subplots(4, 1, figsize=(14, 14), sharex=True)
        abs_ax, rel_ax, dq_abs_ax, dq_rel_ax = axes
    else:
        figure, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
        abs_ax, rel_ax = axes
        dq_abs_ax = None
        dq_rel_ax = None

    for equation in equations:
        abs_ax.plot(
            x_values,
            _column_values(rows, _residual_column(rows, equation)),
            marker="o",
            markersize=3,
            linewidth=1.2,
            label=equation,
        )
        rel_ax.plot(
            x_values,
            _relative_values(rows, equation, relative_mode),
            marker="o",
            markersize=3,
            linewidth=1.2,
            label=equation,
        )

        if dq_abs_ax is not None and dq_rel_ax is not None:
            correction_column = _correction_column(rows, equation)
            if correction_column is None:
                continue

            dq_abs_ax.plot(
                x_values,
                _column_values(rows, correction_column),
                marker="o",
                markersize=3,
                linewidth=1.2,
                label=equation,
            )
            dq_rel_ax.plot(
                x_values,
                _relative_correction_values(rows, equation),
                marker="o",
                markersize=3,
                linewidth=1.2,
                label=equation,
            )

    def has_positive_line_data(axis) -> bool:
        for line in axis.get_lines():
            y_values = np.asarray(line.get_ydata(), dtype=float)
            if np.any(np.isfinite(y_values) & (y_values > 0.0)):
                return True
        return False

    abs_ax.set_yscale("log")
    rel_ax.set_yscale("log")
    abs_ax.set_ylabel("RMS residual")
    rel_ax.set_ylabel(f"Relative residual ({relative_mode})")
    abs_ax.grid(True, which="both", linestyle="--", alpha=0.35)
    rel_ax.grid(True, which="both", linestyle="--", alpha=0.35)
    abs_ax.legend(loc="best")
    rel_ax.legend(loc="best")

    if dq_abs_ax is not None and dq_rel_ax is not None:
        if has_positive_line_data(dq_abs_ax):
            dq_abs_ax.set_yscale("log")
        if has_positive_line_data(dq_rel_ax):
            dq_rel_ax.set_yscale("log")
        dq_abs_ax.set_ylabel("RMS delta Q")
        dq_rel_ax.set_ylabel("Relative delta Q (inner1)")
        dq_abs_ax.grid(True, which="both", linestyle="--", alpha=0.35)
        dq_rel_ax.grid(True, which="both", linestyle="--", alpha=0.35)
        dq_abs_ax.legend(loc="best")
        dq_rel_ax.legend(loc="best")
        dq_rel_ax.set_xlabel("Total step")
    else:
        rel_ax.set_xlabel("Total step")

    rel_ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    if dq_rel_ax is not None:
        dq_rel_ax.ticklabel_format(axis="x", style="plain", useOffset=False)

    figure.suptitle(title or f"Residual history: {csv_path.name}")
    figure.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot Forge residual history from CSV output.")
    parser.add_argument("--input", required=True, help="Path to residual_history.csv")
    parser.add_argument("--output", default="residual_history.png", help="Path to output image")
    parser.add_argument("--equations", help="Comma-separated subset of equations: ro,roUx,roUy,roUz,roe")
    parser.add_argument("--phases", help="Comma-separated subset of phases: outer_begin,inner_begin,inner_iter,outer_end")
    parser.add_argument(
        "--relative-mode",
        choices=["init3", "inner"],
        default="init3",
        help="Relative residual family to draw in the lower panel; computed from RMS values when absent in CSV",
    )
    parser.add_argument("--title", help="Optional plot title")
    args = parser.parse_args()

    equations = _parse_equations(args.equations)
    phases = _parse_phases(args.phases)
    plot_history(
        csv_path=Path(args.input),
        output_path=Path(args.output),
        equations=equations,
        phases=phases,
        relative_mode=args.relative_mode,
        title=args.title,
    )

    print(f"Saved residual plot to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())