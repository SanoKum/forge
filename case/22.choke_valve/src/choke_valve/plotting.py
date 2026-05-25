from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt

from .flow import FlowSolution
from .geometry import AreaProfile


AREA_Y_LIMITS = (0.0, 16000.0)
MACH_Y_LIMITS = (0.0, 2.5)
PRESSURE_Y_LIMITS = (0.0, 6.0e6)
TITLE_FONT_SIZE = 22
LABEL_FONT_SIZE = 18
TICK_FONT_SIZE = 14
LEGEND_FONT_SIZE = 14
ANNOTATION_FONT_SIZE = 13


def plot_results(area_profile: AreaProfile, flow_solution: FlowSolution | None, output_path: str | Path) -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = _create_axes(flow_solution)
    x_values = [sample.x for sample in area_profile.samples]
    wall_values = [sample.wall_radius for sample in area_profile.samples]
    needle_values = [sample.needle_radius for sample in area_profile.samples]
    area_values = [sample.area for sample in area_profile.samples]
    geometry_axis = axes[0]
    geometry_axis.plot(
        x_values,
        wall_values,
        color="#0f766e",
        linewidth=2.0,
        label="r_wall",
    )
    geometry_axis.plot(
        x_values,
        needle_values,
        color="#9a3412",
        linewidth=2.0,
        label="r_needle",
    )
    geometry_axis.scatter(
        [point.x for point in area_profile.wall_points],
        [point.r for point in area_profile.wall_points],
        color="#0f766e",
        s=28,
        marker="o",
        zorder=3,
        label="wall points",
    )
    geometry_axis.scatter(
        [point.x + area_profile.offset for point in area_profile.needle_points],
        [point.r for point in area_profile.needle_points],
        color="#9a3412",
        s=28,
        marker="o",
        zorder=3,
        label="needle points",
    )
    geometry_axis.axvline(area_profile.throat_x, color="#1f2937", linewidth=1.0, linestyle=":")
    if flow_solution is not None and flow_solution.shock_x is not None:
        geometry_axis.axvline(flow_solution.shock_x, color="#dc2626", linewidth=1.0, linestyle="-.", label="Shock")
    geometry_axis.set_ylabel("Radius [mm]", fontsize=LABEL_FONT_SIZE)
    geometry_axis.set_aspect("equal", adjustable="box")
    geometry_axis.set_title(_plot_title(flow_solution), fontsize=TITLE_FONT_SIZE)
    geometry_axis.grid(True, alpha=0.25)
    geometry_axis.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0, fontsize=LEGEND_FONT_SIZE)
    geometry_axis.tick_params(labelsize=TICK_FONT_SIZE)

    area_axis = axes[1]
    area_axis.plot(x_values, area_values, color="#1d4ed8", linewidth=2.0, label="Area")
    area_axis.axvline(area_profile.throat_x, color="#1f2937", linewidth=1.0, linestyle=":")
    if flow_solution is not None and flow_solution.shock_x is not None:
        area_axis.axvline(flow_solution.shock_x, color="#dc2626", linewidth=1.0, linestyle="-.", label="Shock")
    area_axis.set_ylabel("Area [mm^2]", fontsize=LABEL_FONT_SIZE)
    area_axis.set_ylim(*AREA_Y_LIMITS)
    area_axis.grid(True, alpha=0.25)
    area_axis.legend(loc="best", fontsize=LEGEND_FONT_SIZE)
    area_axis.tick_params(labelsize=TICK_FONT_SIZE)

    if flow_solution is not None:
        station_x = [station.x for station in flow_solution.stations]
        mach_values = [station.mach for station in flow_solution.stations]
        pressure_values = [station.pressure for station in flow_solution.stations]
        branch_values = [station.branch for station in flow_solution.stations]

        mach_axis = axes[2]
        mach_axis.plot(station_x, mach_values, color="#1d4ed8", linewidth=2.0)
        mach_axis.axhline(1.0, color="#64748b", linewidth=1.0, linestyle="--")
        mach_axis.axvline(flow_solution.throat_x, color="#1f2937", linewidth=1.0, linestyle=":")
        if flow_solution.shock_x is not None:
            mach_axis.axvline(flow_solution.shock_x, color="#dc2626", linewidth=1.0, linestyle="-.")
        mach_axis.set_ylabel("Mach", fontsize=LABEL_FONT_SIZE)
        mach_axis.set_ylim(*MACH_Y_LIMITS)
        mach_axis.grid(True, alpha=0.25)
        mach_axis.tick_params(labelsize=TICK_FONT_SIZE)
        _annotate_regimes(mach_axis, station_x, mach_values, branch_values)

        pressure_axis = axes[3]
        pressure_axis.plot(station_x, pressure_values, color="#7c3aed", linewidth=2.0)
        pressure_axis.axvline(flow_solution.throat_x, color="#1f2937", linewidth=1.0, linestyle=":")
        if flow_solution.shock_x is not None:
            pressure_axis.axvline(flow_solution.shock_x, color="#dc2626", linewidth=1.0, linestyle="-.")
        pressure_axis.set_ylabel("Pressure [Pa]", fontsize=LABEL_FONT_SIZE)
        pressure_axis.set_ylim(*PRESSURE_Y_LIMITS)
        pressure_axis.set_xlabel("x [mm]", fontsize=LABEL_FONT_SIZE)
        pressure_axis.grid(True, alpha=0.25)
        pressure_axis.tick_params(labelsize=TICK_FONT_SIZE)
    else:
        area_axis.set_xlabel("x [mm]", fontsize=LABEL_FONT_SIZE)

    axes[-1].set_xlim(left=0.0, right=max(x_values))

    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def _create_axes(flow_solution: FlowSolution | None):
    rows = 4 if flow_solution is not None else 2
    figure, axes = plt.subplots(rows, 1, figsize=(13, 4.4 * rows), sharex=True)
    if rows == 1:
        axes = [axes]
    return figure, axes


def _annotate_regimes(axis, x_values: list[float], mach_values: list[float], branch_values: list[str]) -> None:
    if not x_values:
        return
    labels = []
    for x_value, mach, branch in zip(x_values, mach_values, branch_values):
        if branch == "sonic":
            labels.append((x_value, mach, "sonic"))
            break
    if branch_values and branch_values[-1] == "supersonic":
        labels.append((x_values[-1], mach_values[-1], "supersonic"))
    for x_value, mach, branch in zip(x_values, mach_values, branch_values):
        if branch == "shock-downstream":
            labels.append((x_value, mach, "shock"))
            break
    for x_value, mach, label in labels:
        axis.annotate(label, xy=(x_value, mach), xytext=(6, 8), textcoords="offset points", fontsize=ANNOTATION_FONT_SIZE)


def _plot_title(flow_solution: FlowSolution | None) -> str:
    title = "Choke Valve Geometry and Flow"
    if flow_solution is None:
        return title
    regime_label = flow_solution.regime.replace("_", " ")
    if flow_solution.regime == "choked_external_shock":
        regime_label = f"{regime_label}; no internal shock"
    return f"{title} [{regime_label}]"