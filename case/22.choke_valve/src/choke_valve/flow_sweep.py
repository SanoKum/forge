from __future__ import annotations

from dataclasses import dataclass

from .flow import FlowSolution, GasModel, solve_isentropic_flow
from .geometry import GeometryInput, build_area_profile


@dataclass(frozen=True)
class FlowSweepResult:
    offset: float
    regime: str
    choked: bool
    shock_x: float | None
    throat_x: float
    throat_area: float
    mass_flow: float
    exit_mach: float
    exit_pressure: float
    critical_back_pressure: float
    ideal_exit_pressure: float


def build_flow_offset_sweep(
    geometry: GeometryInput,
    offsets: list[float],
    *,
    p0_in: float,
    t0_in: float,
    p_back: float,
    gas: GasModel,
) -> tuple[FlowSweepResult, ...]:
    results: list[FlowSweepResult] = []
    for offset in offsets:
        profile = build_area_profile(geometry, offset=offset)
        solution = solve_isentropic_flow(profile, p0_in=p0_in, t0_in=t0_in, p_back=p_back, gas=gas)
        results.append(_to_flow_sweep_result(offset, solution))
    return tuple(results)


def _to_flow_sweep_result(offset: float, solution: FlowSolution) -> FlowSweepResult:
    return FlowSweepResult(
        offset=offset,
        regime=solution.regime,
        choked=solution.choked,
        shock_x=solution.shock_x,
        throat_x=solution.throat_x,
        throat_area=solution.throat_area,
        mass_flow=solution.mass_flow,
        exit_mach=solution.stations[-1].mach,
        exit_pressure=solution.exit_pressure,
        critical_back_pressure=solution.critical_back_pressure,
        ideal_exit_pressure=solution.ideal_exit_pressure,
    )