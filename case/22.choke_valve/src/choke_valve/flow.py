from __future__ import annotations

import math
from dataclasses import dataclass

from .geometry import AreaProfile


@dataclass(frozen=True)
class GasModel:
    gamma: float
    gas_constant: float
    name: str = "gas"


@dataclass(frozen=True)
class FlowStation:
    x: float
    area: float
    mach: float
    pressure: float
    temperature: float
    density: float
    velocity: float
    branch: str


@dataclass(frozen=True)
class FlowSolution:
    regime: str
    choked: bool
    shock_candidate: bool
    shock_x: float | None
    throat_x: float
    throat_area: float
    mass_flow: float
    critical_back_pressure: float
    ideal_exit_pressure: float
    exit_pressure: float
    exit_shock_pressure: float | None
    stations: tuple[FlowStation, ...]


DEFAULT_METHANE = GasModel(gamma=1.30, gas_constant=518.3, name="CH4")
DEFAULT_AREA_SCALE = 1.0e-6


def solve_isentropic_flow(
    area_profile: AreaProfile,
    p0_in: float,
    t0_in: float,
    p_back: float,
    gas: GasModel = DEFAULT_METHANE,
    area_scale: float = DEFAULT_AREA_SCALE,
) -> FlowSolution:
    if p0_in <= 0.0:
        raise ValueError("p0_in must be positive")
    if t0_in <= 0.0:
        raise ValueError("t0_in must be positive")
    if p_back <= 0.0:
        raise ValueError("p_back must be positive")
    if p_back >= p0_in:
        raise ValueError("p_back must be lower than p0_in")
    if area_profile.blocked:
        raise ValueError("area profile is blocked or in contact")

    gamma = gas.gamma
    area_exit = area_profile.samples[-1].area
    throat_area = area_profile.throat_area
    critical_ratio = pressure_ratio(1.0, gamma)
    area_star_candidate = math.inf
    if p_back / p0_in > critical_ratio:
        exit_mach_subsonic = mach_from_pressure_ratio(p_back / p0_in, gamma)
        if exit_mach_subsonic < 1.0:
            area_star_candidate = area_exit / area_mach_function(exit_mach_subsonic, gamma)

    if area_star_candidate <= throat_area * (1.0 + 1.0e-9):
        stations = _build_station_profile(
            area_profile=area_profile,
            area_star=area_star_candidate,
            p0_in=p0_in,
            t0_in=t0_in,
            gas=gas,
            downstream_branch="subsonic",
        )
        mass_flow = mass_flow_from_station(stations[0], gas, area_scale=area_scale)
        return FlowSolution(
            regime="subsonic",
            choked=False,
            shock_candidate=False,
            shock_x=None,
            throat_x=area_profile.throat_x,
            throat_area=throat_area,
            mass_flow=mass_flow,
            critical_back_pressure=critical_back_pressure(area_profile, p0_in, gas),
            ideal_exit_pressure=stations[-1].pressure,
            exit_pressure=stations[-1].pressure,
            exit_shock_pressure=None,
            stations=stations,
        )

    stations = _build_station_profile(
        area_profile=area_profile,
        area_star=throat_area,
        p0_in=p0_in,
        t0_in=t0_in,
        gas=gas,
        downstream_branch="supersonic",
    )
    ideal_exit_pressure = stations[-1].pressure
    exit_shock_pressure = None
    shock_x = None
    shock_candidate = False
    regime = "choked_supersonic"
    if area_exit > throat_area * (1.0 + 1.0e-9):
        exit_shock_pressure = exit_pressure_for_normal_shock(
            area_profile,
            p0_in,
            t0_in,
            gas,
            area_profile.samples[-1].x - 1.0e-6,
        )

    if area_exit > throat_area * (1.0 + 1.0e-9) and ideal_exit_pressure < p_back < critical_back_pressure(area_profile, p0_in, gas):
        shock_x = find_normal_shock_position(area_profile, p0_in, t0_in, p_back, gas)
        shock_candidate = shock_x is not None
        if shock_x is not None:
            stations = build_normal_shock_profile(area_profile, p0_in, t0_in, gas, shock_x)
            regime = "choked_internal_shock"
        elif exit_shock_pressure is not None and p_back > exit_shock_pressure:
            regime = "choked_internal_shock_unresolved"
        else:
            regime = "choked_external_shock"
    elif exit_shock_pressure is not None and ideal_exit_pressure < p_back <= exit_shock_pressure:
        regime = "choked_external_shock"

    return FlowSolution(
        regime=regime,
        choked=True,
        shock_candidate=shock_candidate,
        shock_x=shock_x,
        throat_x=area_profile.throat_x,
        throat_area=throat_area,
            mass_flow=choked_mass_flow(throat_area, p0_in, t0_in, gas, area_scale=area_scale),
        critical_back_pressure=critical_back_pressure(area_profile, p0_in, gas),
        ideal_exit_pressure=ideal_exit_pressure,
        exit_pressure=stations[-1].pressure,
        exit_shock_pressure=exit_shock_pressure,
        stations=stations,
    )


def critical_back_pressure(area_profile: AreaProfile, p0_in: float, gas: GasModel) -> float:
    exit_area = area_profile.samples[-1].area
    exit_mach_at_choke = mach_from_area_ratio(exit_area / area_profile.throat_area, gas.gamma, branch="subsonic")
    return p0_in * pressure_ratio(exit_mach_at_choke, gas.gamma)


def choked_mass_flow(
    throat_area: float,
    p0_in: float,
    t0_in: float,
    gas: GasModel,
    area_scale: float = DEFAULT_AREA_SCALE,
) -> float:
    gamma = gas.gamma
    return (
        throat_area
        * area_scale
        * p0_in
        * math.sqrt(gamma / (gas.gas_constant * t0_in))
        * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )


def mass_flow_from_station(station: FlowStation, gas: GasModel, area_scale: float = DEFAULT_AREA_SCALE) -> float:
    return station.density * station.velocity * station.area * area_scale


def mach_from_pressure_ratio(pressure_ratio_value: float, gamma: float) -> float:
    if not 0.0 < pressure_ratio_value <= 1.0:
        raise ValueError("pressure ratio must be within (0, 1]")
    exponent = (gamma - 1.0) / gamma
    return math.sqrt((2.0 / (gamma - 1.0)) * (pressure_ratio_value ** (-exponent) - 1.0))


def area_mach_function(mach: float, gamma: float) -> float:
    if mach <= 0.0:
        raise ValueError("mach must be positive")
    factor = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
    return (1.0 / mach) * factor ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))


def mach_from_area_ratio(area_ratio: float, gamma: float, branch: str) -> float:
    if area_ratio < 1.0:
        raise ValueError("area ratio must be >= 1")
    if math.isclose(area_ratio, 1.0, rel_tol=1.0e-10, abs_tol=1.0e-10):
        return 1.0
    if branch not in {"subsonic", "supersonic"}:
        raise ValueError("branch must be 'subsonic' or 'supersonic'")

    low, high = (1.0e-8, 1.0 - 1.0e-8) if branch == "subsonic" else (1.0 + 1.0e-8, 20.0)
    for _ in range(120):
        mid = 0.5 * (low + high)
        value = area_mach_function(mid, gamma)
        if branch == "subsonic":
            if value > area_ratio:
                low = mid
            else:
                high = mid
        else:
            if value > area_ratio:
                high = mid
            else:
                low = mid
    return 0.5 * (low + high)


def pressure_ratio(mach: float, gamma: float) -> float:
    return (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (-gamma / (gamma - 1.0))


def temperature_ratio(mach: float, gamma: float) -> float:
    return 1.0 / (1.0 + 0.5 * (gamma - 1.0) * mach * mach)


def normal_shock_downstream_mach(mach_upstream: float, gamma: float) -> float:
    if mach_upstream <= 1.0:
        raise ValueError("normal shock requires supersonic upstream Mach")
    numerator = 1.0 + 0.5 * (gamma - 1.0) * mach_upstream * mach_upstream
    denominator = gamma * mach_upstream * mach_upstream - 0.5 * (gamma - 1.0)
    return math.sqrt(numerator / denominator)


def normal_shock_pressure_ratio(mach_upstream: float, gamma: float) -> float:
    if mach_upstream <= 1.0:
        raise ValueError("normal shock requires supersonic upstream Mach")
    return 1.0 + (2.0 * gamma / (gamma + 1.0)) * (mach_upstream * mach_upstream - 1.0)


def find_normal_shock_position(
    area_profile: AreaProfile,
    p0_in: float,
    t0_in: float,
    p_back: float,
    gas: GasModel,
) -> float | None:
    throat_x = area_profile.throat_x
    exit_x = area_profile.samples[-1].x
    if exit_x <= throat_x:
        return None

    low = throat_x + 1.0e-6
    high = exit_x - 1.0e-6
    if high <= low:
        return None

    low_residual = exit_pressure_for_normal_shock(area_profile, p0_in, t0_in, gas, low) - p_back
    high_residual = exit_pressure_for_normal_shock(area_profile, p0_in, t0_in, gas, high) - p_back
    if low_residual == 0.0:
        return low
    if high_residual == 0.0:
        return high
    if low_residual * high_residual > 0.0:
        return None

    for _ in range(80):
        mid = 0.5 * (low + high)
        mid_residual = exit_pressure_for_normal_shock(area_profile, p0_in, t0_in, gas, mid) - p_back
        if abs(mid_residual) < 1.0e-6 * p_back:
            return mid
        if low_residual * mid_residual > 0.0:
            low = mid
            low_residual = mid_residual
        else:
            high = mid
    return 0.5 * (low + high)


def exit_pressure_for_normal_shock(
    area_profile: AreaProfile,
    p0_in: float,
    t0_in: float,
    gas: GasModel,
    shock_x: float,
) -> float:
    downstream_state = _normal_shock_state(area_profile, p0_in, t0_in, gas, shock_x)
    exit_area = area_profile.samples[-1].area
    exit_mach = mach_from_area_ratio(exit_area / downstream_state.area_star_downstream, gas.gamma, branch="subsonic")
    return downstream_state.p0_downstream * pressure_ratio(exit_mach, gas.gamma)


@dataclass(frozen=True)
class _ShockState:
    shock_x: float
    shock_area: float
    mach_upstream: float
    mach_downstream: float
    pressure_upstream: float
    pressure_downstream: float
    temperature_upstream: float
    temperature_downstream: float
    p0_downstream: float
    area_star_downstream: float


def _normal_shock_state(
    area_profile: AreaProfile,
    p0_in: float,
    t0_in: float,
    gas: GasModel,
    shock_x: float,
) -> _ShockState:
    shock_area = _interpolate_area(area_profile, shock_x)
    mach_upstream = mach_from_area_ratio(shock_area / area_profile.throat_area, gas.gamma, branch="supersonic")
    temp_upstream = t0_in * temperature_ratio(mach_upstream, gas.gamma)
    pressure_upstream = p0_in * pressure_ratio(mach_upstream, gas.gamma)
    pressure_jump = normal_shock_pressure_ratio(mach_upstream, gas.gamma)
    mach_downstream = normal_shock_downstream_mach(mach_upstream, gas.gamma)
    pressure_downstream = pressure_upstream * pressure_jump
    temp_downstream = t0_in * temperature_ratio(mach_downstream, gas.gamma)
    p0_downstream = pressure_downstream / pressure_ratio(mach_downstream, gas.gamma)
    area_star_downstream = shock_area / area_mach_function(mach_downstream, gas.gamma)
    return _ShockState(
        shock_x=shock_x,
        shock_area=shock_area,
        mach_upstream=mach_upstream,
        mach_downstream=mach_downstream,
        pressure_upstream=pressure_upstream,
        pressure_downstream=pressure_downstream,
        temperature_upstream=temp_upstream,
        temperature_downstream=temp_downstream,
        p0_downstream=p0_downstream,
        area_star_downstream=area_star_downstream,
    )


def build_normal_shock_profile(
    area_profile: AreaProfile,
    p0_in: float,
    t0_in: float,
    gas: GasModel,
    shock_x: float,
) -> tuple[FlowStation, ...]:
    shock_state = _normal_shock_state(area_profile, p0_in, t0_in, gas, shock_x)
    stations: list[FlowStation] = []
    throat_crossed = False
    for sample in area_profile.samples:
        if sample.x > shock_x:
            break
        if math.isclose(sample.area / area_profile.throat_area, 1.0, rel_tol=1.0e-9, abs_tol=1.0e-9):
            mach = 1.0
            branch = "sonic"
            throat_crossed = True
        else:
            branch = "subsonic"
            if throat_crossed or sample.x > area_profile.throat_x:
                branch = "supersonic"
            mach = mach_from_area_ratio(sample.area / area_profile.throat_area, gas.gamma, branch)
        stations.append(_make_station(sample.x, sample.area, mach, p0_in, t0_in, gas, branch))

    stations.append(
        FlowStation(
            x=shock_x,
            area=shock_state.shock_area,
            mach=shock_state.mach_upstream,
            pressure=shock_state.pressure_upstream,
            temperature=shock_state.temperature_upstream,
            density=shock_state.pressure_upstream / (gas.gas_constant * shock_state.temperature_upstream),
            velocity=shock_state.mach_upstream * math.sqrt(gas.gamma * gas.gas_constant * shock_state.temperature_upstream),
            branch="shock-upstream",
        )
    )
    stations.append(
        FlowStation(
            x=shock_x,
            area=shock_state.shock_area,
            mach=shock_state.mach_downstream,
            pressure=shock_state.pressure_downstream,
            temperature=shock_state.temperature_downstream,
            density=shock_state.pressure_downstream / (gas.gas_constant * shock_state.temperature_downstream),
            velocity=shock_state.mach_downstream * math.sqrt(gas.gamma * gas.gas_constant * shock_state.temperature_downstream),
            branch="shock-downstream",
        )
    )

    for sample in area_profile.samples:
        if sample.x <= shock_x:
            continue
        mach = mach_from_area_ratio(sample.area / shock_state.area_star_downstream, gas.gamma, branch="subsonic")
        stations.append(_make_station(sample.x, sample.area, mach, shock_state.p0_downstream, t0_in, gas, "subsonic"))
    return tuple(stations)


def _build_station_profile(
    area_profile: AreaProfile,
    area_star: float,
    p0_in: float,
    t0_in: float,
    gas: GasModel,
    downstream_branch: str,
) -> tuple[FlowStation, ...]:
    stations: list[FlowStation] = []
    throat_crossed = False
    for sample in area_profile.samples:
        area_ratio = sample.area / area_star
        if math.isclose(area_ratio, 1.0, rel_tol=1.0e-9, abs_tol=1.0e-9):
            mach = 1.0
            throat_crossed = True
            branch = "sonic"
        else:
            branch = "subsonic"
            if throat_crossed or sample.x > area_profile.throat_x:
                branch = downstream_branch
            mach = mach_from_area_ratio(area_ratio, gas.gamma, branch)

        temp = t0_in * temperature_ratio(mach, gas.gamma)
        pressure = p0_in * pressure_ratio(mach, gas.gamma)
        density = pressure / (gas.gas_constant * temp)
        velocity = mach * math.sqrt(gas.gamma * gas.gas_constant * temp)
        stations.append(_make_station(sample.x, sample.area, mach, p0_in, t0_in, gas, branch))
    return tuple(stations)


def _make_station(x: float, area: float, mach: float, p0: float, t0: float, gas: GasModel, branch: str) -> FlowStation:
    temp = t0 * temperature_ratio(mach, gas.gamma)
    pressure = p0 * pressure_ratio(mach, gas.gamma)
    density = pressure / (gas.gas_constant * temp)
    velocity = mach * math.sqrt(gas.gamma * gas.gas_constant * temp)
    return FlowStation(
        x=x,
        area=area,
        mach=mach,
        pressure=pressure,
        temperature=temp,
        density=density,
        velocity=velocity,
        branch=branch,
    )


def _interpolate_area(area_profile: AreaProfile, x: float) -> float:
    samples = area_profile.samples
    if x <= samples[0].x:
        return samples[0].area
    if x >= samples[-1].x:
        return samples[-1].area
    for left, right in zip(samples, samples[1:]):
        if left.x <= x <= right.x:
            if math.isclose(left.x, right.x):
                return left.area
            ratio = (x - left.x) / (right.x - left.x)
            return left.area + ratio * (right.area - left.area)
    return samples[-1].area