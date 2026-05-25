from __future__ import annotations

from dataclasses import dataclass

from .geometry import GeometryInput, build_area_profile


@dataclass(frozen=True)
class SweepResult:
    offset: float
    throat_x: float
    throat_area: float
    minimum_gap: float
    blocked: bool


def build_offset_sweep(geometry: GeometryInput, offsets: list[float]) -> tuple[SweepResult, ...]:
    results: list[SweepResult] = []
    for offset in offsets:
        profile = build_area_profile(geometry, offset=offset)
        results.append(
            SweepResult(
                offset=offset,
                throat_x=profile.throat_x,
                throat_area=profile.throat_area,
                minimum_gap=profile.minimum_gap,
                blocked=profile.blocked,
            )
        )
    return tuple(results)