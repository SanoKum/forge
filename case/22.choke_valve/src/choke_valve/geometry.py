from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, order=True)
class Point:
    x: float
    r: float


@dataclass(frozen=True)
class GeometryInput:
    wall: tuple[Point, ...]
    needle: tuple[Point, ...]
    x_min: float | None = None
    x_max: float | None = None
    dx: float = 0.1


@dataclass(frozen=True)
class AreaSample:
    x: float
    wall_radius: float
    needle_radius: float
    gap: float
    area: float


@dataclass(frozen=True)
class AreaProfile:
    offset: float
    throat_x: float
    throat_area: float
    minimum_gap: float
    blocked: bool
    wall_points: tuple[Point, ...]
    needle_points: tuple[Point, ...]
    samples: tuple[AreaSample, ...]


def load_geometry_input(
    path: str | Path,
    *,
    x_min: float | None = None,
    x_max: float | None = None,
    dx: float | None = None,
) -> GeometryInput:
    source = Path(path)
    if source.suffix.lower() == ".csv":
        return _load_geometry_from_csv(source, x_min=x_min, x_max=x_max, dx=dx)
    raw = json.loads(source.read_text(encoding="utf-8"))
    wall = _parse_points(raw["wall_points"], label="wall_points")
    needle = _parse_points(raw["needle_points"], label="needle_points")
    sample_dx = float(raw.get("dx", 0.1) if dx is None else dx)
    if sample_dx <= 0.0:
        raise ValueError("dx must be positive")

    x_min_value = raw.get("x_min") if x_min is None else x_min
    x_max_value = raw.get("x_max") if x_max is None else x_max
    if x_min_value is not None:
        x_min_value = float(x_min_value)
    if x_max_value is not None:
        x_max_value = float(x_max_value)

    return GeometryInput(wall=wall, needle=needle, x_min=x_min_value, x_max=x_max_value, dx=sample_dx)


def _load_geometry_from_csv(path: Path, *, x_min: float | None, x_max: float | None, dx: float | None) -> GeometryInput:
    wall_points: list[list[float]] = []
    needle_points: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(row for row in handle if not row.lstrip().startswith("#"))
        required_columns = {"surface", "x", "r"}
        if reader.fieldnames is None or set(reader.fieldnames) < required_columns:
            raise ValueError("CSV geometry must contain headers: surface,x,r")
        for index, row in enumerate(reader, start=2):
            surface = (row.get("surface") or "").strip().lower()
            x_value = row.get("x")
            r_value = row.get("r")
            if surface not in {"wall", "needle"}:
                raise ValueError(f"CSV row {index} has invalid surface '{surface}'; use 'wall' or 'needle'")
            if x_value is None or r_value is None:
                raise ValueError(f"CSV row {index} must define x and r")
            point = [float(x_value), float(r_value)]
            if surface == "wall":
                wall_points.append(point)
            else:
                needle_points.append(point)

    sample_dx = 0.1 if dx is None else float(dx)
    if sample_dx <= 0.0:
        raise ValueError("dx must be positive")
    return GeometryInput(
        wall=_parse_points(wall_points, label="wall_points"),
        needle=_parse_points(needle_points, label="needle_points"),
        x_min=x_min,
        x_max=x_max,
        dx=sample_dx,
    )


def build_area_profile(geometry: GeometryInput, offset: float) -> AreaProfile:
    wall = _sorted_points(geometry.wall, label="wall")
    needle = _sorted_points(geometry.needle, label="needle")
    x_min = wall[0].x if geometry.x_min is None else float(geometry.x_min)
    x_max = wall[-1].x if geometry.x_max is None else float(geometry.x_max)
    if x_max <= x_min:
        raise ValueError("x_max must be greater than x_min")

    samples: list[AreaSample] = []
    x = x_min
    while x <= x_max + geometry.dx * 0.5:
        wall_radius = _interpolate(wall, x, outside="clamp")
        needle_radius = _interpolate(needle, x - offset, outside="zero")
        gap = wall_radius - needle_radius
        area = math.pi * max(wall_radius * wall_radius - needle_radius * needle_radius, 0.0)
        samples.append(
            AreaSample(
                x=round(x, 12),
                wall_radius=wall_radius,
                needle_radius=needle_radius,
                gap=gap,
                area=area,
            )
        )
        x += geometry.dx

    throat = min(samples, key=lambda sample: sample.area)
    minimum_gap = min(sample.gap for sample in samples)
    blocked = minimum_gap <= 0.0
    return AreaProfile(
        offset=offset,
        throat_x=throat.x,
        throat_area=throat.area,
        minimum_gap=minimum_gap,
        blocked=blocked,
        wall_points=wall,
        needle_points=needle,
        samples=tuple(samples),
    )


def _parse_points(raw_points: object, label: str) -> tuple[Point, ...]:
    if not isinstance(raw_points, list) or len(raw_points) < 2:
        raise ValueError(f"{label} must contain at least two points")

    parsed: list[Point] = []
    for index, item in enumerate(raw_points):
        if not isinstance(item, list | tuple) or len(item) != 2:
            raise ValueError(f"{label}[{index}] must be a two-value list")
        x = float(item[0])
        r = float(item[1])
        if r < 0.0:
            raise ValueError(f"{label}[{index}] radius must be non-negative")
        parsed.append(Point(x=x, r=r))
    return tuple(parsed)


def _sorted_points(points: tuple[Point, ...], label: str) -> tuple[Point, ...]:
    ordered = tuple(sorted(points, key=lambda point: point.x))
    for left, right in zip(ordered, ordered[1:]):
        if math.isclose(left.x, right.x):
            raise ValueError(f"{label} contains duplicate x values at {left.x}")
    return ordered


def _interpolate(points: tuple[Point, ...], x: float, outside: str) -> float:
    first = points[0]
    last = points[-1]
    if x < first.x:
        return 0.0 if outside == "zero" else first.r
    if x > last.x:
        return 0.0 if outside == "zero" else last.r

    for left, right in zip(points, points[1:]):
        if left.x <= x <= right.x:
            if math.isclose(left.x, right.x):
                return left.r
            ratio = (x - left.x) / (right.x - left.x)
            return left.r + ratio * (right.r - left.r)

    return last.r