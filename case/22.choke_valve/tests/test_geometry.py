from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from choke_valve.flow import DEFAULT_METHANE, solve_isentropic_flow
from choke_valve.geometry import GeometryInput, Point, build_area_profile, load_geometry_input


class GeometryTests(unittest.TestCase):
    def test_area_profile_finds_throat(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 5.0), Point(5.0, 4.0), Point(10.0, 5.0)),
            needle=(Point(0.0, 0.0), Point(5.0, 2.0), Point(10.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=1.0,
        )

        profile = build_area_profile(geometry, offset=0.0)

        self.assertAlmostEqual(profile.throat_x, 5.0)
        self.assertAlmostEqual(profile.minimum_gap, 2.0)
        self.assertFalse(profile.blocked)

    def test_offset_moves_needle_into_restriction(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 3.0), Point(10.0, 3.0)),
            needle=(Point(0.0, 0.0), Point(2.0, 2.5), Point(4.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=1.0,
        )

        centered = build_area_profile(geometry, offset=3.0)
        retracted = build_area_profile(geometry, offset=6.0)

        self.assertAlmostEqual(centered.throat_area, retracted.throat_area)
        self.assertAlmostEqual(centered.throat_x, 5.0)
        self.assertAlmostEqual(retracted.throat_x, 8.0)
        self.assertGreater(centered.minimum_gap, 0.0)

    def test_load_geometry_from_csv(self) -> None:
        csv_text = """surface,x,r
wall,0.0,5.0
wall,5.0,4.0
wall,10.0,5.0
needle,0.0,0.0
needle,5.0,2.0
needle,10.0,0.0
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "points.csv"
            path.write_text(csv_text, encoding="utf-8")
            geometry = load_geometry_input(path, x_min=0.0, x_max=10.0, dx=1.0)

        self.assertEqual(len(geometry.wall), 3)
        self.assertEqual(len(geometry.needle), 3)
        self.assertEqual(geometry.dx, 1.0)

    def test_constant_area_stays_subsonic(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 3.0), Point(10.0, 3.0)),
            needle=(Point(0.0, 0.0), Point(10.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=2.0,
        )
        area_profile = build_area_profile(geometry, offset=0.0)

        solution = solve_isentropic_flow(
            area_profile=area_profile,
            p0_in=1_000_000.0,
            t0_in=300.0,
            p_back=950_000.0,
            gas=DEFAULT_METHANE,
        )

        self.assertFalse(solution.choked)
        self.assertEqual(solution.regime, "subsonic")
        self.assertFalse(solution.shock_candidate)
        self.assertTrue(all(station.branch == "subsonic" for station in solution.stations))
        self.assertAlmostEqual(solution.stations[0].mach, solution.stations[-1].mach)

    def test_converging_duct_chokes_below_critical_pressure(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 5.0), Point(10.0, 3.0)),
            needle=(Point(0.0, 0.0), Point(10.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=1.0,
        )
        area_profile = build_area_profile(geometry, offset=0.0)

        solution = solve_isentropic_flow(
            area_profile=area_profile,
            p0_in=1_000_000.0,
            t0_in=300.0,
            p_back=500_000.0,
            gas=DEFAULT_METHANE,
        )

        self.assertTrue(solution.choked)
        self.assertEqual(solution.regime, "choked_supersonic")
        self.assertGreater(solution.critical_back_pressure, 500_000.0)
        self.assertAlmostEqual(solution.stations[-1].mach, 1.0, places=6)

    def test_diverging_nozzle_can_place_normal_shock(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 5.0), Point(5.0, 3.0), Point(10.0, 5.0)),
            needle=(Point(0.0, 0.0), Point(10.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=0.5,
        )
        area_profile = build_area_profile(geometry, offset=0.0)

        solution = solve_isentropic_flow(
            area_profile=area_profile,
            p0_in=1_000_000.0,
            t0_in=300.0,
            p_back=500_000.0,
            gas=DEFAULT_METHANE,
        )

        self.assertTrue(solution.choked)
        self.assertEqual(solution.regime, "choked_internal_shock")
        self.assertTrue(solution.shock_candidate)
        self.assertIsNotNone(solution.shock_x)
        self.assertGreater(solution.shock_x, solution.throat_x)
        self.assertAlmostEqual(solution.exit_pressure, 500_000.0, delta=2_500.0)
        shock_branches = {station.branch for station in solution.stations}
        self.assertIn("shock-upstream", shock_branches)
        self.assertIn("shock-downstream", shock_branches)

    def test_diverging_nozzle_can_classify_external_shock(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 5.0), Point(5.0, 3.0), Point(10.0, 5.0)),
            needle=(Point(0.0, 0.0), Point(10.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=0.5,
        )
        area_profile = build_area_profile(geometry, offset=0.0)

        solution = solve_isentropic_flow(
            area_profile=area_profile,
            p0_in=1_000_000.0,
            t0_in=300.0,
            p_back=350_000.0,
            gas=DEFAULT_METHANE,
        )

        self.assertTrue(solution.choked)
        self.assertEqual(solution.regime, "choked_external_shock")
        self.assertFalse(solution.shock_candidate)
        self.assertIsNone(solution.shock_x)
        self.assertIsNotNone(solution.exit_shock_pressure)
        self.assertGreater(solution.exit_shock_pressure, solution.exit_pressure)
        self.assertGreater(350_000.0, solution.ideal_exit_pressure)
        self.assertLessEqual(350_000.0, solution.exit_shock_pressure)


if __name__ == "__main__":
    unittest.main()