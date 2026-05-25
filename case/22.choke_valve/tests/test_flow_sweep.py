from __future__ import annotations

import unittest

from choke_valve.flow import DEFAULT_METHANE, solve_isentropic_flow
from choke_valve.flow_sweep import build_flow_offset_sweep
from choke_valve.geometry import GeometryInput, Point, build_area_profile


class FlowSweepTests(unittest.TestCase):
    def test_flow_sweep_matches_individual_solutions(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 5.0), Point(5.0, 3.0), Point(10.0, 5.0)),
            needle=(Point(0.0, 0.0), Point(10.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=0.5,
        )
        offsets = [0.0, -1.0]
        results = build_flow_offset_sweep(
            geometry,
            offsets=offsets,
            p0_in=1_000_000.0,
            t0_in=300.0,
            p_back=500_000.0,
            gas=DEFAULT_METHANE,
        )

        self.assertEqual(len(results), 2)
        for result in results:
            profile = build_area_profile(geometry, offset=result.offset)
            solution = solve_isentropic_flow(profile, 1_000_000.0, 300.0, 500_000.0, DEFAULT_METHANE)
            self.assertEqual(result.regime, solution.regime)
            self.assertEqual(result.choked, solution.choked)
            self.assertEqual(result.shock_x, solution.shock_x)
            self.assertAlmostEqual(result.mass_flow, solution.mass_flow)
            self.assertAlmostEqual(result.exit_mach, solution.stations[-1].mach)


if __name__ == "__main__":
    unittest.main()