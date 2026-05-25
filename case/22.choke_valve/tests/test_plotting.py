from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from choke_valve.flow import DEFAULT_METHANE, solve_isentropic_flow
from choke_valve.geometry import GeometryInput, Point, build_area_profile
from choke_valve.plotting import plot_results


class PlottingTests(unittest.TestCase):
    def test_plot_results_writes_png(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 5.0), Point(5.0, 4.0), Point(10.0, 5.0)),
            needle=(Point(0.0, 0.0), Point(5.0, 2.0), Point(10.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=1.0,
        )
        area_profile = build_area_profile(geometry, offset=0.0)
        solution = solve_isentropic_flow(
            area_profile=area_profile,
            p0_in=1_000_000.0,
            t0_in=300.0,
            p_back=700_000.0,
            gas=DEFAULT_METHANE,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "profile.png"
            plot_results(area_profile, solution, output)
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()