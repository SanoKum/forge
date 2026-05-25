from __future__ import annotations

import unittest

from choke_valve.geometry import GeometryInput, Point, build_area_profile
from choke_valve.sweep import build_offset_sweep


class SweepTests(unittest.TestCase):
    def test_offset_sweep_matches_individual_profiles(self) -> None:
        geometry = GeometryInput(
            wall=(Point(0.0, 5.0), Point(10.0, 3.0)),
            needle=(Point(0.0, 0.0), Point(2.0, 2.0), Point(4.0, 0.0)),
            x_min=0.0,
            x_max=10.0,
            dx=1.0,
        )
        offsets = [0.0, -1.0, -2.0]

        results = build_offset_sweep(geometry, offsets=offsets)

        self.assertEqual(len(results), 3)
        self.assertEqual([result.offset for result in results], offsets)

        for result in results:
            profile = build_area_profile(geometry, offset=result.offset)
            self.assertAlmostEqual(result.throat_x, profile.throat_x)
            self.assertAlmostEqual(result.throat_area, profile.throat_area)
            self.assertAlmostEqual(result.minimum_gap, profile.minimum_gap)
            self.assertEqual(result.blocked, profile.blocked)


if __name__ == "__main__":
    unittest.main()