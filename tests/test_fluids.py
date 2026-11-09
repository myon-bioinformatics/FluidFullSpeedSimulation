import math
import random
import unittest

from fixedq import random_scales
from fluids import (
    FLUIDS, Section, channel_flow, fluid_table, get_fluid, normal_depth, poiseuille_tau, torricelli_empty_time,
    torricelli_head,
)


class Properties(unittest.TestCase):
    def test_kinematic_viscosity(self):
        self.assertAlmostEqual(get_fluid("water").nu * 1e6, 1.0034, places=3)
        self.assertAlmostEqual(get_fluid("ethanol").nu * 1e6, 1.520, places=2)
        self.assertAlmostEqual(get_fluid("Methanol").nu * 1e6, 0.742, places=3)

    def test_ranges_and_sources_are_recorded(self):
        for f in FLUIDS.values():
            lo, hi = f.mu_range
            self.assertTrue(lo <= f.mu <= hi, f.key)
            self.assertTrue(f.sources)
        self.assertEqual({row["key"] for row in fluid_table()}, set(FLUIDS))

    def test_unknown_fluid(self):
        with self.assertRaises(KeyError):
            get_fluid("mercury")


class ClosedForms(unittest.TestCase):
    def test_torricelli_empties_at_t_empty(self):
        a = math.pi * 0.002 ** 2
        t = torricelli_empty_time(0.3, 0.01, a, 0.61)
        self.assertEqual(torricelli_head(0.3, t, 0.01, a, 0.61), 0.0)
        self.assertAlmostEqual(torricelli_head(0.3, 0.0, 0.01, a, 0.61), 0.3)

    def test_poiseuille_tau_scales_with_nu(self):
        w, e = get_fluid("water"), get_fluid("ethanol")
        ratio = poiseuille_tau(e.nu, 0.01, 1e-3, 0.5) / poiseuille_tau(w.nu, 0.01, 1e-3, 0.5)
        self.assertAlmostEqual(ratio, e.nu / w.nu, places=12)


class Sections(unittest.TestCase):
    def test_depth_inverts_area(self):
        rng = random.Random(1)
        for sec, cap in ((Section("rect", 0.05), 0.05), (Section("semicircle", 0.04), 0.04)):
            for _ in range(200):
                h = rng.uniform(0, cap)
                self.assertAlmostEqual(sec.depth(sec.area(h)), h, places=10)
        self.assertAlmostEqual(Section("semicircle", 0.04).area(0.04), math.pi * 0.04 ** 2 / 2)

    def test_bad_section(self):
        with self.assertRaises(ValueError):
            Section("triangle", 1)
        with self.assertRaises(ValueError):
            Section("rect", 0.1).max_area()

    def test_manning_ignores_viscosity_laminar_does_not(self):
        sec = Section("rect", 0.05)
        for item in random_scales(9, 5, n_min=3, n_max=4):
            h = item["value"]
            m1 = channel_flow(sec, h, 0.05, "manning", nu=1e-6)
            m2 = channel_flow(sec, h, 0.05, "manning", nu=2e-6)
            self.assertEqual(m1, m2)
            l1 = channel_flow(sec, h, 0.05, "laminar", nu=1e-6)
            l2 = channel_flow(sec, h, 0.05, "laminar", nu=2e-6)
            self.assertAlmostEqual(l1 / l2, 2.0, places=12)
        with self.assertRaises(ValueError):
            channel_flow(sec, 0.01, 0.05, "chezy", nu=1e-6)

    def test_normal_depth_solves_flow(self):
        sec = Section("semicircle", 0.04)
        nu = get_fluid("water").nu
        for item in random_scales(4, 5, n_min=5, n_max=5):
            q = item["value"]
            h = normal_depth(sec, q, 0.05, "manning", nu=nu, h_max=0.04)
            self.assertAlmostEqual(channel_flow(sec, h, 0.05, "manning", nu=nu) / q, 1.0, places=9)


if __name__ == "__main__":
    unittest.main()
