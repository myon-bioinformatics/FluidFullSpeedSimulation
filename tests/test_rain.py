"""🌧️ 斜面に降る雨: 運動波の閉じた式、Network.rain、spec の rains、区間数での収束、E18 の場面。"""

import contextlib
import io
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import random_scales
from fluids import (
    Section, channel_flow, get_fluid, overland_equilibrium, overland_power_law, overland_recession, overland_rise,
    overland_t_e_power_law, plan_width, wave_celerity,
)
from network import Network, build

SEC = Section("rect", 1.0)
NU = get_fluid("water").nu
KW = {"nu": NU, "h_max": 0.01}


class ClosedForms(unittest.TestCase):
    def test_equilibrium_and_power_law_agree_on_a_wide_plane(self):
        for law, i in (("manning", 1e-4), ("laminar", 3e-6)):
            eq = overland_equilibrium(SEC, i, 30.0, 0.05, law, **KW)
            self.assertAlmostEqual(eq["Q_e"], i * 30.0)
            self.assertAlmostEqual(channel_flow(SEC, eq["h_e"], 0.05, law, nu=NU) / eq["Q_e"], 1.0, places=9)
            self.assertAlmostEqual(eq["t_e"], eq["h_e"] / i)
            alpha, m = overland_power_law(0.05, law, nu=NU)
            t_wide = overland_t_e_power_law(i, 30.0, alpha, m)
            self.assertLess(abs(t_wide - eq["t_e"]) / eq["t_e"], 2 * eq["h_e"] / plan_width(SEC))  # R ≈ h の分
            # 立ち上がりは t_e で定常に届き、減衰は定常から単調に減る
            self.assertAlmostEqual(overland_rise(SEC, i, eq["t_e"], 0.05, law, nu=NU) / eq["Q_e"], 1.0, places=9)
            rec = [overland_recession(SEC, i, 30.0, 0.05, law, t, **KW) for t in (0.0, 0.5 * eq["t_e"], eq["t_e"], 3 * eq["t_e"])]
            self.assertAlmostEqual(rec[0] / eq["Q_e"], 1.0, places=6)
            self.assertEqual(rec, sorted(rec, reverse=True))

    def test_celerity_is_m_times_velocity_for_random_depths(self):
        for law in ("manning", "laminar"):
            _alpha, m = overland_power_law(0.05, law, nu=NU)
            for item in random_scales(181, 4, n_min=3, n_max=3):
                h = item["value"]
                v = channel_flow(SEC, h, 0.05, law, nu=NU) / SEC.area(h)
                self.assertAlmostEqual(wave_celerity(SEC, h, 0.05, law, nu=NU) / (m * v), 1.0, delta=3 * h)
        self.assertEqual(wave_celerity(SEC, 0.0, 0.05, "manning", nu=NU), 0.0)
        self.assertEqual(plan_width(Section("semicircle", 0.04)), 0.08)


class RainHelper(unittest.TestCase):
    def test_rain_adds_one_pump_per_cell_with_plan_area_rate(self):
        net = Network(0.1).tank("sky", 1.0, volume=10.0).tank("basin", 1.0)
        net.channel("plane", "basin", length=30.0, cells=6, slope=0.05, section=SEC, max_depth=0.01)
        pumps = net.rain("plane", 2e-5, "sky", until=10.0)
        self.assertEqual(len(pumps), 6)
        self.assertTrue(all(p["kind"] == "pump" and p["rain"] and p["channel"] == "plane" for p in pumps))
        self.assertAlmostEqual(pumps[0]["rate"], 2e-5 * 1.0 * 5.0)
        res = net.run(20.0)
        self.assertEqual(res["channels"]["plane"]["rained_q"], sum(p["flow_sum"] for p in pumps))
        self.assertEqual(res["conservation_error_q"], 0)
        with self.assertRaises(ValueError):
            net.rain("nowhere", 1e-5, "sky")
        with self.assertRaises(ValueError):
            net.rain("plane", -1e-5, "sky")

    def test_spec_rains_matches_hand_built(self):
        spec = {"dt": 0.1, "tanks": [{"name": "sky", "area": 1, "volume": 10}, {"name": "basin", "area": 1}],
                "channels": [{"prefix": "plane", "dst": "basin", "length": 30.0, "cells": 6, "slope": 0.05,
                              "section": {"shape": "rect", "size": 1.0}, "max_depth": 0.01}],
                "rains": [{"prefix": "plane", "intensity": 2e-5, "src": "sky", "until": 10.0}]}
        net = Network(0.1).tank("sky", 1.0, volume=10.0).tank("basin", 1.0)
        net.channel("plane", "basin", length=30.0, cells=6, slope=0.05, section=SEC, max_depth=0.01)
        net.rain("plane", 2e-5, "sky", until=10.0)
        self.assertEqual(build(spec).run(20.0)["volumes_q"], net.run(20.0)["volumes_q"])


class Accuracy(unittest.TestCase):
    def test_steady_state_for_random_intensities(self):
        for item in random_scales(182, 3, n_min=5, n_max=5):
            a = ffs.rain_accuracy("manning", item["value"], 10)
            with self.subTest(i=item["label"]):
                self.assertLessEqual(ffs.rel_err(a["outflow_steady"], a["Q_e"]), 1e-3)
                self.assertLessEqual(ffs.rel_err(a["depth_steady"], a["h_e"]), 1e-3)
                self.assertEqual(a["rained_q"], a["drained_q"])

    def test_t95_error_shrinks_with_cells_and_manning_is_fluid_free(self):
        coarse = ffs.rain_accuracy("manning", ffs.RAIN["manning"], 10)
        fine = ffs.rain_accuracy("manning", ffs.RAIN["manning"], 40)
        self.assertLess(fine["t95_rel_err"], 0.7 * coarse["t95_rel_err"])
        self.assertGreater(coarse["t95"], coarse["t95_exact"])  # 区間に分けた分、遅れる側にずれる
        ethanol = ffs.rain_accuracy("manning", ffs.RAIN["manning"], 10, fluid="ethanol")
        self.assertEqual(ethanol["volumes_q"], coarse["volumes_q"])

    def test_too_much_rain_is_refused(self):
        with self.assertRaises(ValueError):
            ffs.rain_plane(5e-4, 10)


class Cli(unittest.TestCase):
    def test_draw_rain_scene(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(["draw", "🌧️", "--frames", "2", "--rows", "4"])
        self.assertEqual(code, 0)
        self.assertEqual(ffs.SCENE_ICONS["🌧"], "rain-runoff")
        self.assertIn("rain-runoff", ffs.RIG_SCENES)
        self.assertEqual(ffs.resolve("🌧️"), "E18")


if __name__ == "__main__":
    unittest.main()
