"""🏺 断面が高さで変わる水槽（profile）と水時計（E22）。"""

import contextlib
import io
import math
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import from_q, random_scales, to_q
from network import Network, build, simulate


def bowl(h: float) -> float:
    return 0.01 + 0.2 * h  # 上に広がる器


class ProfileTank(unittest.TestCase):
    def test_head_round_trip_table_and_callable(self):
        table = [[0.0, 0.01], [0.1, 0.03], [0.3, 0.07]]
        for item in random_scales(61, 6, n_min=1, n_max=2):
            h = min(item["value"], 0.3)
            net = Network(1.0).tank("t", profile=table, head=h)
            self.assertAlmostEqual(net.head("t"), h, delta=1e-9)
            net2 = Network(1.0).tank("c", 0.3, head=h, profile=bowl)
            self.assertAlmostEqual(net2.head("c"), h, delta=1e-6)  # 台形積分 + 表の補間の分

    def test_volume_matches_integral(self):
        net = Network(1.0).tank("c", 0.3, head=0.3, profile=bowl, samples=4096)
        exact = 0.01 * 0.3 + 0.1 * 0.3 ** 2  # ∫(0.01 + 0.2h) dh
        self.assertAlmostEqual(from_q(net.nodes["c"]["v"]), exact, delta=1e-9)

    def test_above_table_uses_top_area(self):
        net = Network(1.0).tank("t", profile=[[0.0, 0.02], [0.1, 0.02]], head=0.25)
        self.assertAlmostEqual(net.head("t"), 0.25, delta=1e-9)
        net.nodes["t"]["v"] = 0
        self.assertEqual(net.head("t"), 0.0)

    def test_errors(self):
        net = Network(1.0)
        with self.assertRaises(ValueError):
            net.tank("a", profile=[[0.1, 0.02], [0.2, 0.03]])  # 0 から始まらない
        with self.assertRaises(ValueError):
            net.tank("b", profile=[[0.0, 0.0], [0.2, 0.03]])  # 断面 0
        with self.assertRaises(ValueError):
            net.tank("c", profile=bowl)  # callable には area=h_max が要る
        with self.assertRaises(ValueError):
            net.tank("d", 0.3, profile=bowl, samples=1)
        with self.assertRaises(ValueError):
            net.tank("e")
        net.tank("p", profile=[[0.0, 0.01], [0.2, 0.02]], head=0.1).tank("sea", 1.0)
        with self.assertRaises(ValueError):
            net.valve("p", "sea", tau=10.0)  # 弁は断面一定の水槽だけ

    def test_spec_profile_matches_hand_built(self):
        table = [[0.0, 0.01], [0.1, 0.03], [0.3, 0.07]]
        spec = {"dt": 0.5, "tanks": [{"name": "t", "profile": table, "head": 0.25}, {"name": "sea", "area": 1}],
                "links": [{"law": "orifice", "src": "t", "dst": "sea", "area": 1e-5}], "t_end": 300.0}
        net = Network(0.5).tank("t", profile=table, head=0.25).tank("sea", 1.0)
        net.orifice("t", "sea", area=1e-5)
        self.assertEqual(net.run(300.0)["volumes_q"], build(spec).run(300.0)["volumes_q"])
        res = simulate(spec)
        self.assertEqual(res["conservation_error_q"], 0)
        self.assertIn("t", res["profile_samples"])

    def test_pipe_from_profile_tank(self):
        """Q = c·h（前進 Euler）。断面一定の筒なら、既存の厳密離散化と 1 ステップ目の流出が一致する。"""
        net = Network(0.1).tank("p", profile=[[0.0, 0.02], [0.5, 0.02]], head=0.3).tank("sea", 1.0)
        ln = net.pipe("p", "sea", radius=1e-3, length=0.5)
        self.assertIsNone(ln["tau"])
        first = net.step(0)[0]
        self.assertAlmostEqual(from_q(first), ln["conductance"] * 0.3 * 0.1, delta=1e-12)
        res = net.run(50.0)
        self.assertEqual(res["conservation_error_q"], 0)
        self.assertLess(net.head("p"), 0.3)
        self.assertTrue(res["links"][0]["laminar"])


class WaterClock(unittest.TestCase):
    def test_clock_rate_and_linear_descent(self):
        k = ffs.clock_rate(0.3, 0.1, 2e-3, 0.61)
        self.assertAlmostEqual(k, 0.61 * math.pi * 1e-6 * math.sqrt(2 * 9.80665 * 0.3) / (math.pi * 0.01), places=15)
        net = ffs.water_clock("clepsydra", dt=0.25)
        tr = net.run(0.5 * 0.3 / k, trace=["clock"], sample_every=1)["trace"]["clock"]
        self.assertLessEqual(ffs._linear_dev(tr, 0.25, 0.3, k, 0.5 * 0.3 / k), 1e-4)

    def test_random_geometry_stays_linear(self):
        for r_item, d_item in zip(random_scales(62, 3, n_min=1, n_max=1, digits=2), random_scales(63, 3, n_min=3, n_max=3, digits=2)):
            r_top, d = r_item["value"], d_item["value"]
            k = ffs.clock_rate(0.3, r_top, d, 0.61)
            net = ffs.water_clock("clepsydra", r_top=r_top, orifice_d=d)
            tr = net.run(0.9 * 0.3 / k, trace=["clock"], sample_every=1)["trace"]["clock"]
            with self.subTest(r_top=r_item["label"], d=d_item["label"]):
                self.assertLessEqual(ffs._linear_dev(tr, net.dt, 0.3, k, 0.9 * 0.3 / k), 2e-4)

    def test_bad_shape(self):
        with self.assertRaises(ValueError):
            ffs.water_clock("bucket")

    def test_draw_scene(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(["draw", "🏺", "--frames", "3", "--rows", "4"])
        self.assertEqual(code, 0)
        lines = [line for line in out.getvalue().splitlines() if not line.startswith("⏱️")]
        self.assertEqual(len(lines), 3 * 4)
        self.assertTrue(all(len(line) == 3 for line in lines))  # 🏺 🥫 🔺 の 3 列


if __name__ == "__main__":
    unittest.main()
