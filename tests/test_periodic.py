"""🎍🍷🫗 溜まる → 閾値 → 一気に: tip / siphon / mariotte の法則と E19〜E21。"""

import contextlib
import io
import math
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import from_q, qmul, random_scales, to_q
from fluids import (
    G, mariotte_constant_time, mariotte_flow, siphon_drain_time, siphon_drain_time_with_inflow, siphon_flow,
    tip_period, torricelli_empty_time,
)
from network import Network, build, simulate


class ClosedForms(unittest.TestCase):
    def test_tip_period(self):
        self.assertAlmostEqual(tip_period(1e-4, 2e-6), 50.0)
        self.assertAlmostEqual(tip_period(1e-4, 2e-6, 1.5), 51.5)
        with self.assertRaises(ValueError):
            tip_period(0.0, 1e-6)

    def test_siphon_drain_is_torricelli_with_shifted_head(self):
        # h_drop = 0 で吸い口が底なら、小孔の空になる時間と同じ式
        self.assertAlmostEqual(siphon_drain_time(3e-3, 0.6, 2e-5, 0.08, 0.0, 0.0),
                               torricelli_empty_time(0.08, 3e-3, 2e-5, 0.6), places=12)
        self.assertEqual(siphon_flow(0.6, 2e-5, 0.0), 0.0)
        with self.assertRaises(ValueError):
            siphon_drain_time(3e-3, 0.6, 2e-5, 0.0, 0.08, 0.0)

    def test_siphon_reference_converges_and_rejects_strong_inflow(self):
        a = siphon_drain_time_with_inflow(3e-3, 0.6, 2e-5, 0.08, 0.0, 0.05, 5e-6, steps=4000)
        b = siphon_drain_time_with_inflow(3e-3, 0.6, 2e-5, 0.08, 0.0, 0.05, 5e-6, steps=40000)
        self.assertAlmostEqual(a / b, 1.0, places=7)
        self.assertAlmostEqual(siphon_drain_time_with_inflow(3e-3, 0.6, 2e-5, 0.08, 0.0, 0.05, 0.0),
                               siphon_drain_time(3e-3, 0.6, 2e-5, 0.08, 0.0, 0.05), places=8)
        with self.assertRaises(ValueError):
            siphon_drain_time_with_inflow(3e-3, 0.6, 2e-5, 0.08, 0.0, 0.05, 2e-5)

    def test_mariotte(self):
        self.assertAlmostEqual(mariotte_flow(0.61, 1e-5, 0.1), 0.61 * 1e-5 * math.sqrt(2 * G * 0.1))
        self.assertEqual(mariotte_constant_time(0.05, 0.1, 5e-3, 0.61, 1e-5), 0.0)
        for item in random_scales(211, 4, n_min=2, n_max=2):  # 1〜10 cm（h0 = 30 cm より下）
            h = item["value"]
            self.assertAlmostEqual(mariotte_constant_time(0.3, h, 5e-3, 0.61, 1e-5) * mariotte_flow(0.61, 1e-5, h),
                                   (0.3 - h) * 5e-3, places=12)


class TipLaw(unittest.TestCase):
    def test_dump_is_exact_and_periodic_for_random_rates(self):
        for item in random_scales(191, 3, n_min=5, n_max=7):
            q = item["value"]
            dt = 1e-4 / q / 500.0
            with self.subTest(q=item["label"]):
                net = ffs.shishi_layout(q, dt=dt, t_dump=3 * dt, dump_rate=q * 400.0)
                amount = qmul(to_q(q), net.dt_q)
                res = net.run(3.5 * tip_period(1e-4, q, 3 * dt))
                tips = res["links"][1]["tips"]
                closed = [x for x in tips if x["dumped_q"] is not None]
                self.assertEqual(len(tips), 3)
                self.assertGreaterEqual(len(closed), 2)
                for x in closed:
                    self.assertEqual(x["dumped_q"], x["volume_q"] + 2 * amount)
                    self.assertEqual(x["left_q"], amount)
                gaps = [b["t"] - a["t"] for a, b in zip(tips, tips[1:])]
                self.assertLessEqual(max(gaps) - min(gaps), 1e-9)
                self.assertEqual(res["conservation_error_q"], 0)

    def test_dump_rate_too_small_leaves_water_behind(self):
        """落とす流量が足りなければ、筒は空にならずに戻る（残りが 1 ステップ分より多い）。"""
        net = ffs.shishi_layout(2e-6, dt=0.05, t_dump=1.0, dump_rate=5e-5)  # 1 s × 5e-5 = 5e-5 < V_tip
        res = net.run(120.0)
        tips = res["links"][1]["tips"]
        self.assertGreaterEqual(len(tips), 2)
        self.assertGreater(tips[0]["left_q"], qmul(to_q(2e-6), net.dt_q))
        self.assertEqual(res["conservation_error_q"], 0)

    def test_arguments(self):
        net = Network(0.1).tank("a", 1.0, volume=1.0).tank("b", 1.0)
        with self.assertRaises(ValueError):
            net.tip("a", "b", v_tip=0.0)
        with self.assertRaises(ValueError):
            net.tip("a", "b", v_tip=1e-3, dump_rate=0.0)


class SiphonLaw(unittest.TestCase):
    def test_spec_builds_and_matches_hand_built(self):
        spec = {"dt": 0.02, "tanks": [{"name": "bucket", "area": 0.05, "volume": 0.5}, {"name": "cup", "area": 3e-3},
                                      {"name": "bowl", "area": 0.05}],
                "links": [{"law": "pump", "src": "bucket", "dst": "cup", "rate": 5e-6},
                          {"law": "siphon", "src": "cup", "dst": "bowl", "h_crest": 0.08, "h_drop": 0.05,
                           "area": 2e-5, "cd": 0.6}],
                "t_end": 150.0}
        hand = ffs.cup_layout(5e-6, dt=0.02).run(150.0)
        self.assertEqual(build(spec).run(150.0)["volumes_q"], hand["volumes_q"])
        res = simulate(spec)
        self.assertEqual(len(res["links"][1]["primes"]), 2)
        self.assertEqual(res["conservation_error_q"], 0)

    def test_intake_above_bottom_keeps_that_much(self):
        net = Network(0.005).tank("cup", 3e-3, head=0.08).tank("bowl", 0.05)
        net.siphon("cup", "bowl", h_crest=0.08, h_drop=0.05, area=2e-5, h_intake=0.02)
        res = net.run(40.0)
        self.assertEqual(res["volumes_q"]["cup"], qmul(to_q(0.02), to_q(3e-3)))
        self.assertEqual(res["links"][0]["breaks"][0]["left_q"], res["volumes_q"]["cup"])
        self.assertFalse(res["links"][0]["primed"])

    def test_never_primes_below_crest(self):
        net = Network(0.02).tank("cup", 3e-3, head=0.079).tank("bowl", 0.05)
        net.siphon("cup", "bowl", h_crest=0.08, h_drop=0.05, area=2e-5)
        res = net.run(10.0)
        self.assertEqual(res["links"][0]["primes"], [])
        self.assertEqual(res["links"][0]["volume_moved_q"], 0)

    def test_arguments(self):
        net = Network(0.1).tank("a", 1.0, volume=1.0).tank("b", 1.0)
        for bad in ({"h_crest": 0.0, "area": 1e-5}, {"h_crest": 0.1, "h_intake": 0.2, "area": 1e-5},
                    {"h_crest": 0.1, "area": 0.0}, {"h_crest": 0.1, "h_drop": -0.1, "area": 1e-5}):
            with self.assertRaises(ValueError):
                net.siphon("a", "b", **bad)


class MariotteLaw(unittest.TestCase):
    def test_constant_then_torricelli_for_random_h_air(self):
        for item in random_scales(212, 3, n_min=2, n_max=2):  # 1〜10 cm
            h_air = item["value"]
            with self.subTest(h_air=item["label"]):
                net = ffs.bottle_layout(h_air, dt=0.1)
                res = net.run(1.2 * mariotte_constant_time(0.3, h_air, 5e-3, 0.61, 1e-5))
                ln = res["links"][0]
                self.assertAlmostEqual(ln["constant_flow"] / mariotte_flow(0.61, 1e-5, h_air), 1.0, places=6)
                self.assertIsNotNone(ln["falling_from"])
                self.assertLess(net.head("bottle"), h_air)
                self.assertEqual(res["conservation_error_q"], 0)

    def test_level_below_tube_is_plain_orifice(self):
        """空気管の下端より水面が低ければ、最初からふつうの小孔とビット単位で同じ。"""
        a = Network(0.05).tank("t", 5e-3, head=0.08).tank("b", 1.0)
        a.mariotte("t", "b", area=1e-5, h_air=0.1)
        b = Network(0.05).tank("t", 5e-3, head=0.08).tank("b", 1.0)
        b.orifice("t", "b", area=1e-5)
        self.assertEqual(a.run(100.0)["volumes_q"], b.run(100.0)["volumes_q"])

    def test_arguments(self):
        net = Network(0.1).tank("a", 1.0, volume=1.0).tank("b", 1.0)
        with self.assertRaises(ValueError):
            net.mariotte("a", "b", area=1e-5, h_air=0.0)
        with self.assertRaises(ValueError):
            net.mariotte("a", "b", area=0.0, h_air=0.1)


class Cli(unittest.TestCase):
    def call(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(list(argv))
        return code, out.getvalue()

    def test_run_and_draw_by_icon(self):
        for icon, eid in (("🎍", "E19"), ("🍷", "E20"), ("🫗", "E21")):
            self.assertEqual(ffs.resolve(icon), eid)
        for icon in ("🎍", "🍷", "🫗"):
            code, text = self.call("draw", icon, "--frames", "3", "--rows", "4")
            self.assertEqual(code, 0, icon)
            rows = [line for line in text.splitlines() if not line.startswith("⏱️")]
            self.assertEqual(len(rows), 3 * 4)
            self.assertTrue(all(len(line) == 2 for line in rows))
        code, text = self.call("draw", "🫗", "🟩", "--frames", "2")
        self.assertEqual(code, 0)
        self.assertIn("🟩", text)


if __name__ == "__main__":
    unittest.main()
