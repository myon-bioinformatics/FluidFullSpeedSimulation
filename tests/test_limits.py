"""🚧 E23 inertia-limits: U 字管の閉じた式、エンジンの単調な緩和、Ritter の解、ダム決壊の外れ方、Joukowsky。"""

import contextlib
import io
import math
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import random_scales
from fluids import G, get_fluid


class UTubeTheory(unittest.TestCase):
    def test_uniform_tube_period(self):
        nu = get_fluid("water").nu
        for item in random_scales(31, 4, n_min=3, n_max=3):  # 半径 1〜10 mm
            r = item["value"]
            th = ffs.utube_theory(math.pi * r * r, r, 0.3, 0.3, nu)
            self.assertAlmostEqual(th["L_eff"], 2 * 0.3 + 0.3, places=12)  # 一様な管: 水の全長
            self.assertAlmostEqual(th["T0"], 2 * math.pi * math.sqrt(0.9 / (2 * G)), places=12)
            self.assertAlmostEqual(th["tau_engine"], th["b"] / (2 * G), places=12)  # エンジンの τ = b/2g

    def test_zeta_falls_with_radius_and_crossover(self):
        nu = get_fluid("water").nu
        r_star = ffs.utube_crossover_radius(0.3, 0.3, nu)
        self.assertAlmostEqual(ffs.utube_theory(math.pi * r_star ** 2, r_star, 0.3, 0.3, nu)["zeta"], 1.0, places=6)
        zetas = [ffs.utube_theory(math.pi * r * r, r, 0.3, 0.3, nu)["zeta"] for r in (0.3e-3, 0.5e-3, 1e-3, 3e-3)]
        self.assertEqual(zetas, sorted(zetas, reverse=True))
        self.assertTrue(ffs.utube_theory(math.pi * 1e-6, 1e-3, 0.3, 0.3, nu)["oscillates"])
        self.assertFalse(ffs.utube_theory(0.01, 0.75e-3, 0.5, 0.1, nu)["oscillates"])

    def test_overdamped_roots(self):
        th = ffs.utube_theory(0.01, 0.75e-3, 0.5, 0.1, get_fluid("water").nu)
        self.assertGreater(th["zeta"], 1.0)
        self.assertGreater(th["tau_slow"], th["tau_fast"])
        self.assertAlmostEqual(th["tau_slow"] * th["tau_fast"], th["L_eff"] / (2 * G), places=9)  # 根の積 = L_eff/2g
        self.assertLess(abs(th["tau_slow"] / th["tau_engine"] - 1.0), 1.0 / th["zeta"] ** 2)


class UTubeEngine(unittest.TestCase):
    def test_monotone_relaxation_for_random_heads(self):
        for item in random_scales(32, 3, n_min=2, n_max=3):
            dh = item["value"]
            with self.subTest(dh=item["label"]):
                net, th = ffs.utube(0.01, 0.75e-3, 0.5, 0.2, dh, 20.0)
                res = net.run(10 * th["tau_engine"], trace=["left", "right"])
                diff = [a - b for a, b in zip(res["trace"]["left"], res["trace"]["right"])]
                self.assertTrue(all(b <= a for a, b in zip(diff, diff[1:])))
                self.assertGreaterEqual(min(diff), 0.0)
                self.assertLess(diff[-1], 1e-4 * dh)
                self.assertEqual(res["conservation_error_q"], 0)
                # 厳密離散化: 1 ステップ後の差 = Δh·e^(−dt/τ)
                self.assertAlmostEqual(diff[0] / dh, math.exp(-20.0 / th["tau_engine"]), places=9)

    def test_limb_must_not_run_dry(self):
        with self.assertRaises(ValueError):
            ffs.utube(0.01, 0.75e-3, 0.5, 0.2, 0.5, 20.0)

    def test_mirror_symmetry(self):
        a, _ = ffs.utube(0.01, 0.75e-3, 0.5, 0.2, 0.05, 20.0)
        b, _ = ffs.utube(0.01, 0.75e-3, 0.5, 0.2, -0.05, 20.0)
        ra, rb = a.run(2000.0), b.run(2000.0)
        self.assertEqual(ra["volumes_q"]["left"], rb["volumes_q"]["right"])
        self.assertEqual(ra["volumes_q"]["right"], rb["volumes_q"]["left"])


class Ritter(unittest.TestCase):
    def test_profile_landmarks(self):
        h0, t = 0.1, 1.0
        c0 = math.sqrt(G * h0)
        self.assertAlmostEqual(ffs.ritter_depth(0.0, t, h0), 4 * h0 / 9, places=12)
        self.assertAlmostEqual(ffs.ritter_depth(-c0 * t, t, h0), h0, places=12)
        self.assertEqual(ffs.ritter_depth(2 * c0 * t, t, h0), 0.0)
        self.assertEqual(ffs.ritter_depth(-5.0, t, h0), h0)
        self.assertEqual(ffs.ritter_depth(1.0, 0.0, h0), 0.0)
        self.assertAlmostEqual(ffs.ritter_isoline_speed(0.0, h0), 2 * c0, places=12)
        self.assertAlmostEqual(ffs.ritter_isoline_speed(4 / 9, h0), 0.0, places=12)  # 4/9·h0 の点はダム地点に留まる

    def test_front_fit(self):
        self.assertAlmostEqual(ffs._front_fit([0.0, 1.0, 2.0], [0.0, 0.5, 1.0]), 2.0, places=12)
        self.assertIsNone(ffs._front_fit([0.0, 1.0], [0.0, None]))
        self.assertIsNone(ffs._front_fit([0.0, 1.0], [1.0, 1.0]))


class DamBreak(unittest.TestCase):
    def test_engine_never_reaches_the_tip_speed(self):
        for item in random_scales(33, 2, n_min=2, n_max=2):  # h0 = 1〜10 cm
            h0 = item["value"]
            with self.subTest(h0=item["label"]):
                rep = ffs.dam_break_report(h0, 1e-2, t_end=12.0)
                self.assertIsNotNone(rep["tip_ratio"])
                self.assertLess(rep["tip_ratio"], 1.0)
                self.assertTrue(rep["drawdown_reversed"])
                self.assertEqual(rep["res"]["conservation_error_q"], 0)
                self.assertLessEqual(rep["courant"], 1.0)

    def test_initial_column_is_exact(self):
        net, info = ffs.dam_break(0.1, slope=1e-3)
        total = sum(net.nodes[c]["v"] for c in net.channels["bed"]["cells"])
        self.assertEqual(total, info["n_res"] * net.nodes["bed:0"]["v"])
        self.assertAlmostEqual(net.head("bed:0"), 0.1, places=12)
        self.assertEqual(net.head(f"bed:{info['n_res']}"), 0.0)


class Joukowsky(unittest.TestCase):
    def test_surge(self):
        jk = ffs.joukowsky(1.0, 1000.0, c=1400.0)
        self.assertAlmostEqual(jk["delta_p"], 1.4e6)
        self.assertAlmostEqual(jk["head"], 1400.0 / G, places=9)
        self.assertAlmostEqual(jk["delta_p_bar"], 14.0)


class Cli(unittest.TestCase):
    def call(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(list(argv))
        return code, out.getvalue()

    def test_run_and_draw_by_icon(self):
        code, text = self.call("run", "🚧")
        self.assertEqual(code, 0)
        self.assertTrue(text.startswith("✅ E23"))
        code, text = self.call("draw", "🚧", "--frames", "2", "--rows", "4")
        self.assertEqual(code, 0)
        self.assertIn("🟦", text)
        self.assertEqual(sum(line.startswith("⏱️") for line in text.splitlines()), 2)


if __name__ == "__main__":
    unittest.main()
