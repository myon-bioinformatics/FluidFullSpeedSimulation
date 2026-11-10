"""⚗️ E24 ワインの蒸留（still_network、カット、draw の ⚗️）。"""

import contextlib
import io
import math
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import from_q, random_scales
from fluids import abv, mole_fraction, moles_for_abv, rayleigh_log


class Distillation(unittest.TestCase):
    def test_wine_strength_sweep_follows_rayleigh(self):
        """ワインの度数を 7〜12 % で振っても、釜に残る量と組成の関係は Rayleigh 式どおり。"""
        for wine in (7.0, 9.5, 12.0):
            ne, nw = moles_for_abv(wine, 0.01)
            res = ffs.still_network(ne, nw, dt=0.5).run(1200.0)
            pot = ffs._fraction(res, "pot")
            x = mole_fraction(pot["ethanol"], pot["water"])
            w = pot["ethanol"] + pot["water"]
            w_ray = (ne + nw) * math.exp(-rayleigh_log(x, mole_fraction(ne, nw)))
            with self.subTest(wine=wine):
                self.assertLess(abs(w - w_ray) / w_ray, 5e-4)
                self.assertEqual(res["parts_error_q"], {"ethanol": 0, "water": 0})

    def test_cut_point_moves_brouillis(self):
        """初留をどこで切るか（流れ出る度数）で、初留の度数と量が入れ替わる（遅く切るほど薄くて多い）。"""
        ne, nw = moles_for_abv(12.0, 0.01)
        prev = None
        for item in sorted(random_scales(241, 3, n_min=0, n_max=0, digits=2), key=lambda i: -i["value"]):
            cut = item["value"]
            res = ffs.still_network(ne, nw, receivers=("b", "t"), cuts=[{"to": "b", "abv_below": cut}, {"to": "t"}]).run(3000.0)
            b = ffs._fraction(res, "b")
            if prev is not None:
                self.assertLessEqual(b["abv"], prev["abv"] + 1e-9)
                self.assertGreaterEqual(b["litres"], prev["litres"] - 1e-9)
            prev = b

    def test_no_distillate_above_azeotrope(self):
        ne, nw = moles_for_abv(60.0, 0.005)
        res = ffs.still_network(ne, nw).run(300.0, sample_every=5)
        rows = res["series"]
        self.assertTrue(all(abv(r["parts:receiver"]["ethanol"], r["parts:receiver"]["water"]) < 97.2 for r in rows
                            if r["parts:receiver"]["ethanol"] > 0))

    def test_draw_scene(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(["draw", "⚗️", "--frames", "2", "--rows", "4"])
        self.assertEqual(code, 0)
        rows = [line for line in out.getvalue().splitlines() if not line.startswith("⏱️")]
        self.assertTrue(all(len(line) == 5 for line in rows))  # 釜・前留・ハート・セコンド・後留


if __name__ == "__main__":
    unittest.main()
