"""🍷 エタノール・水の混合物: 密度表・度数の換算・気液平衡・Rayleigh 式と、成分モードの水槽（boil・evaporate・カット）。"""

import math
import unittest

from fixedq import from_q, random_scales
from fluids import (
    ABV_CHECK_ROWS, AZEOTROPE_X, M_ETHANOL, M_WATER, abv, bubble_t, density_20, mix_volume, mole_fraction,
    moles_for_abv, rayleigh_log, vle_y,
)
from network import Network, simulate

VOL = lambda p: mix_volume(p["ethanol"], p["water"])


class MixtureData(unittest.TestCase):
    def test_abv_matches_table_rows(self):
        for wt, vol in ABV_CHECK_ROWS:
            w = wt / 100.0
            got = abv(w / M_ETHANOL, (1 - w) / M_WATER)
            with self.subTest(mass_percent=wt):
                self.assertAlmostEqual(got, vol, delta=0.35)  # 71.9 % の行だけ 0.32 ずれる（表の行どうしの食い違い）

    def test_mixing_shrinks(self):
        ne, nw = moles_for_abv(50.0, 1e-3)
        separate = ne * M_ETHANOL / density_20(1.0) + nw * M_WATER / density_20(0.0)
        self.assertLess(mix_volume(ne, nw), separate)  # 混ぜると縮む
        self.assertGreater((separate - mix_volume(ne, nw)) / separate, 0.03)

    def test_moles_for_abv_round_trip(self):
        for item in random_scales(240, 6, n_min=-1, n_max=-1, digits=2):
            a = min(item["value"], 96.0)
            ne, nw = moles_for_abv(a, 0.02)
            self.assertAlmostEqual(abv(ne, nw), a, places=9)
            self.assertAlmostEqual(mix_volume(ne, nw), 0.02, places=12)
        with self.assertRaises(ValueError):
            moles_for_abv(101.0, 1.0)

    def test_vle_and_azeotrope(self):
        self.assertEqual(vle_y(0.0), 0.0)
        self.assertAlmostEqual(vle_y(AZEOTROPE_X), AZEOTROPE_X)
        self.assertAlmostEqual(bubble_t(0.0), 100.0)
        for x in (0.01, 0.05, 0.2, 0.5, 0.8):
            self.assertGreater(vle_y(x), x)  # 共沸点より下では蒸気の方が濃い

    def test_rayleigh_log_matches_numeric_integral(self):
        x0, x1 = 0.12, 0.01
        n = 200000
        h = (x0 - x1) / n
        num = sum(h / (vle_y(x1 + (i + 0.5) * h) - (x1 + (i + 0.5) * h)) for i in range(n))
        self.assertAlmostEqual(rayleigh_log(x1, x0), num, delta=1e-6)
        with self.assertRaises(ValueError):
            rayleigh_log(0.5, 0.1)


class SpeciesNetwork(unittest.TestCase):
    def still(self, dt=1.0, power=5000.0, cuts=None):
        ne, nw = moles_for_abv(12.0, 0.025)
        net = Network(dt, species=("ethanol", "water"), volume_of=VOL)
        net.tank("pot", 0.1, parts={"ethanol": ne, "water": nw}).tank("a", 0.05).tank("b", 0.05)
        ln = net.boil("pot", None if cuts else "a", power=power, cuts=cuts)
        return net, ln, ne, nw

    def test_species_conserved_exactly(self):
        net, _ln, _ne, _nw = self.still()
        net.tank("barrel", 0.2, parts={"ethanol": 10.0, "water": 30.0}).tank("air", 1.0)
        net.evaporate("barrel", "air", rates={"ethanol": 1e-4, "water": 3e-5})
        net.pump("a", "b", rate=1e-3)
        res = net.run(2000.0)
        self.assertEqual(res["parts_error_q"], {"ethanol": 0, "water": 0})
        self.assertEqual(res["conservation_error_q"], 0)
        for nd in res["parts_q"].values():
            self.assertTrue(all(q >= 0 for q in nd.values()))

    def test_rayleigh_first_order(self):
        errs = []
        for dt in (4.0, 1.0):
            net, _ln, ne, nw = self.still(dt)
            res = net.run(3000.0)
            pe, pw = (from_q(res["parts_q"]["pot"][k]) for k in ("ethanol", "water"))
            w0, x0, x = ne + nw, ne / (ne + nw), pe / (pe + pw)
            errs.append(abs((pe + pw) - w0 * math.exp(-rayleigh_log(x, x0))) / (pe + pw))
        self.assertLess(errs[1], errs[0] / 3)  # 1 次: dt を 1/4 にするとほぼ 1/4

    def test_evaporate_is_exact_exponential(self):
        net = Network(10.0, species=("ethanol", "water"))
        net.tank("barrel", 0.2, parts={"ethanol": 10.0, "water": 30.0}).tank("air", 1.0)
        net.evaporate("barrel", "air", rates={"ethanol": 2e-5})
        res = net.run(10000.0)
        self.assertAlmostEqual(from_q(res["parts_q"]["barrel"]["ethanol"]), 10.0 * math.exp(-0.2), delta=1e-9)
        self.assertEqual(res["parts_q"]["barrel"]["water"], net.nodes["barrel"]["parts"]["water"])
        self.assertAlmostEqual(from_q(res["parts_q"]["barrel"]["water"]), 30.0, delta=1e-12)

    def test_proportional_split_keeps_composition(self):
        net = Network(1.0, species=("ethanol", "water"), volume_of=VOL)
        net.tank("a", 0.1, parts={"ethanol": 3.0, "water": 7.0}).tank("b", 0.1)
        net.pump("a", "b", rate=0.5)
        net.run(4.0)
        b = net.nodes["b"]["parts"]
        self.assertAlmostEqual(b["ethanol"] / (b["ethanol"] + b["water"]), 0.3, places=12)

    def test_cuts_switch_receivers(self):
        net, ln, _ne, _nw = self.still(cuts=[{"to": "a", "abv_below": 5.0}, {"to": "b"}])
        res = net.run(4500.0)
        cut = res["links"][0]["cuts"]
        self.assertIn("ended", cut[0])
        self.assertLess(res["links"][0]["last_abv"], 5.0)
        a, b = res["parts_q"]["a"], res["parts_q"]["b"]
        self.assertGreater(abv(from_q(a["ethanol"]), from_q(a["water"])), abv(from_q(b["ethanol"]), from_q(b["water"])))

    def test_errors(self):
        with self.assertRaises(ValueError):
            Network(1.0).boil("a", "b", power=1.0)
        net = Network(1.0, species=("ethanol", "water"))
        with self.assertRaises(ValueError):
            net.tank("a", 0.1, head=0.2)
        with self.assertRaises(ValueError):
            net.tank("a", 0.1, parts={"sugar": 1.0})
        net.tank("a", 0.1, parts={"water": 1.0}).tank("b", 0.1)
        with self.assertRaises(ValueError):
            net.pump("a", "b", rate=1.0, delay=2.0)
        with self.assertRaises(ValueError):
            net.evaporate("a", "b", rates={"sugar": 1.0})
        with self.assertRaises(ValueError):
            net.boil("a", power=1.0, cuts=[{"to": "zzz"}])
        with self.assertRaises(ValueError):
            Network(1.0).tank("a", 1.0, parts={"water": 1.0})

    def test_spec(self):
        ne, nw = moles_for_abv(12.0, 0.01)
        res = simulate({"dt": 2.0, "species": ["ethanol", "water"],
                        "tanks": [{"name": "pot", "area": 0.05, "parts": {"ethanol": ne, "water": nw}},
                                  {"name": "jug", "area": 0.05}],
                        "links": [{"law": "boil", "src": "pot", "dst": "jug", "power": 2000.0}], "t_end": 600.0})
        self.assertEqual(res["parts_error_q"], {"ethanol": 0, "water": 0})
        self.assertGreater(res["final"]["h:jug"], 0.0)


if __name__ == "__main__":
    unittest.main()
