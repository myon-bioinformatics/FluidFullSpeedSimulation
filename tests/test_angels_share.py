"""👼 天使の取り分（E25）: 較正・閉じた式・RH*・40 % を切る年数・木の吸い込み・描画。"""

import contextlib
import io
import math
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import from_q, random_scales
from fluids import abv, mix_volume, moles_for_abv


class Calibration(unittest.TestCase):
    def setUp(self):
        self.rates = ffs.angel_rates()

    def test_scotland_targets_are_hit_exactly(self):
        e1, w1 = ffs.cask_state(1.0, self.rates, ffs.SCOTLAND["rh"])
        self.assertAlmostEqual(mix_volume(e1, w1), ffs.CASK["volume"] * 0.98, places=12)
        self.assertAlmostEqual(abv(e1, w1), ffs.CASK["fill_abv"] - 0.5, places=9)
        self.assertGreater(self.rates["k_e"], self.rates["k_w"])  # 湿った蔵ではエタノールの方が速く抜ける

    def test_water_rate_and_crossover(self):
        r = self.rates
        self.assertEqual(ffs.water_rate(r["k_w"], r["rh_ref"], r["rh_ref"]), r["k_w"])
        self.assertEqual(ffs.water_rate(r["k_w"], 1.0, r["rh_ref"]), 0.0)
        self.assertEqual(ffs.water_rate(r["k_w"], 1.2, r["rh_ref"]), 0.0)  # 湿度 100 % 超は 0 で止める
        rh = ffs.crossover_rh(r)
        self.assertAlmostEqual(ffs.water_rate(r["k_w"], rh, r["rh_ref"]), r["k_e"], places=15)

    def test_direction_follows_rh_star_for_random_humidity(self):
        rh_star = ffs.crossover_rh(self.rates)
        for item in random_scales(251, 8, n_min=1, n_max=1, digits=3):
            rh = item["value"]
            for f in (0.5, 2.0):  # 温度の倍率は向きを変えない
                up = abv(*ffs.cask_state(5.0, self.rates, rh, f)) > ffs.CASK["fill_abv"]
                with self.subTest(rh=item["label"], f=f):
                    self.assertEqual(up, rh < rh_star)

    def test_fit_temperature_and_years(self):
        f = ffs.fit_temperature(self.rates, 0.45, 0.04)
        e1, w1 = ffs.cask_state(1.0, self.rates, 0.45, f)
        self.assertAlmostEqual(1.0 - mix_volume(e1, w1) / ffs.CASK["volume"], 0.04, places=9)
        t40 = ffs.years_to_abv(40.0, self.rates, ffs.SCOTLAND["rh"])
        self.assertAlmostEqual(abv(*ffs.cask_state(t40, self.rates, ffs.SCOTLAND["rh"])), 40.0, places=6)
        self.assertIsNone(ffs.years_to_abv(40.0, self.rates, 0.3))  # 乾いた蔵では度数が上がるので切らない


class Simulation(unittest.TestCase):
    def test_simulation_matches_closed_form_and_conserves(self):
        rates = ffs.angel_rates()
        net = ffs.angel_cask(rates, {"a": {"rh": 0.7, "f": 1.7}}, steps_per_year=52)
        res = net.run(12 * ffs.YEAR)
        self.assertEqual(res["parts_error_q"], {"ethanol": 0, "water": 0})
        e, w = ffs.cask_state(12.0, rates, 0.7, 1.7)
        p = res["parts_q"]["a"]
        self.assertAlmostEqual(from_q(p["ethanol"]) / e, 1.0, delta=1e-10)
        self.assertAlmostEqual(from_q(p["water"]) / w, 1.0, delta=1e-10)
        air = res["parts_q"]["air"]
        ne0, nw0 = moles_for_abv(ffs.CASK["fill_abv"], ffs.CASK["volume"])
        self.assertAlmostEqual(from_q(air["ethanol"]) + from_q(p["ethanol"]), ne0, places=9)

    def test_soak_first_year_then_settles(self):
        soak = ffs.angel_soak(ffs.angel_rates())
        self.assertAlmostEqual(soak["loss_1y"], ffs.SOAK["first_year_loss"], places=6)
        self.assertLess(abs(soak["loss_2y"] - ffs.SCOTLAND["loss"]), 0.002)
        self.assertEqual(soak["parts_error_q"], {"ethanol": 0, "water": 0})
        self.assertGreater(soak["capacity_liters"], 0.0)

    def test_experiment_passes(self):
        r = ffs.run_experiment("👼")
        self.assertTrue(r["ok"], ffs.render(r))
        s = r["summary"]
        self.assertGreater(s["kentucky-top"]["abv_end"], ffs.CASK["fill_abv"])
        self.assertLess(s["scotland"]["abv_end"], ffs.CASK["fill_abv"])
        self.assertTrue(0 <= s["years_to_40_sim"] - s["years_to_40"] < 1.0)
        self.assertTrue(0.0 < s["rh_star"] < ffs.SCOTLAND["rh"])

    def test_draw_scene(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(["draw", "👼", "--frames", "2", "--rows", "4"])
        self.assertEqual(code, 0)
        lines = out.getvalue().splitlines()
        self.assertEqual([line for line in lines if line.startswith("⏱️")], ["⏱️ 1.00", "⏱️ 30.00"])
        self.assertTrue(all(len(line) == 4 for line in lines if not line.startswith("⏱️")))


if __name__ == "__main__":
    unittest.main()
