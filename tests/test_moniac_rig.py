"""E11 moniac-rig を据え置きの試験台にする。

GOLDEN は試験台の最後の水位（有効 9 桁）。回路・法則・離散化を変えたときに、意図した変化かどうかを
ここで気づくための物差しで、正しさそのものの根拠ではない（正しさは E11 の閉じた式との比較が受け持つ）。
意図して変えたときは、理由をコミットに書いて GOLDEN を更新する。
"""

import unittest

import fluidfullspeedsimulation as ffs
from fixedq import random_scales

GOLDEN = {"h:header": 0.300750625, "h:side": 0.135225125, "h:store": 0.200107155, "h:sump": 0.191188099}


class MoniacRig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = ffs.moniac_rig().run(1200.0)

    def test_golden_levels(self):
        got = ffs.rig_fingerprint(self.res)
        self.assertEqual(set(got), set(GOLDEN))
        for key, want in GOLDEN.items():
            self.assertAlmostEqual(got[key] / want, 1.0, delta=1e-7, msg=key)

    def test_conservation(self):
        self.assertEqual(self.res["conservation_error_q"], 0)

    def test_steady_matches_closed_form_for_random_dt(self):
        for item in random_scales(81, 2, n_min=2, n_max=2):
            with self.subTest(dt=item["label"]):
                r = ffs.run_experiment("E11", dt=item["value"])
                self.assertTrue(r["ok"], ffs.render(r))

    def test_more_lag_stages_still_conserves(self):
        res = ffs.moniac_rig(stages=3).run(300.0)
        self.assertEqual(res["conservation_error_q"], 0)

    def test_valve_arguments(self):
        net = ffs.Network(0.1).tank("a", 1.0, volume=1.0).tank("b", 1.0)
        with self.assertRaises(ValueError):
            net.valve("a", "b", tau=1.0, opening=1.5)
        with self.assertRaises(ValueError):
            net.valve("a", "b", tau=1.0, sensor="a")
        with self.assertRaises(ValueError):
            net.valve("a", "b", tau=1.0, sensor="zzz", cam=[[0, 0], [1, 1]])


if __name__ == "__main__":
    unittest.main()
