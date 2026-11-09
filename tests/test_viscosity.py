"""🍯 粘度のつまみ: get_fluid(mu=, mu_scale=)、Network・spec・draw への通し方、粘い小孔の Cd(Re)。"""

import contextlib
import io
import math
import unittest

import fluidfullspeedsimulation as ffs
import pixelwater as pw
from fixedq import random_scales
from fluids import Section, channel_flow, get_fluid, orifice_cd, viscous_orifice
from network import Network, simulate


class Knob(unittest.TestCase):
    def test_mu_and_scale(self):
        w = get_fluid("water")
        self.assertEqual(get_fluid("water", mu=2e-3).mu, 2e-3)
        self.assertAlmostEqual(get_fluid("water", mu_scale=3.0).mu / w.mu, 3.0, places=12)
        self.assertEqual(get_fluid("water", mu_scale=3.0).rho, w.rho)  # 密度は動かさない
        cold = get_fluid("water", 10.0)
        self.assertAlmostEqual(get_fluid("water", 10.0, mu_scale=2.0).mu, cold.mu * 2.0, places=15)
        self.assertEqual(get_fluid("water", 10.0, mu=1e-3).rho, cold.rho)  # 温度で動くのは密度だけ
        self.assertIs(get_fluid("water"), w)  # 何も渡さなければ今までと同じもの
        for bad in ({"mu": 1e-3, "mu_scale": 2.0}, {"mu": 0.0}, {"mu_scale": -1.0}, {"mu_scale": math.inf}):
            with self.assertRaises(ValueError):
                get_fluid("water", **bad)

    def test_laminar_responds_manning_does_not(self):
        sec = Section("semicircle", 0.04)
        for item in random_scales(91, 5, n_min=0, n_max=1):
            k = item["value"]
            nu = get_fluid("water", mu_scale=k).nu
            h = 1e-3
            lam = channel_flow(sec, h, 0.05, "laminar", nu=nu) / channel_flow(sec, h, 0.05, "laminar", nu=nu / k)
            self.assertAlmostEqual(lam, 1.0 / k, places=12)
            self.assertEqual(channel_flow(sec, h, 0.05, "manning", nu=nu),
                             channel_flow(sec, h, 0.05, "manning", nu=nu / k))

    def test_network_and_spec(self):
        self.assertAlmostEqual(Network(0.1, mu_scale=4.0).fluid.nu / get_fluid("water").nu, 4.0, places=12)
        spec = {"dt": 1.0, "mu_scale": 2.0, "tanks": [{"name": "a", "area": 0.01, "head": 0.1}, {"name": "b", "area": 1}],
                "links": [{"law": "pipe", "src": "a", "dst": "b", "radius": 1e-3, "length": 0.5}], "t_end": 10.0}
        slow = simulate(spec)
        fast = simulate({**spec, "mu_scale": 0.5})
        self.assertGreater(slow["links"][0]["tau"], fast["links"][0]["tau"])
        self.assertAlmostEqual(slow["links"][0]["tau"] / fast["links"][0]["tau"], 4.0, places=9)

    def test_pixel_lateral_rate_follows_viscosity(self):
        self.assertLess(pw.lateral_rate("water", mu_scale=2.0), pw.lateral_rate("water"))
        res = pw.run_pixels("dam-break", "water", mu_scale=3.0, ticks=60)
        self.assertEqual(res["count_error"], 0)


class ViscousOrifice(unittest.TestCase):
    def test_cd_falls_with_viscosity_for_random_scales(self):
        d = 4e-3
        prev = math.inf
        for k in sorted(i["A"] * 10 ** i["N"] for i in random_scales(92, 6, n_min=0, n_max=4)):
            cd, re = viscous_orifice(0.2, d, get_fluid("water").nu * k)
            self.assertLessEqual(cd, 0.61)
            self.assertLessEqual(cd, prev + 1e-15)
            self.assertAlmostEqual(cd, orifice_cd(re), places=9)  # 解いた Cd は Re の式と矛盾しない
            prev = cd
        self.assertEqual(viscous_orifice(0.0, d, 1e-6), (0.0, 0.0))
        self.assertEqual(orifice_cd(0.0), 0.0)

    def test_viscous_off_is_unchanged(self):
        def run(**kw):
            net = Network(0.1, mu_scale=50.0).tank("t", 0.01, head=0.3).tank("sea", 1.0)
            net.orifice("t", "sea", area=1e-5, **kw)
            return net.run(30.0)
        plain = run()
        self.assertEqual(plain["volumes_q"], run(viscous=False)["volumes_q"])
        self.assertGreater(run(viscous=True)["volumes_q"]["t"], plain["volumes_q"]["t"])
        self.assertTrue(plain["links"][0]["viscous_regime"])


class Cli(unittest.TestCase):
    def call(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(list(argv))
        return code, out.getvalue()

    def test_draw_with_viscosity(self):
        self.assertEqual(self.call("draw", "🫙", "--viscosity-scale", "8", "--frames", "1")[0], 0)
        self.assertEqual(self.call("draw", "🥢", "--mu", "5", "--frames", "1")[0], 0)
        self.assertEqual(self.call("draw", "🌊", "🟪", "--viscosity-scale", "2", "--frames", "1")[0], 0)
        for argv in (["draw", "🥢", "--mu", "1", "--viscosity-scale", "2"], ["draw", "🥢", "--mu", "0"],
                     ["draw", "🥢", "--viscosity-scale", "nan"]):
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as cm:
                ffs.main(argv)
            self.assertEqual(cm.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
