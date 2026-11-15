"""🧪 作動流体の物性と閉じた式（fluids.py）: 断面・樋・継ぎ目、水温、粘度のつまみ、エタノール・水の混合物。

- fluids: 物性・閉じた式（Torricelli・Poiseuille）・断面・継ぎ目と曲がり
- temperature: 🧊 水温: 温度表から μ・ρ を引く（get_fluid(key, temp_c)）と、Network・spec への通し方。
- viscosity: 🍯 粘度のつまみ: get_fluid(mu=, mu_scale=)、Network・spec・draw への通し方、粘い小孔の Cd(Re)。
- mixture: 🍷 エタノール・水の混合物: 密度表・度数の換算・気液平衡・Rayleigh 式と、成分モードの水槽（boil・evaporate・カット）。
"""

import math
import random
import unittest

import pixelwater as pw
from fixedq import from_q, random_scales
from fluids import (
    ABV_CHECK_ROWS, AZEOTROPE_X, FLUIDS, M_ETHANOL, M_WATER, Section, TEMP_TABLES, abv, bubble_t, channel_flow,
    density_20, fall_time, fluid_table, get_fluid, lip_flow, lip_head, mix_volume, moles_for_abv, normal_depth,
    orifice_cd, poiseuille_tau, rayleigh_log, side_spill, superelevation, table_value, temperature_range,
    torricelli_empty_time, torricelli_head, viscous_orifice, vle_y,
)
from network import Network, build, simulate


# ──── fluids ──────────────────────────────────────────────────────

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

    def test_placeholders_scale_water_and_stay_out_of_experiments(self):
        from fluids import PLACEHOLDER_MU_SCALE, REAL_FLUIDS
        w = get_fluid("water")
        for key, k in PLACEHOLDER_MU_SCALE.items():
            f = get_fluid(key)
            self.assertTrue(f.placeholder)
            self.assertEqual(f.rho, w.rho)
            self.assertAlmostEqual(f.nu / w.nu, k, places=12)
            self.assertNotIn(key, REAL_FLUIDS)
        self.assertLess(get_fluid("namekian-water").nu, w.nu)  # さらさら
        self.assertGreater(get_fluid("poison-water").nu, w.nu)  # どろっと
        self.assertEqual(REAL_FLUIDS, ("water", "ethanol", "methanol"))

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


class JointsAndBends(unittest.TestCase):
    def test_top_width(self):
        sec = Section("semicircle", 0.04)
        self.assertAlmostEqual(sec.top_width(0.04), 0.08)
        self.assertEqual(sec.top_width(0.0), 0.0)
        self.assertEqual(Section("rect", 0.05).top_width(0.01), 0.05)

    def test_lip_head_inverts_lip_flow(self):
        sec = Section("semicircle", 0.04)
        for item in random_scales(12, 6, n_min=4, n_max=6):
            q = item["value"]
            h = lip_head(sec, q, 3e-3)
            self.assertAlmostEqual(lip_flow(sec, h, 3e-3) / q, 1.0, places=9)
        self.assertEqual(lip_flow(sec, 2e-3, 3e-3), 0.0)

    def test_fall_superelevation_side_spill(self):
        self.assertAlmostEqual(fall_time(0.05), math.sqrt(0.1 / 9.80665))
        self.assertEqual(fall_time(0.0), 0.0)
        self.assertAlmostEqual(superelevation(1.0, 0.08, 0.2), 0.08 / (9.80665 * 0.2))
        self.assertEqual(side_spill(0.1, 0.03, 0.04), 0.0)
        self.assertGreater(side_spill(0.1, 0.041, 0.04), 0.0)


# ──── temperature ─────────────────────────────────────────────────

"""🧊 水温: 温度表から μ・ρ を引く（get_fluid(key, temp_c)）と、Network・spec への通し方。"""


class Temperature(unittest.TestCase):
    def test_base_temperature_returns_the_same_fluid(self):
        for key in FLUIDS:
            base = get_fluid(key)
            self.assertIs(get_fluid(key, base.temp_c), base)
            self.assertIs(get_fluid(key, None), base)

    def test_table_nodes_are_reproduced(self):
        for t, mu, rho in TEMP_TABLES["water"]["points"]:
            got_mu, got_rho = table_value("water", t)
            self.assertAlmostEqual(got_mu / mu, 1.0, places=12)
            self.assertAlmostEqual(got_rho / rho, 1.0, places=12)
        # 水は 20 °C の値が表と FLUIDS で同じなので、表の点の μ がそのまま出る
        self.assertAlmostEqual(get_fluid("water", 10.0).mu / 1.3059e-3, 1.0, places=12)

    def test_viscosity_falls_with_temperature_for_random_temperatures(self):
        temps = sorted({4.0 * item["value"] for item in random_scales(91, 8, n_min=0, n_max=0)})
        for key in ("water", "ethanol", "namekian-water"):
            lo, hi = temperature_range(key)
            ts = [t for t in temps if lo <= t <= hi]
            mus = [get_fluid(key, t).mu for t in ts]
            self.assertEqual(mus, sorted(mus, reverse=True), key)
        lo, hi = temperature_range("methanol")
        self.assertEqual((lo, hi), (0.0, 25.0))
        self.assertGreater(get_fluid("methanol", 5.0).mu, get_fluid("methanol", 25.0).mu)

    def test_outside_the_table(self):
        with self.assertRaises(ValueError):
            get_fluid("water", 41.0)
        with self.assertRaises(ValueError):
            get_fluid("methanol", 30.0)

    def test_network_and_spec_take_temperature(self):
        warm = Network(1.0).tank("a", 0.01, head=0.3).tank("b", 1.0)
        cold = Network(1.0, temp_c=10.0).tank("a", 0.01, head=0.3).tank("b", 1.0)
        tw = warm.pipe("a", "b", radius=1e-3, length=0.5)["tau"]
        tc = cold.pipe("a", "b", radius=1e-3, length=0.5)["tau"]
        self.assertAlmostEqual(tc / tw, get_fluid("water", 10.0).nu / get_fluid("water").nu, places=12)
        res = simulate({"dt": 1.0, "temp_c": 10.0, "tanks": [{"name": "a", "area": 1, "volume": 1},
                                                             {"name": "b", "area": 1}],
                        "links": [{"law": "pump", "src": "a", "dst": "b", "rate": 0.1}], "t_end": 2.0})
        self.assertEqual(res["temp_c"], 10.0)
        self.assertEqual(build({"tanks": []}).fluid.temp_c, 20.0)


# ──── viscosity ───────────────────────────────────────────────────

"""🍯 粘度のつまみ: get_fluid(mu=, mu_scale=)、Network・spec・draw への通し方、粘い小孔の Cd(Re)。"""


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


# ──── mixture ─────────────────────────────────────────────────────

"""🍷 エタノール・水の混合物: 密度表・度数の換算・気液平衡・Rayleigh 式と、成分モードの水槽（boil・evaporate・カット）。"""


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
