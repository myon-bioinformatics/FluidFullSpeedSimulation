"""🔬 実験 E01〜E33（ID + 識別子）: 全検証の通しと、場面ごとの閉じた式・回路・定常値。

- experiments: E01〜E33 の全検証が通る・ID と識別子・乱数で選んだ刻みでの通し
- moniac_rig: E11 moniac-rig を据え置きの試験台にする。
- rain: 🌧️ 斜面に降る雨: 運動波の閉じた式、Network.rain、spec の rains、区間数での収束、E18 の場面。
- periodic: 🎍🍷🫗 溜まる → 閾値 → 一気に: tip / siphon / mariotte の法則と E19〜E21。
- water_clock: 🏺 断面が高さで変わる水槽（profile）と水時計（E22）。
- limits: 🚧 E23 inertia-limits: U 字管の閉じた式、エンジンの単調な緩和、Ritter の解、ダム決壊の外れ方、Joukowsky。
- distillation: ⚗️ E24 ワインの蒸留（still_network、カット、draw の ⚗️）。
- angels_share: 👼 天使の取り分（E25）: 較正・閉じた式・RH*・40 % を切る年数・木の吸い込み・描画。
- hub_flows: 💴 🪙 ⛏️ 同じ配管に流す物を置き換える（E26〜E28）・🚿 出口の数（E29）・🫀 血液の循環（E30）・🔥 火力発電（E31、💧 水・蒸気の循環・📉 部分負荷を含む）・❄️ ヒートポンプ（E32）・🧠 脳脊髄液（E33）: 定常・総量・Bateman・決済・描画。
"""

import contextlib
import io
import json
import math
import unittest

import fluidfullspeedsimulation as ffs
import pixelwater as pw
from fixedq import from_q, qmul, random_scales, to_q
from fluids import (
    G, Section, abv, channel_flow, get_fluid, mariotte_constant_time, mariotte_flow, mix_volume, mole_fraction,
    moles_for_abv, overland_equilibrium, overland_power_law, overland_recession, overland_rise,
    overland_t_e_power_law, plan_width, rayleigh_log, siphon_drain_time, siphon_drain_time_with_inflow, siphon_flow,
    tip_period, torricelli_empty_time, wave_celerity,
)
from network import Network, build, simulate


# ──── experiments ─────────────────────────────────────────────────

class Experiments(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = {eid: ffs.run_experiment(eid) for eid in ffs.EXPERIMENTS}

    def test_every_experiment_passes(self):
        for eid, r in self.results.items():
            with self.subTest(eid):
                failed = [c["id"] for c in r["checks"] if not c["ok"]]
                self.assertEqual(failed, [], ffs.render(r))
                self.assertTrue(r["ok"])

    def test_check_ids_are_unique_per_experiment(self):
        for eid, r in self.results.items():
            ids = [c["id"] for c in r["checks"]]
            self.assertEqual(len(ids), len(set(ids)), eid)

    def test_results_survive_json(self):
        for eid, r in self.results.items():
            self.assertEqual(json.loads(json.dumps(r, allow_nan=False)), r, eid)

    def test_ids_and_slugs_are_unique(self):
        slugs = [e["slug"] for e in ffs.list_experiments()]
        self.assertEqual(len(slugs), len(set(slugs)))
        for e in ffs.list_experiments():
            self.assertEqual(ffs.resolve(e["id"]), e["id"])
            self.assertEqual(ffs.resolve(e["slug"]), e["id"])
            self.assertEqual(ffs.resolve(f"{e['id'].lower()}:{e['slug']}"), e["id"])
        with self.assertRaises(KeyError):
            ffs.resolve("E99")


class RandomisedSweeps(unittest.TestCase):
    def test_integrator_exact_for_random_dt_and_rate(self):
        for d, r in zip(random_scales(31, 4, n_min=1, n_max=2), random_scales(32, 4, n_min=4, n_max=6)):
            with self.subTest(dt=d["label"], rate=r["label"]):
                res = ffs.run_experiment("E01", dt=d["value"], t_end=50 * d["value"], rate=r["value"])
                ok = {c["id"]: c["ok"] for c in res["checks"]}
                self.assertTrue(ok["integrate-exact"])
                self.assertTrue(ok["integrate-float"])
                self.assertTrue(ok["conserve-fill"])

    def test_float_valve_kc_split_holds_for_random_gains(self):
        """弱い弁は A×10^-4、強い弁は K_c×(1.5 + A)。どちらも乱数で選び、K_c をまたぐ結論が保たれるか。"""
        kc = ffs.run_experiment("E05")["summary"]["critical_gain"]
        for low, high in zip(random_scales(33, 2, n_min=4, n_max=4), random_scales(34, 2, n_min=0, n_max=0)):
            with self.subTest(low=low["label"], high=f"{1.5 + high['A']:.3g}·Kc"):
                res = ffs.run_experiment("E05", gain_low=low["value"], gain_high=kc * (1.5 + high["A"]))
                self.assertTrue(res["ok"], ffs.render(res))


# ──── moniac_rig ──────────────────────────────────────────────────

"""E11 moniac-rig を据え置きの試験台にする。

GOLDEN は試験台の最後の水位（有効 9 桁）。回路・法則・離散化を変えたときに、意図した変化かどうかを
ここで気づくための物差しで、正しさそのものの根拠ではない（正しさは E11 の閉じた式との比較が受け持つ）。
意図して変えたときは、理由をコミットに書いて GOLDEN を更新する。
"""


# GOLDEN は試験台の最後の水位（有効 9 桁）。回路・法則・離散化を変えたときに、意図した変化かどうかを
# ここで気づくための物差しで、正しさそのものの根拠ではない（正しさは E11 の閉じた式との比較が受け持つ）。
# 意図して変えたときは、理由をコミットに書いて GOLDEN を更新する。
#
# 入る量を「一定の速さで入る」と見なす既定（inflows="uniform"）での値。頭の槽は閉じた式 Q·τ/(A·Σo) = 0.3 と一致する。
# 入った量をステップの終わりに足す扱い（inflows="end"）では 0.300750625・0.135225125・0.200107155・0.191188099 だった
# 開き具合が 1 でない弁 1 本の槽を法則のまま（V·k·開き具合）で取っていたときは h:store 0.200026964・h:sump 0.191199455 だった。
# 1 − e^(−開き具合·dt/τ) に置き直した今の値は、刻みを細かくした先（dt/16 で h:store 0.2・h:sump 0.1912）に 100 倍ほど近い
GOLDEN = {"h:header": 0.3, "h:side": 0.135000562, "h:store": 0.200000201, "h:sump": 0.19119999}


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


# ──── rain ────────────────────────────────────────────────────────

"""🌧️ 斜面に降る雨: 運動波の閉じた式、Network.rain、spec の rains、区間数での収束、E18 の場面。"""


SEC = Section("rect", 1.0)


NU = get_fluid("water").nu


KW = {"nu": NU, "h_max": 0.01}


class RainClosedForms(unittest.TestCase):
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


# ──── periodic ────────────────────────────────────────────────────

"""🎍🍷🫗 溜まる → 閾値 → 一気に: tip / siphon / mariotte の法則と E19〜E21。"""


class PeriodicClosedForms(unittest.TestCase):
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


# ──── water_clock ─────────────────────────────────────────────────

"""🏺 断面が高さで変わる水槽（profile）と水時計（E22）。"""


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


# ──── limits ──────────────────────────────────────────────────────

"""🚧 E23 inertia-limits: U 字管の閉じた式、エンジンの単調な緩和、Ritter の解、ダム決壊の外れ方、Joukowsky。"""


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


# ──── distillation ────────────────────────────────────────────────

"""⚗️ E24 ワインの蒸留（still_network、カット、draw の ⚗️）。"""


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


# ──── angels_share ────────────────────────────────────────────────

"""👼 天使の取り分（E25）: 較正・閉じた式・RH*・40 % を切る年数・木の吸い込み・描画。"""


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


# ──── hub_flows ───────────────────────────────────────────────────

"""💴 🪙 ⛏️ 同じ配管に流す物を置き換える（E26〜E28）・🚿 出口の数（E29）・🫀 血液の循環（E30）・🔥 火力発電（E31）・❄️ ヒートポンプ（E32）: 定常・総量・Bateman・決済・描画。"""


class HubFlowChecks(unittest.TestCase):
    def test_all_checks_pass(self):
        for eid in ("E26", "E27", "E28", "E29", "E30"):
            slug, _icon, _title, fn = ffs.EXPERIMENTS[eid]
            with self.subTest(experiment=slug):
                bad = [c["id"] for c in fn()["checks"] if not c["ok"]]
                self.assertEqual(bad, [])


class Plumbing(unittest.TestCase):
    def test_same_positions_for_every_kind(self):
        shapes = set()
        for kind in ffs.HUB_FLOWS:
            nm = ffs.hub_names(kind)
            net = ffs.hub_network(kind)
            self.assertEqual(len(nm["hubs"]), 3)
            self.assertEqual(len(net.nodes), 8)
            shapes.add(tuple(sorted((ln["src"] == nm["source"], ln["dst"] == nm["outside"]) for ln in net.links
                                    if ln["src"] != nm["outside"])))
        self.assertEqual(len(shapes), 2)  # 💴 🪙 は同じ形、⛏️ は拠点 → 外 の管が無いだけ

    def test_money_source_is_finite_and_drawn_on_a_separate_scale(self):
        initial = ffs.HUB_FLOWS["money"]["source_volume"]
        res, cols = ffs.rig_section("money-flow")
        self.assertLess(from_q(res["volumes_q"]["central-bank"]), initial)
        self.assertEqual(res["total_q"], to_q(initial))
        view = pw.network_frames(res, cols, frames=8, rows=16)
        caps = dict(zip(view["columns"], view["caps"]))
        self.assertGreater(caps["central-bank"], caps["town"])
        self.assertEqual(caps["bank-a"], caps["town"])
        small = ffs.hub_network("money", source_volume=0.0).run(1.0)
        self.assertEqual(small["volumes_q"]["town"], 0)
        self.assertTrue(all(link["volume_moved_q"] == 0 for link in small["links"]))

    def test_money_scales_with_issue_and_gold_does_not_scale_with_speed(self):
        for item in random_scales(260, 4, n_min=-1, n_max=2):
            q = item["value"]
            with self.subTest(issue=item["label"]):
                a, b = ffs.hub_steady("money", issue=q), ffs.hub_steady("money", issue=2 * q)
                self.assertAlmostEqual(b["town"] / a["town"], 2.0, places=12)
        base = ffs.hub_steady("gold")
        p = ffs.HUB_FLOWS["gold"]
        for item in random_scales(270, 4, n_min=-1, n_max=1):
            c = item["value"]
            over = {k: c * p[k] for k in ("tau_out", "tau_arm", "tau_sat", "tau_back", "delay_arm")}
            slow = ffs.hub_steady("gold", tau_feed=tuple(c * t for t in p["tau_feed"]), **over)
            with self.subTest(slow=item["label"]):
                for name in ("vault", "bank-a", "sec-c", "town"):
                    self.assertAlmostEqual(slow[name], base[name], places=9)

    def test_gold_steady_adds_up_to_the_total(self):
        st = ffs.hub_steady("gold")
        held = sum(v for k, v in st.items() if k not in ("feeds", "cycle_time"))
        self.assertAlmostEqual(held, ffs.HUB_FLOWS["gold"]["source_volume"], places=9)
        with self.assertRaises(ValueError):
            ffs.hub_steady("rare-metal")

    def test_gold_total_is_exact_at_every_sample(self):
        net = ffs.hub_network("gold")
        total = net.total_q()
        for n in range(400):
            net.step(n)
            self.assertEqual(net.total_q(), total)
        self.assertEqual(total, to_q(128.0))

    def test_bateman(self):
        for item in random_scales(280, 4, n_min=1, n_max=2):
            k1 = item["value"]
            k2, k3 = 2.5 * k1, 0.4 * k1
            with self.subTest(k1=item["label"]):
                t = 1.7 / k1
                n = ffs.bateman([k1, k2, k3], t, 3.0)
                self.assertAlmostEqual(n[0], 3.0 * math.exp(-k1 * t), places=12)
                self.assertAlmostEqual(n[1], 3.0 * k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t)), places=12)
                # 3 段目は数値積分と比べる（dN3/dt = k2·N2 − k3·N3）
                steps, acc = 20000, 0.0
                h = t / steps
                for i in range(steps):
                    s = (i + 0.5) * h
                    acc += ffs.bateman([k1, k2], s, 3.0)[1] * k2 * math.exp(-k3 * (t - s)) * h
                self.assertAlmostEqual(n[2], acc, places=6)
                tp = ffs.bateman_peak_time(k1, k2)
                peak = ffs.bateman([k1, k2], tp)[1]
                self.assertGreater(peak, ffs.bateman([k1, k2], tp * 0.99)[1])
                self.assertGreater(peak, ffs.bateman([k1, k2], tp * 1.01)[1])
        self.assertEqual(ffs.bateman([1.0, 2.0], -1.0), [0.0, 0.0])

    def test_settlement_shortfall_is_exact(self):
        for item in random_scales(261, 4, n_min=1, n_max=4):
            x = item["value"]
            res = ffs.settlement(1.0 - x)
            with self.subTest(short=item["label"]):
                self.assertEqual(res["links"][0]["volume_moved_q"], to_q(1.0 - x))
                self.assertEqual(res["conservation_error_q"], 0)
        res = ffs.settlement(0.25 + 0.75 / 64, overlap=True)
        self.assertEqual(res["links"][0]["volume_moved_q"], to_q(1.0))
        res = ffs.settlement(0.25, overlap=True)  # 1 ステップ分が無いと払い切れない
        self.assertLess(res["links"][0]["volume_moved_q"], to_q(1.0))


class Outlets(unittest.TestCase):
    def test_one_outlet_is_exact_and_more_outlets_drain_faster(self):
        for item in random_scales(290, 4, n_min=0, n_max=2):  # dt/τ = 0.001〜9.99 … 粗い刻みは 0.5 まで
            ratio = min(item["value"], 0.5)
            with self.subTest(dt_over_tau=item["label"]):
                factors = [ffs.outlet_factor([1.0 / n] * n, ratio) for n in (1, 2, 4, 8, 64)]
                self.assertAlmostEqual(factors[0], math.exp(-ratio), places=15)
                self.assertEqual(factors, sorted(factors, reverse=True))
                self.assertGreater(factors[-1], 1.0 - ratio)  # n → ∞ でも 1 − K·dt（前進オイラー）より下には行かない

    def test_engine_matches_the_factor_for_uneven_outlets(self):
        for item in random_scales(291, 3, n_min=1, n_max=2):
            a = item["value"]
            shares = [a, 0.75 * (1.0 - a), 0.25 * (1.0 - a)]
            res = ffs.outlet_network([shares], tau=30.0, dt=0.5, outlets="independent").run(60.0)
            want = 128.0 * ffs.outlet_factor([s / 30.0 for s in shares], 0.5) ** 120
            with self.subTest(share=item["label"]):
                self.assertAlmostEqual(from_q(res["volumes_q"]["pool-1"]), want, places=8)
                self.assertEqual(res["conservation_error_q"], 0)

    def test_scene_is_a_staircase(self):
        res, cols = ffs.rig_section("outlet-count")
        view = pw.network_frames(res, cols, frames=None, rows=16)
        self.assertEqual(ffs.scene_arg("🚿"), "outlet-count")
        for fr in view["frames"]:  # 別々に取る 1・2・4・8 本は階段、厳密に解く 8 本（5 列目）は 1 本と同じ
            self.assertEqual(fr["counts"][:4], sorted(fr["counts"][:4], reverse=True))
            self.assertEqual(fr["counts"][4], fr["counts"][0])
        self.assertGreater(view["frames"][1]["counts"][0], view["frames"][1]["counts"][3])

    def test_exact_outlets_for_any_count_split_and_step(self):
        import math as m
        for item in random_scales(292, 4, n_min=0, n_max=2):  # dt/τ = 0.001〜9.99 … 粗い刻みは 0.5 まで
            h = min(item["value"], 0.5)
            rng_shares = [[1.0], [0.5, 0.5], [0.7, 0.2, 0.1], [1 / 16] * 16]
            steps = max(1, round(2.0 / h))
            res = ffs.outlet_network(rng_shares, tau=30.0, dt=30.0 * h).run(steps * 30.0 * h)
            want = 128.0 * m.exp(-steps * h)
            with self.subTest(dt_over_tau=item["label"]):
                for i in range(1, 5):
                    self.assertAlmostEqual(from_q(res["volumes_q"][f"pool-{i}"]) / want, 1.0, places=9)
                moved = [from_q(res["volumes_q"][f"out-3-{j}"]) for j in (1, 2, 3)]
                for got, share in zip(moved, (0.7, 0.2, 0.1)):
                    self.assertAlmostEqual(got / sum(moved), share, places=9)
                self.assertEqual(res["conservation_error_q"], 0)

    def test_exact_outlets_with_a_pump_and_when_the_tank_runs_short(self):
        import math as m
        tau, vol = 30.0, 128.0
        for share, t in ((1 / 8, 60.0), (1 / 64, 15.0)):
            q = vol / tau * share
            want = (vol + q * tau) * m.exp(-t / tau) - q * tau
            for n in (1, 2, 4):
                res = ffs.outlet_network([[1.0 / n] * n], tau=tau, dt=0.5, pump_extra=share).run(t)
                with self.subTest(share=share, n=n):
                    self.assertAlmostEqual(from_q(res["volumes_q"]["pool-1"]) / want, 1.0, places=9)
                    self.assertAlmostEqual(from_q(res["volumes_q"]["pumped-1"]), q * t, places=9)
        res = ffs.outlet_network([[0.5, 0.5]], tau=tau, dt=0.5, pump_extra=2.0).run(60.0)  # 途中で空になる
        self.assertEqual(res["volumes_q"]["pool-1"], 0)
        self.assertEqual(res["conservation_error_q"], 0)

    def test_outlet_mode_switch(self):
        import network
        with self.assertRaises(ValueError):
            network.Network(dt=1.0, outlets="sometimes")
        with self.assertRaises(ValueError):
            network.Network(dt=1.0).tank("a", 1.0, outlets="sometimes")
        shares = [[0.25] * 4, [0.25] * 4]
        mixed = ffs.outlet_network(shares, tau=30.0, dt=0.5, modes=("independent", None)).run(60.0)
        both = ffs.outlet_network(shares, tau=30.0, dt=0.5, outlets="independent", modes=(None, "exact")).run(60.0)
        self.assertEqual(mixed["volumes_q"], both["volumes_q"])  # 槽ごとの指定が全体の指定より先
        self.assertLess(mixed["volumes_q"]["pool-1"], mixed["volumes_q"]["pool-2"])
        spec = {"dt": 0.5, "t_end": 60.0, "tanks": [{"name": "a", "area": 1.0, "volume": 128.0},
                                                      {"name": "b", "area": 1.0}, {"name": "c", "area": 1.0}],
                "links": [{"law": "valve", "src": "a", "dst": "b", "tau": 60.0},
                          {"law": "valve", "src": "a", "dst": "c", "tau": 60.0}]}
        left = {mode: network.simulate({**spec, "outlets": mode})["volumes_q"]["a"] for mode in network.OUTLET_MODES}
        self.assertLess(left["independent"], left["exact"])
        self.assertAlmostEqual(from_q(left["exact"]) / (128.0 * 0.1353352832366127), 1.0, places=9)

    def test_single_valve_tanks_are_untouched_by_the_mode(self):
        a = ffs.outlet_network([[1.0]], tau=30.0, dt=0.5).run(60.0)
        b = ffs.outlet_network([[1.0]], tau=30.0, dt=0.5, outlets="independent").run(60.0)
        self.assertEqual(a["volumes_q"], b["volumes_q"])

    def test_uniform_inflows_pump_fed_and_combined_tank_are_exact(self):
        import math as m
        import network
        for item in random_scales(293, 4, n_min=0, n_max=2):
            h = min(item["value"], 0.5)
            tau, dt = 20.0, 20.0 * h
            steps = max(1, round(2.0 / h))
            t = steps * dt
            got = {}
            for mode in network.INFLOW_MODES:
                net = network.Network(dt=dt, inflows=mode)
                net.tank("sea", 1.0, volume=1e6).tank("b", 1.0).tank("s", 1.0)
                net.pump("sea", "b", rate=3.0)
                net.valve("b", "s", tau=tau)
                res = net.run(t)
                got[mode] = from_q(res["volumes_q"]["b"]) / (3.0 * tau * -m.expm1(-t / tau)) - 1.0
                self.assertEqual(res["conservation_error_q"], 0)
            with self.subTest(dt_over_tau=item["label"]):
                self.assertAlmostEqual(got["uniform"], 0.0, places=9)
                self.assertAlmostEqual(got["end"], h / -m.expm1(-h) - 1.0, places=9)  # 入った量が半ステップ減りそこねる
            # 入るポンプ + 弁 3 本 + 出るポンプ: V∞ + (V0 − V∞)·e^(−K·t)
            net = network.Network(dt=dt)
            net.tank("sea", 1.0, volume=1e6).tank("p", 1.0, volume=100.0).tank("gone", 1.0)
            net.pump("sea", "p", rate=2.0)
            for j, share in enumerate((0.6, 0.3, 0.1)):
                net.tank(f"o{j}", 1.0)
                net.valve("p", f"o{j}", tau=tau / share)
            net.pump("p", "gone", rate=0.5)
            res = net.run(t)
            settle = (2.0 - 0.5) * tau
            with self.subTest(combined=item["label"]):
                self.assertAlmostEqual(from_q(res["volumes_q"]["p"]) / (settle + (100.0 - settle) * m.exp(-t / tau)), 1.0, places=9)
                self.assertEqual(res["conservation_error_q"], 0)

    def test_uniform_inflows_fill_an_empty_tank_and_keep_totals_in_a_ring(self):
        import network
        net = network.Network(dt=4.0)
        net.tank("a", 1.0, volume=128.0).tank("b", 1.0)
        net.valve("a", "b", tau=16.0)
        net.valve("b", "a", tau=8.0)
        total = net.total_q()
        for n in range(200):
            net.step(n)
            self.assertEqual(net.total_q(), total)
            self.assertGreaterEqual(min(net.nodes["a"]["v"], net.nodes["b"]["v"]), 0)
        # 行き渡った後の取り分は滞在時間の比（16 : 8）
        self.assertAlmostEqual(from_q(net.nodes["a"]["v"]) / 128.0, 16.0 / 24.0, places=9)

    def test_inflow_mode_switch(self):
        import network
        with self.assertRaises(ValueError):
            network.Network(dt=1.0, inflows="sometimes")
        spec = {"dt": 2.0, "t_end": 40.0, "tanks": [{"name": "a", "area": 1.0, "volume": 128.0},
                                                      {"name": "b", "area": 1.0}, {"name": "c", "area": 1.0}],
                "links": [{"law": "valve", "src": "a", "dst": "b", "tau": 20.0},
                          {"law": "valve", "src": "b", "dst": "c", "tau": 10.0}]}
        left = {mode: network.simulate({**spec, "inflows": mode})["volumes_q"]["b"] for mode in network.INFLOW_MODES}
        self.assertEqual(network.simulate(spec)["volumes_q"]["b"], left["uniform"])  # 既定は uniform
        self.assertGreater(left["end"], left["uniform"])  # 終わりに足す扱いは下の槽が多く出る

    def test_short_gifs_play_slower(self):
        self.assertEqual(pw._section_delays(40), [80] + [10] * 38 + [200])
        self.assertEqual(pw._section_delays(8), [80] + [30] * 6 + [200])


class Scenes(unittest.TestCase):
    def test_columns_and_shared_caps(self):
        for scene, icon in (("money-flow", "💴"), ("gold-flow", "🪙"), ("rare-metal", "⛏️")):
            with self.subTest(scene=scene):
                self.assertEqual(ffs.scene_arg(icon), scene)
                res, cols = ffs.rig_section(scene)
                self.assertEqual(res["conservation_error_q"], 0)
                self.assertEqual(sum(">" in c for c in cols), 3)
                view = pw.network_frames(res, cols, frames=8, rows=16)
                tanks = [c for c, k in zip(view["caps"], view["kinds"]) if k == "tank"]
                self.assertEqual(len(set(tanks[1:])), 1)  # 源以外は同じ縦いっぱい
                if scene != "money-flow":
                    self.assertEqual(set(tanks), {128.0})
                    self.assertEqual(view["frames"][0]["counts"][0], max(view["frames"][0]["counts"]))

    def test_gold_blocks_show_the_shares(self):
        res, cols = ffs.rig_section("gold-flow")
        view = pw.network_frames(res, cols, frames=8, rows=16)
        last = dict(zip(view["columns"], view["frames"][-1]["counts"]))
        steady = ffs.hub_steady("gold")
        for name in ("vault", "bank-a", "town"):
            self.assertAlmostEqual(last[name], steady[name] / 8.0, delta=0.75)  # 1 段 = 8 t、±半段 + 行き渡る途中


class Circulation(unittest.TestCase):
    """🫀 心臓 → 動脈 → 器官 7 つ → 静脈 → 肺 → 心臓。"""

    def test_steady_adds_up_and_every_tank_is_flow_times_transit(self):
        p, tau = ffs.CIRCULATION, ffs.circulation_taus()
        for dt in (None, p["dt"]):
            st = ffs.circulation_steady(dt=dt)
            with self.subTest(dt=dt):
                self.assertFalse(st["limited"])
                self.assertAlmostEqual(sum(st[n] for n in ffs.circulation_columns()), p["blood"], places=12)
        st = ffs.circulation_steady()
        for organ, share in p["flows"].items():  # 中心容積の原理 V = F·T（どの器官も滞在時間は同じ）
            self.assertAlmostEqual(st[organ], p["output"] * share * tau["organ"], places=12)
        self.assertAlmostEqual(st["heart"], p["blood"] * 0.07, places=12)  # 表の残り 7 % が心臓
        self.assertAlmostEqual(sum(p["flows"].values()), 1.0, places=12)

    def test_heart_sends_only_what_comes_back(self):
        p = ffs.CIRCULATION
        base = ffs.circulation_steady(dt=p["dt"])
        keep = base["heart"] - base["output"] * p["dt"]  # これだけ失っても送る速さは変わらない
        for item in random_scales(300, 6, n_min=1, n_max=3):  # 失う割合 0.1〜99.9 %
            x = item["value"]
            st = ffs.circulation_steady(blood=p["blood"] * (1.0 - x), dt=p["dt"])
            with self.subTest(loss=item["label"]):
                self.assertEqual(st["limited"], p["blood"] * x > keep)
                if st["limited"]:  # 空の心臓: 送る速さは血液量に比例
                    self.assertAlmostEqual(st["output"] / (1.0 - x), p["blood"] / (base["rest_time"] + p["dt"]),
                                           places=12)
                    self.assertAlmostEqual(st["heart"], st["output"] * p["dt"], places=12)
                else:
                    self.assertEqual(st["output"], p["output"])

    def test_total_is_exact_and_brain_balance_settles_to_zero(self):
        net = ffs.circulation_network()
        total = net.total_q()
        for n in range(400):
            net.step(n)
            self.assertEqual(net.total_q(), total)
        self.assertEqual(total, to_q(ffs.CIRCULATION["blood"]))
        res = ffs.circulation_network().run(40.0 * ffs.circulation_taus()["vein"])
        r = ffs.node_rates(res)["brain"]
        self.assertAlmostEqual(r["in"], r["out"], places=12)
        self.assertAlmostEqual(r["in"] * 60.0, 0.15 * 5.0, places=9)  # 脳 750 mL/min

    def test_scene_heart_drains_and_others_fill(self):
        self.assertEqual(ffs.scene_arg("🫀"), "circulation")
        res, cols = ffs.rig_section("circulation")
        self.assertEqual(cols, ffs.circulation_columns())
        view = pw.network_frames(res, cols, frames=None, rows=16)
        heart = [fr["counts"][0] for fr in view["frames"]]
        self.assertEqual(heart, sorted(heart, reverse=True))  # 心臓は減るだけ（行き渡るまで戻りが送る速さに届かない）
        for col in range(1, len(cols)):
            filled = [fr["counts"][col] for fr in view["frames"]]
            self.assertEqual(filled, sorted(filled))
        self.assertEqual(res["conservation_error_q"], 0)


class ThermalPower(unittest.TestCase):
    """🔥 燃料 → 熱が溜まる段 → 電気・煙突・海（エネルギーだけを流す）。"""

    def test_plan_meets_the_efficiency_and_every_stage_splits_to_one(self):
        for kind, p in ffs.THERMAL["plants"].items():
            plan = ffs.thermal_plan(kind)
            st = ffs.thermal_steady(kind)
            with self.subTest(plant=kind):
                for _name, _tau, outs in plan["stages"]:
                    self.assertAlmostEqual(sum(outs.values()), 1.0, places=12)
                    self.assertTrue(all(0.0 < v < 1.0 for v in outs.values()))
                self.assertAlmostEqual(st["gross"] / st["fuel_lhv"], p["efficiency"], places=12)
                out = st["grid"] + st["aux"] + st["stack"] + st["sea"]  # 定常では燃やした HHV が全部どこかへ出る
                self.assertAlmostEqual(out / st["fuel_hhv"], 1.0, places=12)
                self.assertGreater(st["sea"] + st["stack"], st["grid"] * (1.0 - p["efficiency"]) / p["efficiency"] * 0.9)

    def test_combined_cycle_formula(self):
        outs = {name: o for name, _t, o in ffs.thermal_plan("lng")["stages"]}
        gt, rec, st = outs["gt"]["power"], outs["hrsg"]["steam"], outs["steam"]["power"]
        self.assertAlmostEqual(gt + (1.0 - gt) * rec * st, ffs.THERMAL["plants"]["lng"]["efficiency"], places=12)

    def test_curzon_ahlborn_table_and_carnot(self):
        published = {"West Thurrock": (0.64, 0.40), "CANDU": (0.48, 0.28), "Larderello": (0.33, 0.178)}
        for name, tc, th, _seen in ffs.THERMAL["curzon_ahlborn"]:
            cn, ca = published[name]
            with self.subTest(plant=name):
                self.assertAlmostEqual(ffs.carnot(tc, th), cn, delta=0.0051)  # Larderello は 0.325（表により 0.33 / 0.32）
                self.assertAlmostEqual(ffs.curzon_ahlborn(tc, th), ca, delta=0.0005 if name == "Larderello" else 0.005)
        for item in random_scales(310, 6, n_min=-3, n_max=-2):  # 高温 100〜9990 °C
            th = item["value"]
            with self.subTest(t_hot=item["label"]):
                self.assertLess(ffs.curzon_ahlborn(25.0, th), ffs.carnot(25.0, th))

    def test_lag_scales_with_slow_and_shares_do_not(self):
        for kind in ffs.THERMAL["plants"]:
            base = ffs.thermal_steady(kind)
            for item in random_scales(311, 4, n_min=-1, n_max=1):  # 0.1〜99.9 倍
                c = item["value"]
                st = ffs.thermal_steady(kind, slow=c)
                with self.subTest(plant=kind, slow=item["label"]):
                    self.assertAlmostEqual(st["lag"] / base["lag"], c, places=9)
                    self.assertAlmostEqual(st["grid"], base["grid"], places=9)

    def test_scenes(self):
        self.assertEqual((ffs.scene_arg("🪨"), ffs.scene_arg("🌀")), ("coal-power", "combined-cycle"))
        for scene, kind in ffs.THERMAL_SCENES.items():
            res, cols = ffs.rig_section(scene)
            self.assertEqual(cols, ffs.thermal_columns(kind))
            view = pw.network_frames(res, cols, frames=None, rows=16)
            fuel = [fr["counts"][0] for fr in view["frames"]]
            self.assertEqual(fuel, sorted(fuel, reverse=True))
            for col in range(len(cols) - 4, len(cols)):  # 出口は溜まる一方
                got = [fr["counts"][col] for fr in view["frames"]]
                self.assertEqual(got, sorted(got))
            self.assertEqual(res["volumes_q"]["fuel"], 0)  # ちょうど燃やし切る
            self.assertEqual(res["conservation_error_q"], 0)


class SteamLoop(unittest.TestCase):
    """💧 E31 の 2 段目: 水は回り、熱は水に乗って運ばれる（🪨 橘湾 2 号の流量 + 蒸気表、🌀 単圧の排熱回収ボイラ）。"""

    def _balance(self, plan):
        inflow, heat = {}, {}
        for src, w in plan["water"].items():
            for dst, share in w["outs"].items():
                inflow[dst] = inflow.get(dst, 0.0) + w["flow"] * share
                heat[dst] = heat.get(dst, 0.0) + w["flow"] * share * w["h"]
        return inflow, heat

    def test_every_tank_balances_water_and_mixing(self):
        for kind, n in (("coal", None), ("lng", None), ("coal", 1), ("coal", 3), ("coal", 8)):
            plan = ffs.steam_loop_plan(kind, n)
            inflow, heat = self._balance(plan)
            heated = {d for spec in plan["heat"].values() for d in spec["outs"]}
            with self.subTest(kind=kind, stages=n):
                for name, w in plan["water"].items():
                    self.assertAlmostEqual(inflow[name] / w["flow"], 1.0, places=9)  # 入る水 = 出る水
                    if name not in heated and "pull" not in w:  # 混ざるだけの槽（給水加熱器・給水）は混ざった h
                        self.assertAlmostEqual(heat[name] / inflow[name], w["h"], places=6)
                q_in = plan["q_in"]
                self.assertAlmostEqual(q_in, sum(plan["work"].values()) - plan["work"].get("gt", 0.0) + plan["sea"], places=6)

    def test_more_stages_move_the_exhaust_towards_wet(self):
        h_f, h_g = ffs.STEAM_LOOP["saturation"][:2]
        xs = [(ffs.steam_loop_plan("coal", n)["h_ex"] - h_f) / (h_g - h_f) for n in range(1, 9)]
        self.assertEqual(xs, sorted(xs, reverse=True))
        self.assertGreater(xs[0], 1.0)  # まとめた 1 本だと過熱に見える

    def test_water_does_not_move_while_the_plant_heats_up(self):
        for kind in ("coal", "lng"):
            for item in random_scales(312, 2, n_min=1, n_max=1):  # dt 0.1〜0.999 s
                net = ffs.steam_loop_network(kind, dt=item["value"], t_fuel=32.0)
                start = {n: net.nodes[n]["parts"]["water"] for n in ffs.steam_loop_plan(kind)["water"]}
                res = net.run(32.0)
                with self.subTest(kind=kind, dt=item["label"]):
                    for n, w in start.items():
                        self.assertAlmostEqual(res["parts_q"][n]["water"] / w, 1.0, places=9)
                    self.assertEqual(res["parts_error_q"], {"water": 0, "heat": 0})

    def test_scenes(self):
        self.assertEqual((ffs.scene_arg("💧"), ffs.scene_arg("♨️")), ("steam-loop", "hrsg-loop"))
        for scene, hot, cold in (("steam-loop", "reheater", "condenser"), ("hrsg-loop", "hrsg", "condenser")):
            res, cols = ffs.rig_section(scene)
            got = dict(zip(cols, pw.network_frames(res, cols, frames=None, rows=16)["frames"][-1]["counts"]))
            with self.subTest(scene=scene):
                self.assertEqual(max(got.values()), got[hot])
                self.assertEqual(min(got.values()), got[cold])


class PartLoad(unittest.TestCase):
    """📉 Willans の直線: 無負荷の燃料は出力によらないので、出力を下げるほど効率が下がる。"""

    def test_curve(self):
        for kind, p in ffs.THERMAL["plants"].items():
            a = ffs.willans_no_load(kind)
            with self.subTest(kind=kind):
                self.assertTrue(0.0 < a < 0.2)
                self.assertAlmostEqual(ffs.part_load_efficiency(kind, 1.0), p["efficiency"], places=12)
                xs = sorted(item["value"] for item in random_scales(313, 6, n_min=1, n_max=1))
                etas = [ffs.part_load_efficiency(kind, x) for x in xs]
                self.assertEqual(etas, sorted(etas))  # 出力が上がるほど効率も上がる
        self.assertAlmostEqual(ffs.part_load_efficiency("lng", 0.5) / ffs.THERMAL["plants"]["lng"]["efficiency"], 0.9,
                               places=12)
        self.assertGreater(ffs.part_load_efficiency("coal", 0.5) / ffs.THERMAL["plants"]["coal"]["efficiency"],
                           ffs.part_load_efficiency("lng", 0.5) / ffs.THERMAL["plants"]["lng"]["efficiency"])

    def test_load_step_burns_the_fuel_exactly(self):
        self.assertEqual(ffs.scene_arg("📉"), "part-load")
        res, cols = ffs.rig_section("part-load")
        self.assertEqual(res["volumes_q"]["fuel"], 0)
        self.assertEqual(res["conservation_error_q"], 0)
        net = ffs.thermal_network("coal", willans=False)
        self.assertEqual(len([ln for ln in net.links if ln["kind"] == "pump"]), 2)  # 全負荷・線形は今までどおり 2 本


class HeatPump(unittest.TestCase):
    """❄️ 電気 + 外気の熱 → 部屋（🔥 の裏返し）。"""

    def test_cop_balance_and_carnot(self):
        p = ffs.HEAT_PUMP
        st = ffs.heat_pump_steady()
        self.assertAlmostEqual(st["cop"], p["heat"] / p["power"], places=12)
        self.assertAlmostEqual(st["room"], p["heat"], places=12)
        self.assertLess(st["cop"], ffs.carnot_cop(p["t_out"], p["t_room"]))
        for item in random_scales(320, 6, n_min=-2, n_max=0):  # Tc 1〜999 °C、Th = Tc + 0.1〜99.9 K
            tc = item["value"]
            th = tc + item["A"] * 10 ** (-item["N"] - 1)
            with self.subTest(t_cold=item["label"]):
                self.assertAlmostEqual(ffs.carnot_cop(tc, th) * ffs.carnot(tc, th), 1.0, places=12)
                self.assertGreater(ffs.carnot_cop(tc, th), 1.0)

    def test_running_on_plant_power_beats_burning_at_home(self):
        cop = ffs.heat_pump_steady()["cop"]
        for kind in ffs.THERMAL["plants"]:
            ts = ffs.thermal_steady(kind)
            eta = ts["grid"] / ts["fuel_hhv"]
            with self.subTest(plant=kind):
                self.assertGreater(eta * cop, 1.0)
                self.assertLess(eta, 1.0)

    def test_scene(self):
        self.assertEqual(ffs.scene_arg("❄️"), "heat-pump")
        res, cols = ffs.rig_section("heat-pump")
        self.assertEqual(cols, ffs.heat_pump_columns())
        view = pw.network_frames(res, cols, frames=None, rows=16)
        for col in (0, 1):  # 電気と外気は減るだけ、部屋は溜まるだけ
            got = [fr["counts"][col] for fr in view["frames"]]
            self.assertEqual(got, sorted(got, reverse=True))
        room = [fr["counts"][-1] for fr in view["frames"]]
        self.assertEqual(room, sorted(room))
        self.assertEqual((res["volumes_q"]["power"], res["volumes_q"]["outdoor"]), (0, 0))
        self.assertEqual(res["conservation_error_q"], 0)


class CerebrospinalFluid(unittest.TestCase):
    """🧠 血液 → 脳室 → くも膜下腔 → 血液。"""

    def test_steady_turnover_and_opening(self):
        st = ffs.csf_steady()
        self.assertAlmostEqual(st["ventricles"] + st["subarachnoid"], 150.0, places=12)
        self.assertAlmostEqual(st["turnover_per_day"], 1440.0 * 0.35 / 150.0, places=12)
        for item in random_scales(330, 5, n_min=1, n_max=2):  # 吸収の通りやすさ 0.01〜0.999 倍
            o = item["value"]
            with self.subTest(opening=item["label"]):
                self.assertAlmostEqual(ffs.csf_steady(opening=o)["subarachnoid"], 125.0 / o, places=9)
                self.assertAlmostEqual(ffs.csf_steady(opening=o)["ventricles"], 25.0, places=12)

    def test_running_head_does_not_move_and_blocked_head_gains_production(self):
        res = ffs.csf_network((1.0, 0.0), steady=True).run(720.0)
        st = ffs.csf_steady(dt=ffs.CSF["dt"])
        self.assertAlmostEqual(from_q(res["volumes_q"]["subarachnoid-1"]), st["subarachnoid"], places=9)
        gained = from_q(res["volumes_q"]["subarachnoid-2"] + res["volumes_q"]["ventricles-2"])
        self.assertAlmostEqual(gained - (st["subarachnoid"] + st["ventricles"]), 0.35 * 720.0, places=9)
        self.assertEqual(res["conservation_error_q"], 0)

    def test_scene(self):
        self.assertEqual(ffs.scene_arg("🧠"), "csf")
        res, cols = ffs.rig_section("csf")
        view = pw.network_frames(res, cols, frames=None, rows=16)
        normal = {fr["counts"][1] for fr in view["frames"]}
        blocked = [fr["counts"][3] for fr in view["frames"]]
        self.assertEqual(len(normal), 1)  # ふつうの頭は動かない（±0）
        self.assertEqual(blocked, sorted(blocked))
        self.assertEqual(blocked[-1], 16)


if __name__ == "__main__":
    unittest.main()
