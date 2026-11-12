"""💴 🪙 ⛏️ 同じ配管に流す物を置き換える（E26〜E28）: 定常・総量・Bateman・決済・描画。"""
import math
import unittest

import fluidfullspeedsimulation as ffs
import pixelwater as pw
from fixedq import from_q, random_scales, to_q


class Experiments(unittest.TestCase):
    def test_all_checks_pass(self):
        for eid in ("E26", "E27", "E28", "E29"):
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
            res = ffs.outlet_network([shares], tau=30.0, dt=0.5).run(60.0)
            want = 128.0 * ffs.outlet_factor([s / 30.0 for s in shares], 0.5) ** 120
            with self.subTest(share=item["label"]):
                self.assertAlmostEqual(from_q(res["volumes_q"]["pool-1"]), want, places=8)
                self.assertEqual(res["conservation_error_q"], 0)

    def test_scene_is_a_staircase(self):
        res, cols = ffs.rig_section("outlet-count")
        view = pw.network_frames(res, cols, frames=None, rows=16)
        self.assertEqual(ffs.scene_arg("🚿"), "outlet-count")
        for fr in view["frames"]:
            self.assertEqual(fr["counts"], sorted(fr["counts"], reverse=True))
        self.assertGreater(view["frames"][1]["counts"][0], view["frames"][1]["counts"][-1])

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


if __name__ == "__main__":
    unittest.main()
