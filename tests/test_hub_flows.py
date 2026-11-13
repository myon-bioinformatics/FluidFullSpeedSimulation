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


if __name__ == "__main__":
    unittest.main()


# --- 🧮 丸め方 carry64（係数を 2^-64 刻みにして、切り捨ての残りを繰り越す） ---------------------
def _hp_two_stage(dt: float, tau1: float, tau2: float, steps: int, v0: float):
    """2 段の槽（uniform）を 60 桁の小数で進めた、下の槽の量。"""
    from decimal import Decimal, localcontext

    with localcontext() as ctx:
        ctx.prec = 60
        x1, x2 = Decimal(dt) / Decimal(tau1), Decimal(dt) / Decimal(tau2)
        g1, g2 = 1 - (-x1).exp(), 1 - (-x2).exp()
        through = 1 - g2 / x2
        a, b = Decimal(v0), Decimal(0)
        for _ in range(steps):
            out = a * g1
            a, b = a - out, b + out - (b * g2 + out * through)
        return b


def _two_stage(rounding: str, dt: float, v0: float, steps: int):
    from network import Network

    net = Network(dt=dt, rounding=rounding)
    net.tank("a", 1.0, volume=v0)
    net.tank("b", 1.0)
    net.tank("sink", 1.0)
    net.valve("a", "b", tau=32.0)
    net.valve("b", "sink", tau=19.2)
    start = net.total_q()
    for n in range(steps):
        net.step(n)
    return net, start


def test_exp_neg_q_matches_decimal():
    from decimal import Decimal, localcontext
    from fractions import Fraction

    from fixedq import exp_neg_q

    for x in (Fraction(1, 512), Fraction(5, 96), Fraction(3, 4), Fraction(7, 2), Fraction(40), Fraction(0)):
        with localcontext() as ctx:
            ctx.prec = 60
            want = (-(Decimal(x.numerator) / Decimal(x.denominator))).exp() * (1 << 64)
            assert abs(exp_neg_q(x, 64) - want) <= Decimal("0.5000001")


def test_carry64_keeps_total_and_stays_within_a_few_ulps():
    from fixedq import ONE

    dt, steps = 0.25, 384
    worst = {}
    for rounding in ("q48", "carry64"):
        for v0 in (128.0, 2.0 ** -12):
            net, start = _two_stage(rounding, dt, v0, steps)
            assert net.total_q() == start  # 総量は 1 ulp も変わらない
            ulps = abs(net.nodes["b"]["v"] - _hp_two_stage(dt, 32.0, 19.2, steps, v0) * ONE)
            worst[rounding] = max(worst.get(rounding, 0), ulps)
    assert worst["carry64"] <= 4
    assert worst["q48"] > worst["carry64"]


def test_carry64_only_touches_valve_only_tanks_and_default_is_q48():
    from network import ROUNDING_MODES, Network, build

    assert ROUNDING_MODES[0] == "q48" and Network(dt=1.0).rounding == "q48"
    try:
        Network(dt=1.0, rounding="q64")
    except ValueError:
        pass
    else:
        raise AssertionError("知らない丸め方は断る")
    spec = {"dt": 1.0, "rounding": "carry64", "tanks": [{"name": "a", "area": 1.0, "volume": 8.0}, {"name": "b", "area": 1.0},
                                                         {"name": "c", "area": 1.0}],
            "links": [{"law": "valve", "src": "a", "dst": "b", "tau": 4.0}, {"law": "valve", "src": "a", "dst": "c", "tau": 8.0},
                      {"law": "pump", "src": "b", "dst": "c", "rate": 0.01}]}
    net = build(spec)
    start = net.total_q()
    for n in range(64):
        net.step(n)
    assert net.total_q() == start
    assert net._carry_tank("a", net._by_src()["a"]) is not None  # 弁 2 本の槽は carry64
    assert net._carry_tank("b", net._by_src()["b"]) is None  # 弁が無い槽（ポンプだけ）は q48 のまま


def test_carry64_with_a_pump_out_matches_the_exact_answer():
    from decimal import Decimal, localcontext

    from fixedq import ONE
    from network import Network

    dt, steps, v0, tau, rate = 0.25, 128, 128.0, 32.0, 0.5  # dV/dt = −V/τ − q。空になる前に測る
    worst = {}
    for rounding in ("q48", "carry64"):
        net = Network(dt=dt, rounding=rounding)
        net.tank("a", 1.0, volume=v0)
        net.tank("sink", 1.0)
        net.valve("a", "sink", tau=tau)
        net.pump("a", "sink", rate=rate)
        start = net.total_q()
        for n in range(steps):
            net.step(n)
        assert net.total_q() == start
        with localcontext() as ctx:
            ctx.prec = 60
            k, q, t = 1 / Decimal(tau), Decimal(rate), Decimal(dt) * steps
            want = ((Decimal(v0) + q / k) * (-k * t).exp() - q / k) * ONE
            worst[rounding] = abs(net.nodes["a"]["v"] - want)
    assert worst["carry64"] <= 4 and worst["q48"] > worst["carry64"]


def test_carry64_follows_a_moving_valve():
    from network import Network

    net = Network(dt=0.5, rounding="carry64")
    net.tank("sea", 1.0, volume=64.0)
    net.tank("clock", 1.0)
    net.tank("a", 1.0, volume=8.0)
    net.tank("b", 1.0)
    net.pump("sea", "clock", rate=0.25)
    net.valve("a", "b", tau=16.0, sensor="clock", cam=[[0.0, 1.0], [16.0, 0.25]])
    start, seen = net.total_q(), set()
    for n in range(96):
        net.step(n)
        seen.add(net._carry_tank("a", net._by_src()["a"])[1])  # Φ は開き具合が変わるたびに出し直される
    assert net.total_q() == start and len(seen) > 32
