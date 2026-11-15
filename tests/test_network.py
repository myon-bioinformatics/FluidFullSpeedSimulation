"""🌊 水力回路のエンジン（network.py）: 総量の保存・組み立て・樋・運ばれる物、竹の傷み、麺、単弁の厳密解、丸め方 carry64。

- network: 総量の保存・組み立て（build / spec）・樋・運ばれる物・測る道具（critical_gain ほか）
- bamboo_wear: 🎋 竹の傷み: 節の削り残しで竹が分かれる・割れ目からの漏れも総量に入る・おかしな傷みは断る。
- noodles: 🍜 長さ・たまり・箸の仕組み（network の麺）を小さな回路で確かめる。E12〜E14 の通しは test_experiments が回す。
- single_valve_review: Independent analytic oracles; keep the default rounding policy.
"""

import json
import math
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import from_q, qmul, random_scales, to_q
from network import (
    EJECT, LAWS, Network, PACKING, build, crest, critical_gain, jsonable, pile_height, register_law, settling_time,
    simulate,
)


# ──── network ─────────────────────────────────────────────────────

def mixed(dt: float, rate: float, fluid: str = "water") -> Network:
    """管の法則を一通り混ぜた回路（総量の保存を確かめる用）。"""
    net = Network(dt, fluid).tank("sea", 1.0, volume=1.0).tank("a", 0.01, head=0.2).tank("b", 0.02).tank("c", 0.01)
    net.pump("sea", "a", rate=rate)
    net.pipe("a", "b", radius=1e-3, length=0.5, delay=3 * dt)
    net.orifice("a", "c", area=1e-6)
    net.weir("b", "c", shape="v")
    net.cam("c", "sea", table=[[0, 0], [0.1, 1e-4]], delay=dt)
    net.float_valve("sea", "b", sensor="b", setpoint=0.05, gain=1e-3, qmax=1e-4, stages=2, tau=1.0)
    net.channel("gutter", "sea", length=1.0, cells=5, slope=0.05, section={"shape": "semicircle", "size": 0.02})
    net.pump("c", "gutter:0", rate=rate / 2)
    return net


class Conservation(unittest.TestCase):
    def test_total_is_exact_for_random_dt_and_rates(self):
        dts = random_scales(21, 4, n_min=1, n_max=2)
        rates = random_scales(22, 4, n_min=5, n_max=6)
        for d, r in zip(dts, rates):
            with self.subTest(dt=d["label"], rate=r["label"]):
                net = mixed(d["value"], r["value"])
                res = net.run(200 * d["value"])
                self.assertEqual(res["conservation_error_q"], 0)
                self.assertTrue(all(v >= 0 for v in res["volumes_q"].values()))

    def test_starved_tank_is_shared_without_going_negative(self):
        net = Network(1.0).tank("a", 1.0, volume=1e-9).tank("b", 1.0).tank("c", 1.0)
        net.pump("a", "b", rate=1.0)
        net.pump("a", "c", rate=3.0)
        res = net.run(3.0)
        self.assertEqual(res["volumes_q"]["a"], 0)
        self.assertEqual(res["conservation_error_q"], 0)
        self.assertLess(res["volumes_q"]["b"], res["volumes_q"]["c"])


class Building(unittest.TestCase):
    SPEC = {
        "dt": 0.1, "fluid": "ethanol",
        "tanks": [{"name": "A", "area": 0.01, "head": 0.05}, {"name": "B", "area": 0.01}],
        "links": [
            {"law": "pipe", "src": "A", "dst": "B", "radius": 2e-3, "length": 0.5, "delay": 5.0},
            {"law": "pipe", "src": "B", "dst": "A", "radius": 2e-3, "length": 1.0, "delay": 10.0},
        ],
        "t_end": 300.0,
    }

    def test_spec_matches_hand_built_bit_for_bit(self):
        net = Network(0.1, "ethanol").tank("A", 0.01, head=0.05).tank("B", 0.01)
        net.pipe("A", "B", radius=2e-3, length=0.5, delay=5.0)
        net.pipe("B", "A", radius=2e-3, length=1.0, delay=10.0)
        hand = net.run(300.0)
        spec = build(self.SPEC).run(300.0)
        self.assertEqual(hand["volumes_q"], spec["volumes_q"])

    def test_simulate_is_json_safe(self):
        back = json.loads(json.dumps(simulate(self.SPEC)))
        hand = build(self.SPEC).run(300.0)["volumes_q"]
        self.assertEqual({k: int(v) for k, v in back["volumes_q"].items()}, hand)
        self.assertEqual(jsonable(2 ** 60), str(2 ** 60))  # 2^53 超は文字列（JS の Number で桁落ちしない）
        self.assertEqual(jsonable(2 ** 53), 2 ** 53)
        self.assertEqual(jsonable([math.inf, -math.inf]), ["inf", "-inf"])

    def test_spec_errors(self):
        with self.assertRaises(ValueError):
            build({"tanks": [], "bogus": 1})
        with self.assertRaises(ValueError):
            simulate({"tanks": []})
        with self.assertRaises(ValueError):
            build({"tanks": [{"name": "a", "area": 1}], "links": [{"law": "nope", "src": "a", "dst": "a"}]})

    def test_registered_law_is_usable_from_spec(self):
        @register_law("half")
        def _half(net, ln, n):
            return net.nodes[ln["src"]]["v"] // 2

        try:
            res = simulate({"dt": 1.0, "tanks": [{"name": "a", "area": 1, "volume": 1}, {"name": "b", "area": 1}],
                            "links": [{"law": "half", "src": "a", "dst": "b"}], "t_end": 3.0})
            self.assertAlmostEqual(from_q(int(res["volumes_q"]["a"])), 1 / 8)
            self.assertEqual(res["conservation_error_q"], 0)
        finally:
            LAWS.pop("half")

    def test_bad_wiring(self):
        net = Network(1.0).tank("a", 1.0)
        with self.assertRaises(ValueError):
            net.tank("a", 1.0)
        with self.assertRaises(ValueError):
            net.pump("a", "missing", rate=1.0)
        with self.assertRaises(ValueError):
            net.pump("a", "a", rate=1.0)
        with self.assertRaises(ValueError):
            Network(0)


class Flumes(unittest.TestCase):
    BAMBOO = {"shape": "semicircle", "size": 0.04}

    def spec(self):
        return {
            "dt": 0.01,
            "tanks": [{"name": "sea", "area": 1, "volume": 1}, {"name": "floor", "area": 1}],
            "flumes": [{"prefix": "take", "dst": "sea", "spill": "floor", "segments": [
                {"length": 1.0, "cells": 5, "slope": 0.05, "section": self.BAMBOO, "joint": {"drop": 0.02}},
                {"length": 1.0, "cells": 5, "slope": 0.03, "section": self.BAMBOO, "bend_radius": 0.2,
                 "joint": {"lip": 0.002}},
                {"length": 0.5, "cells": 5, "slope": 0.05, "section": self.BAMBOO},
            ]}],
            "links": [{"law": "pump", "src": "sea", "dst": "take/0:0", "rate": 5e-5}],
            "t_end": 10.0,
        }

    def test_spec_flume_matches_hand_built(self):
        spec = self.spec()
        net = Network(0.01).tank("sea", 1, volume=1).tank("floor", 1)
        net.flume("take", "sea", spill="floor", segments=spec["flumes"][0]["segments"])
        net.pump("sea", "take/0:0", rate=5e-5)
        self.assertEqual(net.run(10.0)["volumes_q"], build(spec).run(10.0)["volumes_q"])
        res = simulate(spec)
        self.assertEqual(res["conservation_error_q"], 0)
        self.assertEqual(set(res["channels"]), {"take/0", "take/1", "take/2"})

    def test_flume_errors(self):
        net = Network(0.01).tank("sea", 1, volume=1)
        with self.assertRaises(ValueError):
            net.flume("f", "sea", [])
        with self.assertRaises(ValueError):
            net.flume("f", "sea", [{"length": 1, "cells": 2, "slope": 0.05, "section": self.BAMBOO,
                                    "joint": {"gap": 0.01}}])
        with self.assertRaises(ValueError):
            net.channel("g", "sea", length=1, cells=2, slope=0.05, section=self.BAMBOO, bend_radius=0.2)

    def test_overtopping_goes_to_spill_node(self):
        """縁の容量を超えて入れた水は、そのまま spill の節点へ落ちる（総量は保存）。"""
        net = Network(0.01).tank("sea", 1, volume=1).tank("floor", 1)
        net.channel("c", "sea", length=0.1, cells=1, slope=0.05, section=self.BAMBOO, lip=0.039, spill="floor")
        net.pump("sea", "c:0", rate=1e-3)
        res = net.run(2.0)
        self.assertGreater(res["volumes_q"]["floor"], 0)
        self.assertEqual(res["conservation_error_q"], 0)

    def test_thick_tracer_waits_at_lip(self):
        net = Network(0.05).tank("sea", 1, volume=1)
        net.channel("c", "sea", length=0.5, cells=5, slope=0.05, section=self.BAMBOO, law="laminar", beta=3.0,
                    lip=0.003)
        net.pump("sea", "c:0", rate=3e-6)
        net.tracer("thin", "c:0", release=20.0, stop="sea", size=1e-3)
        net.tracer("thick", "c:0", release=20.0, stop="sea", size=2e-3)
        res = net.run(60.0)
        tr = {t["name"]: t for t in res["tracers"]}
        self.assertEqual(tr["thin"]["state"], "done")
        self.assertEqual(tr["thick"]["state"], "node")
        self.assertGreater(tr["thick"]["stuck_time"], 0)
        self.assertEqual(tr["thick"]["stuck_at"], "c:4")


class Tracers(unittest.TestCase):
    def test_tank_chain_residence_and_pipe_delay(self):
        """一定流量の水槽 2 つ＋遅延管: 通過時間 = V1/Q + 遅れ + V2/Q（先入れ先出し）。

        t2 は管が満ちるまでの遅れの間だけ出す一方なので、定常の V2 = 3e-4 − Q·遅れ。
        """
        q, dt = 1e-5, 0.05
        net = Network(dt).tank("sea", 1.0, volume=1.0).tank("t1", 0.01, volume=2e-4).tank("t2", 0.01, volume=3e-4)
        net.pump("sea", "t1", rate=q)
        net.pump("t1", "t2", rate=q, delay=4.0)
        net.pump("t2", "sea", rate=q)
        net.tracer("leaf", "t1", release=1.0, stop="sea")
        res = net.run(200.0)
        tr = res["tracers"][0]
        self.assertEqual(tr["state"], "done")
        self.assertEqual([p["node"] for p in tr["path"]], ["t1", "t2", "sea"])
        took = tr["arrival"] - tr["path"][0]["t"]
        expected = 2e-4 / q + 4.0 + (3e-4 - q * 4.0) / q
        self.assertLessEqual(abs(took - expected), 2 * dt)

    def test_unknown_tracer_node(self):
        net = Network(1.0).tank("a", 1.0)
        with self.assertRaises(ValueError):
            net.tracer("x", "b")
        with self.assertRaises(ValueError):
            net.tracer("x", "a", stop="b")


class Measures(unittest.TestCase):
    def test_critical_gain(self):
        self.assertEqual(critical_gain(1.0, 1.0, 1.0, 1), math.inf)
        # 3 段の等しい遅れだけ（leak → 大）なら古典値 K_c = 8（位相 −π で各段 60°）
        kc = critical_gain(1e-12, 1.0, 1.0, 3)
        self.assertAlmostEqual(kc, 8.0, places=6)

    def test_settling_time(self):
        self.assertEqual(settling_time([0.0, 0.5, 0.99, 1.0, 1.0], 1.0, 1.0), 2.0)
        self.assertEqual(settling_time([1.0, 1.0], 1.0, 1.0), 0.0)
        self.assertEqual(settling_time([1.0, 0.0], 1.0, 1.0), math.inf)

    def test_pump_start_and_until(self):
        net = Network(0.5).tank("a", 1.0, volume=1.0).tank("b", 1.0)
        net.pump("a", "b", rate=0.1, until=2.0)
        net.pump("a", "b", rate=0.1, start=9.0)
        res = net.run(10.0)
        self.assertEqual(res["volumes_q"]["b"], (4 + 2) * qmul(to_q(0.1), to_q(0.5)))


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


class Rounding(unittest.TestCase):
    def test_exp_neg_q_matches_decimal(self):
        from decimal import Decimal, localcontext
        from fractions import Fraction

        from fixedq import exp_neg_q

        for x in (Fraction(1, 512), Fraction(5, 96), Fraction(3, 4), Fraction(7, 2), Fraction(40), Fraction(0)):
            with localcontext() as ctx:
                ctx.prec = 60
                want = (-(Decimal(x.numerator) / Decimal(x.denominator))).exp() * (1 << 64)
                assert abs(exp_neg_q(x, 64) - want) <= Decimal("0.5000001")

    def test_carry64_keeps_total_and_stays_within_a_few_ulps(self):
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

    def test_carry64_only_touches_valve_only_tanks_and_default_is_q48(self):
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

    def test_carry64_with_a_pump_out_matches_the_exact_answer(self):
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

    def test_carry64_follows_a_moving_valve(self):
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

    def test_single_partly_open_valve_is_exact_by_default(self):
        """🚿 開き具合が 1 でない弁 1 本の槽: 既定（槽ごとに厳密）は V0·e^(−開き具合·t/τ) と合う。法則のまま（independent）は刻みの 1 次ずれる。"""
        import math

        from fixedq import ONE
        from network import Network

        got = {}
        for mode in ("exact", "independent"):
            net = Network(dt=4.0, outlets=mode)
            net.tank("a", 1.0, volume=128.0)
            net.tank("sink", 1.0)
            net.valve("a", "sink", tau=32.0, opening=0.5)
            start = net.total_q()
            for n in range(16):  # t = 2τ
                net.step(n)
            assert net.total_q() == start
            got[mode] = net.nodes["a"]["v"] / ONE / (128.0 * math.exp(-1.0)) - 1.0
        assert abs(got["exact"]) < 1e-12
        assert 0.03 < got["independent"] < 0.033  # 開き具合·(1 − 開き具合)·t·dt/(2τ²) = 1/32 のあたり


# ──── bamboo_wear ─────────────────────────────────────────────────

"""🎋 竹の傷み: 節の削り残しで竹が分かれる・割れ目からの漏れも総量に入る・おかしな傷みは断る。"""


BAMBOO = {"shape": "semicircle", "size": 0.04}


class BambooWear(unittest.TestCase):
    def flume(self, **seg):
        net = Network(0.01).tank("sea", 1, volume=1).tank("floor", 1)
        base = {"length": 2.0, "cells": 20, "slope": 0.05, "section": BAMBOO}
        return net, net.flume("take", "sea", [{**base, **seg}], spill="floor")

    def test_node_remnant_splits_the_bamboo(self):
        net, infos = self.flume(node_remnant={"at": 0.3, "height": 0.002}, cracks=[{"cell": 2, "area": 1e-6},
                                                                                  {"cell": 15, "area": 1e-6}])
        a, b = infos[0], infos[0]["after_node"]
        self.assertEqual((len(a["cells"]), len(b["cells"])), (6, 14))
        self.assertAlmostEqual(a["length"] + b["length"], 2.0)
        self.assertEqual(a["lip"], 0.002)
        self.assertEqual([c["cell"] for c in a["cracks"]], [2])
        self.assertEqual([c["cell"] for c in b["cracks"]], [9])
        cracks = [ln for ln in net.links if ln.get("crack")]
        self.assertEqual([ln["src"] for ln in cracks], ["take/0b:9", "take/0:2"])

    def test_crack_leak_is_conserved_for_random_areas(self):
        for item in random_scales(92, 3, n_min=6, n_max=7):
            with self.subTest(area=item["label"]):
                net, _ = self.flume(cracks=[{"cell": 5, "area": item["value"]}])
                net.pump("sea", "take/0:0", rate=5e-5)
                res = net.run(10.0)
                crack = next(ln for ln in res["links"] if ln["kind"] == "orifice")
                self.assertEqual(res["conservation_error_q"], 0)
                self.assertEqual(res["volumes_q"]["floor"], crack["volume_moved_q"])
                self.assertGreater(crack["volume_moved_q"], 0)

    def test_bad_wear(self):
        for seg in ({"node_remnant": {"at": 0.0, "height": 0.002}}, {"node_remnant": {"at": 0.5, "height": 0}},
                    {"node_remnant": {"at": 0.5, "height": 0.002, "gap": 1}}, {"cracks": [{"cell": 20, "area": 1e-6}]},
                    {"cracks": [{"cell": 1, "area": 1e-6, "shape": "slit"}]}):
            with self.subTest(seg=seg), self.assertRaises(ValueError):
                self.flume(**seg)
        net = Network(0.01).tank("sea", 1, volume=1)
        with self.assertRaises(ValueError):
            net.channel("c", "sea", length=1, cells=2, slope=0.05, section=BAMBOO, cracks=[{"cell": 0, "area": 1e-6}])


# ──── noodles ─────────────────────────────────────────────────────

"""🍜 長さ・たまり・箸の仕組み（network の麺）を小さな回路で確かめる。E12〜E14 の通しは test_experiments が回す。"""


def lip_chain(q: float = 3e-6, dt: float = 0.05, cells: int = 5) -> Network:
    net = Network(dt).tank("sea", 1.0, volume=1.0)
    net.channel("c", "sea", length=0.5, cells=cells, slope=0.05, section=BAMBOO, law="laminar", beta=3.0,
                lip=0.003)
    net.pump("sea", "c:0", rate=q)
    return net


class Bodies(unittest.TestCase):
    def test_body_length_in_cells(self):
        for item in random_scales(91, 4, n_min=1, n_max=1):
            length = item["value"] * 4  # 0.04〜0.4 m（区間 0.1 m の何個分か）
            with self.subTest(length=item["label"]):
                net = lip_chain(q=5e-5, dt=0.01)
                tr = net.tracer("n", "c:0", release=5.0, stop="sea", length=length)
                net.run(5.4)
                body = net.body(tr)
                self.assertEqual(body[0], tr["node"])
                self.assertLessEqual(len(body), math.ceil(length / 0.1 - 1e-9))
                self.assertGreaterEqual(sum(net.nodes[c]["length"] for c in body), min(length, 0.1 * len(body)) - 1e-9)

    def test_bad_noodle(self):
        net = lip_chain()
        with self.assertRaises(ValueError):
            net.tracer("x", "c:0", length=-1.0)
        with self.assertRaises(ValueError):
            net.tracer("x", "c:0", strands=0)


class Piles(unittest.TestCase):
    def test_pile_raises_the_crest(self):
        net = lip_chain()
        a = net.tracer("a", "c:0", release=20.0, size=1.5e-3, length=0.2, strands=10, friction=10.0)
        b = net.tracer("b", "c:0", release=25.0, size=1.5e-3, length=0.2, strands=10, friction=10.0)
        net.run(60.0)
        lip = next(ln for ln in net.links if ln["kind"] == "lip")
        self.assertTrue(a["held"] and b["held"])
        want = 2 * 10 * math.pi * (1.5e-3) ** 2 / 4.0 / (PACKING * lip["width"])
        self.assertAlmostEqual(pile_height(net, lip), want, places=12)
        self.assertAlmostEqual(crest(net, lip), 0.003 + want, places=12)
        self.assertEqual(lip["releases"], [])

    def test_floating_noodles_never_pile(self):
        """水より軽い麺は支える力が 0 以下なので、引っかかってもすぐ押し流される。"""
        net = lip_chain()
        tr = net.tracer("light", "c:0", release=20.0, stop="sea", size=1.5e-3, length=0.2, density=900.0)
        res = net.run(60.0)
        self.assertEqual(tr["state"], "done")
        lip = res["lips"][0]
        self.assertEqual([r["members"] for r in lip["releases"]], [["light"]])
        self.assertEqual(res["conservation_error_q"], 0)

    def test_point_noodles_keep_the_old_rule(self):
        """長さ 0 の麺は積もらない（縁は上がらず、越流水深 < 太さなら待つだけ）。"""
        net = lip_chain()
        tr = net.tracer("p", "c:0", release=20.0, stop="sea", size=1.5e-3)
        net.run(60.0)
        lip = next(ln for ln in net.links if ln["kind"] == "lip")
        self.assertTrue(tr["held"])
        self.assertEqual(pile_height(net, lip), 0.0)

    def test_release_count_matches_prediction_for_random_friction(self):
        for item in random_scales(92, 2, n_min=1, n_max=1):
            mu = 0.1 + 0.2 * (item["A"] - 1.0) / 9.0  # 0.1〜0.3 の中で乱数（A を写す）
            pred = ffs.clog_prediction(3e-6, 10, ffs.STRAND["size"], 0.2, ffs.STRAND["density"], mu)
            with self.subTest(friction=round(mu, 4), k_star=pred["k_star"]):
                r = ffs.run_experiment("E13", friction=mu, bunches=pred["k_star"] + 2)
                ok = {c["id"]: c["ok"] for c in r["checks"]}
                self.assertTrue(ok["release-count"], ffs.render(r))
                self.assertTrue(ok["release-force"])
                self.assertTrue(ok["conserve"])


class Eject(unittest.TestCase):
    def test_rule(self):
        ln = {"radius": 0.12, "rim": 0.04, "last_outer": 0.039}
        self.assertTrue(EJECT["spill"](None, ln, {"size": 1.5e-3}))
        self.assertFalse(EJECT["spill"](None, ln, {"size": 0.5e-3}))
        self.assertFalse(EJECT["spill"](None, dict(ln, radius=None), {"size": 1.5e-3}))

    def test_eject_window_is_inside_the_rim(self):
        q = ffs.eject_flow(1.5e-3)
        self.assertAlmostEqual(ffs.bend_outer(q), 0.04 - 0.75e-3, places=9)


class Catchers(unittest.TestCase):
    def test_fixed_grab_time(self):
        """t_grab を固定すると、区間にいる時間がそれより長い麺は全部取れ、短ければ全部取り逃す。"""
        for t_grab, caught in ((0.1, 5), (5.0, 0)):
            net = lip_chain(q=5e-5, dt=0.01)
            for i in range(5):
                net.tracer(f"s{i}", "c:0", release=5.0 + i, stop="sea")
            net.catcher("hand", "c:2", t_grab=t_grab)
            res = net.run(15.0)
            c = res["catchers"][0]
            self.assertEqual(len(c["caught"]), caught)
            self.assertEqual(len(c["caught"]) + len(c["missed"]), 5)
            self.assertEqual(res["conservation_error_q"], 0)
            self.assertTrue(all(t["state"] in ("caught", "done") for t in res["tracers"]))

    def test_bad_catcher(self):
        net = lip_chain()
        with self.assertRaises(ValueError):
            net.catcher("x", "nowhere")
        with self.assertRaises(ValueError):
            net.catcher("x", "c:0", grab=(0.5, 0.1))


class Draw(unittest.TestCase):
    def test_scenes(self):
        for scene in ("noodle-length", "chopsticks"):
            with self.subTest(scene=scene):
                res, cols = ffs.rig_section(scene)
                self.assertEqual(res["conservation_error_q"], 0)
                self.assertTrue(all(row["t"] >= 20.0 for row in res["series"]))
        self.assertEqual(ffs.scene_arg("🍽️"), "chopsticks")
        self.assertEqual(ffs.scene_arg("🧱"), "noodle-clog")
        res, cols = ffs.rig_section("noodle-length")
        view = ffs.pw.network_frames(res, cols, frames=None)
        widest = max(sum(line.count("🍜") for line in fr["lines"]) for fr in view["frames"])
        self.assertGreaterEqual(widest, 3)


# ──── single_valve_review ─────────────────────────────────────────

"""Independent analytic oracles; keep the default rounding policy."""


class SingleValveReview(unittest.TestCase):
    def circuit(self, dt, a, count=1, incoming=0, outgoing=0, mode='exact', inflows='uniform', rounding='q48'):
        net = Network(dt, outlets=mode, inflows=inflows, rounding=rounding)
        net.tank('a', 1, volume=8).tank('sink', 1).tank('feed', 1, volume=128)
        for _ in range(count):
            net.valve('a', 'sink', tau=1, opening=a/count)
        if incoming:
            net.pump('feed', 'a', rate=incoming)
        if outgoing:
            net.pump('a', 'sink', rate=outgoing)
        return net

    def test_analytic_matrix(self):
        for a in (0, .25, .5, .75, 1):
            for dt in (1/8, 1/64):
                for count in (1, 2):
                    for incoming, outgoing in ((0, 0), (.25, 0), (0, .125), (.25, .125)):
                        with self.subTest(a=a, dt=dt, count=count, incoming=incoming, outgoing=outgoing):
                            net = self.circuit(dt, a, count, incoming, outgoing)
                            self.assertEqual(net.rounding, 'q48')
                            r = net.run(2)
                            want = 8+(incoming-outgoing)*2 if not a else 8*math.exp(-2*a)+(incoming-outgoing)*(-math.expm1(-2*a))/a
                            self.assertAlmostEqual(from_q(r['volumes_q']['a']), want, delta=2e-10)
                            self.assertEqual(r['conservation_error_q'], 0)
                            self.assertTrue(all(v >= 0 for v in r['volumes_q'].values()))

    def test_legacy_is_not_euler(self):
        dt, a, n = 1/8, .25, 16
        r = self.circuit(dt, a, mode='independent').run(n*dt)
        want = 8*(1-a*(-math.expm1(-dt)))**n
        self.assertAlmostEqual(from_q(r['volumes_q']['a']), want, delta=2e-10)
        self.assertGreater(want, 8*math.exp(-a*n*dt))
        self.assertLess(8*(1-a*dt)**n, 8*math.exp(-a*n*dt))

    def test_end_inflow(self):
        dt, a, incoming, n = 1/8, .25, .25, 16
        r = self.circuit(dt, a, incoming=incoming, inflows='end').run(n*dt)
        decay = math.exp(-a*dt)
        want = 8*decay**n+incoming*dt*(1-decay**n)/(1-decay)
        self.assertAlmostEqual(from_q(r['volumes_q']['a']), want, delta=2e-10)
        self.assertEqual(r['conservation_error_q'], 0)

    def test_rounding_modes(self):
        for rounding in ('q48', 'carry64'):
            r = self.circuit(1/8, .25, rounding=rounding).run(2)
            self.assertAlmostEqual(from_q(r['volumes_q']['a']), 8*math.exp(-.5), delta=2e-10)
            self.assertEqual(r['conservation_error_q'], 0)

    def test_starvation(self):
        r = self.circuit(1/8, .5, outgoing=100).run(2)
        self.assertEqual(r['conservation_error_q'], 0)
        self.assertTrue(all(v >= 0 for v in r['volumes_q'].values()))


if __name__ == "__main__":
    unittest.main()
