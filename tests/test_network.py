import json
import math
import unittest

from fixedq import from_q, qmul, random_scales, to_q
from network import LAWS, Network, build, critical_gain, jsonable, register_law, settling_time, simulate


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


if __name__ == "__main__":
    unittest.main()
