"""🍜 長さ・たまり・箸の仕組み（network の麺）を小さな回路で確かめる。E12〜E14 の通しは test_experiments が回す。"""

import math
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import random_scales
from network import EJECT, PACKING, Network, crest, pile_height

BAMBOO = {"shape": "semicircle", "size": 0.04}


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
        res, cols = ffs.rig_section("noodle-length")
        view = ffs.pw.network_frames(res, cols, frames=None)
        widest = max(sum(line.count("🍜") for line in fr["lines"]) for fr in view["frames"])
        self.assertGreaterEqual(widest, 3)


if __name__ == "__main__":
    unittest.main()
