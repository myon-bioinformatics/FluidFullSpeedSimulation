import contextlib
import io
import unittest

import fluidfullspeedsimulation as ffs
from fixedq import random_scales
from fluids import FLUIDS, TEMP_TABLES, get_fluid, table_value, temperature_range
from network import Network, build, simulate

BAMBOO = {"shape": "semicircle", "size": 0.04}


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


class Draw(unittest.TestCase):
    def test_draw_scenes(self):
        for scene in ("cold-water", "bamboo-wear"):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = ffs.main(["draw", scene, "--frames", "2"])
            self.assertEqual(code, 0, scene)
            self.assertIn("🍜", out.getvalue(), scene)


if __name__ == "__main__":
    unittest.main()
