import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import fluidfullspeedsimulation as ffs
from fixedq import random_scales


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


class Cli(unittest.TestCase):
    def call(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(list(argv))
        return code, out.getvalue()

    def test_list_fluids_run(self):
        code, text = self.call("list")
        self.assertEqual(code, 0)
        self.assertIn("E08", text)
        code, text = self.call("fluids")
        self.assertEqual(code, 0)
        self.assertIn("ethanol", text)
        code, text = self.call("run", "tank-integrator", "-v")
        self.assertEqual(code, 0)
        self.assertTrue(text.startswith("✅ E01"))

    def test_json_and_jsonl(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "plot.jsonl"
            code, text = self.call("run", "E03", "--json", "--jsonl", str(path))
            self.assertEqual(code, 0)
            self.assertEqual(json.loads(text)["id"], "E03")
            rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
            self.assertTrue(rows)
            self.assertEqual({r["id"] for r in rows}, {"E03"})

    def test_bad_choice_exits_2(self):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as cm:
            ffs.main(["run", "E99"])
        self.assertEqual(cm.exception.code, 2)
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as cm:
            ffs.main([])
        self.assertEqual(cm.exception.code, 2)

    def test_failed_check_exits_1(self):
        original = ffs.EXPERIMENTS["E01"]
        ffs.EXPERIMENTS["E01"] = original[:3] + (lambda: {"checks": [ffs.check("x", "❌", False)], "runs": {}},)
        try:
            code, text = self.call("run", "E01")
        finally:
            ffs.EXPERIMENTS["E01"] = original
        self.assertEqual(code, 1)
        self.assertIn("❌", text)


if __name__ == "__main__":
    unittest.main()
