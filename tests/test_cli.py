"""⌨️ コマンド（list・fluids・run・draw）: 終了コード、絵文字の引数、draw の出力、場面ごとの run / draw、vendor の固定ファイル。

- CliCase: ``call``（標準出力ごと）と ``usage_error``（argparse の終了コード 2）を 1 か所に
- Commands: list・fluids・run、--json・--jsonl、終了コード 1 / 2
- EmojiArguments: 絵文字をそのままコマンドに書く（run 🥢 / draw 🏛️ 🟩 / draw 🥢 🟥）
- Draw: draw のファイル書き出し・コマ・仮置きの流体・試験台、断る引数
- 場面ごと: 🍯 粘度・🧊🎋 水温と竹の傷み・🌧️ 雨・🎍🍷🫗 溜まる → 一気に・🚧 慣性
- Vendor: vendor.lock.json と vendor/ のバイト列、cli_args は vendor から読む
"""

import contextlib
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import fluidfullspeedsimulation as ffs

ROOT = Path(__file__).resolve().parents[1]


class CliCase(unittest.TestCase):
    def call(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(list(argv))
        return code, out.getvalue()

    def usage_error(self, *argv):
        with contextlib.redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit) as cm:
            ffs.main(list(argv))
        self.assertEqual(cm.exception.code, 2)
        return err.getvalue()


class Commands(CliCase):
    def test_list_fluids_run(self):
        code, text = self.call("list")
        self.assertEqual(code, 0)
        self.assertIn("E08", text)
        self.assertIn("E10", text)
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
        self.usage_error("run", "E99")
        self.usage_error()

    def test_failed_check_exits_1(self):
        original = ffs.EXPERIMENTS["E01"]
        ffs.EXPERIMENTS["E01"] = original[:3] + (lambda: {"checks": [ffs.check("x", "❌", False)], "runs": {}},)
        try:
            code, text = self.call("run", "E01")
        finally:
            ffs.EXPERIMENTS["E01"] = original
        self.assertEqual(code, 1)
        self.assertIn("❌", text)


class EmojiArguments(CliCase):
    def test_experiment_icons(self):
        for e in ffs.list_experiments():
            self.assertEqual(ffs.resolve(e["icon"]), e["id"])
        self.assertEqual(ffs.resolve("🏛"), "E11")  # 異体字セレクタなしでも同じ
        code, text = self.call("run", "🧮")
        self.assertEqual(code, 0)
        self.assertTrue(text.startswith("✅ E01"))

    def test_fluid_blocks_and_icons(self):
        aliases = ffs.fluid_aliases()
        self.assertEqual(aliases["🟩"], "namekian-water")
        self.assertEqual(aliases["☠"], "poison-water")
        self.assertEqual(aliases["🍶"], "ethanol")
        code, text = self.call("draw", "🏛️", "🟩", "--frames", "1")
        self.assertEqual(code, 0)
        self.assertIn("🟩", text)
        code, text = self.call("draw", "flume", "--fluid", "🧴", "--frames", "1")
        self.assertEqual(code, 0)
        self.assertIn("🟧", text)

    def test_scene_icons_and_temperature(self):
        code, text = self.call("draw", "🧊", "--frames", "1", "--temp-c", "5")
        self.assertEqual(code, 0)
        code, text = self.call("draw", "🪵", "🟪", "--frames", "1")
        self.assertEqual(code, 0)
        self.assertIn("🟪", text)
        for icon, scene in ffs.SCENE_ICONS.items():
            self.assertIn(scene, ffs.RIG_SCENES + ffs.pw.SECTION_CASES + ffs.pw.SCENES, icon)

    def test_any_emoji_paints(self):
        code, text = self.call("draw", "🥢", "🟥", "--frames", "1")
        self.assertEqual(code, 0)
        self.assertIn("🟥", text)
        self.assertNotIn("🟦", text)
        code, text = self.call("draw", "🫙", "🍵", "--frames", "1")
        self.assertEqual(code, 0)
        self.assertIn("🍵", text)

    def test_rejections(self):
        self.usage_error("draw", "🥢", "hoge")
        self.usage_error("draw", "🍣")
        self.usage_error("run", "🍣")
        self.assertIn("流体が 2 つ", self.usage_error("draw", "🥢", "🟩", "--fluid", "ethanol"))
        self.assertIn("色が 2 つ", self.usage_error("draw", "🥢", "🟥", "--fill", "🟫"))
        self.usage_error("draw", "🥢", "--fram", "2")  # 省略形は受け付けない（cli_args の方針）


class Draw(CliCase):
    def test_draw_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, text = self.call("draw", "moniac-drain", "--view", "files", "--fluid", "ethanol", "--ticks", "30",
                                   "--scale", "4", "--out-dir", tmp, "--file-prefix", "cli-", "--format", "gif")
            self.assertEqual(code, 0)
            self.assertIn("🎞️", text)
            self.assertIn("🧾", text)
            names = sorted(p.name for p in Path(tmp).iterdir())
            self.assertEqual(names, ["cli-moniac-drain_ethanol_s0_t0030_r01.gif",
                                     "cli-moniac-drain_ethanol_s0_t0030_r01.json"])

    def test_draw_rejects_unknown_scene(self):
        for argv in (["draw", "waterfall"], ["draw", "flume", "--view", "movie"],
                     ["draw", "flume", "--cap-mode", "median"], ["draw", "flume", "--rows", "x"]):
            with self.subTest(argv=argv):
                self.usage_error(*argv)

    def test_draw_frames(self):
        code, text = self.call("draw", "flume", "--frames", "3", "--rows", "8", "--fill", "🟩")
        self.assertEqual(code, 0)
        lines = text.splitlines()
        self.assertEqual(len(lines), 3 * (1 + 8))  # （⏱️ + 8 段）× 3 時刻
        self.assertIn("🟩", text)
        self.assertNotIn("🟦", text)
        code, text = self.call("draw", "dam-break", "--frames", "2", "--fill", "🟧")
        self.assertEqual(code, 0)
        self.assertIn("🟧", text)
        self.assertEqual(sum(line.startswith("⏱️") for line in text.splitlines()), 2)

    def test_placeholder_fluids(self):
        code, text = self.call("draw", "dam-break", "--fluid", "namekian-water", "--frames", "1")
        self.assertEqual(code, 0)
        self.assertIn("🟩", text)
        code, text = self.call("draw", "flume", "--fluid", "poison-water", "--frames", "2")
        self.assertEqual(code, 0)
        self.assertIn("🟪", text)
        code, text = self.call("fluids")
        self.assertIn("🏷️", text)
        self.assertNotIn("namekian-water", ffs.FLUID_KEYS)

    def test_draw_rigs(self):
        code, text = self.call("draw", "moniac-rig", "--frames", "2")
        self.assertEqual(code, 0)
        rows = [line for line in text.splitlines() if not line.startswith("⏱️")]
        self.assertEqual(len(rows), 2 * 8)
        self.assertTrue(all(len(line) == 4 for line in rows))  # 頭・貯める・戻りの管・横
        code, text = self.call("draw", "somen-flume", "--frames", "2", "--cap-mode", "rim")
        self.assertEqual(code, 0)
        self.assertIn("🍜", text)


class ViscosityFlags(CliCase):
    def test_draw_with_viscosity(self):
        self.assertEqual(self.call("draw", "🫙", "--viscosity-scale", "8", "--frames", "1")[0], 0)
        self.assertEqual(self.call("draw", "🥢", "--mu", "5", "--frames", "1")[0], 0)
        self.assertEqual(self.call("draw", "🌊", "🟪", "--viscosity-scale", "2", "--frames", "1")[0], 0)
        for argv in (["draw", "🥢", "--mu", "1", "--viscosity-scale", "2"], ["draw", "🥢", "--mu", "0"],
                     ["draw", "🥢", "--viscosity-scale", "nan"]):
            with self.subTest(argv=argv):
                self.usage_error(*argv)


class TemperatureWearScenes(CliCase):
    def test_draw_scenes(self):
        for scene in ("cold-water", "bamboo-wear"):
            code, text = self.call("draw", scene, "--frames", "2")
            self.assertEqual(code, 0, scene)
            self.assertIn("🍜", text, scene)


class RainScene(CliCase):
    def test_draw_rain_scene(self):
        code, _ = self.call("draw", "🌧️", "--frames", "2", "--rows", "4")
        self.assertEqual(code, 0)
        self.assertEqual(ffs.SCENE_ICONS["🌧"], "rain-runoff")
        self.assertIn("rain-runoff", ffs.RIG_SCENES)
        self.assertEqual(ffs.resolve("🌧️"), "E18")


class PeriodicScenes(CliCase):
    def test_run_and_draw_by_icon(self):
        for icon, eid in (("🎍", "E19"), ("🍷", "E20"), ("🫗", "E21")):
            self.assertEqual(ffs.resolve(icon), eid)
        for icon in ("🎍", "🍷", "🫗"):
            code, text = self.call("draw", icon, "--frames", "3", "--rows", "4")
            self.assertEqual(code, 0, icon)
            rows = [line for line in text.splitlines() if not line.startswith("⏱️")]
            self.assertEqual(len(rows), 3 * 4)
            self.assertTrue(all(len(line) == 2 for line in rows))
        code, text = self.call("draw", "🫗", "🟩", "--frames", "2")
        self.assertEqual(code, 0)
        self.assertIn("🟩", text)


class InertiaLimitsScene(CliCase):
    def test_run_and_draw_by_icon(self):
        code, text = self.call("run", "🚧")
        self.assertEqual(code, 0)
        self.assertTrue(text.startswith("✅ E23"))
        code, text = self.call("draw", "🚧", "--frames", "2", "--rows", "4")
        self.assertEqual(code, 0)
        self.assertIn("🟦", text)
        self.assertEqual(sum(line.startswith("⏱️") for line in text.splitlines()), 2)


class Vendor(unittest.TestCase):
    def test_vendored_files_match_the_lock(self):
        lock = json.loads((ROOT / "vendor.lock.json").read_text(encoding="utf-8"))
        self.assertEqual(lock["schema"], "vendor-lock/1")
        for item in lock["files"]:
            data = (ROOT / item["destination"]).read_bytes()
            self.assertEqual(hashlib.sha256(data).hexdigest(), item["sha256"], item["destination"])
            blob = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
            self.assertEqual(blob, item["blob_sha"], item["destination"])
            self.assertEqual(item["ref"], item["commit"])
            self.assertRegex(item["commit"], r"^[0-9a-f]{40}$")

    def test_cli_args_comes_from_vendor(self):
        cli = ffs.vendored("cli_args")
        self.assertEqual(Path(cli.__file__).resolve(), (ROOT / "vendor" / "cli_args.py").resolve())
        self.assertEqual((cli.EXIT_OK, cli.EXIT_FAILURE, cli.EXIT_USAGE), (0, 1, 2))


if __name__ == "__main__":
    unittest.main()
