"""絵文字をそのままコマンドに書く（run 🥢 / draw 🏛️ 🟩 / draw 🥢 🟥）と、vendor の固定ファイル。"""

import contextlib
import hashlib
import io
import json
import unittest
from pathlib import Path

import fluidfullspeedsimulation as ffs

ROOT = Path(__file__).resolve().parents[1]


class EmojiArguments(unittest.TestCase):
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
