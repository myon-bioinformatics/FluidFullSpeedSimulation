"""🎞️ 断面図の GIF・PNG と、README に貼るデモ（場面ごとに GIF 1 枚）。"""
import contextlib
import hashlib
import io
import json
import re
import tempfile
import unittest
from pathlib import Path

import fluidfullspeedsimulation as ffs
import pixelwater as pw

ROOT = Path(__file__).resolve().parent.parent
DEMO_DIR = ROOT / "docs" / "demos"


def flat(rows) -> bytes:
    return b"".join(bytes(r) for r in rows)


class SectionImages(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res, cls.cols = pw.section_case("flume-hopper")
        cls.view = pw.network_frames(cls.res, cls.cols, frames=12, rows=16)

    def test_every_gif_frame_is_the_section_pixels(self):
        for scale, gap in ((8, 1), (4, 0), (16, 2)):
            with self.subTest(scale=scale, gap=gap):
                gif = pw.read_gif(pw.section_gif_bytes(self.view, scale=scale, gap=gap))
                self.assertEqual(len(gif["frames"]), len(self.view["frames"]))
                for i, frame in enumerate(gif["frames"]):
                    self.assertEqual(frame, flat(pw.section_pixels(self.view, i, scale=scale, gap=gap)), i)

    def test_png_is_the_last_time(self):
        png = pw.read_png(pw.section_png_bytes(self.view))
        self.assertEqual(flat(png["rows"]), flat(pw.section_pixels(self.view, -1)))

    def test_eight_colour_table_and_holds(self):
        data = pw.section_gif_bytes(self.view)
        self.assertEqual(data[10] & 7, 2)  # 全体の色表 2^(2+1) = 8 色（7 色 + 埋め草）
        self.assertEqual(bytes(data[13:16]), bytes(pw.section_palette("water")[0]))
        delays = [int.from_bytes(data[m.start() + 4:m.start() + 6], "little")
                  for m in re.finditer(rb"\x21\xf9\x04", data)]
        self.assertEqual(delays, pw._section_delays(len(self.view["frames"])))
        self.assertEqual((delays[0], delays[1], delays[-1]), pw.SECTION_DELAYS)

    def test_layout_is_tiles_of_eight(self):
        lay = pw.section_layout(self.view)
        self.assertEqual((lay["width"] % 8, lay["height"] % 8), (0, 0))
        # じょうご（4 マス）| 1 マス | 管 | 1 マス | 樋 20 区間 + 周りの余白 1 マス
        self.assertEqual(lay["width"], 8 * (1 + 4 + 1 + 1 + 1 + 20 + 1))
        self.assertEqual(lay["height"], 8 * (16 + 4))
        self.assertEqual([w for _x, w in lay["spans"]][:3], [32, 8, 8])
        with self.assertRaises(ValueError):
            pw.section_layout(self.view, scale=1)
        with self.assertRaises(ValueError):
            pw.section_pixels(self.view, scale=8, gap=8)

    def test_blocks_follow_counts_and_the_pipe_joins_its_tanks(self):
        lay = pw.section_layout(self.view)
        for i, fr in enumerate(self.view["frames"]):
            px = pw.section_pixels(self.view, i, layout=lay)
            for c, (x, _w) in enumerate(lay["spans"]):
                if self.view["kinds"][c] == "pipe":
                    continue
                colours = [px[8 + r * 8][x] for r in range(16)]  # 各段の左上
                n = fr["counts"][c]
                body = [k for k in colours if k in (pw.S_BODY, pw.S_SURFACE, pw.S_CARRIED)]
                self.assertEqual(len(body), n, (i, c))
                if n and fr["cells"][16 - n][c] == pw.CELL_FILL:
                    self.assertEqual(colours[16 - n], pw.S_SURFACE)  # いちばん上は水面の色
        px = pw.section_pixels(self.view, -1, layout=lay)
        floor = 8 + 15 * 8 + 6  # 底の段の下の方
        (hx, hw), (px_, _pw), (fx, _fw) = lay["spans"][:3]
        self.assertTrue(all(px[floor][x] in (pw.S_PIPE, pw.S_CARRIED) for x in range(hx + hw - 1, fx)))

    def test_time_bar_runs_from_empty_to_full(self):
        lay = pw.section_layout(self.view)
        y = lay["bar_y"] + 3
        first, last = pw.section_pixels(self.view, 0, layout=lay), pw.section_pixels(self.view, -1, layout=lay)
        self.assertNotIn(pw.S_BAR, first[y])
        self.assertEqual(last[y][8:lay["width"] - 9], bytes([pw.S_BAR]) * (lay["width"] - 17))

    def test_tracer_colour_differs_from_the_fill(self):
        self.assertEqual(pw.block_color("water", "🟩"), pw.IRONMATE_BASE["C01"])
        self.assertEqual(pw.block_color("water", "🟥"), pw.EXTRA_BLOCKS["🟥"])
        self.assertEqual(pw.block_color("ethanol", "🍣"), pw.IRONMATE_BASE["C09"])  # 好きな絵文字は流体の色
        for fluid, fill in (("water", None), ("ethanol", None), ("water", "🟨"), ("methanol", "🟧")):
            pal = pw.section_palette(fluid, fill)
            with self.subTest(fluid=fluid, fill=fill):
                self.assertEqual(len(set(pal)), len(pal))
                self.assertNotEqual(pal[pw.S_CARRIED], pal[pw.S_BODY])

    def test_noodle_shows_in_the_pipe(self):
        view = pw.network_frames(self.res, self.cols, frames=None, rows=8)
        piped = [i for i, fr in enumerate(view["frames"]) if any(r[1] == pw.CELL_IN_PIPE for r in fr["cells"])]
        self.assertTrue(piped)
        lay = pw.section_layout(view)
        x = lay["spans"][1][0]
        floor = 8 + 7 * 8 + 6
        self.assertEqual(pw.section_pixels(view, piped[0], layout=lay)[floor][x], pw.S_CARRIED)
        self.assertEqual(pw.section_pixels(view, piped[-1] + 1, layout=lay)[floor][x], pw.S_PIPE)


class DrawSection(unittest.TestCase):
    def test_files_numbering_and_receipt(self):
        res, cols = pw.section_case("flume")
        with tempfile.TemporaryDirectory() as tmp:
            a = pw.draw_section(res, cols, "flume", tmp, frames=8, rows=8, file_prefix="t-", file_suffix="_x")
            b = pw.draw_section(res, cols, "flume", tmp, frames=8, rows=8, file_prefix="t-", file_suffix="_x",
                                formats=("gif",))
            self.assertEqual(Path(a["paths"]["gif"]).name, "t-flume_water_n008_h08_r01_x.gif")
            self.assertEqual(Path(b["paths"]["gif"]).name, "t-flume_water_n008_h08_r02_x.gif")
            self.assertNotIn("png", b["paths"])
            receipt = json.loads(Path(a["paths"]["json"]).read_text(encoding="utf-8"))
            self.assertEqual(receipt["schema"], pw.SECTION_RECEIPT_SCHEMA)
            self.assertEqual(receipt["plan"]["rows"], 8)
            self.assertEqual(receipt["conservation_error_q"], 0)
            for kind in ("gif", "png"):
                data = Path(a["paths"][kind]).read_bytes()
                self.assertEqual(receipt["files"][kind]["sha256"], hashlib.sha256(data).hexdigest())
            self.assertEqual(Path(a["paths"]["gif"]).read_bytes(), Path(b["paths"]["gif"]).read_bytes())

    def test_cli_files_for_section_scenes(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = ffs.main(["draw", "🍷", "🟪", "--view", "files", "--frames", "6", "--rows", "8",
                                 "--out-dir", tmp, "--format", "gif"])
            self.assertEqual(code, 0)
            self.assertIn("🎞️", out.getvalue())
            names = sorted(p.name for p in Path(tmp).iterdir())
            self.assertEqual(names, ["pythagorean-cup_poison-water_n006_h08_r01.gif",
                                     "pythagorean-cup_poison-water_n006_h08_r01.json"])
            receipt = json.loads((Path(tmp) / names[1]).read_text(encoding="utf-8"))
            self.assertEqual(receipt["plan"]["fill"], "🟪")

    def test_frames_view_keeps_its_defaults(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(ffs.main(["draw", "flume"]), 0)
        self.assertEqual(len(out.getvalue().splitlines()), 6 * (1 + pw.SECTION_ROWS))


class Demos(unittest.TestCase):
    def test_one_demo_per_draw_scene(self):
        scenes = [scene for scene, _f, _fill in ffs.DEMOS]
        self.assertEqual(sorted(scenes), sorted(ffs.RIG_SCENES + pw.SECTION_CASES + pw.SCENES))
        self.assertEqual(ffs.demo_ids()[0], "D01")
        for text, want in (("D03", "D03"), ("d3", "D03"), ("🏛️", "D01"), ("flume-hopper", "D04")):
            self.assertEqual(ffs.demo_arg(text), want)

    def test_write_then_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = ffs.write_demos(tmp, only=["D01", "🍷"])
            self.assertTrue(a["ok"])
            self.assertEqual([d["id"] for d in a["demos"]], ["D01", "D13"])
            b = ffs.write_demos(tmp, only=["D03"])
            manifest = json.loads((Path(tmp) / "demos.json").read_text(encoding="utf-8"))
            self.assertEqual([d["id"] for d in manifest["demos"]], ["D01", "D03", "D13"])  # 差し替えて並べ直す
            ok = ffs.write_demos(tmp, check=True, only=["D01", "D03", "D13"])
            self.assertEqual([d["status"] for d in ok["demos"]], ["same"] * 3)
            path = Path(tmp) / b["demos"][0]["file"]
            path.write_bytes(path.read_bytes()[:-1])
            path13 = Path(tmp) / a["demos"][1]["file"]
            path13.unlink()
            bad = ffs.write_demos(tmp, check=True, only=["D01", "D03", "D13"])
            self.assertFalse(bad["ok"])
            self.assertEqual([d["status"] for d in bad["demos"]], ["same", "differs", "missing"])

    def test_cli_demos_check_exit_codes(self):
        with tempfile.TemporaryDirectory() as tmp:
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(ffs.main(["demos", "--out-dir", tmp, "--only", "D01"]), 0)
                self.assertEqual(ffs.main(["demos", "--out-dir", tmp, "--only", "D01", "--check"]), 0)
                self.assertEqual(ffs.main(["demos", "--out-dir", tmp, "--only", "D02", "--check"]), 1)
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as cm:
                ffs.main(["demos", "--only", "D99"])
            self.assertEqual(cm.exception.code, 2)

    def test_committed_demos_match_their_manifest_and_the_readme(self):
        """コミットした GIF と demos.json と README の 3 つが揃っているか（作り直しは CI の demos --check）。"""
        manifest = json.loads((DEMO_DIR / "demos.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema"], ffs.DEMOS_SCHEMA)
        self.assertEqual([d["id"] for d in manifest["demos"]], ffs.demo_ids())
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        for d, (scene, fluid, fill) in zip(manifest["demos"], ffs.DEMOS):
            with self.subTest(demo=d["id"]):
                self.assertEqual((d["scene"], d["fluid"], d["fill"]), (scene, fluid, fill))
                data = (DEMO_DIR / d["file"]).read_bytes()
                self.assertEqual(hashlib.sha256(data).hexdigest(), d["sha256"])
                gif = pw.read_gif(data)
                self.assertEqual(len(gif["frames"]), d["frames"])
                self.assertEqual(f"{gif['width']}x{gif['height']}", d["canvas"])
                self.assertLess(len(data), 128 * 1024)
                self.assertIn(f"docs/demos/{d['file']}", readme)
        self.assertEqual(sorted(p.name for p in DEMO_DIR.glob("*.gif")), sorted(d["file"] for d in manifest["demos"]))


if __name__ == "__main__":
    unittest.main()
