"""🟦 描画（pixelwater.py）: ドット絵の水、断面図、PNG・GIF と受領記録、README のデモ GIF。

- pixelwater: ビット盤の水（落下・横移動）・PNG / GIF・ファイル書き出し・断面図
- draw_viscosity: File rendering must retain viscosity controls and reproducible receipts.
- demos: 🎞️ 断面図の GIF・PNG と、README に貼るデモ（場面ごとに GIF 1 枚）。
"""

import contextlib
import hashlib
import io
import json
import random
import re
import tempfile
import unittest
from pathlib import Path

import fluidfullspeedsimulation as ffs
import pixelwater as pw
from fixedq import ONE, random_scales


# ──── pixelwater ──────────────────────────────────────────────────

def column(ticks: int, height: int = 10) -> dict:
    """幅 3 の筒のてっぺんに水 1 粒。"""
    lines = ["#~#"] + ["#.#"] * (height - 2) + ["###"]
    return pw.run_pixels("custom", "water", ticks=ticks, frame_every=1, map_lines=lines)


class Maps(unittest.TestCase):
    def test_scene_maps_parse(self):
        for name in pw.SCENES:
            g = pw.parse_map(pw.scene_map(name))
            self.assertGreater(g.count(), 0)
            self.assertEqual(g.overlap(), 0)
        with self.assertRaises(KeyError):
            pw.scene_map("waterfall")

    def test_bad_maps(self):
        for bad in ([], ["#.#", "##"], ["#x#"]):
            with self.assertRaises(ValueError):
                pw.parse_map(bad)


class Motion(unittest.TestCase):
    def test_falls_one_row_per_tick(self):
        res = column(12)
        rows = [next(y for y, r in enumerate(water) if r) for water, _f in res["frames"]]
        self.assertEqual(rows[:9], list(range(9)))  # 1 ティックで 1 行
        self.assertEqual(rows[-1], 8)  # 床の上で止まる
        self.assertEqual(res["first_floor_tick"], 8)

    def test_single_pixel_bounces_between_walls(self):
        lines = ["#......#", "#~.....#", "########"]
        res = pw.run_pixels("custom", "methanol", ticks=80, frame_every=1, map_lines=lines, seed=None)
        xs = {(water[1] & -water[1]).bit_length() - 1 for water, _f in res["frames"][2:]}
        self.assertEqual(xs, set(range(1, 7)))  # 壁の間を全部たどる
        self.assertEqual(res["count_error"], 0)

    def test_count_conserved_for_random_rates_and_seeds(self):
        rates = random_scales(51, 4, n_min=1, n_max=1)
        seeds = random.Random(52).sample(range(10 ** 6), 4)
        for item, seed in zip(rates, seeds):
            for scene in pw.SCENES:
                with self.subTest(rate=item["label"], seed=seed, scene=scene):
                    res = pw.run_pixels(scene, "water", ticks=150, seed=seed, base_rate=item["value"])
                    self.assertEqual(res["count_error"], 0)
                    self.assertEqual(res["lateral_passes"], 150 * pw.lateral_rate("water", item["value"]) // ONE)

    def test_viscosity_orders_settling(self):
        for seed in (0, 7):
            for scene in pw.SCENES:
                t = {k: pw.run_pixels(scene, k, seed=seed)["settle_tick"] for k in ("methanol", "water", "ethanol")}
                with self.subTest(seed=seed, scene=scene, ticks=t):
                    self.assertTrue(t["methanol"] < t["water"] < t["ethanol"])

    def test_deterministic(self):
        a, b = pw.run_pixels("dam-break", seed=None), pw.run_pixels("dam-break", seed=None)
        self.assertEqual(a["hash"], b["hash"])
        self.assertNotEqual(a["hash"], pw.run_pixels("dam-break", seed=3)["hash"])


class Images(unittest.TestCase):
    def test_palette_follows_ironmate_rule(self):
        pal = pw.palette("water")
        self.assertEqual(pal[2][:3], (0x63, 0xB3, 0xFF))  # Ironmate C02 main
        self.assertEqual(pal[3][:3], (161, 209, 255))  # highlight = round((3·main + 2·255)/5)
        art = pw.to_ironmate_artifact(pw.run_pixels("dam-break", ticks=5))
        self.assertEqual(art["schema"], "ironmate-pixel-grid/1")
        self.assertEqual(art["colors"][0], "#00000000")
        self.assertEqual(len(art["pixels"]), art["height"])

    def test_png_round_trip_and_crc(self):
        res = pw.run_pixels("moniac-drain", "ethanol", ticks=30)
        for scale in (1, 3, 8):
            png = pw.png_bytes(res, -1, scale=scale)
            got = pw.read_png(png)
            self.assertEqual((got["width"], got["height"]), (res["width"] * scale, res["height"] * scale))
            idx = pw.indices(res, -1)
            self.assertTrue(all(got["rows"][y * scale][x * scale] == idx[y][x]
                                for y in range(res["height"]) for x in range(res["width"])))
        broken = bytearray(png)
        broken[40] ^= 0xFF
        with self.assertRaises(ValueError):
            pw.read_png(bytes(broken))
        with self.assertRaises(ValueError):
            pw.png_bytes(res, scale=0)

    def test_gif_every_frame_round_trips(self):
        res = pw.run_pixels("dam-break", "methanol", ticks=60, frame_every=3)
        for scale in (1, 2, 8):
            gif = pw.read_gif(pw.gif_bytes(res, scale=scale))
            self.assertEqual(len(gif["frames"]), len(res["frames"]))
            for i, frame in enumerate(gif["frames"]):
                idx = pw.indices(res, i)
                self.assertTrue(all(frame[(y * scale) * gif["width"] + x * scale] == idx[y][x]
                                    for y in range(res["height"]) for x in range(res["width"])), (scale, i))

    def test_lzw_round_trip_through_table_resets(self):
        rng = random.Random(9)
        data = bytes(rng.randrange(4) for _ in range(30000))
        self.assertEqual(pw._unlzw(pw._lzw(data, 2), 2, len(data)), data)
        flat = bytes(40000)
        self.assertEqual(pw._unlzw(pw._lzw(flat, 2), 2, len(flat)), flat)
        with self.assertRaises(ValueError):
            pw.read_gif(b"GIF87a" + bytes(10))


class Files(unittest.TestCase):
    def test_semantic_numbering_and_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pw.draw("dam-break", "water", tmp, ticks=20, file_prefix="demo-", file_suffix="_x")
            b = pw.draw("dam-break", "water", tmp, ticks=20, file_prefix="demo-", file_suffix="_x")
            self.assertEqual(Path(a["paths"]["gif"]).name, "demo-dam-break_water_s0_t0020_r01_x.gif")
            self.assertEqual(Path(b["paths"]["gif"]).name, "demo-dam-break_water_s0_t0020_r02_x.gif")
            receipt = json.loads(Path(a["paths"]["json"]).read_text(encoding="utf-8"))
            self.assertEqual(receipt["schema"], pw.RECEIPT_SCHEMA)
            for kind in ("gif", "png"):
                data = Path(a["paths"][kind]).read_bytes()
                self.assertEqual(receipt["files"][kind]["sha256"], hashlib.sha256(data).hexdigest())
            again = pw.draw(out_dir=tmp, **{k: v for k, v in receipt["plan"].items() if k != "ticks"},
                            ticks=receipt["plan"]["ticks"])
            self.assertEqual(Path(again["paths"]["gif"]).read_bytes(), Path(a["paths"]["gif"]).read_bytes())

    def test_summary_is_json_safe(self):
        out = pw.simulate_pixels("moniac-drain", "ethanol", ticks=40)
        json.dumps(out, allow_nan=False)
        self.assertIn("🟨", out["final"])


class Sections(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.flume, cls.cols = pw.section_case("flume", "water")
        cls.hopper, cls.hcols = pw.section_case("flume-hopper", "water")

    def test_blocks_round_and_read_back(self):
        self.assertEqual([pw.blocks(h, 1.0, 8) for h in (0, 0.06, 0.07, 0.5, 1.0, 2.0, -1)], [0, 0, 1, 4, 8, 8, 0])
        self.assertEqual(pw.blocks(1.0, 0.0, 8), 0)
        for item in random_scales(71, 8, n_min=1, n_max=3):
            h, cap = item["value"], 1.0
            for rows in (8, 16, 24):
                n = pw.blocks(h, cap, rows)
                self.assertLessEqual(abs(pw.readback(n, cap, rows) - h), cap / (2 * rows) + 1e-15)

    def test_frames_stack_from_the_bottom(self):
        view = pw.network_frames(self.hopper, self.hcols, frames=None, rows=8)
        for fr, cnt in zip(view["frames"], view["counts"]):
            col = [line[0] for line in fr["lines"]]
            water = [c != pw.EMPTY_BLOCK for c in col]
            self.assertEqual(water, sorted(water))  # 上は空、下に積む（間に空きがない）
            self.assertTrue(cnt[0] <= sum(water) <= max(cnt[0], 1))  # 🍜 は水面の段（空なら底）に置く
        self.assertTrue(all(len(line) == len(view["columns"]) for fr in view["frames"] for line in fr["lines"]))

    def test_front_moves_downstream_for_random_rows(self):
        for item in random_scales(72, 3, n_min=0, n_max=0, digits=1):
            rows = 8 * int(item["A"])
            with self.subTest(rows=rows):
                view = pw.network_frames(self.flume, self.cols, rows=rows)
                arr = view["arrivals"]
                self.assertNotIn(None, arr)
                self.assertTrue(all(a < b for a, b in zip(arr, arr[1:])))
                self.assertLessEqual(view["readback_error"], view["half_block"] + 1e-15)

    def test_caps_and_fill(self):
        rim = pw.network_frames(self.flume, self.cols, cap_mode="rim", frames=1)
        self.assertEqual(rim["caps"][0], 0.005)
        fixed = pw.network_frames(self.flume, self.cols, caps={"flume": 0.001}, frames=1, fill="🟨")
        self.assertEqual(fixed["caps"][0], 0.001)
        self.assertIn("🟨", "".join(fixed["frames"][0]["lines"]))
        with self.assertRaises(ValueError):
            pw.network_frames(self.flume, self.cols, cap_mode="median")
        with self.assertRaises(ValueError):
            pw.network_frames(self.flume, ["nowhere"])
        with self.assertRaises(ValueError):
            pw.network_frames(self.flume, self.cols, rows=0)

    def test_tracer_moves_downstream_and_pipe_shows_delay(self):
        view = pw.network_frames(self.flume, self.cols, frames=None)
        cols = [c[0] for c in view["tracer_columns"] if c[0] is not None]
        self.assertEqual(cols, sorted(cols))
        hv = pw.network_frames(self.hopper, self.hcols, frames=None)
        span = [(a, b) for a, b, place in hv["tracer_spans"][0] if place == "hopper>flume:0"]
        self.assertEqual(len(span), 1)
        self.assertAlmostEqual(span[0][1] - span[0][0], 3.0)
        fallback = pw.network_frames(self.hopper, ["hopper", "flume"], frames=None)
        self.assertTrue(any(pw.IN_PIPE in "".join(fr["lines"]) for fr in fallback["frames"]))

    def test_pixel_frames_and_arrivals(self):
        res = pw.run_pixels("dam-break", "water", frame_every=1)
        view = pw.pixel_frames(res, frames=3, fill="🟩")
        self.assertEqual(len(view["frames"]), 3)
        self.assertIn("🟩", view["frames"][0]["lines"][-2])
        arr = view["arrivals"]
        self.assertEqual(len(arr), 30)
        self.assertNotIn(None, arr)


# ──── draw_viscosity ──────────────────────────────────────────────

"""File rendering must retain viscosity controls and reproducible receipts."""


class DrawViscosity(unittest.TestCase):
    def test_cli_files_viscosity_and_receipt_replay(self):
        hashes = []
        for flags, knobs in (([], {}), (["--viscosity-scale", "100"], {"mu_scale": 100.0}),
                             (["--mu", "100"], {"mu": 0.1})):
            with self.subTest(flags=flags), tempfile.TemporaryDirectory() as tmp:
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(ffs.main(["draw", "dam-break", "--view", "files", "--ticks", "60",
                                               "--out-dir", tmp, "--file-prefix", "check_",
                                               "--file-suffix", "_visc", *flags]), 0)
                receipt_path = next(Path(tmp).glob("*.json"))
                self.assertEqual(receipt_path.name, "check_dam-break_water_s0_t0060_r01_visc.json")
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                res = pw.run_pixels("dam-break", "water", ticks=60, **knobs)
                self.assertEqual(receipt["hash"], res["hash"])
                for key, value in knobs.items():
                    self.assertEqual(receipt["plan"][key], value)
                expected = {"gif": pw.gif_bytes(res, scale=8), "png": pw.png_bytes(res, -1, scale=8)}
                replay = pw.draw(out_dir=tmp, **receipt["plan"])
                for kind, blob in expected.items():
                    actual = (Path(tmp) / receipt["files"][kind]["name"]).read_bytes()
                    self.assertEqual(actual, blob)
                    self.assertEqual(Path(replay["paths"][kind]).read_bytes(), blob)
                    self.assertEqual(receipt["files"][kind]["sha256"], hashlib.sha256(blob).hexdigest())
                hashes.append(receipt["hash"])
        self.assertEqual(len(hashes), 3)
        self.assertNotEqual(hashes[0], hashes[1])
        self.assertNotEqual(hashes[0], hashes[2])

    def test_default_bytes_and_receipt_plan_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = pw.draw(ticks=60, out_dir=tmp)
            receipt = json.loads(Path(result["paths"]["json"]).read_text(encoding="utf-8"))
            self.assertEqual(receipt["plan"], {"scene": "dam-break", "fluid": "water", "ticks": 60,
                "seed": 0, "scale": 8, "frame_every": 2, "base_rate": 0.5})
            self.assertEqual(receipt["files"]["gif"]["sha256"],
                "6e7c35b98e1ee9f58fac531ebe62a5496d9d88b7b98e135cf4913301bc0835b4")
            # zlib streams may differ across platforms; pin decoded pixels instead.
            # The test above also checks exact PNG bytes and receipt replay locally.
            png = pw.read_png(Path(result["paths"]["png"]).read_bytes())
            self.assertEqual((png["width"], png["height"]), (256, 192))
            self.assertEqual(hashlib.sha256(b"".join(png["rows"])).hexdigest(),
                "1e93b498f0b7b4aae7e0170ee39a876d3883c138805bb3e1ba5e9bfe7eb8871b")


# ──── demos ───────────────────────────────────────────────────────

"""🎞️ 断面図の GIF・PNG と、README に貼るデモ（場面ごとに GIF 1 枚）。"""


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
        self.assertEqual((delays[0], delays[1], delays[-1]), (pw.SECTION_DELAYS[0], 240 // 12, pw.SECTION_DELAYS[2]))
        self.assertEqual(tuple(pw._section_delays(40)[i] for i in (0, 1, -1)), pw.SECTION_DELAYS)

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


class KymoImages(unittest.TestCase):
    """🕰️ 時空図: 横 = 断面図と同じ列、縦 = 時刻（1 時刻 = scale/2 px）、色 = 水深の段、🍜 の通り道は gold。"""

    @classmethod
    def setUpClass(cls):
        cls.res, cls.cols = pw.section_case("flume")
        cls.view = pw.network_frames(cls.res, cls.cols, frames=24, rows=16)

    def test_every_gif_frame_is_the_kymo_pixels(self):
        for scale, gap in ((8, 0), (4, 1), (16, 2)):
            with self.subTest(scale=scale, gap=gap):
                gif = pw.read_gif(pw.kymo_gif_bytes(self.view, scale=scale, gap=gap))
                self.assertEqual(len(gif["frames"]), len(self.view["frames"]))
                for i, frame in enumerate(gif["frames"]):
                    self.assertEqual(frame, flat(pw.kymo_pixels(self.view, i, scale=scale, gap=gap)), i)
        png = pw.read_png(pw.kymo_png_bytes(self.view))
        self.assertEqual(flat(png["rows"]), flat(pw.kymo_pixels(self.view, -1)))

    def test_rows_are_times_and_colours_are_depths(self):
        lay = pw.kymo_layout(self.view)
        self.assertEqual((lay["width"] % 8, lay["height"] % 8), (0, 0))
        self.assertEqual(lay["height"], 8 * 2 + 24 * 4)
        self.assertEqual(len(pw.kymo_palette("water", rows=16)), 7 + 16)
        px = pw.kymo_pixels(self.view, -1, layout=lay)
        half = pw.kymo_pixels(self.view, 11, layout=lay)
        gold = 0
        for i, fr in enumerate(self.view["frames"]):
            y = 8 + i * lay["row_h"]
            for c, (x, w) in enumerate(lay["spans"]):
                n = fr["counts"][c]
                with self.subTest(t=i, column=c):
                    self.assertEqual(px[y][x], pw.S_BAR + n if n else pw.S_EMPTY)  # 端は水深の色
                    cells = {fr["cells"][r][c] for r in range(16)}
                    carried = bool(cells & {pw.CELL_CARRIED, pw.CELL_IN_PIPE})  # 管の列が無いと ⏳ は行き先の底
                    self.assertEqual(px[y][x + w // 2] == pw.S_CARRIED, carried)  # 真ん中は 🍜 ⏳ があれば gold
                    self.assertEqual(half[y][x], px[y][x] if i <= 11 else pw.S_BG)  # まだの時刻は背景
                gold += carried
        self.assertGreater(gold, 0)
        arrivals = self.view["arrivals"]
        self.assertEqual(arrivals, sorted(arrivals))  # 水の先端は下流ほど遅く届く（時空図の斜めの境目）
        self.assertLess(arrivals[0], arrivals[-1])

    def test_cli_kymo_files_and_pixel_scenes(self):
        with tempfile.TemporaryDirectory() as tmp:
            with contextlib.redirect_stdout(io.StringIO()):
                code = ffs.main(["draw", "🎋", "--view", "kymo", "--frames", "8", "--rows", "8", "--out-dir", tmp])
            self.assertEqual(code, 0)
            names = sorted(p.name for p in Path(tmp).iterdir())
            self.assertEqual(names, [f"flume_water_kymo_n008_h08_r01.{e}" for e in ("gif", "json", "png")])
            receipt = json.loads((Path(tmp) / names[1]).read_text(encoding="utf-8"))
            self.assertEqual((receipt["plan"]["view"], receipt["canvas"]), ("kymo", "176x48"))
        with contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(ffs.main(["draw", "🌊", "--view", "kymo"]), 2)  # ドット絵の場面には列がない
        self.assertIn("🕰️", err.getvalue())


class Demos(unittest.TestCase):
    def test_one_demo_per_draw_scene(self):
        scenes = [scene for scene, _f, _fill, *how in ffs.DEMOS if not how]
        self.assertEqual(sorted(scenes), sorted(ffs.RIG_SCENES + pw.SECTION_CASES + pw.SCENES))
        kymo = [scene for scene, _f, _fill, *how in ffs.DEMOS if how]
        self.assertTrue(kymo and set(kymo) <= set(ffs.RIG_SCENES + pw.SECTION_CASES))  # 🕰️ は断面図の場面の見せ方
        self.assertEqual({tuple(d[3:]) for d in ffs.DEMOS if len(d) > 3}, {("kymo",)})
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
        self.assertEqual([d["id"] for d in manifest["demos"]], ffs.demo_ids() + ffs.figure_ids())
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        figures = manifest["demos"][len(ffs.DEMOS):]
        for d, (name, icon, _build, _args) in zip(figures, ffs.FIGURES):  # 📊 SVG の図
            with self.subTest(figure=d["id"]):
                data = (DEMO_DIR / d["file"]).read_bytes()
                self.assertEqual((d["scene"], d["icon"], d["kind"]), (name, icon, "figure"))
                self.assertEqual(hashlib.sha256(data).hexdigest(), d["sha256"])
                self.assertIn(f"docs/demos/{d['file']}", readme)
        self.assertEqual(sorted(p.name for p in DEMO_DIR.glob("*.svg")), sorted(d["file"] for d in figures))
        for d, (scene, fluid, fill, *how) in zip(manifest["demos"], ffs.DEMOS):
            with self.subTest(demo=d["id"]):
                self.assertEqual((d["scene"], d["fluid"], d["fill"], d["kind"] == "kymo"), (scene, fluid, fill, bool(how)))
                data = (DEMO_DIR / d["file"]).read_bytes()
                self.assertEqual(hashlib.sha256(data).hexdigest(), d["sha256"])
                gif = pw.read_gif(data)
                self.assertEqual(len(gif["frames"]), d["frames"])
                self.assertEqual(f"{gif['width']}x{gif['height']}", d["canvas"])
                self.assertLess(len(data), 128 * 1024)
                self.assertIn(f"docs/demos/{d['file']}", readme)
        self.assertEqual(sorted(p.name for p in DEMO_DIR.glob("*.gif")),
                         sorted(d["file"] for d in manifest["demos"][:len(ffs.DEMOS)]))


if __name__ == "__main__":
    unittest.main()
