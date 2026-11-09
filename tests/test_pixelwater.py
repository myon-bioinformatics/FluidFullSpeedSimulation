import contextlib
import hashlib
import io
import json
import random
import tempfile
import unittest
from pathlib import Path

import fluidfullspeedsimulation as ffs
import pixelwater as pw
from fixedq import ONE, random_scales


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
        self.assertIn("🟪", out["final"])


class Cli(unittest.TestCase):
    def call(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = ffs.main(list(argv))
        return code, out.getvalue()

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
        for argv in (["draw", "waterfall"], ["draw", "flume", "--view", "files"],
                     ["draw", "flume", "--cap-mode", "median"], ["draw", "flume", "--rows", "x"]):
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as cm:
                ffs.main(argv)
            self.assertEqual(cm.exception.code, 2)

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


if __name__ == "__main__":
    unittest.main()
