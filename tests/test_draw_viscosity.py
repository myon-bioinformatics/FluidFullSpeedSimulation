"""File rendering must retain viscosity controls and reproducible receipts."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import fluidfullspeedsimulation as ffs
import pixelwater as pw


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
            self.assertEqual(receipt["files"]["png"]["sha256"],
                "0957eed2eda4fb051cdeae8642a83ceaad24edcd81f4d86b729c4ee731769226")
