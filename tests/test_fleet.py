"""🔒 共通の検査（親の vendor/fleet_check.py）: vendor.lock.json と vendor/ のバイト列、親 workflow の固定 commit。"""

import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_PIN = "8a930eae6987a142851a94fa59cdd3bd740d16a3"  # Workflows only; vendored files keep the commits in ACCEPTED.
ACCEPTED = {
    "vendor/cli_args.py": "9ff06c4326fc7ba934112abecbc59451ab07be6a",
    "vendor/py_launch.cjs": "2e450118fde453272de9997608e25c31a984463c",
    "vendor/fleet_check.py": "64db25eb0fd9025dbe3cbcbced8b88a39ecb9266",
}


class Fleet(unittest.TestCase):
    def test_vendored_files_match_the_lock(self):
        rows = fleet.check_lock(ROOT, commits=ACCEPTED)
        self.assertEqual({w for w in rows if not w.endswith("-LICENSE")}, set(ACCEPTED))

    def test_parent_workflows_use_the_reviewed_commit(self):
        self.assertEqual(set(fleet.parent_refs(ROOT)), {PARENT_PIN})


if __name__ == "__main__":
    unittest.main()


# 🧰 Scripts: 📝 each one describes itself and ⌨️ one that reads sys.argv declares its arguments (argparse /
# cli_args). Read with ast, never run; the portfolio's script catalog reads the same facts.

class ScriptsDescribeThemselves(unittest.TestCase):
    def test_every_script_has_a_description_and_declared_arguments(self):
        fleet.check_scripts(fleet.repo_root(__file__))
