"""🔒 共通の検査（親の vendor/fleet_check.py）: vendor.lock.json と vendor/ のバイト列、親 workflow の固定 commit。"""

import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_PIN = "98a0e611cc879fde2d047da806de5a2121a65524"  # Workflows only; vendored files keep the commits in ACCEPTED.
ACCEPTED = {
    "vendor/cli_args.py": "e5d4bf4f39d6fc3d5f7890a7c8e9cffb723c1ece",
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
