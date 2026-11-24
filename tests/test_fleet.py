"""🔒 共通の検査（親の vendor/fleet_check.py）: vendor.lock.json と vendor/ のバイト列、親 workflow の固定 commit。"""

import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_PIN = "1368174821fb03be258ecf3312ce0ea809f6ea22"  # Workflows only; vendored files keep the commits in ACCEPTED.
ACCEPTED = {
    "vendor/cli_args.py": "5ae22b1a1b808eca93e3c385a8db867fbeb8ee32",
    "vendor/py_launch.cjs": "de63762d8485e4d398a651715d34ebf4a04c23f2",
    "vendor/fleet_check.py": "bc57f08a590ee1cb33a8bb37c1bb3b8f28a2f3c9",
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

    def test_root_scripts_import_only_the_standard_library(self):
        root = fleet.repo_root(__file__)
        for path in sorted(root.glob("*.py")):
            if path.name != "conftest.py":
                fleet.check_imports(root, path.name)
