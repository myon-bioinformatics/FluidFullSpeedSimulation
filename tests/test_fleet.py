"""🔒 共通の検査（親の vendor/fleet_check.py）: vendor.lock.json と vendor/ のバイト列、親 workflow の固定 commit。"""

import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_PIN = "2e450118fde453272de9997608e25c31a984463c"
ACCEPTED = {
    "vendor/cli_args.py": "9ff06c4326fc7ba934112abecbc59451ab07be6a",
    "vendor/py_launch.cjs": PARENT_PIN,
    "vendor/fleet_check.py": "6c29503cf17eda856766ddfdf4d66a5b0965a8b1",
}


class Fleet(unittest.TestCase):
    def test_vendored_files_match_the_lock(self):
        rows = fleet.check_lock(ROOT, commits=ACCEPTED)
        self.assertEqual({w for w in rows if not w.endswith("-LICENSE")}, set(ACCEPTED))

    def test_parent_workflows_use_the_reviewed_commit(self):
        self.assertEqual(set(fleet.parent_refs(ROOT)), {PARENT_PIN})


if __name__ == "__main__":
    unittest.main()
