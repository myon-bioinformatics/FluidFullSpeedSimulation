"""🔒 共通の検査（親の vendor/fleet_check.py）: vendor.lock.json と vendor/ のバイト列、親 workflow の固定 commit。"""

import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_PIN = "43cfe1df924ea96ca44f47a6771c18102c38f7f3"  # Workflows only; vendored files keep the commits in ACCEPTED.
LEAN_PARENT_PIN = "538293c3a4cb3b84ae22ff5a57a79a4abfec715d"  # New formal lane; existing workflow pins stay unchanged.
ACCEPTED = {
    "vendor/cli_args.py": "5ae22b1a1b808eca93e3c385a8db867fbeb8ee32",
    "vendor/py_launch.cjs": "de63762d8485e4d398a651715d34ebf4a04c23f2",
    "vendor/fleet_check.py": "20de01921e50e9339aa9271b0547a5be7649b3b3",
}


class Fleet(unittest.TestCase):
    def test_vendored_files_match_the_lock(self):
        rows = fleet.check_lock(ROOT, commits=ACCEPTED)
        self.assertEqual({w for w in rows if not w.endswith("-LICENSE")}, set(ACCEPTED))

    def test_parent_workflows_use_the_reviewed_commit(self):
        refs = fleet.parent_refs(ROOT)
        self.assertEqual(set(refs), {PARENT_PIN, LEAN_PARENT_PIN})
        # Only this new reusable workflow may use the Lean verifier revision.
        # A changed pin anywhere else, including a checkout, still fails.
        self.assertEqual(refs.count(LEAN_PARENT_PIN), 1)
        self.assertEqual(
            fleet.check_actions(ROOT, workflow="lean-verification.yml"),
            [("myon-bioinformatics/myon-bioinformatics/.github/workflows/reusable-lean-verify.yml",
              LEAN_PARENT_PIN, None)],
        )


if __name__ == "__main__":
    unittest.main()


# 🧰 Scripts: 📝 each one describes itself and ⌨️ one that reads sys.argv declares its arguments (argparse /
# cli_args). Read with ast, never run; the portfolio's script catalog reads the same facts.

class ScriptsDescribeThemselves(unittest.TestCase):
    def test_every_script_has_a_description_and_declared_arguments(self):
        fleet.check_scripts(fleet.repo_root(__file__))

    def test_pull_request_ci_stays_quick(self):
        # ⏱ push limited to branches, superseded PR runs cancelled, every job timed (fleet_check.check_ci)
        fleet.check_ci(fleet.repo_root(__file__))

    def test_root_scripts_import_only_the_standard_library(self):
        root = fleet.repo_root(__file__)
        for path in sorted(root.glob("*.py")):
            if path.name != "conftest.py":
                fleet.check_imports(root, path.name)

