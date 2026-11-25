import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch, Mock

import get_smiles_InChI_IUPAC_props as lookup


class LookupTests(unittest.TestCase):
    def test_cas_lookup_applies_timeout_and_preserves_legacy_smiles(self):
        replies = [
            Mock(json=lambda: {"IdentifierList": {"CID": [702]}}),
            Mock(json=lambda: {"PropertyTable": {"Properties": [{"ConnectivitySMILES": "CCO", "SMILES": "CCO"}]}}),
            Mock(json=lambda: {}),
            Mock(json=lambda: {}),
        ]
        with patch.object(lookup, "pubchem_get", side_effect=replies) as get:
            result = lookup.get_compound_info("64-17-5")
        self.assertEqual(result["IsomericSMILES"], "CCO")
        self.assertEqual(result["CanonicalSMILES"], "CCO")
        self.assertEqual(get.call_count, 4)
        self.assertTrue(all(call.kwargs["timeout"] == 30 for call in get.call_args_list))

    def test_missing_executable_is_reported_before_importing_optional_dependencies(self):
        script = Path(__file__).resolve().parents[1] / "HSPiP_CLI_v7.py"
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "smiles.txt"
            source.write_text("CCO\n")
            env = {**os.environ, "HSPIP_PATH": str(Path(tmp) / "missing")}
            proc = subprocess.run([sys.executable, str(script), str(source)], env=env,
                                  capture_output=True, text=True, timeout=5)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("HSPIP_PATH", proc.stderr)

    def test_invalid_workdir_is_reported_after_validating_executable(self):
        script = Path(__file__).resolve().parents[1] / "HSPiP_CLI_v7.py"
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            exe = folder / "HSPiP.exe"
            exe.touch()
            exe.chmod(0o755)
            source = folder / "smiles.txt"
            source.write_text("CCO\n")
            env = {**os.environ, "HSPIP_PATH": str(exe), "HSPIP_WORKDIR": str(folder / "missing")}
            proc = subprocess.run([sys.executable, str(script), str(source)], env=env,
                                  capture_output=True, text=True, timeout=5)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("HSPIP_WORKDIR", proc.stderr)
