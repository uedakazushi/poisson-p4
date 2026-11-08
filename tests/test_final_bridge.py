"""Smoke and mutation tests for the independent final-bridge verifier."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'code' / 'verify_final_bridge.py'


class FinalBridgeTests(unittest.TestCase):
    def test_exact_family(self):
        run = subprocess.run([sys.executable, str(SCRIPT)], text=True,
                             capture_output=True, timeout=60, check=True)
        report = json.loads(run.stdout)
        self.assertEqual(report['status'], 'PASS')
        self.assertEqual(report['check_count'], 81)
        self.assertIn('Ga9_potential_full_six', report['checks'])

    def test_wrong_sign_is_detected(self):
        source = SCRIPT.read_text(encoding='utf-8')
        old = '(-1+2*A)*x*y+y/6-r/3'
        self.assertEqual(source.count(old), 1)
        with tempfile.TemporaryDirectory() as tmp:
            mutant = Path(tmp) / 'mutant.py'
            mutant.write_text(source.replace(old, '(-1+2*A)*x*y+y/6+r/3'), encoding='utf-8')
            run = subprocess.run([sys.executable, str(mutant)], text=True,
                                 capture_output=True, timeout=60)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('six_brackets_from_product', run.stderr)


if __name__ == '__main__':
    unittest.main()
