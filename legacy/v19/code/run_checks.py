#!/usr/bin/env python3
"""Run all nine v19 exact verification groups in disposable copies."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
JOBS = [
    ('m3_rank_one', 'm3', 'verify_rank1_exclusion.py', 'rank1_exact_report.json'),
    ('m3_squarefree', 'm3', 'verify_squarefree.py', 'verification.json'),
    ('m3_double_root', 'm3', 'verify_double_root_exclusion.py', 'double_root_verification.json'),
    ('m4_fixed_image', 'm4', 'verify_m4_fixed_image.py', 'm4_fixed_image_report.json'),
    ('m4_constant_regular', 'm4', 'verify_constant_regular.py', 'constant_regular_certificate.json'),
    ('m4_independent_pencil', 'm4', 'verify_varying_independent.py', 'varying_independent_certificate.json'),
    ('m4_two_varying_pencils', 'm4', 'verify_varying_pencils.py', 'varying_pencil_certificate.json'),
    ('independent_moving_rank_one', 'audits', 'verify_m4_moving_rank1.py', 'm4_moving_rank1_report.json'),
    ('independent_regular_audit', 'audits', 'independent_regular_audit.py', 'independent_regular_audit_report.json'),
]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    out=(args.output_dir or Path(tempfile.mkdtemp(prefix='p4_v19_reports_'))).resolve()
    (out/'logs').mkdir(parents=True,exist_ok=True)
    (out/'reports').mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='p4_v19_checks_') as tmp:
        base=Path(tmp)
        def run(job):
            label,source,script,report=job
            work=base/label
            shutil.copytree(ROOT/'proofs'/source,work)
            (work/report).unlink(missing_ok=True)
            result=subprocess.run([sys.executable,str(work/script)],cwd=work,
                                  text=True,capture_output=True,check=False)
            (out/'logs'/(label+'.txt')).write_text(result.stdout+result.stderr)
            if result.returncode!=0:
                raise RuntimeError(f'{label} failed; see log')
            data=json.loads((work/report).read_text())
            if data.get('all_passed') is not True:
                raise RuntimeError(f'{label} did not report all_passed')
            (out/'reports'/(label+'.json')).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
            return {'name':label,'all_passed':True}
        with ThreadPoolExecutor(max_workers=3) as pool:
            results=list(pool.map(run,JOBS))
    summary={'all_passed':True,'version':19,'check_count':len(results),
             'checks':results,'archived_originals_modified':False}
    (out/'verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'output_directory':str(out),**summary},indent=2))

if __name__=='__main__':
    main()
