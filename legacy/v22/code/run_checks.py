"""Run v22 exact certificate audits with explicit, separate output files.

Default: four independent Python-stdlib certificate audits.  --symbolic adds
six SymPy checks, including a separate coordinate-change reconstruction of all
25 first jets.  No original model, certificate, or previous edition is written.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output-directory', type=Path, required=True)
    ap.add_argument('--symbolic', action='store_true')
    ap.add_argument('--jobs', type=int, default=2)
    args = ap.parse_args()
    out = args.output_directory.resolve()
    out.mkdir(parents=True, exist_ok=True)
    model = BASE / 'previous_v21/data/tree_jacobi_models'
    jet = BASE / 'data/e8_jets'
    code = BASE / 'code'
    jobs = [
        ('22_stdlib', 'audit_certificates_stdlib.py', [
            '--model', model/'22_node.json', '--data-directory', BASE/'data/22',
            '--output', out/'22_stdlib.json']),
        ('211_stdlib', 'audit_211_certificate_stdlib.py', [
            '--model', model/'211_triple.json', '--jets', jet/'211_triple.json',
            '--certificate', BASE/'data/211_triple_certificate.json',
            '--output', out/'211_stdlib.json']),
        ('31_ss_stdlib', 'audit_31_ss_stdlib.py', [
            '--model', model/'31_external.json', '--jets', jet/'31_external.json',
            '--certificate-directory', BASE/'data/31_ss',
            '--output', out/'31_ss_stdlib.json']),
        ('31_nn_stdlib', 'audit_31_nn_stdlib.py', [
            '--model', model/'31_external.json', '--jets', jet/'31_external.json',
            '--certificate-directory', BASE/'data/31_nn',
            '--output', out/'31_nn_stdlib.json']),
    ]
    if args.symbolic:
        jobs += [
            ('e8_taylor', 'audit_e8_taylor.py', [model, jet, out/'e8_taylor.json']),
            ('211_short', 'audit_211_short.py', [model/'211_triple.json',
                jet/'211_triple.json', BASE/'data/211_triple_certificate.json',
                out/'211_short.json']),
            ('31_nn_family', 'audit_31_nn_rational.py', [model/'31_external.json',
                jet/'31_external.json', BASE/'data/31_nn/31_NN_k1_IB_rational_certificate.json',
                out/'31_nn_family.json']),
            ('fourfold', 'check_fourfold_slice.py', [model/'1111_concurrent.json', out/'fourfold.json']),
            ('line_line', 'verify_line_line.py', ['--output', out/'line_line.json']),
            ('local_weights', 'check_local_weights.py', ['--output', out/'local_weights.json']),
        ]

    def run(job):
        name, script, argv = job
        command = [sys.executable, str(code/script)] + [str(v) for v in argv]
        process = subprocess.run(command, text=True, capture_output=True)
        (out/(name+'.log')).write_text(process.stdout + process.stderr)
        if process.returncode:
            raise RuntimeError(f'{name} failed; see {out/(name+".log")}')
        report = json.loads((out/(name+'.json')).read_text())
        passed = report.get('all_passed', report.get('all_checks_passed', False))
        # Some smaller scripts record the individual Boolean checks only.
        if not passed and isinstance(report.get('checks'), dict):
            passed = all(report['checks'].values())
        if not passed:
            raise AssertionError(f'{name}: report does not certify all checks')
        count = report.get('total_checks', report.get('check_count', report.get('count', report.get('checks'))))
        if isinstance(count, (list, dict)):
            count = len(count)
        result = {'name': name, 'all_passed': True, 'checks': count}
        print(json.dumps(result), flush=True)
        return result

    reports = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for future in as_completed([pool.submit(run, job) for job in jobs]):
            reports.append(future.result())
    reports.sort(key=lambda r: r['name'])
    total = sum(r['checks'] for r in reports if isinstance(r['checks'], int))
    summary = {'all_passed': True, 'group_count': len(reports), 'checks': total,
               'symbolic_checks_included': args.symbolic, 'groups': reports,
               'scope': 'Exact identities and certificate coverage; geometric completeness is proved in the manuscript.'}
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
