#!/usr/bin/env python3
"""Independent integer-arithmetic audit of the after33 scalar certificate."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import resource
import shlex
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
PARENT = 'S61_BINARY_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v62'
STEM = 'S61_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v63_REVIEW'
PD = ROOT / 'reports' / PARENT

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def stats(lines):
    total = active = degree = bits = 0
    for line in lines:
        a = line.split(); n = int(a[0]); assert len(a) == 1 + 4*n
        active += bool(n); total += n; seen = set()
        for j in range(n):
            q = Fraction(a[4*j+1]); ex = tuple(map(int, a[4*j+2:4*j+5]))
            assert q and all(x >= 0 for x in ex) and ex not in seen
            seen.add(ex); degree = max(degree, sum(ex))
            bits = max(bits, abs(q.numerator).bit_length(), q.denominator.bit_length())
    return dict(terms=total, active=active, degree=degree, coefficient_bits=bits)

def cap():
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))

def run(output, work):
    begin = time.monotonic(); output.mkdir(parents=True, exist_ok=False); work.mkdir(parents=True, exist_ok=False)
    parent = json.loads((PD/'report.json').read_text())
    assert parent['all_passed'] and parent['executed_native_unit_pivots'] == 33
    assert not parent['native34th_pivot_executed']
    for name, digest in parent['input_sha256'].items():
        assert sha(ROOT/name) == digest, name
    info = json.loads((PD/'source.json').read_text()); source = ROOT/info['source']
    prior = ROOT/'reports/S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62/report.json'
    old = json.loads(prior.read_text()); digest = sha(source)
    assert old['all_passed'] and old['generated_output']['sha256'] == digest == info['source_sha256']
    with gzip.open(PD/'source_row.txt.gz', 'rt') as f:
        assert list(map(int, f.readline().split())) == [1,58,0]
        row = list(f); assert len(row) == 58
    with source.open() as f:
        assert list(map(int, f.readline().split())) == [24,58,1]
        metadata = [f.readline().strip(), f.readline().strip()]
        for i in range(8):
            current = [f.readline() for _ in range(58)]
            if i == 7:
                assert all(x.split() == y.split() for x,y in zip(current,row))
    row_stats = stats(row); assert row_stats['terms'] == 50256
    with gzip.open(PD/'cofactor.txt.gz', 'rt') as f:
        assert int(f.readline()) == 58
        cof = list(f); assert len(cof) == 58
    cof_stats = stats(cof)
    assert cof_stats == dict(terms=975, active=33, degree=7, coefficient_bits=3299)
    cpp = ROOT/'code/polynomial_scalar_integer_review_20261103_v63.cpp'
    binary = work/'integer-check'
    flags = shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'], text=True))
    subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(binary),*flags],check=True,capture_output=True,text=True)
    rp = work/'row'; lp = work/'cofactor'; tp = work/'target'
    rp.write_text('1 58 0\n'+''.join(row)); lp.write_text('58\n'+''.join(cof))
    # P is the third serialized variable, not the original d before P=ad.
    tp.write_text('3 5 1 0 3 -6 1 1 2 -4 1 0 2\n')
    p = subprocess.run([str(binary),str(rp),str(lp),str(tp)],check=True,capture_output=True,text=True,preexec_fn=cap)
    integer = json.loads(p.stdout); assert integer['all_passed'] and integer['components'] == 58
    assert integer['cofactor_terms'] == 975 and integer['active'] == 33 and integer['integer_result_terms'] == 3
    (output/'integer-check.log').write_text(p.stdout)
    # Independent inhomogeneous two-component fixture: (e,ae-1).(a,-1)=1.
    rp.write_text('1 2 0\n1 1 0 1 0\n2 1 1 1 0 -1 0 0 0\n')
    lp.write_text('2\n1 1 1 0 0\n1 -1 0 0 0\n')
    tp.write_text('1 1 0 0 0\n')
    f = subprocess.run([str(binary),str(rp),str(lp),str(tp)],check=True,capture_output=True,text=True,preexec_fn=cap)
    fixture = json.loads(f.stdout); assert fixture['all_passed']
    tp.write_text('1 2 0 0 0\n')
    bad = subprocess.run([str(binary),str(rp),str(lp),str(tp)],capture_output=True,text=True,preexec_fn=cap)
    assert bad.returncode == 3 and not json.loads(bad.stdout)['all_passed']
    files = [Path(__file__),cpp,ROOT/'research/notes'/f'{STEM}.md',PD/'report.json',prior]
    files += [ROOT/n for n in parent['input_sha256']]
    report = dict(id=STEM,all_passed=True,method='global denominator clearing; integer multivariate polynomial dot product',
        actual_source=dict(path=info['source'],sha256=digest,bytes=source.stat().st_size,shape=[24,58],row=7,all58_entries_bound=True,metadata=metadata),
        row_stats=row_stats,cofactor_stats=cof_stats,integer_check=integer,
        fixture_passed=True,wrong_target_rejected=True,
        exact_target='a*P^2*(5*P-6*e-4)',known_unit_exponents_a_F_B0_C_J_P=[1,0,0,0,1,2],
        executed_expanded_native_pivots=33,native34_rerun=False,
        conclusions=['an exact next localized unit row is available','append this image column, pivot, and scale by the known unit preserves the localized cokernel'],
        not_claimed=['a replacement column without a separate invertible-column argument','full product41 or full cubic57 or quartic coverage','any QQ nonmembership from the bounded modular search refusal'],
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(files)},seconds=time.monotonic()-begin)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',STEM,report['seconds'])

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    a = p.parse_args();run(a.output,a.work)
