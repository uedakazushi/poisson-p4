#!/usr/bin/env python3
"""Audit the transition of the complete test curve and replay the zero-pair certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
MAIN='S61_BINARY_ZERO_PAIR_STRATUM_REMOVAL_20261102'
STEM='S61_BINARY_ZERO_PAIR_STRATUM_REMOVAL_REVIEW_20261102'


def verify(replay_path=None):
    source=ROOT/'code'/(MAIN.lower()+'.py')
    primary=ROOT/'reports'/(MAIN+'.json')
    original=json.loads(primary.read_text())
    assert original['all_passed']
    s=sp.Symbol('s',nonzero=True)
    x,z=s*s,1/s
    transformed=[sp.cancel(value*character) for value,character in
                 zip([s**3,6,12,16],[x**-2/z,1/x,x*z**2,1/x],strict=True)]
    assert transformed==[1,6/s**2,12,16/s**2]
    assert sum(original['infinity_hilbert'])==original['infinity_length_mod_32003']==72
    assert all(int(sp.Rational(v).p)%32003 and int(sp.Rational(v).q)%32003
               for v in original['constant_pivot_determinants'])
    with tempfile.TemporaryDirectory(prefix='zero-pair-independent-review-')as tmp:
        if replay_path is None:
            replay_path=Path(tmp)/'audit.json'
            subprocess.run([sys.executable,str(source),'--output',str(replay_path)],
                           cwd=ROOT,check=True,timeout=120,capture_output=True,text=True)
        replay=json.loads(Path(replay_path).read_text())
        seconds=replay.pop('seconds')
        original.pop('seconds')
        assert replay==original
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
               for path,h in original['input_sha256'].items())
    files=[Path(__file__),source,primary,ROOT/'research/notes'/(MAIN+'.md'),
           ROOT/'research/notes'/(STEM+'.md')]
    return dict(id=STEM,all_passed=True,reviewed=MAIN,primary_replay_seconds=seconds,
                original_input_hash_count=len(original['input_sha256']),
                independent_complete_curve_transition=list(map(str,transformed)),
                scope='Producer replay, independent chart arithmetic and geometry/source audit; no second full matrix implementation.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in files})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x')as stream:
        json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['primary_replay_seconds'])
