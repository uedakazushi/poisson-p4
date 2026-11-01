#!/usr/bin/env python3
"""Exact coordinate-factor audit and explicitly labelled producer replays."""
import argparse, hashlib, json, time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def strip_times(value):
    if isinstance(value,dict):return {k:strip_times(v) for k,v in value.items() if k!='seconds'}
    if isinstance(value,list):return list(map(strip_times,value))
    return value

def verify(seven_replay=None,triangular_replay=None):
    start=time.time()
    stems=['S61_BINARY_DEEPER_SEVEN_SLICE','S61_BINARY_DEEPER_TRIANGULAR_CHART']
    paths=[ROOT/'reports'/(n+'_20261101.json') for n in stems]
    old=[json.loads(p.read_text()) for p in paths]
    for r in old:
        assert r['all_passed']
        for p,digest in r['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==digest,p
    if seven_replay is None:
        from s61_binary_deeper_seven_slice_20261101 import verify as run
        seven=run()
    else:seven=json.loads(Path(seven_replay).read_text())
    if triangular_replay is None:
        from s61_binary_deeper_triangular_chart_20261101 import verify as run
        triangular=run()
    else:triangular=json.loads(Path(triangular_replay).read_text())
    assert strip_times(seven)==strip_times(old[0])
    assert strip_times(triangular)==strip_times(old[1])
    b,e,d=s.symbols('b e d');sub={b:(e+1)/5}
    for lhs,rhs in [(20*b-37*e,4-33*e),(40*b+7*e,15*e+8),
                    (35*b-31*e,7-24*e),(5*b-3*e,1-2*e),(5*b-11*e,1-10*e)]:
        assert s.expand(lhs.subs(sub)-rhs)==0
    R=d*(100*b*b-505*b*e+187*e*e)-1000*(2*b-e)
    assert s.expand(R.subs(sub)-(90*d*e*e-93*d*e+4*d+600*e-400))==0
    assert s.sympify(seven['slice_normal_determinant'])==(33*e-4)/100
    for case in triangular['removed_artifact_divisors']:
        assert case['joint_dimension']==16 and case['entire_general_beta_stratum_removed']
        assert s.Rational(case['beta_minor'])!=0
    paths += [Path(__file__),ROOT/'research/notes/S61_BINARY_DEEPER_COORDINATE_REVIEW_20261101.md']
    return dict(id='S61_BINARY_DEEPER_COORDINATE_REVIEW_20261101',all_passed=True,accepted=True,
        verification_kind='Independent coordinate-factor and proof audit; exact supplied-producer replays',
        seven_slice_replay_seconds=seven['seconds'],triangular_replay_seconds=triangular['seconds'],
        all_six_coordinate_factor_identifications_passed=True,removed_divisors=3,
        entire_remaining_deeper_chart_free65_asserted=False,
        limitations=['No independent implementation of high-order point numerics.','No full affine exhaustion.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path)
    p.add_argument('--seven-replay',type=Path);p.add_argument('--triangular-replay',type=Path);a=p.parse_args()
    r=verify(a.seven_replay,a.triangular_replay)
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
