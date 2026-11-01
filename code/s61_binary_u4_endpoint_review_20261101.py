#!/usr/bin/env python3
"""Parent independent local-data review of the nonresonant u4 endpoint."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_binary_endpoints_review_20261101 import review_case
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def review():
    path=ROOT/'reports/S61_BINARY_U4_ENDPOINT_20261101.json';r=json.loads(path.read_text());assert r['all_passed'];case=review_case(r['endpoint'])
    m=S61();u=s.symbols('u1:62');one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];beta={i:u[i] for i in one}
    for i in [3,11,12,34,49]:beta[i]=0
    two,F=specialize(m,beta,True);a=s.Symbol('a');q=s.Matrix(r['exceptional_crossing']['projective_q_coordinates']).applyfunc(s.sympify)
    values=[s.factor(sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2)) for f in F]
    assert [(i,v) for i,v in enumerate(values) if v]!=[]
    assert [i for i,v in enumerate(values) if v]==[70] and s.expand(values[70]+3*a*(5*u[13]-2*u[48])/10)==0
    paths=[Path(__file__),path,ROOT/'code/s61_binary_endpoints_review_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'research/notes/S61_BINARY_U4_ENDPOINT_20261101.md']
    return dict(id='S61_BINARY_U4_ENDPOINT_REVIEW_20261101',all_passed=True,local_case=case,universal_crossing_equation=str(values[70]),accepted=['The nonresonant u49=16 point has free rank65, I^4=I*t=t^2=0, t nonzero and quotient free64.','The same relative14D and smooth equivariant descent as the two earlier endpoints removes the entire generic-beta16D stratum.','The old resonant sample has a genuine transverse R0 line, so its length69 truncation is not a finite transverse algebra.','The extra resonance is preserved by the projective stabilizer and has joint-image dimension<=15.'],limitations=['No saturation statement at the exceptional crossing itself.','Not all binary rank-drop strata.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=review()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,local_case=r['local_case'],crossing=r['universal_crossing_equation'])))
