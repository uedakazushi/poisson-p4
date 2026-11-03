#!/usr/bin/env python3
"""Independent exact DAG substitution and two-sided six-generator audit."""
import argparse,gzip,hashlib,json,subprocess,sys
from pathlib import Path
from fractions import Fraction as F
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_GLOBAL_WEIGHT6_TWO_AXES_20261103_v63'

def plus(a,b):
    r=a.copy()
    for m,c in b.items():
        z=r.get(m,F(0))+c
        if z:r[m]=z
        elif m in r:del r[m]
    return r
def times(a,b):
    r={}
    for m,c in a.items():
        for n,d in b.items():
            t=tuple(x+y for x,y in zip(m,n));v=r.get(t,F(0))+c*d
            if v:r[t]=v
            elif t in r:del r[t]
    return r
def values(nodes,inputs):
    out=[]
    for op,*data in nodes:
        if op=='q':
            c=F(data[0]);v={(0,0,0,0):c}if c else{}
        elif op=='v':v=inputs[data[0]]
        elif op=='+':
            v={}
            for j in data[0]:v=plus(v,out[j])
        elif op=='*':
            v={(0,0,0,0):F(1)}
            for j in data[0]:v=times(v,out[j])
        else:raise ValueError(op)
        out.append(v)
    return out

def main(out):
    out.mkdir(parents=True,exist_ok=False);folder=ROOT/'reports'/STEM
    d=json.loads((folder/'report.json').read_text());c=json.load(gzip.open(folder/'certificate.json.gz','rt'))
    pp=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/model37.json.gz';p=json.load(gzip.open(pp,'rt'))
    inputs=[{}for _ in p['input_names']];inputs[0]={(1,0,0,0):F(1)}
    for col,name in enumerate(['z31','z32','z33'],1):
        i=p['input_names'].index(name);assert p['input_weights'][i]==6
        mon=tuple(int(j==col)for j in range(4));inputs[i]={mon:F(1)}
    ev=values(p['nodes'],inputs)
    R,a,x,y,z=ring('a,x,y,z',s.QQ);loc=dict(zip(['a','x','y','z'],s.symbols('a x y z')))
    def parse(t):return R.from_expr(s.sympify(t,locals=loc))
    def poly(dic):return R.from_dict({m:s.QQ(c.numerator,c.denominator)for m,c in dic.items()})
    raw=[poly(ev[i])for i in p['raw150_reconstruction_nodes']]
    fs=[poly(ev[i])for i in p['all350_residual_nodes']]
    assert raw==list(map(parse,c['raw150'])) and fs==list(map(parse,c['source_all350']))
    S,b,W,T,eta=ring('a,W1,W2,eta',s.QQ);sloc=dict(zip(['a','W1','W2','eta'],s.symbols('a W1 W2 eta')))
    def sp(t):return S.from_expr(s.sympify(t,locals=sloc))
    alpha=90*(9*b-8)/41;beta=-360*(9*b-8)/41
    changed=[]
    for f in fs:
        v=S.zero
        for (ia,ix,iy,iz),co in f.items():v+=S(co)*b**ia*(W-alpha*eta)**ix*(T-beta*eta)**iy*eta**iz
        changed.append(v)
    assert changed==list(map(sp,c['transformed_all350']))
    G=[W*W,W*T,T*T,W*eta,T*eta,(9*b-8)*(9*b+8)*eta*eta]
    assert G==list(map(sp,c['canonical_generators']))
    for f,row in zip(changed,c['forward350_by6']):assert f==sum((sp(t)*g for t,g in zip(row,G)),S.zero)
    for g,row in zip(G,c['reverse6_by350_numerator_over_a']):assert b*g==sum((sp(t)*f for t,f in zip(row,changed)),S.zero)
    # Separate small exact ideal/intersection and depth checks, no parent GB.
    m2='''A=QQ[a,ainv,W1,W2,eta]/(a*ainv-1);
m=ideal(W1,W2,eta); g=(9*a-8)*(9*a+8);
I=ideal(W1^2,W1*W2,W2^2,W1*eta,W2*eta,g*eta^2);
assert(I==intersect{m^2,ideal(9*a-8,W1,W2),ideal(9*a+8,W1,W2)});
N=ideal(W1,W2,g*eta); assert(isSubset(N^2,I)); assert(dim(A/I)==1);
R=QQ[u,x,y,z]; J=ideal(x^2,x*y,y^2,x*z,y*z,u*z^2);
assert(dim(R/J)==1); assert(pdim(coker gens J)==3); assert(degree(R/J)==5);
S=QQ[x,y,z]; K=ideal(x^2,x*y,y^2,x*z,y*z);
assert(dim(S/K)==1); assert(pdim(coker gens K)==3);
print "PRIMARY_AND_LOCAL_MULTIPLICITY_PASS"; exit 0;
'''
    script=out/'review.m2';script.write_text(m2)
    rr=subprocess.run(['M2','--script',str(script)],capture_output=True,text=True,timeout=30)
    (out/'review.log').write_text(rr.stdout+rr.stderr);assert rr.returncode==0 and'PRIMARY_AND_LOCAL_MULTIPLICITY_PASS'in rr.stdout
    # One saved fresh producer covers the unchanged all150 T40-axis comparison.
    fresh=ROOT/'work/old-plane-v62/weight6-review-fresh/report.json';fd=json.loads(fresh.read_text())
    assert {k:v for k,v in fd.items()if k!='seconds'}=={k:v for k,v in d.items()if k!='seconds'}
    paths=[Path(__file__),folder/'report.json',folder/'certificate.json.gz',pp,
        ROOT/'research/notes'/f'{STEM}.md',ROOT/'research/notes'/f'{STEM}_REVIEW.md']
    for name,want in d['input_sha256'].items():
        q=ROOT/name;assert hashlib.sha256(q.read_bytes()).hexdigest()==want;paths.append(q)
    report=dict(id=STEM+'_REVIEW',all_passed=True,
        independent_Fraction_all_parent_DAG_substitution_and150_350_coefficients=True,
        all350_forward_and6_reverse_identities_in_polynomial_ring_checked=True,
        only_a_inverse_needed_for_ideal_equality=True,
        independent_small_M2_primary_intersection_and_nilradical_square_zero=True,
        total_closed_scheme_dimension=1,total_closed_scheme_Cohen_Macaulay=True,
        total_closed_crossing_multiplicity=5,special_fiber_depth=0,
        fresh_negative_axis_all150_T40_comparison_unchanged=True,
        positive_axis_component_identification_or_full_R37_nonflatness_claimed=False,
        input_sha256={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest()for q in paths},
        evidence_sha256={q.name:hashlib.sha256(q.read_bytes()).hexdigest()for q in [script,out/'review.log']})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('Weight6 closed two-axis independent review PASS')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
