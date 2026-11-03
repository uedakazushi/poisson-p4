#!/usr/bin/env python3
"""Independent Fraction-monomial arithmetic review of the uniform Qa shear."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
N=6;ZERO={};ONE={(0,)*N:Q(1)}
def add(*polys):
    out={}
    for p in polys:
        for m,c in p.items():out[m]=out.get(m,Q(0))+c
    return {m:c for m,c in out.items()if c}
def scale(p,c):return {m:v*Q(c)for m,v in p.items()if v*Q(c)}
def mul(*polys):
    out=ONE
    for p in polys:
        tmp={}
        for m,c in out.items():
            for n,d in p.items():
                k=tuple(i+j for i,j in zip(m,n));tmp[k]=tmp.get(k,Q(0))+c*d
        out={m:c for m,c in tmp.items()if c}
    return out
def power(p,n):
    out=ONE
    for _ in range(n):out=mul(out,p)
    return out
def derivative(p,i):
    out={}
    for m,c in p.items():
        if m[i]:n=list(m);n[i]-=1;out[tuple(n)]=c*m[i]
    return out
def substitute(p,values):
    terms=[]
    for m,c in p.items():terms.append(scale(mul(*[power(v,k)for v,k in zip(values,m)]),c))
    return add(*terms)
def variable(i):m=[0]*N;m[i]=1;return {tuple(m):Q(1)}
V=[variable(i)for i in range(N)];a,h,x,y,z,w=V
def bracket(f,g,P):
    return add(*[mul(v,add(mul(derivative(f,i),derivative(g,j)),scale(mul(derivative(f,j),derivative(g,i)),-1)))for(i,j),v in P.items()])
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def parse(text,names):
    polynomial=s.Poly(s.sympify(text),*names,domain=s.QQ)
    return {m:Q(int(c.p),int(c.q))for m,c in polynomial.terms()if c}

def build(args):
    start=time.monotonic();folder=ROOT/'reports/Q_QA_UNIFORM_DARBOUX_20261103_v63'
    report_path=folder/'report.json';cert_path=folder/'certificate.json';report=json.loads(report_path.read_text());cert=json.loads(cert_path.read_text());assert report['all_passed']
    replay=json.loads((args.replay/'report.json').read_text());assert {k:v for k,v in replay.items()if k!='seconds'}=={k:v for k,v in report.items()if k!='seconds'}
    assert cert_path.read_bytes()==(args.replay/'certificate.json').read_bytes()
    inputs=[Path(__file__),ROOT/'research/notes/Q_QA_UNIFORM_DARBOUX_20261103_v63_REVIEW.md',report_path,cert_path]
    for name,digest in report['input_sha256'].items():path=ROOT/name;assert sha(path)==digest,name;inputs.append(path)
    source={(2,5):mul(a,h,h),(3,4):scale(mul(h,h),Q(-2,9)),(3,5):scale(mul(h,x),Q(4,3)),(4,5):add(scale(mul(h,y),Q(2,3)),mul(x,x))}
    constant={(2,5):mul(a,h,h),(3,4):scale(mul(h,h),Q(-2,9))}
    correction=add(scale(mul(h,x,z),6),scale(mul(h,y,y),Q(-3,2)),scale(mul(x,x,y),Q(-9,2)))
    F=[power(h,3),mul(h,h,x),mul(h,h,y),mul(h,h,z),add(mul(h,h,w),correction)]
    G=[*F[:4],add(mul(h,h,w),scale(correction,-1))]
    source_names=s.symbols('a h x y z w');target_names=s.symbols('a H X Y Z W')
    assert [parse(v,source_names)for v in cert['homogeneous_forward']]==F
    assert [parse(v,target_names)for v in cert['homogeneous_inverse']]==G
    valuesF=[a,*F];valuesG=[a,*G]
    for f,v in zip(F,V[1:]):assert substitute(f,valuesG)==mul(power(h,8),v)
    for g,v in zip(G,V[1:]):assert substitute(g,valuesF)==mul(power(h,8),v)
    for i,j in itertools.combinations(range(5),2):
        assert bracket(F[i],F[j],source)==substitute(constant.get((i+1,j+1),ZERO),valuesF)
        assert bracket(G[i],G[j],constant)==substitute(source.get((i+1,j+1),ZERO),valuesG)
    affine_values=[a,ONE,x,y,z,w]
    af=[substitute(f,affine_values)for f in F[1:]];ag=[substitute(g,affine_values)for g in G[1:]]
    assert [parse(v,source_names)for v in cert['forward']]==af
    assert [parse(v,source_names)for v in cert['inverse']]==ag
    assert af[:3]==ag[:3]==[x,y,z]
    assert derivative(af[3],5)==derivative(ag[3],5)==ONE
    assert derivative(af[3],0)==derivative(ag[3],0)==ZERO
    for f,v in zip(af,V[2:]):assert substitute(f,[a,ONE,*ag])==v
    for g,v in zip(ag,V[2:]):assert substitute(g,[a,ONE,*af])==v
    affineP={ij:substitute(v,affine_values)for ij,v in source.items()}
    for i,j in itertools.combinations(range(4),2):
        expected={ (0,3):a,(1,2):scale(ONE,Q(-2,9))}.get((i,j),ZERO)
        assert bracket(af[i],af[j],affineP)==expected
    at0=[ZERO,ONE,x,y,z,w];P0={ij:substitute(v,at0)for ij,v in affineP.items()}
    for casimir in [x,af[3]]:
        for v in [x,y,z,w]:assert bracket(casimir,v,P0)==ZERO
    # Uniform gcds and degrees: only h could divide every coordinate because
    # F0=G0=h^3; both last coordinates have fixed +/-9x^2y/2 mod h.
    assert substitute(F[4],[a,ZERO,x,y,z,w])==scale(mul(x,x,y),Q(-9,2))
    assert substitute(G[4],[a,ZERO,x,y,z,w])==scale(mul(x,x,y),Q(9,2))
    assert all(max(sum(m[1:])for m in f)==3 for f in [*F,*G])
    result=dict(id='Q_QA_UNIFORM_DARBOUX_20261103_v63_REVIEW',all_passed=True,
      independent_arithmetic='Dictionary monomials with fractions.Fraction over QQ[a,h,x,y,z,w]; SymPy only parses frozen display strings',
      fresh_producer_report_and_all_certificate_bytes_match=True,
      all_twenty_raw_homogeneous_Poisson_brackets_independently_rebuilt=True,
      both_five_coordinate_homogeneous_compositions_exact_h8_identity=True,
      affine_forward_inverse_all_six_brackets_and_Jacobian1_rebuilt=True,
      no_a_inverse_or_a_dependent_coordinate_map=True,
      bidegree_3_3_uniform_gcd_one_fixed_monomial_certificate=True,
      a_nonzero_entire_affine_open_one_symplectic_leaf=True,
      a_zero_polynomial_Casimirs_x_W_and_all_A2_fibers_leaves=True,
      no_projective_PGL_or_known_component_or_boundary_exhaustion_inferred=True,
      input_sha256={str(path.relative_to(ROOT)):sha(path)for path in dict.fromkeys(inputs)},seconds=time.monotonic()-start)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('INDEPENDENT Fraction ALL-Qa CUBIC DARBOUX/CREMONA REVIEW PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--replay',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args())
