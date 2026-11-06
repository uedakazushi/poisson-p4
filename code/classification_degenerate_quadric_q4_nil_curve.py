#!/usr/bin/env python3
"""Check a specified quadratic curve, with no iterative jet search."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def decode(data):
    P=s.zeros(5)
    for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(data.get(str(ij),'0'));P[ij[::-1]]=-P[ij]
    return P
def enc(P):return {str(ij):str(s.expand(P[ij])) for ij in itertools.combinations(range(5),2) if P[ij]!=0}
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();xx=s.symbols('x0:5');t=s.Symbol('t')
    paths=['reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_COMPONENTS_20261107.json','reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_NORMAL_20261107.json']
    src=json.loads((ROOT/paths[0]).read_text());normal=json.loads((ROOT/paths[1]).read_text())
    P=decode(next(c['point_brackets'] for c in src['cases'] if c['source']=='nil31_long_II'));R=decode(normal['literal_normal'])
    T=s.zeros(5);T[0,1]=xx[2]**2/5;T[0,3]=2*xx[2]**2/5
    for ij in itertools.combinations(range(5),2):T[ij[::-1]]=-T[ij]
    Q=P+t*R+t*t*T
    J={ijk:s.expand(sum(Q[i,l]*s.diff(Q[j,k],xx[l])+Q[j,l]*s.diff(Q[k,i],xx[l])+Q[k,l]*s.diff(Q[i,j],xx[l]) for l in range(5))) for ijk in itertools.combinations(range(5),3) for i,j,k in [ijk]}
    full={ijkl:s.expand(sum((-1)**at*xx[ell]*J[ijkl[:at]+ijkl[at+1:]] for at,ell in enumerate(ijkl))) for ijkl in itertools.combinations(range(5),4)}
    pf=[]
    for missing in range(5):
        i,j,k,l=[v for v in range(5) if v!=missing]
        pf.append(s.expand((-1)**missing*(Q[i,j]*Q[k,l]-Q[i,k]*Q[j,l]+Q[i,l]*Q[j,k])))
    potential=s.factor(sum(x*f for x,f in zip(xx,pf))/5)
    out={'task_id':'E_reduced_nonnormal','specified_candidate_checked':True,'family':enc(Q),'second_order_correction':enc(T),'all_parameter_full_projective_Jacobi_zero':all(v==0 for v in full.values()),'all_parameter_strict_Jacobi_zero':all(v==0 for v in J.values()),'nonzero_projective_Jacobi':{str(k):str(s.factor(v)) for k,v in full.items() if v},'Euler_contracted_Pfaffian':str(potential),'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start}
    with output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({k:out[k] for k in ['all_parameter_full_projective_Jacobi_zero','all_parameter_strict_Jacobi_zero','nonzero_projective_Jacobi','Euler_contracted_Pfaffian','seconds']},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
