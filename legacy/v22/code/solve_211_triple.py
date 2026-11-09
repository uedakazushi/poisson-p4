"""Explore the exhaustive nilpotent-residue / fixed-center 211c slice.

Geometric prerequisites are recorded separately. The output is a rational
linear-combination certificate, never a claim from numerical sampling.
"""
import argparse
import json
import time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix


def linear_consequences(equations,variables):
    polys=[s.Poly(e,*variables,domain=s.QQ) for e in equations]
    mons=sorted({m for p in polys for m,c in p.terms()},key=lambda m:(sum(m),m),reverse=True)
    pos={m:i for i,m in enumerate(mons)}
    rows={i:({pos[m]:s.QQ.from_sympy(c) for m,c in p.terms()}|{len(mons)+i:s.QQ.one}) for i,p in enumerate(polys)}
    rr,piv=DomainMatrix(rows,(len(polys),len(mons)+len(polys)),s.QQ).rref()
    data=rr.to_dod();cert=[];linear=[]
    for i,col in enumerate(piv):
        if col<len(mons) and sum(mons[col])<=1:
            target=sum(s.QQ.to_sympy(c)*s.prod(v**j for v,j in zip(variables,mons[k])) for k,c in data[i].items() if k<len(mons))
            linear.append(target)
            cert.append({'target':str(target),'combination':{str(k-len(mons)):str(s.QQ.to_sympy(c)) for k,c in data[i].items() if k>=len(mons)}})
    return linear,cert

def main(model_path, jet_path, output_path):
    start=time.monotonic();z=s.symbols('z0:43')
    model=json.loads(model_path.read_text());jets=json.loads(jet_path.read_text())
    original=[s.sympify(e) for e in model['jacobi_equations']+jets['all_e8_equations']]
    sub={z[i]:s.S.Zero for i in [37,39,41,26,28,30,22,16,32,23,24,25,19,33]}
    sub.update({z[i]:s.S.One for i in [38,40,42,31]})
    sub.update({z[11]:2*z[20],z[8]:2*z[18],z[12]:2*z[21],z[36]:-s.Rational(3,2)*z[29],
                z[13]:z[34]**2,z[9]:-z[35]**2,z[14]:z[34]*z[35]})
    initial={str(k):str(v) for k,v in sub.items()}
    active=[v for v in z if v not in sub];rounds=[];unit=False
    equations=original
    for stage in range(12):
        equations=sorted(set(s.expand(e.subs(sub,simultaneous=True)) for e in equations)-{s.S.Zero},key=str)
        print('stage',stage,'variables',len(active),'equations',len(equations),flush=True)
        linear,cert=linear_consequences(equations,active)
        if not linear:break
        M,rhs=s.linear_eq_to_matrix(linear,active)
        rr,piv=M.row_join(rhs).rref()
        new={active[j]:s.expand(rr[i,-1]-sum(rr[i,k]*active[k] for k in range(len(active)) if k!=j)) for i,j in enumerate(piv) if j<len(active)}
        rounds.append({'input_equations':[str(e) for e in equations], 'certificates':cert,
                       'linear_equations':[str(e) for e in linear], 'new_substitution':{str(k):str(v) for k,v in new.items()}})
        if len(active) in piv:
            unit=True;break
        sub={k:s.expand(v.subs(new,simultaneous=True)) for k,v in sub.items()}|new
        active=[v for v in active if v not in new]
    output={'model':'211_triple','scope':'Assumes geometric reductions stated in proof; no parameter division in this slice',
            'initial_substitution':initial,'unit':unit,'rounds':rounds,'remaining_parameters':[str(v) for v in active],
            'remaining_equations':[str(e) for e in equations],'final_substitution':{str(k):str(v) for k,v in sub.items()},
            'seconds':round(time.monotonic()-start,3)}
    output_path.write_text(json.dumps(output,indent=2)+'\n')
    print('result',unit,'parameters',len(active),'seconds',output['seconds'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('model',type=Path);ap.add_argument('jets',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args();main(args.model,args.jets,args.output)
