#!/usr/bin/env python3
"""Uniform generic-beta exclusion of regular nilpotent and minimal open strata."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X
from s61_nil22_infinity_20261101 import tensor,beta_matrix,polynomial_minor,sm
from weight1_structure import J,minor_witness,sp_basis
ROOT=Path(__file__).resolve().parents[1]
def run():
    m=S61();x,y,z,w=X[1:];v=s.Matrix([x,y,z,w]);c,p,q,r,B=s.symbols('c p q r B')
    Q=tensor([(s.Rational(1,2)-c)*x*x-x*z,0,-x*x/3,(-3*c*c+3*c/4)*x*x+(-6*c+s.Rational(1,2))*x*z-3*z*z,-x*y/3+p*x*x+q*x*z+(4*c-s.Rational(1,2))*x*w+r*z*z+5*z*w,c*x*x+x*z])
    T=s.Matrix([[1,0,0,0],[0,1,0,-B],[B,0,1,0],[0,0,0,1]]);assert T*J*T.T==J
    pulled=(T.inv()*Q.subs(dict(zip(v,T*v,strict=True)),simultaneous=True)*T.inv().T).applyfunc(s.expand)
    transformed={c:c+4*B/3,p:p+B*q+B*B*r,q:q+2*B*r,r:r}
    assert (pulled-Q.subs(transformed,simultaneous=True)).applyfunc(s.expand)==s.zeros(4)
    Q0=Q.subs(c,0);M=beta_matrix(m,Q0);rr,cc,dd=polynomial_minor(M,dict.fromkeys([p,q,r],0),[p,q,r])
    assert len(cc)==12 and dd==-s.Rational(130977,327680)
    A=s.Matrix([sum(s.diff(Q[i,j],v[j]) for j in range(4)) for i in range(4)]).jacobian(v)
    assert A.rank()==1 and A*A==s.zeros(4)
    trace_orbit=s.Matrix.hstack(*[s.Matrix(aa*A-A*aa).reshape(16,1) for aa in sp_basis()]);trace_witness=minor_witness(trace_orbit);assert trace_witness['rank']==4
    regular=[]
    for label,values in [('A',[z*z/4,x*x/2,-x*z/2,-x*y,y*z-w*w/2,x*w+z*z/4]),('B',[3*x*w/2,-3*x*x/2,3*x*z/2,3*x*y/2-z*w/2,-3*y*z/2+w*w/2,-x*w/2+z*z/2])]:
        R=tensor(values);mat=beta_matrix(m,R);witness=minor_witness(mat);assert witness['rank']==12
        aa=s.Matrix([sum(s.diff(R[i,j],v[j]) for j in range(4)) for i in range(4)]).jacobian(v);assert [((aa**i).rank()) for i in range(1,5)]==[3,2,1,0]
        regular.append({'label':label,'Q':sm(R),'beta_matrix':sm(mat),'rank_minor':witness})
    sources=[ROOT/'legacy/v05/data/minimal_nilpotent_certificate.json',ROOT/'legacy/v05/data/regular_nilpotent_certificate.json']
    inputs=[Path(__file__),ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_infinity_flags.py',ROOT/'code/s61_model.py',m.source/'symplectic_slice_reduced_equations.json']+sources
    return {'id':'S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101','all_passed':True,'field':'QQ','minimal_open':{'Q_after_SL2_normalization':sm(Q),'translation_matrix':sm(T),'translation_parameter_action':{str(a):str(b) for a,b in transformed.items()},'set_B':'-3*c/4','normalized_Q':sm(Q0),'beta_matrix':sm(M),'minor_rows':rr,'minor_columns':cc,'constant_minor':str(dd),'uniform_rank_lower_bound':12,'modular_orbit_minor':trace_witness,'fixed_trace_dimension':6,'affine_stratum_dimension':10,'projective_stratum_dimension':9,'beta_incidence_dimension_upper_bound':13},'regular_nilpotent':regular,'regular_beta_image_dimension_upper_bound':14,'scope':'All regular-nilpotent modular points and the open (a,b)!=0 minimal-nilpotent stratum are absent over sufficiently general beta. Boundary of the latter lies in M0/M1 and is not excluded here. No full infinity or affine-fiber exhaustion.', 'historical_dependencies':[str(p.relative_to(ROOT)) for p in sources],'proof_note':'research/notes/S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101.md','input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'minimal_open_uniform_rank':12,'regular_ranks':[12,12],'output':str(args.output)}))
if __name__=='__main__':main()
