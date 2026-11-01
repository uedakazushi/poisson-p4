#!/usr/bin/env python3
"""A zero-first-jet toric boundary point outside the finite algebra model."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X,jacobi_projective
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    a,b,t=s.symbols('a b t');x0,x1,x2,x3,x4=X
    B=s.zeros(5)
    B[1,2]=(x0+a*x1)*(x0-a*x1+t*x2);B[2,1]=-B[1,2]
    B[3,4]=(x0+b*x3)*(x0-b*x3+t*x4);B[4,3]=-B[3,4]
    assert all(v==0 for v in jacobi_projective(B))
    G=s.Matrix([[1,0,0,0,0],[1,a,0,0,0],[1,-a,t,0,0],[1,0,0,b,0],[1,0,0,-b,t]])
    assert s.expand(G.det()-a*b*t*t)==0
    z=G*s.Matrix(X);target=s.zeros(5);target[1,2]=a*t*z[1]*z[2];target[2,1]=-target[1,2]
    target[3,4]=b*t*z[3]*z[4];target[4,3]=-target[3,4]
    assert (G*B*G.T-target).applyfunc(s.expand)==s.zeros(5)
    model=S61();boundary=model.normalize(B.subs(t,0));verification=model.verify_point(boundary)
    weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1]
    assert all(boundary['parameters'][i]==0 for i in one)
    assert any(boundary['parameters'][i]!=0 for i,w in enumerate(weights) if w==2)
    assert all(v==0 for v in boundary['perspective'])
    family=model.normalize(B);model.verify_equations(family['parameters'])
    assert all(s.expand(f.subs(t,0)-g)==0 for f,g in zip(family['parameters'],boundary['parameters'],strict=True))
    report=dict(id='TORIC_MISSING_BOUNDARY_20261101',all_passed=True,field='QQ[a,b,t]',
        homogeneous_nonzero_pairs={'12':str(B[1,2]),'34':str(B[3,4])},
        projective_coordinate_change=[[str(v) for v in row] for row in G.tolist()],determinant=str(G.det()),
        transformed_toric_coefficients={'12':'a*t','34':'b*t','others':'0'},
        toric_open='a*b*t != 0',universal_projective_jacobi=True,coordinate_change_identity=True,
        family_all134_verified=True,family_S61_coordinates=list(map(str,family['parameters'])),
        boundary_S61_coordinates=list(map(str,boundary['parameters'])),boundary_first_jet_zero=True,
        boundary_nonzero_q=True,boundary_verification=verification,
        boundary_pfaffian='(1-a^2*y1^2)*(1-b^2*y3^2)',
        consequence='The reduced T01 closure contains nonzero q over beta=0. The finite trace-zero algebra model has q=qbar(beta)+N(beta)c, so all of its beta=0 points have q=0; it does not cover the whole T01 boundary.',
        limitations=['No full T01 boundary classification.','No new Poisson component is counted.'],
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),ROOT/'code/s61_model.py',ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json')})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,explicit_toric_boundary_outside_finite_algebra_model=True)))

if __name__=='__main__':main()
