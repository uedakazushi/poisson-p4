#!/usr/bin/env python3
"""Independent exact checks for the q(5,2) first-jet reconstruction proof."""
import argparse,hashlib,json
from pathlib import Path
from itertools import combinations
import sympy as s
from elliptic_six_variable_fiber import affine_family

ROOT=Path(__file__).resolve().parents[1]
PAIRS4=list(combinations(range(4),2));PAIRS5=list(combinations(range(5),2))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    lam,y,P,_=affine_family(2);r=s.symbols('r0:5')
    Phi=s.Matrix([[0,r[3]/lam,lam*r[4],-lam*r[0],-r[1]/lam],
                  [-r[3]/lam,0,r[0],r[1],-lam*r[2]],
                  [-lam*r[4],-r[0],0,r[2],r[3]],
                  [lam*r[0],-r[1],-r[2],0,r[4]],
                  [r[1]/lam,lam*r[2],-r[3],-r[4],0]])
    A=s.Matrix.hstack(-s.Matrix(y),s.eye(4))
    N=s.Matrix([[0,P[2,3],-P[1,3],P[1,2]],[-P[2,3],0,P[0,3],-P[0,2]],
                 [P[1,3],-P[0,3],0,P[0,1]],[-P[1,2],P[0,2],-P[0,1],0]])
    pf=s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
    assert (P*N+pf*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
    projected=(A*Phi*A.T).applyfunc(s.expand)
    identities=[]
    for rr in r:
        residual=s.expand(sum(N[i,j]*s.diff(projected[i,j],rr) for i,j in PAIRS4))
        assert residual==0
        identities.append(str(residual))
    # The displayed five bivectors span the full constant kernel, not a subset.
    H=(A.T*N*A).applyfunc(s.expand)
    polys=[s.Poly(H[i,j],*y) for i,j in PAIRS5]
    mons=sorted(set().union(*(set(p.monoms()) for p in polys)))
    mat=s.Matrix([[p.coeff_monomial(m) for p in polys] for m in mons])
    kernel=s.Matrix([[s.diff(Phi[i,j],rr) for rr in r] for i,j in PAIRS5])
    assert mat*kernel==s.zeros(len(mons),5)
    assert kernel.rank()==5
    # A sample rank5 minor gives the needed upper bound over QQ(lambda).
    sample=mat.subs(lam,2)
    cols=list(sample.rref()[1]);rows=list(sample[:,cols].T.rref()[1])
    assert len(cols)==len(rows)==5
    determinant=sample.extract(rows,cols).det();assert determinant!=0
    projection=s.Matrix([[s.diff(Phi[i+1,j+1],rr) for rr in r] for i,j in PAIRS4])
    projrows=list(projection.T.rref()[1]);assert len(projrows)==5
    projdet=s.factor(projection[projrows,:].det());assert projdet!=0
    # Closed first derivatives map isomorphically to the graph's20 entries.
    J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    b=s.symbols('b0:24');B=[s.zeros(4) for _ in range(4)]
    for h,(i,j) in enumerate(PAIRS4):
        for a in range(4):B[a][i,j]=b[4*h+a];B[a][j,i]=-b[4*h+a]
    closed=s.Matrix([B[i][j,k]-B[j][i,k]+B[k][i,j] for i,j,k in combinations(range(4),3)]).jacobian(b)
    basis=s.Matrix.hstack(*closed.nullspace());assert basis.cols==20
    pl=s.symbols('z0:5');Z=s.zeros(4)
    for (i,j),value in zip(PAIRS4,[*pl,-pl[0]]):Z[i,j]=value;Z[j,i]=-value
    ell=s.Matrix([sum(B[a][i,j]*Z[i,j] for i,j in PAIRS4) for a in range(4)])
    graph=s.Matrix([s.diff(v,zz) for v in -J*ell for zz in pl]).jacobian(b)
    assert graph.rank()==20 and (graph*basis).det()!=0
    # Verify the derivative sign using formal y,w,Z and Omega0=J.
    w=s.Matrix(s.symbols('w0:4'));yy=s.Matrix(s.symbols('u0:4'))
    perspective=yy*w.T-w*yy.T
    assert s.expand(sum(J[i,j]*perspective[i,j] for i,j in PAIRS4)-(yy.T*J*w)[0])==0
    report={'id':'ELLIPTIC_Q52_RECONSTRUCTION_REVIEW_20261101','all_passed':True,
        'field':'QQ(lambda)','universal_annihilation_identities':5,
        'constant_kernel_dimension':5,'kernel_coefficient_matrix_shape':[len(mons),10],
        'sample_lambda':2,'kernel_rank_minor_rows':rows,'kernel_rank_minor_columns':cols,
        'kernel_rank_minor_determinant':str(determinant),
        'spatial_projection_minor_rows':projrows,'spatial_projection_minor_determinant':str(projdet),
        'closed_derivative_dimension':20,'closed_derivative_to_graph_determinant':str((graph*basis).det()),
        'graph_reconstruction_sign':'w=-J*ell_beta',
        'accepted_component_degree':'d_T40=1, using the accompanying geometric argument and established dominance',
        'scope':['Generic q(5,2) component only.','No q(5,1) degree or full S61 exhaustion.'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),ROOT/'code/elliptic_six_variable_fiber.py',ROOT/'code/elliptic_weight1_dominance.py')}}
    with args.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({'all_passed':True,'kernel_dimension':5,'accepted_generic_q52_degree':1}))

if __name__=='__main__':main()
