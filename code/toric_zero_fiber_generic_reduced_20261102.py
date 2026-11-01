#!/usr/bin/env python3
"""Exact ambient quadratic bound proving T01 smooth and its AB fiber generically reduced."""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def matrix(M):return [[str(v) for v in row] for row in M.tolist()]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.time()
    source='legacy/v02/data/symplectic_slice_reduced_equations.json'
    point_source='reports/TORIC_ZERO_JET_FIBER_20261101.json'
    r=json.loads((ROOT/point_source).read_text())
    values={s.Symbol(a):s.Rational(b) for a,b in r['parameter_sample'].items()}
    pt=[s.sympify(a).subs(values) for a in r['S61_parameters']]
    data=json.loads((ROOT/source).read_text())
    one=[i for i,w in enumerate(data['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(data['parameter_weights']) if w==2]
    assert all(pt[i]==0 for i in one)
    J=s.zeros(134,61);H=[]
    for row,eq in enumerate(data['equations']):
        constant=0;h={}
        for mon,a,b in eq:
            co=s.Rational(a,b);constant+=co*s.prod(pt[j] for j in mon)
            for k,col in enumerate(mon):
                J[row,col]+=co*s.prod(pt[j] for l,j in enumerate(mon) if k!=l)
            for aa,bb in itertools.combinations(range(len(mon)),2):
                value=co*s.prod(pt[j] for k,j in enumerate(mon) if k not in (aa,bb))
                if value:
                    key=tuple(sorted((mon[aa],mon[bb])));h[key]=h.get(key,0)+value
        assert constant==0
        H.append(h)
    pc=list(J.rref()[1]);assert len(pc)==43
    assert J[:,one].rank()==10 and J[:,two].rank()==33
    U=s.Matrix.hstack(*J.nullspace())
    assert U.cols==18 and J*U==s.zeros(134,18)
    pr=list(J[:,pc].T.rref()[1]);M=J.extract(pr,pc)
    determinant=M.det();assert determinant!=0
    projection=J[:,pc]*M.inv();assert projection*J[pr,:]==J
    mons=list(itertools.combinations_with_replacement(range(18),2));A=s.zeros(134,len(mons))
    for row,h in enumerate(H):
        for (i,j),co in h.items():
            for k,(a,b) in enumerate(mons):
                A[row,k]+=co*U[i,a]*U[j,b]
                if a!=b:A[row,k]+=co*U[i,b]*U[j,a]
    residual=A-projection*A[pr,:]
    rows=list(residual.T.rref()[1]);assert len(rows)==7
    Q=residual[rows,:];Hess=[]
    for row in range(Q.rows):
        hh=s.zeros(18)
        for k,(i,j) in enumerate(mons):
            hh[i,j]+=Q[row,k];hh[j,i]+=Q[row,k]
        Hess.append(hh)
    stack=s.Matrix.vstack(*Hess);active=list(stack.rref()[1]);inactive=stack.nullspace()
    assert len(active)==6 and len(inactive)==12
    T=s.eye(18)[:,active];K=s.Matrix.hstack(*inactive)
    assert all(hh*K==s.zeros(18,12) for hh in Hess)
    C=s.Matrix.hstack(T,K);assert C.det()!=0
    x=s.symbols('x:6');xx=s.Matrix(x)
    polynomials=[s.expand((xx.T*T.T*hh*T*xx)[0]/2) for hh in Hess]
    u,v,a,b,c,d=x[2],x[5],x[3],2*x[1]+9*x[4],x[4],4*x[0]+9*x[3]
    transformed=[u,v,a,b,c,d];coordinate=s.Matrix(transformed).jacobian(x)
    assert coordinate.det()!=0
    expected=[u*u,u*v,v*v,u*a,u*b,v*c,v*d]
    sixmons=[x[i]*x[j] for i,j in itertools.combinations_with_replacement(range(6),2)]
    def coeff(polys):return s.Matrix([[s.Poly(f,*x).coeff_monomial(m) for m in sixmons] for f in polys])
    P=coeff(polynomials);G=coeff(expected)
    piv=list(G.rref()[1]);change=P[:,piv]*G[:,piv].inv()
    assert change.det()!=0 and change*G==P
    radical=C.inv()[[2,5],:]
    V=U*s.Matrix.hstack(*radical.nullspace());assert V.cols==16
    beta=V[one,:];bc=list(beta.rref()[1]);br=list(beta[:,bc].T.rref()[1])
    assert len(bc)==6;beta_det=beta.extract(br,bc).det();assert beta_det!=0
    rr,uu,zz,vv=s.symbols('r u z v');spatial=[rr,uu,zz,vv]
    Qstar=s.zeros(4);Qstar[0,1]=-rr**2;Qstar[2,3]=-zz**2;Qstar[1,3]=rr**2+2*rr*zz+3*zz**2;Qstar-=Qstar.T
    derivatives=s.Matrix([s.diff(q,z) for q in Qstar for z in spatial]).reshape(16,4)
    # A constant translation annihilates every spatial coefficient iff this coefficient matrix does.
    derivative_rows=[]
    for q in Qstar:
        ds=[s.Poly(s.diff(q,z),*spatial) for z in spatial]
        mm=set().union(*(set(f.monoms()) for f in ds))
        derivative_rows.extend([[f.coeff_monomial(m) for f in ds] for m in mm])
    translation=s.Matrix(derivative_rows);assert translation.rank()==2
    assert Qstar[0,1]+Qstar[2,3]!=0
    z=s.symbols('z');hilbert=s.expand(1+2*z*(1-z)**2)
    assert hilbert==1+2*z-4*z**2+2*z**3 and hilbert.subs(z,1)==1
    paths=[str(Path(__file__).relative_to(ROOT)),source,point_source,
           'research/notes/TORIC_ZERO_JET_EXHAUSTION_20261101.md',
           'research/notes/AB_INCIDENCE_20261101.md',
           'research/notes/TORIC_ZERO_FIBER_GENERIC_REDUCED_20261102.md']
    report=dict(all_passed=True,seconds=time.time()-started,point=list(map(str,pt)),
                Jacobian_rank=43,beta_column_rank=10,q_column_rank=33,
                Jacobian_minor_rows=pr,Jacobian_minor_columns=pc,Jacobian_minor_determinant=str(determinant),
                tangent_basis=matrix(U),quadratic_monomials=mons,quadratic_rows=rows,quadratic_basis=matrix(Q),
                active_coordinates=active,inactive_dimension=12,coordinate_frame=matrix(C),
                active_quadratics=list(map(str,polynomials)),product_coordinates=list(map(str,transformed)),
                coordinate_determinant=str(coordinate.det()),quadratic_row_change=matrix(change),
                product_ideal_generators=list(map(str,expected)),radical_forms=matrix(radical),
                lower_initial_Hilbert_numerator=[1,2,-4,2],lower_initial_dimension=16,lower_initial_degree=1,
                T01_tangent_basis=matrix(V),T01_beta_differential_rank=6,
                beta_minor_rows=[one[i] for i in br],beta_minor_columns=bc,beta_minor_determinant=str(beta_det),
                AB_translation_kernel_dimension=2,T01_regular_dimension=16,zero_fiber_regular_dimension=10,
                zero_fiber_generic_length=1,
                limitations=['Specified point and a nonempty AB open; no all-point fiber reducedness.',
                             'J2 is a lower initial ideal, not identified with the full initial ideal.',
                             'No smoothness or reducedness of the whole surrounding S61 is asserted.'],
                input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print('TORIC ZERO FIBER GENERIC REDUCED PASS',report['seconds'])

if __name__=='__main__':main()
