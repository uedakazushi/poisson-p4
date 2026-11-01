#!/usr/bin/env python3
"""Uniform rank bounds on all rank-two isotropic square-zero tangent fields.

The proof note uses a proper tangent-field incidence to include every L0
boundary point. The calculation deliberately keeps the three free quadratic
coefficients, and every nonzero t, rather than just a dense L0 sample.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_infinity_flags import leading_equations_matrix
from s61_model import S61,X
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]


def sm(matrix):return [[str(v) for v in row] for row in matrix.tolist()]


def tangent_kernel(Q,q,z,U):
    equations=s.Poly(Q[0,1]+Q[2,3],*z).coeffs()
    for i,j,k in combinations(range(4),3):
        expression=s.expand(U[i]*Q[j,k]-U[j]*Q[i,k]+U[k]*Q[i,j])
        equations+=s.Poly(expression,*z).coeffs()
    M=s.Matrix(equations).jacobian(q)
    dr,pivots=M.to_DM().rref();R=dr.to_Matrix()
    free=[i for i in range(45) if i not in pivots]
    N=s.zeros(45,len(free))
    for j,column in enumerate(free):
        N[column,j]=1
        for i,pivot in enumerate(pivots):N[pivot,j]=-R[i,column]
    assert len(pivots)==41 and len(free)==4 and (M*N).applyfunc(s.cancel)==s.zeros(M.rows,4)
    rows=list(M[:,list(pivots)].T.to_DM().rref()[1])
    det=s.factor(M.extract(rows,list(pivots)).det())
    return M,N,free,rows,list(pivots),det


def run():
    model=S61();two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    q=s.symbols('q0:45');parameters=[s.S.Zero]*61
    for i,value in zip(two,q,strict=True):parameters[i]=value
    bivector=model.reconstruct(parameters).subs(X[0],1)
    assert bivector[0,:]==s.zeros(1,5)
    Q=bivector[1:,1:]-J;z=list(X[1:]);x,y,zz,w=z
    t=s.Symbol('t');f=s.symbols('f0:3');results=[]
    for label,U in [('off_diagonal',s.Matrix([0,zz,0,t*x])),
                    ('exceptional_form',s.Matrix([0,x+zz,0,-x]))]:
        M,N,free,rows,columns,det=tangent_kernel(Q,q,z,U)
        if label=='off_diagonal':assert det==t**13
        else:assert det!=0 and not det.free_symbols
        p=[s.S.Zero]*61
        for i,value in zip(two,N*s.Matrix([*f,1]),strict=True):p[i]=value
        beta_matrix=leading_equations_matrix(model,p)[:,1:]
        parametric_Q=Q.subs(dict(zip(q,N*s.Matrix([*f,s.Symbol('c3')]))),simultaneous=True).applyfunc(s.factor)
        # The c3=0 subspace has only the constant-plane component Q_yw.
        for i,j in combinations(range(4),2):
            if (i,j)!=(1,3):assert parametric_Q[i,j].subs(s.Symbol('c3'),0)==0
        result={'label':label,'U':list(map(str,U)),'constraints':sm(M),'kernel':sm(N),
                'free_q_columns':free,'constraint_minor_rows':rows,'constraint_minor_columns':columns,
                'constraint_minor_determinant':str(det),'parametric_Q':sm(parametric_Q),
                'beta_leading_matrix':sm(beta_matrix)}
        if label=='off_diagonal':
            blocks=[('b0',[88,98,119,131],[0],1),('b1',[90,102,115,122],[1],1),
                ('b7',[87,114],[7],1),('b2',[97,130],[2],1),
                ('b6_15',[85,92,116],[6,15],2),('b4_13',[89,101,121],[4,13],2),
                ('b5',[84],[5],1),('b3',[100],[3],1),
                ('b11_12',[86],[11,12],1),('b9_14',[91],[9,14],1)]
            used=set();certificates=[];critical=[]
            for name,rr,cc,rank in blocks:
                # Rows have no entries in later blocks or unused columns.
                assert all(beta_matrix[i,j]==0 for i in rr for j in range(16) if j not in used|set(cc))
                diagonal=beta_matrix.extract(rr,cc)
                assert not any(value.has(*f) for value in diagonal)
                minors=[]
                for ri in combinations(range(len(rr)),rank):
                    for ci in combinations(range(len(cc)),rank):
                        minors.append(s.factor(diagonal.extract(ri,ci).det()))
                numerators=[s.fraction(s.cancel(value))[0] for value in minors]
                gcd=s.S.Zero
                for numerator in numerators:gcd=s.gcd(gcd,numerator)
                while gcd!=0 and s.rem(gcd,t,t)==0:gcd=s.cancel(gcd/t)
                gcd=s.Poly(gcd,t).monic().as_expr();assert gcd!=0
                roots=s.solve(gcd,t)
                assert all(value.is_Rational and value!=0 for value in roots)
                for value in roots:assert diagonal.subs(t,value).rank()>=rank-1
                if roots:critical.append(gcd)
                certificates.append({'label':name,'rows':rr,'columns':cc,'generic_rank':rank,
                    'diagonal_block':sm(diagonal),'maximal_minors':list(map(str,minors)),
                    'critical_polynomial_away_from_t_zero':str(gcd),
                    'critical_roots':list(map(str,roots)),'rank_loss_at_most_one':True})
                used.update(cc)
            assert sum(b[3] for b in blocks)==12
            assert all(s.degree(s.gcd(a,b),t)==0 for a,b in combinations(critical,2))
            assert sum(s.degree(a,t) for a in critical)==10
            result.update({'triangular_blocks':certificates,'uniform_beta_rank_lower_bound':11,
                           'exception_polynomials_pairwise_coprime':True})
        else:
            zero=beta_matrix.subs(dict.fromkeys(f,0));cols=list(zero.rref()[1]);rr=list(zero[:,cols].T.rref()[1])
            determinant=s.factor(beta_matrix.extract(rr,cols).to_DM().det().as_expr())
            assert len(cols)==12 and determinant==s.Rational(4886335125,64)
            result.update({'uniform_beta_rank_lower_bound':12,'beta_minor_rows':rr,'beta_minor_columns':cols,
                           'beta_minor_determinant':str(determinant)})
        results.append(result)
    inputs=[Path(__file__),ROOT/'code/s61_infinity_flags.py',ROOT/'code/s61_model.py',
            model.source/'symplectic_slice_reduced_equations.json']
    return {'id':'S61_L0_INFINITY_20261101','all_passed':True,'field':'QQ(t,f0,f1,f2), with all t!=0 certified',
        'weight_two_u_indices_zero_based':two,'compatible_quadratic_matrix':sm(Q),'cases':results,
        'projective_isotropic_squarezero_ranktwo_field_dimension':6,
        'projective_compatible_tangent_fiber_dimension':3,
        'beta_incidence_dimension_bound':14,
        'scope':'For sufficiently general beta, every L0 point of original134 unsaturated infinity '
                'lies in R0. The proper tangent-field incidence and congruence classification are '
                'proved in the note; no affine-fiber or nonzero-modular-stratum exhaustion.',
        'proof_note':'research/notes/S61_L0_INFINITY_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();report=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'uniform_beta_ranks':[11,12],'incidence_bound':14,'output':str(args.output)}))


if __name__=='__main__':main()
