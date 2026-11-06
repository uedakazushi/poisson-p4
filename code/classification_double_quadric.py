#!/usr/bin/env python3
"""Exact inverse-form spaces and one rank-four double-quadric family.

All linear-space calculations keep all 90 two-form coefficients. The family
mode tests a full 150-coefficient tangent space, not a restricted ansatz.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import sys
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
V = s.symbols('u v z w')
PAIRS = list(itertools.combinations(range(4),2))
TRIPLES = list(itertools.combinations(range(4),3))


def monomials(variables, degree):
    return [s.prod(variables[i] for i in ids)
            for ids in itertools.combinations_with_replacement(range(len(variables)),degree)]


def exact_minor(M, rank):
    dm = DomainMatrix.from_Matrix(M).convert_to(s.QQ)
    cols = list(dm.rref()[1])
    assert len(cols) == rank
    selected = dm.extract(range(M.rows),cols)
    rows = list(selected.transpose().rref()[1])
    assert len(rows) == rank
    minor = dm.extract(rows,cols).det()
    assert minor
    return dict(rows=rows,columns=cols,determinant=str(minor))


def closed_space(Q, closure_rank=71, extra_beta=(), extra_alpha=()):
    mons = sum([monomials(V,d) for d in range(3)],[])
    targets = sum([monomials(V,d) for d in range(4)],[])
    alpha = s.symbols('a:90')
    entries = {pair:sum(alpha[15*i+j]*m for j,m in enumerate(mons))
               for i,pair in enumerate(PAIRS)}
    def ent(i,j):
        return entries[(i,j)] if i<j else -entries[(j,i)]
    expressions=[]
    for i,j,k in TRIPLES:
        da=s.diff(ent(j,k),V[i])-s.diff(ent(i,k),V[j])+s.diff(ent(i,j),V[k])
        qa=s.diff(Q,V[i])*ent(j,k)-s.diff(Q,V[j])*ent(i,k)+s.diff(Q,V[k])*ent(i,j)
        poly=s.Poly(s.expand(Q*da-2*qa),*V)
        expressions.extend(poly.coeff_monomial(m) for m in targets)
    D=s.linear_eq_to_matrix(expressions,alpha)[0]
    betas=s.symbols('b:20')
    beta=[sum(betas[5*i+j]*m for j,m in enumerate([s.S.One,*V])) for i in range(4)]
    expressions=[]
    for i,j in PAIRS:
        poly=s.Poly(s.expand(Q*(s.diff(beta[j],V[i])-s.diff(beta[i],V[j]))
                            -s.diff(Q,V[i])*beta[j]+s.diff(Q,V[j])*beta[i]),*V)
        expressions.extend(poly.coeff_monomial(m) for m in mons)
    B=s.linear_eq_to_matrix(expressions,betas)[0]
    assert D*B==s.zeros(D.rows,20)
    dminor=exact_minor(D,closure_rank);bminor=exact_minor(B,19)
    extras = []
    for beta_extra in extra_beta:
        extras.append([s.expand(Q*(s.diff(beta_extra[j],V[i])-s.diff(beta_extra[i],V[j]))
                              -s.diff(Q,V[i])*beta_extra[j]+s.diff(Q,V[j])*beta_extra[i])
                       for i,j in PAIRS])
    extras.extend(extra_alpha)
    extra_columns = []
    for form in extras:
        polys = [s.Poly(f,*V) for f in form]
        assert all(poly.total_degree() <= 2 for poly in polys)
        extra_columns.append(s.Matrix([poly.coeff_monomial(m) for poly in polys for m in mons]))
    enlarged = B.row_join(s.Matrix.hstack(*extra_columns)) if extra_columns else B
    assert D*enlarged == s.zeros(D.rows,enlarged.cols)
    enlarged_minor = exact_minor(enlarged,90-closure_rank)
    return dict(Q=str(Q),closure_shape=list(D.shape),primitive_shape=list(B.shape),
                closure_rank=closure_rank,primitive_rank=19,kernel_dimension=90-closure_rank,
                composition_zero=True,kernel_equals_primitive_image=(closure_rank==71),
                closure_minor=dminor,primitive_minor=bminor,
                extra_affine_primitive_classes=[[str(f) for f in form] for form in extras],
                extra_quadratic_one_forms=[[str(f) for f in beta] for beta in extra_beta],
                enlarged_image_dimension=90-closure_rank,enlarged_image_minor=enlarged_minor,
                kernel_equals_enlarged_image=True,
                matrix_sha256={name:hashlib.sha256(json.dumps(matrix.tolist(),default=str).encode()).hexdigest()
                               for name,matrix in [('D',D),('B',B),('enlarged',enlarged)]})


def linear_spaces():
    u,v,z,w=V
    cases={'rank5_transverse':1+u*u+v*v+z*z+w*w,
           'rank5_tangent':w+u*u+v*v+z*z,
           'rank4_H_avoids_vertex':u*u+v*v+z*z+w*w,
           'rank4_H_contains_vertex_transverse':1+u*u+v*v+z*z,
           'rank4_H_contains_vertex_tangent':w+u*u+v*v}
    return {name:closed_space(Q) for name,Q in cases.items()}


def rank_three_spaces():
    u,v,z,w = V
    # The transverse extra class is Q du^dv, hence (du^dv)/Q as a rational form.
    transverse = 1+u*u+v*v
    extra_alpha = [[transverse,0,0,0,0,0]]
    # For Q=w+u^2 and any quadratic r(v,z,w), the cubic part of
    # Q*d(u*dr-r*du)-dQ^(u*dr-r*du) cancels identically.
    extras_beta = [[u*s.diff(r,x)-(r if i==0 else 0) for i,x in enumerate(V)]
                  for r in monomials((v,z,w),2)]
    return {
        'rank3_H_not_containing_vertex':closed_space(u*u+v*v+z*z),
        'rank3_H_contains_vertex_transverse':closed_space(transverse,70,extra_alpha=extra_alpha),
        'rank3_H_contains_vertex_tangent':closed_space(w+u*u,65,extra_beta=extras_beta),
    }


def full_tangent_ranks(point, xx, parameter_directions=()):
    import numpy as np
    sys.path.insert(0,str(ROOT/'legacy/v03/code'))
    import poisson_p4_components_check as p
    def sparse(R):
        return {(i,j):{e:int(a) for e,a in s.Poly(R[i,j],*xx).terms() if a}
                for i,j in itertools.combinations(range(5),2) if R[i,j]!=0}
    ps=sparse(point)
    D=p.linearized_jacobi_matrix(ps);vertical=p.vertical_bivector_matrix()
    orbit=[]
    for a,b in itertools.product(range(5),repeat=2):
        A=s.zeros(5);A[a,b]=1;X=A*s.Matrix(xx)
        action=sum([X[k]*point.diff(xx[k]) for k in range(5)],s.zeros(5))-A*point-point*A.T
        orbit.append(p.bivector_vector(sparse(action)))
    M=np.column_stack([vertical,*orbit,p.bivector_vector(ps)])
    F=np.column_stack([M,*[p.bivector_vector(sparse(direction)) for direction in parameter_directions]])
    assert not np.any(D@F)
    certificates={name:p.modular_rank_certificate(mat)for name,mat in [('Jacobi',D),('vertical',vertical),('orbit_scale',M),('family',F)]}
    ranks={name:entry['rank_mod_prime']for name,entry in certificates.items()}
    return dict(matrix_ranks_mod_prime=ranks,certificates=certificates,
                exact_Jacobi_times_all_actual_directions_zero=True,
                projective_tangent_dimension_upper=125-ranks['Jacobi'],
                projective_orbit_dimension_lower=ranks['orbit_scale']-25,
                projective_family_dimension_lower=ranks['family']-25,
                componenthood_proved=(ranks['family']-25==125-ranks['Jacobi']))


def check_affine_poisson(P, xx):
    jac=[s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])
                      +P[k,l]*s.diff(P[i,j],xx[l]) for l in range(5)))
         for i,j,k in itertools.combinations(range(5),3)]
    assert jac==[0]*10
    return s.factor(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])


def tangent_family():
    h=s.symbols('h');c=s.symbols('c');u,v,z,w=V
    xx=(h,*V)
    lam=2*c-1;mu=1-c;rho=2-2*c
    expressions=(lam*h*u,-mu*h*w-c*u*v,-lam*u*u,
                 -mu*v*v-rho*h*z,-h*w-rho*u*v,c*v*w-rho*u*z)
    P=s.zeros(5)
    for (i,j),f in zip(PAIRS,expressions):P[i+1,j+1]=f;P[j+1,i+1]=-f
    pf=check_affine_poisson(P,xx)
    assert s.expand(pf-(c-1)*(h*w+u*v)**2)==0
    result=full_tangent_ranks(P.subs(c,2),xx,(P.diff(c).subs(c,2),))
    result.update(parameter='c',parameter_open='c != 1',
                  universal_brackets={f'{i+1},{j+1}':str(f)for(i,j),f in zip(PAIRS,expressions)},
                  universal_affine_jacobi_zero=True,universal_projective_divisor=str(h*pf),witness_c=2)
    return result


def anisotropic_tangent():
    h=s.symbols('h');u,v,z,w=V;xx=(h,*V)
    Q=w+u*u+v*v
    Omega=s.zeros(4);Omega[0,2]=Omega[1,3]=1;Omega[2,0]=Omega[3,1]=-1
    Z=s.Matrix([u,v,2*z,2*w]);beta=-Omega*Z
    alpha={(i,j):s.expand(3*Q*Omega[i,j]-s.diff(Q,V[i])*beta[j]+s.diff(Q,V[j])*beta[i])
           for i,j in PAIRS}
    coeffs=[alpha[2,3],-alpha[1,3],alpha[1,2],alpha[0,3],-alpha[0,2],alpha[0,1]]
    P=s.zeros(5)
    for (i,j),f in zip(PAIRS,coeffs):
        f=s.Poly(f,*V)
        P[i+1,j+1]=sum(a*h**(2-sum(e))*s.prod(x**n for x,n in zip(V,e))for e,a in f.terms())
        P[j+1,i+1]=-P[i+1,j+1]
    pf=check_affine_poisson(P,xx)
    assert s.expand(pf+3*(h*w+u*u+v*v)**2)==0
    result=full_tangent_ranks(P,xx)
    result.update(normal_form_Q=str(Q),normalized_c='2/3',
                  universal_brackets={f'{i+1},{j+1}':str(P[i+1,j+1])for i,j in PAIRS},
                  affine_jacobi_zero=True,projective_divisor=str(h*pf))
    return result


def rank_three_extra_sources():
    h=s.symbols('h');u,v,z,w=V;xx=(h,*V)
    bv,bz,bw,A,B,C=s.symbols('bv bz bw A B C');parameters=(bv,bz,bw,A,B,C)
    result={}
    for name,r in [('square',v*v),('mixed',w*v)]:
        if name=='square':
            av=(16*bz*A-16*C-bv)/3;az=-3*bz;aw=-2*bw-8*bz*bz
            sigmav=bv-av;sigmaz=4*bz;aconst=-sigmav*C+4*bz*B
            eta0=s.Matrix([A*v-s.Rational(3,4)*z-2*bz*w+B,v/4+C,-2*bz*v])
            k=bw+8*bz*bz;scale=4
        else:
            av=(bv-4*bz*A)/3;az=3*bz;aw=-2*bw-2*C
            sigmav=bv-av;sigmaz=-2*bz;aconst=-sigmav*C-2*bz*B
            eta0=s.Matrix([A*v-s.Rational(3,2)*z+B,-v/2+C,2*bz*w])
            k=bw+2*C;scale=2
        au=sigmaz*bv-sigmav*bz
        beta=s.Matrix([au*u+av*v+az*z+aw*w+aconst,*list(u*s.Matrix([bv,bz,bw])+eta0)])
        gamma=beta+s.Matrix([-r,*[u*s.diff(r,x)for x in (v,z,w)]])
        Q=w+u*u
        alpha={(i,j):s.expand(Q*(s.diff(gamma[j],V[i])-s.diff(gamma[i],V[j]))
                              -s.diff(Q,V[i])*gamma[j]+s.diff(Q,V[j])*gamma[i])for i,j in PAIRS}
        pf=s.expand(alpha[0,1]*alpha[2,3]-alpha[0,2]*alpha[1,3]+alpha[0,3]*alpha[1,2])
        assert s.expand(pf-k*Q**2)==0
        coeffs=[alpha[2,3],-alpha[1,3],alpha[1,2],alpha[0,3],-alpha[0,2],alpha[0,1]]
        P=s.zeros(5)
        for (i,j),f in zip(PAIRS,coeffs):
            poly=s.Poly(f,*V)
            assert poly.total_degree()<=2
            P[i+1,j+1]=s.expand(scale*sum(a*h**(2-sum(e))*s.prod(x**n for x,n in zip(V,e))for e,a in poly.terms()))
            P[j+1,i+1]=-P[i+1,j+1]
        assert s.expand(check_affine_poisson(P,xx)-scale**2*k*(h*w+u*u)**2)==0
        point_sub={bv:0,bz:0,bw:1,A:0,B:0,C:0}
        point=P.subs(point_sub)
        directions=[(3*P.diff(t).subs(point_sub)).applyfunc(s.expand)for t in parameters]
        ranks=full_tangent_ranks(point,xx,directions)
        ranks.update(parameter_names=list(map(str,parameters)),nondegeneracy_polynomial=str(k),
                     primitive_beta=list(map(str,beta)),extra_r=str(r),bivector_scale=scale,
                     all_parameter_Jacobi_and_Pfaffian_exact=True,
                     universal_brackets={f'{i+1},{j+1}':str(P[i+1,j+1])for i,j in PAIRS},
                     witness_substitution={str(a):str(b)for a,b in point_sub.items()},
                     witness_brackets={f'{i+1},{j+1}':str(point[i+1,j+1])for i,j in PAIRS})
        result[name]=ranks
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['linear','rank3','family','anisotropic','rank3_extra']);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));resource.setrlimit(resource.RLIMIT_CPU,(55,56))
    start=time.monotonic()
    value={'linear':linear_spaces,'rank3':rank_three_spaces,'family':tangent_family,'anisotropic':anisotropic_tangent,'rank3_extra':rank_three_extra_sources}[args.mode]()
    report=dict(date='2026-11-07',mode=args.mode,result=value,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-start,global_exhaustion_proved=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x')as stream:json.dump(report,stream,indent=2,default=int);stream.write('\n')
    summary = ({name:{k:case[k]for k in ('matrix_ranks_mod_prime','projective_tangent_dimension_upper','projective_family_dimension_lower')}
                for name,case in value.items()} if args.mode=='rank3_extra' else
               {k:value[k] for k in ('matrix_ranks_mod_prime','projective_tangent_dimension_upper',
                                     'projective_family_dimension_lower','componenthood_proved')}
               if args.mode in ('family','anisotropic') else {k:v['kernel_dimension']for k,v in value.items()})
    print(json.dumps({'mode':args.mode,'seconds':report['elapsed_seconds'],'output':str(args.output),'result_summary':summary},default=int))

if __name__=='__main__':main()
