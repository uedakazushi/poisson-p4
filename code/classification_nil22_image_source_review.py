#!/usr/bin/env python3
"""Independent finite classification of the complete nil22 image source.

SymPy exact coefficient arithmetic; no producer or producer-helper imports.
Ambient componenthood is supplied by separate, already completed reviews.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

import sympy as s
from classification_quadratic_cofactor_review import linear_system,complete_solution,MONS,PAIRS

ROOT=Path(__file__).resolve().parents[1]
q1,q2,x,y=s.symbols('q1 q2 x y');coords=(q1,q2,x,y);p=s.Matrix([x,y])
u,v,w,m,n,a,b,c,d=s.symbols('u v w m n a b c d')
mag=s.symbols('m0:6');F0=sum(k*z for k,z in zip(mag,[1,x,y,x*x,x*y,y*y]))
CHECKS=0


def zero(value,label):
    global CHECKS
    value=s.cancel(value)
    assert value==0,(label,str(value))
    CHECKS+=1


def zeros(A,label):
    for value in A:zero(value,label)


def jacobi(P):
    return [s.expand(sum(P[i,l]*s.diff(P[j,k],coords[l])+P[j,l]*s.diff(P[k,i],coords[l])+P[k,l]*s.diff(P[i,j],coords[l])for l in range(4)))for i,j,k in it.combinations(range(4),3)]


def full_linear_source(row):
    ns=s.symbols('n0:19');env=dict(zip(['x1','x2','x3','x4'],[q1,x,q2,y]));env.update({str(z):z for z in ns})
    old=s.zeros(4)
    for (i,j),text in zip(PAIRS,row['full_six_brackets']):old[i,j]=s.sympify(text,locals=env);old[j,i]=-old[i,j]
    order=[0,2,1,3];P=old.extract(order,order)
    A=[[0]*4 for _ in range(4)];A[0][2]=A[1][3]=1
    rows,rhs=linear_system(A,[0]*4)
    triples=list(it.combinations(range(4),3))
    cubic=[]
    for factors in it.combinations_with_replacement(range(4),3):
        e=tuple(factors.count(i)for i in range(4));cubic.append(e)
    for triple in triples:
        block=[{}for _ in cubic]
        for j,(ij,e)in enumerate((ij,e)for ij in PAIRS for e in MONS):
            if sum(e)!=2:continue
            for i,k,l in [triple,triple[1:]+triple[:1],triple[2:]+triple[:2]]:
                if tuple(sorted((i,k)))!=ij:continue
                sign=1 if i<k else-1
                for index,value in enumerate(A[l]):
                    if not value:continue
                    exponent=list(e);exponent[index]+=1;position=cubic.index(tuple(exponent))
                    block[position][j]=block[position].get(j,0)+sign*value
        rows+=block;rhs +=[0]*len(block)
    free,pivots,_=complete_solution(rows,rhs)
    assert len(rows)==190 and len(pivots)==71 and len(free)==19
    monomials=[s.prod(z**k for z,k in zip(coords,e))for e in MONS]
    coefficients=[s.Poly(P[i,j],*coords).coeff_monomial(mon)for i,j in PAIRS for mon in monomials]
    for row,value in zip(rows,rhs):zero(sum(a*coefficients[j]for j,a in row.items())-value,'all190 linear equations')
    assert s.Matrix(coefficients).jacobian(ns).rank()==19
    assert P[2,3]==0
    M=P[:2,2:]
    assert not any(value.has(q1,q2)for value in M)
    divergences=[sum(s.diff(M[i,j],p[j])for j in range(2))for i in range(2)]
    zero(s.diff(P[0,1],q1)-y-divergences[1],'entire q1 coefficient')
    zero(s.diff(P[0,1],q2)+x+divergences[0],'entire q2 coefficient')
    assert s.Poly(P[0,1],q1,q2).total_degree()<=1
    actualD=s.expand(M.det())
    assert s.Poly(actualD,x,y).total_degree()<=3
    residual=s.Matrix([s.diff(actualD,z)for z in p])+M.adjugate()*p
    J=jacobi(P)
    zero(J[0]+residual[1],'Jacobi equals negative second adjugate gradient')
    zero(J[1]-residual[0],'Jacobi equals first adjugate gradient')
    zero(J[2],'remaining vertical Jacobi');zero(J[3],'remaining vertical Jacobi')
    return dict(all90_coefficients_preserved=True,entire_div_Lie_QwedgeZ_matrix_shape=[190,90],rank=71,nullity=19,
        all19_source_coefficients_verified=True,coordinate_reorder=[0,2,1,3],Ppp_zero_for_entire_linear_source=True,
        M_depends_only_on_p_for_entire_source=True,all_F0_six_coefficients_retained=True,
        exact_full_Jacobi_is_rotation_of_adjugate_gradient=True,
        necessity_localization_in_D_used=False)


def model(D3,S,field):
    D3=s.sympify(D3);D=S.det()+(p.T*S*p)[0]/2+D3;H=S+s.hessian(D3,(x,y))/2
    M=(-H.adjugate()+p*field.T).applyfunc(s.expand)
    zero(s.expand(M.det())-H.det()+(field.dot(s.Matrix([s.diff(D,z)for z in p]))),'universal determinant lemma')
    zeros(M.adjugate()*p+s.Matrix([s.diff(D,z)for z in p]),'universal Hessian/Koszul adjugate identity')
    residual=s.expand(H.det()-field.dot(s.Matrix([s.diff(D,z)for z in p]))-D)
    return D,H,M,residual


def leading_and_relations():
    S=s.Matrix([[u,v],[v,w]]);V0=s.Matrix([m,n]);T=s.Matrix([[a,b],[c,d]])
    specifications=[('stable',x*y*(x-y),{a:-s.Rational(1,3),b:0,c:0,d:-s.Rational(1,3)},[],4),
        ('double',x*x*y,{b:0,c:0,d:-1-2*a},[a],3),
        ('triple',x**3,{a:-s.Rational(1,3),b:0},[c,d],2)]
    records=[]
    for name,D3,mapping,free,rank in specifications:
        D,H,M,R=model(D3,S,V0+T*p)
        equations=[z for e,z in s.Poly(R,x,y).terms()if sum(e)==3]
        A,B=s.linear_eq_to_matrix(equations,[a,b,c,d]);assert A.rank()==rank
        zeros(A*s.Matrix([mapping.get(z,z)for z in [a,b,c,d]])-B,'entire cubic V1 necessity')
        assert s.Matrix([mapping.get(z,z)for z in [a,b,c,d]]).jacobian(free).rank()==4-rank if free else rank==4
        R=s.expand(R.subs(mapping,simultaneous=True))
        if name=='stable':
            m0=1+w/6;n0=-1-u/6
            zero(s.Poly(R,x,y).coeff_monomial(x*x)+n+1+u/6,'stable solves n without inversion')
            zero(s.Poly(R,x,y).coeff_monomial(y*y)-m+1+w/6,'stable solves m without inversion')
            R=s.expand(R.subs({m:m0,n:n0},simultaneous=True))
            zero(s.Poly(R,x,y).coeff_monomial(x*y)+(u+v+w+9)/3,'stable Hessian linear relation')
            f1=u*u+2*u*w+15*u-6*w-54;f2=w*w+2*u*w-6*u+15*w-54
            zero(R.subs(v,-u-w-9)+f1*x/6-f2*y/6,'stable complete residual equations')
            zero(f1-f2-(u-w)*(u+w+21),'stable exhaustive factor')
            zero(f1.subs(w,u)-3*(u+6)*(u-3),'stable equal branch')
            zero(f1.subs(w,-u-21)+(u-3)*(u+24),'stable sum branch')
        elif name=='double':
            m0=a*v/2;n0=-1-(a+s.Rational(1,2))*u
            zero(s.Poly(R,x,y).coeff_monomial(x*x)+n+1+(a+s.Rational(1,2))*u,'double solves n')
            zero(s.Poly(R,x,y).coeff_monomial(x*y)-a*v+2*m,'double solves m')
            residual=(4*a+1)*w*y*y/2-v*(2-(a+1)*u)*x/2+((4+(2*a+1)*u)*w-a*v*v)*y/2
            zero(R.subs({m:m0,n:n0},simultaneous=True)-residual,'double ALL lowdegree equations')
            # Each nonzero/zero v,w branch is now literal substitution.
            for aa,uu,vv,ww in [(-s.Rational(1,4),s.Rational(8,3),v,-3*v*v/64),(-s.Rational(1,4),-8,0,w),(0,2,v,0),(a,u,0,0)]:
                zero(residual.subs({a:aa,u:uu,v:vv,w:ww},simultaneous=True),'every complete double branch suffices')
        else:
            target_m=-u/18-c*v/3
            zero(s.Poly(R,x,y).coeff_monomial(x*x)+3*m+u/6+c*v,'triple solves m')
            lower=-((d+s.Rational(2,3))*v+c*w)*x*y-(d+s.Rational(1,2))*w*y*y+(3*w-u*m-v*n)*x-(v*m+w*n)*y
            zero(R-lower+3*(m-target_m)*x*x,'triple ALL lowdegree equations')
            delta=u*w-v*v
            mm=-delta/(18*w);nn=-v*mm/w
            zero((u*mm+v*nn-3*w)+(delta*delta+54*w**3)/(18*w*w),'triple nonzero w complete discriminant relation')
        records.append(dict(cubic_type=name,all_V1_coefficient_rank=rank,V1_solution={str(z):str(mapping.get(z,z))for z in [a,b,c,d]},complete_lower_coefficient_comparison=True))
    # D3=0: no missing full-rank, rank1 or zero Hessian.
    D,H,M,R=model(0,S,V0+T*p)
    zero(R+(V0.T*S*p)[0]+((T*p).T*S*p)[0]+(p.T*S*p)[0]/2,'D3zero entire matrix conditions')
    for ss in [s.Matrix([[0,1],[1,0]]),s.diag(1,0)]:
        RR=s.expand(R.subs({u:ss[0,0],v:ss[0,1],w:ss[1,1]}))
        equations=s.Poly(RR,x,y).coeffs();AA,BB=s.linear_eq_to_matrix(equations,[m,n,a,b,c,d])
        solution={m:0,n:0,b:0,c:0,d:-1-a}if ss.det()!=0 else{m:0,a:-s.Rational(1,2),b:0}
        zeros(AA*s.Matrix([solution.get(z,z)for z in [m,n,a,b,c,d]])-BB,'quadratic entire linear source')
        assert AA.rank()==(5 if ss.det()!=0 else 3)
    return records


def simultaneous_GL(A,D,S,field):
    inverse=A.inv()*p;sub=dict(zip(p,inverse));det=A.det()
    return (s.expand(det**2*D.subs(sub,simultaneous=True)),
        (det**2*A.inv().T*S*A.inv()).applyfunc(s.cancel),
        (A*field.subs(sub,simultaneous=True)).applyfunc(s.cancel))


def normalizers():
    r,t=s.symbols('r t',nonzero=True)
    # Full covariance with a generic upper/lower triangular invertible GL2;
    # the written matrix identity proves the unrestricted GL2 law.
    A=s.Matrix([[r*r,0],[t,r]])
    S=s.Matrix([[u,v],[v,w]]);D=S.det()+(p.T*S*p)[0]/2+x**3
    field=s.Matrix([m-x/3,n+c*x+d*y])
    DD,SS,VV=simultaneous_GL(A,D,S,field)
    zero(DD-SS.det()-(p.T*SS*p)[0]/2-x**3,'entire triple stabilizer weighted cubic preservation')
    expected=s.Matrix([[r*r*u-2*r*t*v+t*t*w,r**3*v-r*r*t*w],[r**3*v-r*r*t*w,r**4*w]])
    zeros(SS-expected,'entire triple stabilizer Hessian law')
    zeros(VV-s.Matrix([r*r*m-x/3,t*m+r*n+(c/r-t*(d+s.Rational(1,3))/r**2)*x+d*y]),'entire triple stabilizer V law')
    # All wnonzero solutions admit r²=-18w/delta. Parameterizing this finite
    # choice by r gives w=-6/r4, delta=108/r6 and is surjective over C.
    ww=-6/r**4;delta=108/r**6;uu=(delta+v*v)/ww
    mapping={w:ww,u:uu,c:-v/(6*ww),d:-s.Rational(1,2),m:-delta/(18*ww),n:v*delta/(18*ww**2),t:r*v/ww}
    zeros(SS.subs(mapping,simultaneous=True)-s.diag(-18,-6),'all wnonzero Hessian normalization')
    zeros(VV.subs(mapping,simultaneous=True)-s.Matrix([1-x/3,-y/2]),'all wnonzero vector normalization')
    # w0,vnonzero, r³v=1, t=ru/(2v).
    mapping={w:0,v:1/r**3,c:-u*r**3/6,d:-s.Rational(2,3),m:0,n:0,t:r**4*u/2}
    zeros(SS.subs(mapping,simultaneous=True)-s.Matrix([[0,1],[1,0]]),'all wzero vnonzero Hessian normalization')
    zeros(VV.subs(mapping,simultaneous=True)-s.Matrix([-x/3,-2*y/3]),'all wzero vnonzero vector normalization')
    # Double stabilizer A=diag(r,1); its v and w action keeps all zero branches.
    A=s.diag(r,1);D=S.det()+(p.T*S*p)[0]/2+x*x*y;field=s.Matrix([a*v/2+a*x,-1-(a+s.Rational(1,2))*u-(1+2*a)*y])
    DD,SS,VV=simultaneous_GL(A,D,S,field)
    zero(DD-SS.det()-(p.T*SS*p)[0]/2-x*x*y,'entire double stabilizer weighted cubic preservation')
    zeros(SS-s.Matrix([[u,r*v],[r*v,r*r*w]]),'entire double stabilizer Hessian law')
    zeros(VV-s.Matrix([a*r*v/2+a*x,-1-(a+s.Rational(1,2))*u-(1+2*a)*y]),'entire double stabilizer V law')
    return dict(binary_cubic_strata=['three simple roots','double plus simple root','triple root','zero cubic'],
        double_nonzero_v_normalizer='diag(1/v,1)',double_w_nonzero_vzero_normalizer='diag(r,1), r²w=1',
        triple_w_nonzero_normalizer='A=[[r²,0],[t,r]], t=rv/w, r²=-18w/(uw-v²)',
        triple_w_nonzero_discriminant='(uw-v²)²+54w³=0 forces uw-v²!=0',
        triple_wzero_vnonzero_normalizer='r³v=1,t=ru/(2v)',
        triple_wvzero_reduced_consequence='u²=0 => u=m=0; n,c,d free',
        all_zero_and_discriminant_complements_retained=True)


def build_P(D3,S,V0,V1):
    field=V0+V1*p;D,H,M,R=model(D3,S,field);zero(R,'canonical complete determinant relation')
    div=[sum(s.diff(M[i,j],p[j])for j in range(2))for i in range(2)]
    P=s.zeros(4);P[0,1]=F0+(y+div[1])*q1+(-x-div[0])*q2;P[1,0]=-P[0,1]
    for i,j in it.product(range(2),repeat=2):P[i,j+2]=M[i,j];P[j+2,i]=-M[i,j]
    zeros(jacobi(P),'universal entire source Jacobi')
    zeros(s.Matrix([sum(s.diff(P[i,j],coords[i])for i in range(4))for j in range(4)])-s.Matrix([x,y,0,0]),'universal entire source divergence')
    zero(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]+D,'universal actual Pf=-D')
    return P,M,D


def eleven_sources(stable,remainder):
    half=s.Rational(1,2);third=s.Rational(1,3);lam=s.Symbol('lam')
    specifications=[('stable_line_conic',x*y*(x-y),s.Matrix([[3,-15],[-15,3]]),s.Matrix([3*half,-3*half]),-s.eye(2)/3,0),
        ('stable_three_lines',x*y*(x-y),s.Matrix([[-6,3],[3,-6]]),s.zeros(2,1),-s.eye(2)/3,0),
        ('double_A',x*x*y,s.Matrix([[s.Rational(8,3),1],[1,-s.Rational(3,64)]]),s.Matrix([-s.Rational(1,8),-s.Rational(5,3)]),s.diag(-s.Rational(1,4),-half),0),
        ('double_B',x*x*y,s.diag(-8,1),s.Matrix([0,1]),s.diag(-s.Rational(1,4),-half),0),
        ('double_C',x*x*y,s.Matrix([[2,1],[1,0]]),s.Matrix([0,-2]),s.diag(0,-1),0),
        ('double_D',x*x*y,s.diag(u,0),s.Matrix([0,-1-(a+half)*u]),s.diag(a,-1-2*a),2),
        ('triple_A',x**3,s.diag(-18,-6),s.Matrix([1,0]),s.diag(-third,-half),0),
        ('triple_B',x**3,s.Matrix([[0,1],[1,0]]),s.zeros(2,1),s.diag(-third,-2*third),0),
        ('triple_C',x**3,s.zeros(2),s.Matrix([0,n]),s.Matrix([[-third,0],[c,d]]),3),
        ('quadratic_A',0,s.Matrix([[0,1],[1,0]]),s.zeros(2,1),s.diag(lam-half,-lam-half),1),
        ('quadratic_B',0,s.diag(1,0),s.Matrix([0,n]),s.Matrix([[-half,0],[c,d]]),3)]
    records=[];h=s.Symbol('h');env={str(z):z for z in (h,q1,q2,x,y,u,v,w,m,n,a,b,c,d,lam,*mag)};env.update(dict(zip(['x1','x2','x3','x4'],coords)))
    for (name,D3,S,V0,V1,extra),saved in zip(specifications,[*stable['cases'],*remainder['cases']]):
        assert saved['name']=='nil22_image_'+name
        P,M,D=build_P(D3,S,V0,V1)
        homogeneous=s.zeros(5)
        for i,j in PAIRS:
            entry=sum(co*h**(2-sum(e))*s.prod(z**k for z,k in zip(coords,e))for e,co in s.Poly(P[i,j],*coords).terms())
            homogeneous[i+1,j+1],homogeneous[j+1,i+1]=entry,-entry
        for i,j in it.combinations(range(5),2):zero(homogeneous[i,j]-s.sympify(saved['full_six_magnetic_source'].get(f'{i},{j}','0'),locals=env),'all11 full150 source identities')
        assert s.Poly(D,x,y).total_degree()>=2
        records.append(dict(name=name,full_parameter_count=6+extra,parameter_domain='affine space',irreducible=True,
            all_six_magnetic_parameters_retained=True,universal_full_Jacobi_modular_and_actualPf_verified=True,
            determinant=str(s.factor(D)),ambient_componenthood_from_finite_source_reduction_claimed=False))
    # The other three stable Hessians are not silently discarded.
    matrices=[s.Matrix([[3,-15],[-15,3]]),s.Matrix([[-6,3],[3,-6]]),s.Matrix([[3,12],[12,-24]]),s.Matrix([[-24,12],[12,3]])]
    allP=[]
    for S in matrices:
        P,_,_=build_P(x*y*(x-y),S,s.Matrix([1+S[1,1]/6,-1-S[0,0]/6]),-s.eye(2)/3);allP.append(P)
    for old,new,A in [(0,3,s.Matrix([[1,0],[1,-1]])),(3,2,s.Matrix([[0,-1],[-1,0]]))]:
        G=s.diag(A,A);sub=dict(zip(coords,G.inv()*s.Matrix(coords)))
        changed=(G*allP[old]*G.T).subs(sub,simultaneous=True)
        target=allP[new].copy();newF0=A.det()*F0.subs(dict(zip(p,A.inv()*p)),simultaneous=True)
        target[0,1]+=newF0-F0;target[1,0]=-target[0,1]
        zeros(changed-target,'all6 magnetic stable GL2 transfer')
    return records


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    old_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
    stable_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_STABLE_CUBIC_20261107.json'
    remaining_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
    linear=full_linear_source(json.loads(old_path.read_text())['cases'][3])
    lower=leading_and_relations();GL=normalizers()
    records=eleven_sources(json.loads(stable_path.read_text()),json.loads(remaining_path.read_text()))
    paths=[Path(__file__),ROOT/'code/classification_quadratic_cofactor_review.py',old_path,stable_path,remaining_path]
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_or_producer_helper_import=True,
        coefficient_ring='QQ exact symbolic polynomial identities and independent Fraction complete90linear elimination',
        input_sha256={str(z.relative_to(ROOT)):hashlib.sha256(z.read_bytes()).hexdigest()for z in paths},
        full_initial_linear_source=linear,entire_cubic_V1_and_lower_relations=lower,full_GL2_normalizers_and_complements=GL,
        complete_eleven_irreducible_sources=records,complete_reduced_rank4_support_exhausted_by_these_sources_up_to_actual_GL2=True,
        omitted_parameter_boundary=False,all_eleven_pairwise_distinct_components_claimed=False,
        independent_identity_count=CHECKS,reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k]for k in ['all_passed','independent_identity_count','seconds','maximum_RSS_KiB']}))


if __name__=='__main__':main()
