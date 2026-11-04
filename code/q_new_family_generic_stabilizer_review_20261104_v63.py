#!/usr/bin/env python3
"""Independent general-affine-flag exhaustion and exact generic stabilizer audit."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63_REVIEW'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read(p):return json.loads(p.read_text())


def mm(A,B,zero):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),zero)
             for j in range(len(B[0]))] for i in range(len(A))]


def det(A,zero,one):
    r=zero
    for p in itertools.permutations(range(len(A))):
        c=one*((-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p))))
        for i,j in enumerate(p):c*=A[i][j]
        r+=c
    return r


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    if out.exists():
        assert not any(out.iterdir()), 'refuse to overwrite evidence'
    else:out.mkdir(parents=True)
    base=ROOT/'reports/Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63'
    author=read(base/'report.json');cert=read(base/'certificate.json')
    assert author['all_passed']
    for rel,hsh in author['input_sha256'].items():assert sha(ROOT/rel)==hsh,rel
    for rel,hsh in author['evidence_sha256'].items():assert sha(base/rel)==hsh,rel
    geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    sig=ROOT/'reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63/certificate.json'
    sigreview=ROOT/'reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63_REVIEW/report.json'
    assert read(sigreview)['all_passed']
    geometry=read(geo)
    u,v=s.symbols('u v');n0,n1,n2=s.symbols('n0 n1 n2')
    x,y,z,w=s.symbols('x y z w');coords=(x,y,z,w)
    P=s.Matrix(geometry['affine_P']).subs({n0:u,n1:1,n2:v})
    div=s.Matrix([s.expand(sum(s.diff(P[i,j],coords[i]) for i in range(4)))for j in range(4)])
    assert div==s.Matrix([0,4*u/3,u*x+u*v*y+s.Rational(1,48),0])
    # Begin with the ENTIRE affine flag forced by H, vertex and compatible
    # constant pencil. In particular b,d,f,kappa are initially independent.
    a,b,c,d,e,f,g,H,I,J,K,tx,T,Z,tw,kap=s.symbols('a b c d e f g H I J K tx T Z tw kap')
    A=s.Matrix([[a,0,0,b],[e,f,0,g],[H,I,J,K],[c,0,0,d]])
    image=A*s.Matrix(coords)+s.Matrix([tx,T,Z,tw])
    left=kap*A*P*A.T
    assert s.expand(left[0,3]-s.Rational(8,9)-s.Rational(8,9)*(kap*(a*d-b*c)-1))==0
    xy=s.Poly(s.expand(left[0,1]-(s.Rational(1,36)-256*image[3]/3)),*coords)
    yw=s.Poly(s.expand(left[1,3]-4*image[0]/3),*coords)
    assert s.expand(xy.coeff_monomial(x)-s.Rational(4,3)*(64*c-kap*f*b))==0
    assert s.expand(xy.coeff_monomial(w)-s.Rational(256,3)*(d-kap*f*a))==0
    assert s.expand(yw.coeff_monomial(x)-s.Rational(4,3)*(kap*f*d-a))==0
    assert s.expand(yw.coeff_monomial(w)-s.Rational(4,3)*(64*kap*f*c-b))==0
    cov=kap*A*div-div.subs(dict(zip(coords,image)),simultaneous=True)
    assert s.expand(cov[1]-4*u*(kap*f-1)/3)==0
    # Thus eta=kap*f=1, b=64c,d=a. Divergence z coefficients give j,e,g
    # without selecting the a!=0 or c!=0 chart.
    D=a*a-64*c*c
    initial={b:64*c,d:a,f:D,kap:1/D}
    c2=s.Poly(s.expand(cov[2].subs(initial,simultaneous=True)),*coords)
    assert s.cancel(c2.coeff_monomial(y)-u*v*(J-D*D)/D)==0
    assert s.cancel(c2.coeff_monomial(x)-u*(J/D-a-v*e))==0
    assert s.cancel(c2.coeff_monomial(w)+u*(64*c+v*g))==0
    reduction={**initial,J:D*D,e:(D-a)/v,g:-64*c/v}
    tx0=-c/48+2*(a-1)/(3*v)
    tw0=(1-a)/s.Integer(3072)+2*c/(3*v)
    i0=3*D*(tx0+v*T)/4+D*(1-D)/(64*u)
    for q in (xy.coeff_monomial(1),yw.coeff_monomial(1)):
        assert s.cancel(q.subs(reduction,simultaneous=True).subs({tx:tx0,tw:tw0},simultaneous=True))==0
    assert s.cancel(c2.coeff_monomial(1).subs({J:D*D,I:i0,tx:tx0,tw:tw0},simultaneous=True))==0
    # Independent sparse coefficient-ring extraction of ALL six leading
    # brackets after reduction. Only the three Pyz coefficients can remain.
    F=s.QQ.frac_field(u,v,a,c);R,XX,YY,WW=ring('x,y,w',F)
    cv=lambda q:F.from_sympy(s.sympify(q))
    Q=[[R.zero for _ in range(4)]for _ in range(4)]
    for i,j in itertools.combinations(range(4),2):
        for mon,co in s.Poly(P[i,j],*coords).terms():
            if sum(mon)==2:
                assert mon[2]==0
                Q[i][j][(mon[0],mon[1],mon[3])]=cv(co)
        Q[j][i]=-Q[i][j]
    reduced=[[cv(a),F.zero,F.zero,cv(64*c)],
             [cv((D-a)/v),cv(D),F.zero,cv(-64*c/v)],
             [F.zero,F.zero,cv(D*D),F.zero],
             [cv(c),F.zero,F.zero,cv(a)]]
    inputs=[XX,YY,R.zero,WW]
    images=[sum((R(reduced[i][k])*inputs[k]for k in range(4)),R.zero)for i in range(4)]
    def compose(q,im):
        ans=R.zero
        for mon,co in q.items():
            zz=R(co)
            for r0,power in zip(im,mon):
                if power:zz*=r0**power
            ans+=zz
        return ans
    differences={}
    for i,j in itertools.combinations(range(4),2):
        push=sum((R(reduced[i][k]*reduced[j][l]/cv(D))*Q[k][l]
                  for k in range(4)for l in range(4)),R.zero)
        target=compose(Q[i][j],[images[0],images[1],images[3]])
        difference=push-target
        if (i,j)!=(1,2):assert not difference
        else:differences=difference
    p=32*u+v
    q1=-p*D*D+p*a*a+4096*a*c-64*v*c*c
    q2=-D*D+a*a+u*a*c+64*c*c
    q3=v*D*D-v*a*a+4096*a*c+64*p*c*c
    assert differences==R.from_dict({(2,0,0):cv(q1/(64*v)),(1,0,1):cv(64*q2/v),(0,0,2):cv(q3/v)})
    assert s.expand(q1+q3+(v-p)*q2+32*a*c*(u*u-256))==0
    assert s.expand(q3+v*q2-a*c*(4096+u*v)-64*c*c*(p+v))==0
    assert s.expand(q2.subs(c,0)-a*a*(1-a*a))==0
    # General-field matrix arithmetic for N,E and the two full remaining
    # flag families. No author equation evaluator or solver is used.
    F=s.QQ.frac_field(u,v,T);cv=lambda q:F.from_sympy(s.sympify(q))
    zero,one=F.zero,F.one
    eye=[[one if i==j else zero for j in range(5)]for i in range(5)]
    N=[[zero for _ in range(5)]for _ in range(5)]
    N[2][0]=cv(4*u/3);N[3][0]=cv(s.Rational(1,48));N[3][1]=cv(u);N[3][2]=cv(u*v)
    N2=mm(N,N,zero)
    assert mm(N2,N,zero)==[[zero]*5 for _ in range(5)]
    def ga(t):return [[eye[i][j]+t*N[i][j]+t*t*N2[i][j]/2 for j in range(5)]for i in range(5)]
    E=[[cv(q) for q in row]for row in [[1,0,0,0,0],[-4/(3*v),-1,0,0,0],
       [4/(3*v*v),2/v,1,0,0],[-2/(3*v**3),-1/v**2,0,1,0],[s.Rational(1,1536),0,0,0,-1]]]
    assert mm(E,E,zero)==eye and mm(E,N,zero)==mm(N,E,zero)
    assert E==[[cv(q)for q in row]for row in cert['E']]
    sigma=read(sig)['branches'];Mminus=[[cv(q)for q in row]for row in sigma[1]['homogeneous_matrix']]
    Mplus=[[cv(q)for q in row]for row in sigma[0]['homogeneous_matrix']]
    RR=[[cv(q)for q in row]for row in cert['R']]
    assert mm(Mminus,RR,zero)==Mplus
    assert mm(RR,RR,zero)==ga(cv(-1/(128*u*v)))
    assert mm(ga(cv(1/(256*u*v))),RR,zero)==E
    # Full raw homogeneous source, now over QQ(u,v,T,H,K,Z).
    Bfield=s.QQ.frac_field(u,v,T,H,K,Z);bc=lambda q:Bfield.from_sympy(s.sympify(q))
    S,HH,SX,SY,SZ,SW=ring('h,x,y,z,w',Bfield);vs=[HH,SX,SY,SZ,SW]
    raw=[[S.zero for _ in range(5)]for _ in range(5)]
    for i,j in itertools.combinations(range(4),2):
        for mon,co in s.Poly(P[i,j],*coords).terms():raw[i+1][j+1][(2-sum(mon),*mon)]=bc(co)
        raw[j+1][i+1]=-raw[i+1][j+1]
    def subst(q,im):
        ans=S.zero
        for mon,co in q.items():
            zz=S(co)
            for r0,power in zip(im,mon):
                if power:zz*=r0**power
            ans+=zz
        return ans
    branch_evidence=[]
    for branch in cert['all_full_flag_branches']:
        aa=branch['a'];tx1=2*(aa-1)/(3*v);tw1=(1-aa)/s.Integer(3072);ii=3*(tx1+v*T)/4
        Ar=[[aa,0,0,0],[(1-aa)/v,1,0,0],[H,ii,1,K],[0,0,0,aa]]
        tr=[tx1,T,Z,tw1]
        M=[[Bfield.one if i==j else Bfield.zero for j in range(5)]for i in range(5)]
        for i in range(4):
            M[i+1][0]=bc(tr[i])
            for j in range(4):M[i+1][j+1]=bc(Ar[i][j])
        im=[sum((S(M[i][j])*vs[j]for j in range(5)),S.zero)for i in range(5)]
        actual={}
        for i,j in itertools.combinations(range(5),2):
            push=sum((S(M[i][k]*M[j][l])*raw[k][l]for k in range(5)for l in range(5)),S.zero)
            difference=push-subst(raw[i][j],im)
            for mon,co in difference.items():
                assert i>0 and j>0
                actual[(i-1,j-1,*mon[1:])]=co
        expected={tuple(lab):bc(q)for lab,q in zip(branch['coefficient_labels'],branch['full_linear_equations'])}
        assert actual==expected
        basefield=s.QQ.frac_field(u,v,T)
        L=[];const=[]
        for label in branch['coefficient_labels']:
            pp=s.Poly(Bfield.to_sympy(actual[tuple(label)]),H,K,Z,domain=basefield)
            assert pp.total_degree()<=1
            L.append([basefield.from_sympy(pp.coeff_monomial(q))for q in (H,K,Z)])
            const.append(basefield.from_sympy(pp.coeff_monomial(1)))
        chosen=None
        for rows in itertools.combinations(range(len(L)),3):
            dd=det([L[i]for i in rows],basefield.zero,basefield.one)
            if dd:chosen=rows;break
        assert chosen is not None
        timepar=cv(3*T/(4*u) if aa==1 else 3*(T-4/(3*v*v))/(4*u))
        desired=ga(timepar)if aa==1 else mm(ga(timepar),E,zero)
        hs=[desired[3][1],desired[3][4],desired[3][0]]
        assert all(sum((row[j]*hs[j]for j in range(3)),co)==0 for row,co in zip(L,const))
        replacement={H:F.to_sympy(hs[0]),K:F.to_sympy(hs[1]),Z:F.to_sympy(hs[2])}
        assert [[cv(Bfield.to_sympy(q).subs(replacement,simultaneous=True))for q in row]for row in M]==desired
        branch_evidence.append(dict(a=aa,all150_homogeneous_coefficients_checked=True,
            nonzero_coefficient_count=len(actual),own_rank3_rows=list(chosen),
            own_rank3_minor=str(basefield.to_sympy(dd)),Ga_time=str(F.to_sympy(timepar))))
    evidence=dict(general_flag_started_without_block_or_scalar_restriction=True,
        divergence_y='kappa*f=1',linear_brackets='b=64c,d=a',
        divergence_z='j=D^2,e=(D-a)/v,g=-64c/v',
        quadratic_coefficients=[str(q1),str(q2),str(q3)],
        generic_unit_factors=['u','v','32*u+v','u^2-256','32*u+2*v'],
        exhaustion='ac=0; second combination forces c=0; D!=0 and Q2 then give a=+1 or -1 over every geometric extension',
        branches=branch_evidence,ga_embedding_injective='y translation=4u*t/3',
        group='Ga x mu2',E_outside_Ga='E has x,w linear block -I whereas Ga has I',
        group_scheme='The displayed families and their inverse y-translation coordinate are regular; characteristic zero algebraic stabilizers are smooth.',
        intrinsic_coverage='unique triple H, unique vertex of Q|H, unique compatible constant bivector line; affine transformations preserve compatibility and constancy',
        no_special_parameter_exhaustion=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',
           ROOT/'code/q_new_family_generic_stabilizer_20261104_v63.py',
           ROOT/'research/notes/Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63.md',geo,sig,sigreview]
    report=dict(id=ID,all_passed=True,
        independent_unrestricted_affine_flag_scalar_and_divergence_reduction=True,
        no_swapped_line_or_nonrational_quadratic_case_omitted=True,
        independent_all_six_leading_brackets_and_both_quadratic_combinations=True,
        independent_full150_homogeneous_coefficients_for_both_flag_families=True,
        own_rank3_remaining_linear_minors_and_unique_Ga_parameterizations=True,
        sigma_difference_E_square_commutation_and_Ga_exclusion_verified=True,
        actual_generic_full_PGL5_stabilizer='Ga x mu2 over QQ(u,v), C(u,v), and geometric generic fields',
        no_special_parameter_or_boundary_stabilizer_classification=True,
        no_point_prime_GB_solver_search_or_author_producer_replay=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in ('all_passed','actual_generic_full_PGL5_stabilizer','seconds')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
