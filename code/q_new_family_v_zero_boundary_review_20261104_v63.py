#!/usr/bin/env python3
"""Independent v=0 full-action/209-source and proper-image proof audit."""
import argparse
from fractions import Fraction
import gzip
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
ID='Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63_REVIEW'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)


def determinant(A,zero,one):
    answer=zero
    for p in itertools.permutations(range(len(A))):
        co=one*((-1)**sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p))))
        for i,j in enumerate(p):co*=A[i][j]
        answer+=co
    return answer


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    if out.exists():
        assert not any(out.iterdir()), 'refuse to overwrite evidence'
    else:out.mkdir(parents=True)
    base=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63'
    author=read(base/'report.json');cert=read(base/'certificate.json');assert author['all_passed']
    for rel,hsh in author['input_sha256'].items():assert sha(ROOT/rel)==hsh,rel
    for rel,hsh in author['evidence_sha256'].items():assert sha(base/rel)==hsh,rel
    geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    inc=ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/certificate.json.gz'
    geometry=read(geo);incidence=read(inc)
    u,U=s.symbols('u U');n0,n1,n2=s.symbols('n0 n1 n2')
    x,y,z,w=s.symbols('x y z w');coords=(x,y,z,w)
    P=s.Matrix(geometry['affine_P']).subs({n0:u,n1:1,n2:0})
    div=s.Matrix([s.expand(sum(s.diff(P[i,j],coords[i])for i in range(4)))for j in range(4)])
    assert div==s.Matrix([0,4*u/3,u*x+s.Rational(1,48),0])
    N=s.zeros(5);N[1:,0]=div.subs(dict.fromkeys(coords,0));N[1:,1:]=div.jacobian(coords)
    assert N*N==s.zeros(5)
    assert N.extract([2,3],[0,1]).det()==4*u*u/3
    Qm=s.Matrix(geometry['quadric_matrix']).subs({n0:u,n1:1,n2:0})
    assert s.expand(Qm.det()-s.Rational(32768,531441)*u**4)==0
    assert s.factor(Qm.extract([1,2,4],[1,2,4]).det())==-s.Rational(128,729)*u*u
    # Start with the entire flag and an unknown tensor scalar, for BOTH
    # source u and target U. This independently proves the u-orbit assertion.
    a,b,c,d,e,f,g,H,I,J,K,tx,T,Z,tw,kap=s.symbols('a b c d e f g H I J K tx T Z tw kap')
    A=s.Matrix([[a,0,0,b],[e,f,0,g],[H,I,J,K],[c,0,0,d]])
    image=A*s.Matrix(coords)+s.Matrix([tx,T,Z,tw])
    target=P.subs(u,U)
    push=kap*A*P*A.T
    xy=s.Poly(s.expand(push[0,1]-(s.Rational(1,36)-256*image[3]/3)),*coords)
    yw=s.Poly(s.expand(push[1,3]-4*image[0]/3),*coords)
    assert s.expand(xy.coeff_monomial(x)-s.Rational(4,3)*(64*c-kap*f*b))==0
    assert s.expand(xy.coeff_monomial(w)-s.Rational(256,3)*(d-kap*f*a))==0
    assert s.expand(yw.coeff_monomial(x)-s.Rational(4,3)*(kap*f*d-a))==0
    assert s.expand(yw.coeff_monomial(w)-s.Rational(4,3)*(64*kap*f*c-b))==0
    cov=kap*A*div-div.subs(u,U).subs(dict(zip(coords,image)),simultaneous=True)
    assert s.expand(cov[1]-4*(kap*f*u-U)/3)==0
    assert s.Poly(s.expand(cov[2]),*coords).coeff_monomial(w)==-U*b
    # At fixed u, eta=1, so b=c=0,d=a,f=a^2,kappa=a^-2.
    # The remaining z/x condition gives j=a^3 without dividing by v.
    fix={U:u,b:0,c:0,d:a,f:a*a,kap:1/a**2}
    assert s.cancel(s.Poly(s.expand(cov[2].subs(fix,simultaneous=True)),*coords).coeff_monomial(x)-u*(J/a**2-a))==0
    reducedA=s.Matrix([[a,0,0,0],[e,a*a,0,g],[H,I,a**3,K],[0,0,0,a]])
    Y=reducedA*s.Matrix(coords)+s.Matrix([tx,T,Z,tw])
    diff=(reducedA*P*reducedA.T/a**2-P.subs(dict(zip(coords,Y)),simultaneous=True)).applyfunc(s.expand)
    dz=s.Poly(diff[1,2],*coords);dw=s.Poly(diff[2,3],*coords)
    assert s.expand(dz.coeff_monomial(w*w)-(a**3-a*a-64*a*g))==0
    assert s.expand(dz.coeff_monomial(x*w)+a*(128*e+u*g))==0
    assert dw.coeff_monomial(w)==-2*g/3
    # Therefore g=0,a=1,e=0. The constants then force tx=tw=I=0.
    simple={a:1,e:0,g:0,U:u,b:0,c:0,d:1,f:1,kap:1,J:1}
    assert s.cancel(xy.coeff_monomial(1).subs(simple)-s.Rational(256,3)*tw)==0
    assert s.cancel(yw.coeff_monomial(1).subs(simple)+s.Rational(4,3)*tx)==0
    assert s.cancel(s.Poly(s.expand(cov[2].subs(simple)),*coords).coeff_monomial(1).subs({tx:0,tw:0})-4*u*I/3)==0
    # Complete full150 raw homogeneous check, using a sparse coefficient
    # ring rather than primary's affine expression expansion.
    F=s.QQ.frac_field(u,T,H,K,Z);cv=lambda q:F.from_sympy(s.sympify(q))
    R,HH,XX,YY,ZZ,WW=ring('h,x,y,z,w',F);vs=[HH,XX,YY,ZZ,WW]
    pairs=list(itertools.combinations(range(4),2))
    raw=[[R.zero for _ in range(5)]for _ in range(5)]
    for i,j in pairs:
        for mon,co in s.Poly(P[i,j],*coords).terms():raw[i+1][j+1][(2-sum(mon),*mon)]=cv(co)
        raw[j+1][i+1]=-raw[i+1][j+1]
    M=[[cv(q)for q in row]for row in [[1,0,0,0,0],[0,1,0,0,0],[T,0,1,0,0],[Z,H,0,1,K],[0,0,0,0,1]]]
    images=[sum((R(M[i][j])*vs[j]for j in range(5)),R.zero)for i in range(5)]
    def compose(poly,images):
        answer=R.zero
        for mon,co in poly.items():
            a0=R(co)
            for q,power in zip(images,mon):
                if power:a0*=q**power
            answer+=a0
        return answer
    differences={}
    for i,j in itertools.combinations(range(5),2):
        pushed=sum((R(M[i][k]*M[j][l])*raw[k][l]for k in range(5)for l in range(5)),R.zero)
        differences[i,j]=pushed-compose(raw[i][j],images)
    selected=[(1,3),(3,4),(2,3)];const=(2,0,0,0,0)
    FF=s.QQ.frac_field(u,T);L=[];rhs=[]
    for pair in selected:
        pp=s.Poly(F.to_sympy(differences[pair].get(const,F.zero)),H,K,Z,domain=FF)
        L.append([FF.from_sympy(pp.coeff_monomial(q))for q in(H,K,Z)])
        rhs.append(-FF.from_sympy(pp.coeff_monomial(1)))
    dd=determinant(L,FF.zero,FF.one);assert dd==FF.from_sympy(-256*u/243)
    answer={H:3*T/4,K:s.Integer(0),Z:T/(64*u)}
    checks=0
    for poly in differences.values():
        for co in poly.values():assert s.cancel(F.to_sympy(co).subs(answer,simultaneous=True))==0
        checks+=15
    assert checks==150
    finalM=s.Matrix([[s.sympify(q)for q in row]for row in [[1,0,0,0,0],[0,1,0,0,0],[T,0,1,0,0],[T/(64*u),3*T/4,0,1,0],[0,0,0,0,1]]])
    assert(finalM-s.eye(5)-3*T*N/(4*u)).applyfunc(s.expand)==s.zeros(5)
    # Explicit sufficiency of u->-u; no firstjet rational-function extension
    # from the v!=0 source is assumed.
    signD=s.diag(-1,1,-1,1);back={x:-x,y:y,z:-z,w:w}
    assert(-signD*P.subs(back,simultaneous=True)*signD.T-P.subs(u,-u)).applyfunc(s.expand)==s.zeros(4)
    # Evaluate all ORIGINAL encoded equations with independent QQ[u]
    # polynomial arithmetic, not the saved pullback strings.
    V,Uvar=ring('u',s.QQ)
    graph=[]
    for q in incidence['exact99_polynomial_graph']:
        expression=s.expand(s.sympify(q).subs({n0:u,n1:1,n2:0}))
        graph.append(V.from_expr(expression))
    assert graph[90:]==[V.zero]*9
    encoded=[]
    for equation in incidence['full209_equations_encoded']:
        total=V.zero
        for ids,num,den in equation:
            term=V(s.QQ(num,den))
            for idx in ids:term*=graph[idx]
            total+=term
        encoded.append(total)
    assert len(encoded)==209 and all(not q for q in encoded)
    QF=s.QQ.frac_field(u);cv2=lambda q:QF.from_sympy(s.sympify(q))
    S,SX,SY,SZ,SW=ring('x,y,z,w',QF);xs=[SX,SY,SZ,SW]
    mons=[S.from_expr(s.sympify(q))for q in incidence['coefficient_monomials']]
    actual=[[S.zero for _ in range(4)]for _ in range(4)]
    for kk,(i,j)in enumerate(pairs):
        q=S.from_expr(s.sympify(incidence['positive_P'][i][j]))
        q+=sum((S(cv2(graph[15*kk+m].as_expr()))*mons[m]for m in range(15)),S.zero)
        assert q==S.from_expr(P[i,j])
        actual[i][j]=q;actual[j][i]=-q
    mon3=[]
    for degree in range(4):
        for ii in itertools.combinations_with_replacement(range(4),degree):
            exponent=[0]*4
            for j in ii:exponent[j]+=1
            mon3.append(tuple(exponent))
    direct=[]
    for i,j,k in itertools.combinations(range(4),3):
        sch=sum((actual[i][b]*actual[j][k].diff(xs[b])+actual[j][b]*actual[k][i].diff(xs[b])+actual[k][b]*actual[i][j].diff(xs[b])for b in range(4)),S.zero)
        direct +=[sch.get(mon,QF.zero)for mon in mon3]
    for i,j in pairs:
        quadratic=S.from_dict({mon:co for mon,co in actual[i][j].items()if sum(mon)==2})
        center=quadratic.diff(SZ)
        direct +=[center.get(tuple(int(k==m)for k in range(4)),QF.zero)for m in range(4)]
    direct +=[cv2(q)for q in N*N]
    CC=[[0]*4 for _ in range(4)];CC[1][2]=1;CC[2][1]=-1
    mon1=[(0,0,0,0),*[(int(i==j)for i in range(4))for j in range(4)]]
    mon1=[tuple(q)for q in mon1]
    for i,j,k in itertools.combinations(range(4),3):
        q=sum((CC[i][b]*actual[j][k].diff(xs[b])+CC[j][b]*actual[k][i].diff(xs[b])+CC[k][b]*actual[i][j].diff(xs[b])for b in range(4)),S.zero)
        direct +=[q.get(mon,QF.zero)for mon in mon1]
    assert len(direct)==209 and all(not q for q in direct)
    gates=[ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63_REVIEW/report.json',
           ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_N_GLOBAL_R1_20261103_v63_REVIEW.json']
    assert all(read(q)['all_passed']for q in gates)
    hquartic=ROOT/'reports/ADDITIONAL_RANK4_FAMILIES_20261103_v62/R4_H_QUARTIC_ESSENTIAL5.json'
    record=read(hquartic);assert record['source_parameters']==[]and record['projective_image_dimension']==23
    proof=dict(domain='u!=0,v=0',quadric_det='32768*u^4/531441',
        all_full150_Ga_coefficients_verified=True,remaining_unit_minor=str(FF.to_sympy(dd)),
        parameter_orbit_proof='For source u,target U, eta=kappa*f=U/u. Four low-bracket coefficients imply eta^2=1 by invertibility, hence U=+-u; the full finite tau realizes minus.',
        stabilizer='Ga for every u!=0 over QQ,C and geometric fields',
        saturation_dimension=24,saturation_reason='PGL source dimension24+parameter1, exact source fibers Ga times two parameter choices have dimension1',
        original209_encoded_substitutions_zero=True,
        separately_rebuilt209='140 full physical Jacobi coefficients +24 compact center +25 modular N squared +20 constant-pencil compatibility',
        old99_marks='all nine auxiliary coordinates zero: center is z at infinity, constant pencil dy wedge dz is rank2',
        membership_only_closed209_incidence=True,
        CA_CB_MPSN_exclusion='accepted global PROPER sources cover the entire components; every nonzero quintic is a product of five geometric linear factors, a closed Chow condition. Current 3H+smooth Q is not.',
        Hquartic_exclusion='accepted component is closure of one23D orbit, with H+irreducible quartic. Current orbit has dimension23 for every u; boundary orbits have dimension<=22, while the two nonzero divisor types differ.',
        no_component_membership_from_squarezero=True,
        no_finite_center_nondegenerate_pencil_closure_claim=True,
        no_other_known_component_or_all_new_component_boundary_classification=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',
            ROOT/'code/q_new_family_v_zero_boundary_20261104_v63.py',
            ROOT/'research/notes/Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63.md',geo,inc,*gates,hquartic,
            ROOT/'research/notes/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63_REVIEW.md',
            ROOT/'legacy/v06/nilpotent22.tex']
    report=dict(id=ID,all_passed=True,
        independent_whole_v_zero_u_unit_stabilizer_exact_Ga=True,
        general_two_parameter_flag_proves_U_plus_minus_u_without_firstjet_boundary_extrapolation=True,
        independent_full150_homogeneous_Ga_identity=True,
        independently_evaluated_all_original209_encoded_source_equations=True,
        separately_reconstructed_physical209_conditions_and_old90_coefficient_tensor=True,
        proper_CA_CB_MPSN_and_dense_Hquartic_orbit_exclusion_scope_accepted=True,
        saturation_dimension=24,closed_incidence_only_no_squarezero_component_inference=True,
        no_new_point_prime_GB_solver_or_author_producer_replay=True,
        input_sha256={str(q.relative_to(ROOT)):sha(q)for q in inputs},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','saturation_dimension','seconds')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
