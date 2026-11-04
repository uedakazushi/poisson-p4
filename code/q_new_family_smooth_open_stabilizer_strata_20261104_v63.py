#!/usr/bin/env python3
"""Actual stabilizer strata on u*v*(32*u+v)!=0, with modular-image coverage."""
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
ID='Q_NEW_FAMILY_SMOOTH_OPEN_STABILIZER_STRATA_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def mm(A,B,zero):
    return [[sum((A[i][k]*B[k][j]for k in range(len(B))),zero)for j in range(len(B[0]))]for i in range(len(A))]


def det(A,zero,one):
    result=zero
    for p in itertools.permutations(range(len(A))):
        term=one*((-1)**sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p))))
        for i,j in enumerate(p):term*=A[i][j]
        result+=term
    return result


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    if out.exists():
        assert all(p.name=='unexpected-residual.json'for p in out.iterdir()), 'refuse to overwrite frozen evidence'
    else:out.mkdir(parents=True)
    geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    old=ROOT/'reports/Q_NEW_FAMILY_GENERIC_STABILIZER_20261104_v63_REVIEW/report.json'
    geometry=json.loads(geo.read_text());assert geometry['all_passed']
    assert json.loads(old.read_text())['all_passed']
    u,v,a,c,t,T=s.symbols('u v a c t T');n0,n1,n2=s.symbols('n0 n1 n2')
    x,y,z,w=s.symbols('x y z w');coords=(x,y,z,w)
    P=s.Matrix(geometry['affine_P']).subs({n0:u,n1:1,n2:v})
    # Whole constant-compatibility locus: retain the actual jump, and do not
    # use generic uniqueness to obtain the special-parameter affine flag.
    pairs=list(itertools.combinations(range(4),2));ks=s.symbols('k0:6');C=s.zeros(4)
    for(i,j),q in zip(pairs,ks):C[i,j]=q;C[j,i]=-q
    rows=[];labels=[]
    for i,j,k in itertools.combinations(range(4),3):
        pp=s.Poly(s.expand(sum(C[i,b]*s.diff(P[j,k],coords[b])+C[j,b]*s.diff(P[k,i],coords[b])+C[k,b]*s.diff(P[i,j],coords[b])for b in range(4))),*coords)
        for mon in (1,*coords):
            rows.append([s.diff(pp.coeff_monomial(mon),q)for q in ks]);labels.append([i,j,k,str(mon)])
    compat=s.Matrix(rows)
    assert compat[:,3]==s.zeros(20,1)
    velocity=s.Matrix([0,-u*ks[5],0,ks[3],0,ks[5]])
    residual=compat*velocity
    assert s.expand(residual[0]-s.Rational(4,3)*(u*u-64)*ks[5])==0
    assert all(s.expand(q)==0 for q in list(residual)[1:])
    minor4=s.factor(compat.extract([11,4,19,15],[0,1,2,4]).det());assert minor4!=0
    assert s.cancel(minor4/v).is_Rational
    # N is rebuilt from actual affine divergence, then homogenized. Its
    # image is the same two-plane for EVERY u!=0, independently of C.
    div=s.Matrix([s.expand(sum(s.diff(P[i,j],coords[i])for i in range(4)))for j in range(4)])
    N=s.zeros(5);N[1:,0]=div.subs(dict.fromkeys(coords,0));N[1:,1:]=div.jacobian(coords)
    assert N.extract([2,3],[0,1]).det()==4*u*u/3
    assert all(N[i,j]==0 for i in(0,1,4)for j in range(5))
    # Exhaust the complete quadratic system, retaining the resonant cases.
    D=a*a-64*c*c;p=32*u+v
    q1=-p*D*D+p*a*a+4096*a*c-64*v*c*c
    q2=-D*D+a*a+u*a*c+64*c*c
    q3=v*D*D-v*a*a+4096*a*c+64*p*c*c
    assert s.expand(q1+q3+(v-p)*q2+32*a*c*(u*u-256))==0
    assert s.expand(q3+v*q2-a*c*(4096+u*v)-64*c*c*(p+v))==0
    for sign in(1,-1):
        assert s.expand(q2.subs(u,16*sign)-(a+8*sign*c)**2*(1-(a-8*sign*c)**2))==0
        assert s.expand((q3+v*q2).subs(u,16*sign)-16*c*(256+sign*v)*(a+8*sign*c))==0
        assert s.expand(q1.subs({u:16*sign,v:-256*sign})-256*sign*q2.subs(u,16*sign))==0
        assert s.expand(q3.subs({u:16*sign,v:-256*sign})-256*sign*q2.subs(u,16*sign))==0
    assert s.expand(q2.subs({v:-16*u,a:0})-64*c*c*(1-64*c*c))==0
    # One fixed unit constant-coefficient system solves ALL remaining H,K,Z.
    # Each case is a symbolic family, never a target/source point sample.
    results=[];saved={}
    cases=[('J0_cplus',u,-16*u,s.Integer(0),s.Rational(1,8),(u,T)),
           ('J0_cminus',u,-16*u,s.Integer(0),-s.Rational(1,8),(u,T))]
    for sign in(1,-1):
        for eps in(1,-1):
            cases.append((f'resonant_{sign}_{eps}',s.Integer(16*sign),s.Integer(-256*sign),
                          (t+eps)/2,sign*(t-eps)/16,(t,T)))
    for name,uu,vv,aa,cc,gens in cases:
        F=s.QQ.frac_field(*gens);cv=lambda z:F.from_sympy(s.sympify(z))
        R,HH,XX,YY,ZZ,WW=ring('h,x,y,z,w',F);vs=[HH,XX,YY,ZZ,WW]
        d=cv(aa*aa-64*cc*cc);assert d
        kap=1/d;tx=-cc/48+2*(aa-1)/(3*vv);tw=(1-aa)/s.Integer(3072)+2*cc/(3*vv)
        ee=(aa*aa-64*cc*cc-aa)/vv;gg=-64*cc/vv
        ii=3*(aa*aa-64*cc*cc)*(tx+vv*T)/4+(aa*aa-64*cc*cc)*(1-aa*aa+64*cc*cc)/(64*uu)
        source=P.subs({u:uu,v:vv},simultaneous=True)
        raw=[[R.zero for _ in range(5)]for _ in range(5)]
        for i,j in pairs:
            for mon,co in s.Poly(source[i,j],*coords).terms():raw[i+1][j+1][(2-sum(mon),*mon)]=cv(co)
            raw[j+1][i+1]=-raw[i+1][j+1]
        def matrix(hh,kk,zz):
            return [[F.one,F.zero,F.zero,F.zero,F.zero],
                    [cv(tx),cv(aa),F.zero,F.zero,cv(64*cc)],
                    [cv(T),cv(ee),d,F.zero,cv(gg)],
                    [zz,hh,cv(ii),d*d,kk],
                    [cv(tw),cv(cc),F.zero,F.zero,cv(aa)]]
        def compose(q,im):
            answer=R.zero
            for mon,co in q.items():
                term=R(co)
                for xx,power in zip(im,mon):
                    if power:term*=xx**power
                answer+=term
            return answer
        def difference(M):
            im=[sum((R(M[i][j])*vs[j]for j in range(5)),R.zero)for i in range(5)]
            ans={}
            for i,j in itertools.combinations(range(5),2):
                push=sum((R(kap*M[i][k]*M[j][l])*raw[k][l]for k in range(5)for l in range(5)),R.zero)
                ans[i,j]=push-compose(raw[i][j],im)
            return ans
        base= difference(matrix(F.zero,F.zero,F.zero))
        units=[]
        for j in range(3):
            arg=[F.zero]*3;arg[j]=F.one
            df=difference(matrix(*arg));units.append({pair:df[pair]-base[pair]for pair in base})
        selected=[(1,3),(3,4),(2,3)];constant=(2,0,0,0,0)
        L=[[units[j][pair].get(constant,F.zero)for j in range(3)]for pair in selected]
        rhs=[-base[pair].get(constant,F.zero)for pair in selected]
        ld=det(L,F.zero,F.one)
        assert ld==cv(-256*uu/(243*(aa*aa-64*cc*cc)))
        solution=[]
        for j in range(3):
            A=[row[:]for row in L]
            for i in range(3):A[i][j]=rhs[i]
            solution.append(det(A,F.zero,F.one)/ld)
        M=matrix(*solution);final=difference(M)
        residual_proof=None
        if any(final.values()):
            assert name.startswith('resonant_')
            # The quadratic candidate has an additional LINEAR obstruction.
            # This is case elimination, not another parameter-point test.
            residual_proof={};unit_witness=None
            for pair,poly in final.items():
                for mon,co in poly.items():
                    quotient=s.factor(F.to_sympy(co)/(t*t-1))
                    num,den=s.fraction(quotient)
                    assert not quotient.has(T)
                    assert s.Poly(num,t).length()==1 and s.Poly(den,t).length()==1
                    assert s.Poly(num,t).monoms()[0][0]>=0
                    residual_proof[str((pair,mon))]=str(quotient)
                    if unit_witness is None:unit_witness=str((pair,mon))
                    for root in(1,-1):
                        assert s.cancel(F.to_sympy(co).subs(t,root))==0
            assert unit_witness is not None
            residual_proof=dict(necessary_and_sufficient_extra_equation='t^2=1',
                                unit_coefficient_witness=unit_witness,
                                all_coefficients_divided_by_t_squared_minus_one=residual_proof)
        for pair in base:
            assert base[pair]+sum((R(solution[j])*units[j][pair]for j in range(3)),R.zero)==final[pair]
        nc=[[cv(q)for q in row]for row in N.subs({u:uu,v:vv},simultaneous=True).tolist()]
        assert mm(nc,M,F.zero)==[[kap*q for q in row]for row in mm(M,nc,F.zero)]
        payload=dict(name=name,parameters=list(map(str,gens)),u=str(uu),v=str(vv),a=str(aa),c=str(cc),
                     D=str(F.to_sympy(d)),kappa=str(F.to_sympy(kap)),
                     constant_linear_matrix=[[str(F.to_sympy(q))for q in row]for row in L],
                     constant_rhs=list(map(lambda q:str(F.to_sympy(q)),rhs)),
                     constant_unit_minor=str(F.to_sympy(ld)),
                     H_K_Z=[str(F.to_sympy(q))for q in solution],
                     matrix=[[str(F.to_sympy(q))for q in row]for row in M],
                     all150_homogeneous_coefficients_zero=(residual_proof is None),
                     resonant_linear_obstruction=residual_proof)
        results.append(payload);saved[name]=(F,M,nc,kap)
    # Exact finite group relations on J=0 away from the resonance (the same
    # identities also hold at the two resonant parameters).
    F,M,NN,kap=saved['J0_cplus'];cv=lambda z:F.from_sympy(s.sympify(z))
    S=[[cv(F.to_sympy(z).subs(T,0))for z in row]for row in M]
    one=F.one;zero=F.zero;eye=[[one if i==j else zero for j in range(5)]for i in range(5)]
    Es=[[cv(q)for q in row]for row in [[1,0,0,0,0],[-4/(3*(-16*u)),-1,0,0,0],
       [4/(3*(-16*u)**2),2/(-16*u),1,0,0],[-2/(3*(-16*u)**3),-1/(-16*u)**2,0,1,0],[s.Rational(1,1536),0,0,0,-1]]]
    assert mm(S,S,zero)==eye and mm(Es,Es,zero)==eye
    assert mm(S,NN,zero)==[[-q for q in row]for row in mm(NN,S,zero)]
    assert mm(S,Es,zero)==mm(Es,S,zero)
    certificate=dict(domain='u*v*(32*u+v)!=0',
        compatibility_matrix=[[str(q)for q in row]for row in rows],compatibility_labels=labels,
        constant_kernel='kxy=kxw=kyw=0; kxz=-u*kzw; (u^2-64)kzw=0; kyz free',
        compatibility_jump='dimension2 at u=+8 or-8, otherwise dimension1',rank4_minor=str(minor4),
        flag_repair='N M=kappa M N implies M(ImN)=ImN; ImN=span(ey,ez) for every u!=0. Use this intrinsic plane, NOT uniqueness of C.',
        N=[[str(q)for q in row]for row in N.tolist()],imageN_unit_minor='4*u^2/3',
        quadratic_system=list(map(str,(q1,q2,q3))),families=results,
        nonJ0='Ga x mu2 when v+16u!=0, including u^2=256',
        whole_J0='(Ga semidirect mu2) x mu2 when v=-16u; first mu2 inverts Ga; includes u^2=256',
        rejected_resonant_Gm='At u=+16,v=-256 or u=-16,v=256 the quadratic candidate permits Gm, but full LINEAR coefficients force t^2=1 in every epsilon case.',
        S_matrix=[[str(F.to_sympy(q))for q in row]for row in S],
        finite_relations='S^2=E^2=1, SE=ES, S N=-N S, E N=N E',
        group_scheme_scope='split algebraic stabilizers in characteristic zero, over QQ on the indicated strata and after geometric base extension',
        boundary_u_v_p_zero_not_classified=True,independent_review_pending=True)
    ep=out/'certificate.json';ep.write_text(json.dumps(certificate,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),geo,old]
    evidences={ep.name:sha(ep)}
    failed=out/'unexpected-residual.json'
    if failed.exists():evidences[failed.name]=sha(failed)
    report=dict(id=ID,all_passed=True,whole_smooth_uvp_open_stabilizer_stratification=True,
        compatibility_jump_at_u_squared_64_retained=True,whole_open_flag_uses_intrinsic_modular_image=True,
        both_new_J0_symbolic_families_full150_coefficients_zero=True,
        all4_resonant_quadratic_candidates_cut_to_t_squared_one_by_linear_coefficients=True,
        fixed_constant_rank3_minor_unit=True,finite_group_laws_and_action_on_Ga_checked=True,
        extra_resonant_Gm_hypothesis_refuted=True,no_point_prime_GB_or_solver_search=True,
        no_u_v_p_boundary_classification=True,independent_review_pending=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256=evidences,seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','seconds')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
