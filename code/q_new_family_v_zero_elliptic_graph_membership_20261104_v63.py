#!/usr/bin/env python3
"""One exact QQ[u] full-family T39/T40 graph-membership calculation."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

from s61_model import S61,X,J,coefficients,pushforward

ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_ELLIPTIC_GRAPH_MEMBERSHIP_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read(p):return json.loads(p.read_text())


def encode(poly):
    return [[int(mon[0]),int(co.p),int(co.q)]for mon,co in poly.terms()if co]


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    generic=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_GENERIC_KNOWN_EXCLUSIONS_20261104_v63/report.json'
    boundary=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63_REVIEW/report.json'
    globalgate=ROOT/'reports/Q_EXACT_KNOWN_MPS_MEMBERSHIP_20261103_v63.json'
    G=read(geo);assert G['all_passed']and read(generic)['all_passed']and read(boundary)['all_passed']
    assert read(globalgate)['global_reduced_T39_T40_S61_graph_identification_from_existing_universal_proofs']
    inputpaths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),geo,generic,boundary,globalgate,
                ROOT/'code/s61_model.py',
                ROOT/'legacy/v02/data/symplectic_slice_certificate.json',
                ROOT/'legacy/v02/data/symplectic_slice_equations.json',
                ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    graphs=[]
    for component,stem in(('T39','Q51'),('T40','Q52')):
        path=ROOT/f'reports/{stem}_UNIVERSAL_SECTION_20261101.json'
        review=ROOT/f'reports/{stem}_UNIVERSAL_SECTION_REVIEW_20261101.json'
        data=read(path);assert data['all_passed']and read(review)['all_passed']
        graphs.append((component,data,path));inputpaths +=[path,review]
    u=s.Symbol('u');n0,n1,n2=s.symbols('n0 n1 n2');x,y,z,w=s.symbols('x y z w')
    affine=s.Matrix(G['affine_P']).subs({n0:u,n1:1,n2:0});spatial=(x,y,z,w)
    P=s.zeros(5)
    for i in range(4):
        for j in range(i+1,4):
            q=0
            for mon,co in s.Poly(affine[i,j],*spatial).terms():
                q+=co*X[0]**(2-sum(mon))*s.prod(X[k+1]**mon[k]for k in range(4))
            P[i+1,j+1]=s.expand(q);P[j+1,i+1]=-P[i+1,j+1]
    R=s.Matrix([[1,0,0,0],[0,0,0,s.Rational(9,8)],
                [0,1,0,-s.Rational(1,32)],[0,0,-s.Rational(9,2),0]])
    value=P.subs(dict(zip(X,[1,0,0,0,0])))[1:,1:]
    assert R*value*R.T==J[1:,1:]
    assert R.det()==-s.Rational(81,16)
    M=s.eye(5);M[1:,1:]=R
    preflight=dict(field='QQ[u], no source/target point specialization',scope='u!=0,v=0',
        constant_Darboux_R=[[str(q)for q in row]for row in R.tolist()],det_R='-81/16',
        fixed_affine_origin_value=[[str(q)for q in row]for row in value.tolist()],
        historical_normalizer='frozen S61.normalize after this fixed-point value normalization; perspective response matrix constant',
        exact_membership='entire closed T39/T40 global45 graphs, not a discriminant or chosen-frame open',
        single_budget='60s/2GB/32MiB, no new numeric point/prime/high invariant/GB')
    (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    model=S61();normal=model.normalize(pushforward(P,M));values=list(map(s.cancel,normal['parameters']))
    checks=model.verify_point(normal)
    assert checks==dict(reconstruction_coefficients=150,reduced_equations=134,projective_jacobi=True)
    assert all(s.denom(q).is_Rational for q in [*values,*normal['perspective']])
    stages=dict(full61=list(map(str,values)),full150=list(map(str,normal['coefficients'])),
                perspective4=list(map(str,normal['perspective'])),physical_checks=checks)
    (out/'normalized-source.json').write_text(json.dumps(stages,indent=2)+'\n')
    us=s.symbols('u1:62');records=[]
    for component,data,path in graphs:
        one=data['weight1_u_indices_zero_based'];two=data['weight2_u_indices_zero_based']
        assert len(one)==16 and len(two)==45 and sorted(one+two)==list(range(61))
        substitution={us[i]:values[i]for i in one}
        residual=[]
        for index,formula in zip(two,data['quadratic_section_q']):
            expr=s.cancel(values[index]-s.sympify(formula).subs(substitution,simultaneous=True))
            assert s.denom(expr).is_Rational
            residual.append(s.Poly(expr,u,domain=s.QQ))
        nonzero=[p for p in residual if not p.is_zero];assert nonzero
        gcdpoly=nonzero[0].monic()
        prefixes=[]
        for p in nonzero[1:]:
            gcdpoly=s.gcd(gcdpoly,p).monic();prefixes.append(str(gcdpoly.as_expr()))
        valuation=0;away=gcdpoly
        while away.degree()>0 and away.eval(0)==0:
            away=s.exquo(away,s.Poly(u,u,domain=s.QQ));valuation+=1
        assert away.eval(0)!=0
        records.append(dict(component=component,graph_source=str(path.relative_to(ROOT)),
            all45_polynomial_residuals=[encode(p)for p in residual],
            nonzero_ordinals=[i for i,p in enumerate(residual)if not p.is_zero],
            gcd_monic=str(gcdpoly.as_expr()),gcd_encoded=encode(gcdpoly),
            gcd_u_valuation=valuation,gcd_on_u_nonzero=str(away.monic().as_expr()),
            all_nonzero_u_excluded=(away.degree()==0),
            exact_membership_for_nonzero_u='gcd_on_u_nonzero(u)=0',
            gcd_prefixes=prefixes))
    ep=out/'graph-residuals.json';ep.write_text(json.dumps(records,indent=2)+'\n')
    excluded=all(q['all_nonzero_u_excluded']for q in records)
    report=dict(id=ID,all_passed=True,constant_value_frame_and_historical_S61_gauge_bound=True,
        full150_reconstruction_original134_and_projective_Jacobi_exact_QQ_u=True,
        all61_and_all90_graph_residuals_polynomial_without_variable_denominators=True,
        gcd_results={q['component']:q['gcd_monic']for q in records},
        gcd_results_after_excluding_u_zero={q['component']:q['gcd_on_u_nonzero']for q in records},
        all_u_nonzero_outside_T39_and_T40=excluded,
        actual_membership_equivalence_from_entire_closed_global_graphs=True,
        u_zero_type_and_membership_not_claimed=True,
        no_new_numeric_point_prime_GB_Inv26_or_automatic_budget_extension=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputpaths},
        evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},
        seconds=time.monotonic()-start,frozen=True,independent_review_pending=True)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','gcd_results','gcd_results_after_excluding_u_zero','all_u_nonzero_outside_T39_and_T40','seconds')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
