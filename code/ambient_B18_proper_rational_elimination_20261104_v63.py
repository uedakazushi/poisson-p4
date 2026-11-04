#!/usr/bin/env python3
"""ONE explicit unit-rational five-pivot elimination and twelve flag quotients."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encpoly(p,vars):
    if p==0:return []
    return [[list(m),int(s.numer(c)),int(s.denom(c))] for m,c in s.Poly(p,*vars).terms() if c]
def expression(terms,vars):return sum((s.Rational(n,d)*s.prod(v**k for v,k in zip(vars,m)) for m,n,d in terms),s.Integer(0))
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    pp=ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63/actual-proper-source14-equations17.json';pr=pp.parent/'report.json'
    P=read(pp);assert read(pr)['all_passed']
    v=s.symbols('v0:9');a,b,c,d,e=s.symbols('a b c d e');vars=tuple(v)+(a,b,c,d,e);R=(a,d)+tuple(v[2:]);D=1-d*d
    eq=[expression(q,vars) for q in P['equations14variable_sparse']]
    assert len(eq)==17 and P['five_unit_pivot_rows']==[0,1,2,3,4]
    V0=2*d-a;B=d*(a+d-v[2])/2
    EN=s.expand(B*B+d*d*B/2-a*d*d*(d-v[2])-d*d*v[4]+d*v[6])
    VN=s.expand(B*D/2+2*EN);CN=s.expand(EN*(a+d+v[2])/2-d*v[3]*D)
    numer={v[0]:V0,v[1]:VN,b:B,c:CN,e:EN};dpow={v[0]:0,v[1]:1,b:0,c:1,e:1}
    sub={q:numer[q]/D**dpow[q] for q in numer}
    # Guard BEFORE powers of substitution numerators in the source12 terms.
    rootcounts={q:len(s.Poly(numer[q],*R).terms()) for q in numer}
    clearing=[];bounds=[]
    for idx in range(5,17):
        terms=P['equations14variable_sparse'][idx]
        m=max((sum(exp[j]*dpow.get(q,0) for j,q in enumerate(vars)) for exp,n,de in terms),default=0);clearing.append(m)
        bound=0
        for exp,n,de in terms:
            dexp=sum(exp[j]*dpow.get(q,0) for j,q in enumerate(vars))
            count=2**(m-dexp)
            for j,q in enumerate(vars):
                if q in numer:count*=rootcounts[q]**exp[j]
            bound+=count
        bounds.append(bound)
    totalbound=sum(bounds)
    guard={'root_numerator_term_counts':{str(q):rootcounts[q] for q in numer},
      'D_clearing_powers12':clearing,'conservative_remaining12_expansion_term_bounds':bounds,
      'total_conservative_remaining12_expansion_terms':totalbound,'term_cap':150000,
      'only_new_denominators':['2','1-d^2'],'D_origin':'1'}
    (out/'expansion-guard.json').write_text(json.dumps(guard,indent=2)+'\n')
    assert totalbound<=150000,('STOP before expansion: certified bound exceeds budget',totalbound)
    # A triangular quotient proof, not merely a point/substitution check.
    stages=[(2,v[0],V0,s.Rational(1,3)),(0,b,B,s.Integer(2)),
      (4,v[1],B/2+2*e,s.Integer(2)),(3,e,EN/D,-D),(1,c,CN/D,s.Integer(2))]
    prior={};identities=[]
    for idx,q,rhs,lead in stages:
        reduced=eq[idx].subs(prior,simultaneous=True)
        assert s.cancel(reduced-lead*(q-rhs))==0
        identities.append({'source_equation_row':idx,'pivot':str(q),'root':str(rhs),'leading_unit':str(lead),
          'prior_substitutions':{str(k):str(z) for k,z in prior.items()},'exact_reduced_identity_verified':True})
        prior[q]=rhs
    assert all(s.cancel(q.subs(sub,simultaneous=True))==0 for q in eq[:5])
    records=[];quotients=[];allterms=0
    for k,idx in enumerate(range(5,17)):
        m=clearing[k];expanded=0
        for ex,num,den in P['equations14variable_sparse'][idx]:
            term=s.Rational(num,den);de=0
            for j,q in enumerate(vars):
                power=ex[j]
                if not power:continue
                if q in numer:term*=numer[q]**power;de+=dpow[q]*power
                else:term*=q**power
            expanded+=term*D**(m-de)
        expanded=s.expand(expanded)
        assert s.cancel(expanded-D**m*eq[idx].subs(sub,simultaneous=True))==0
        poly=s.Poly(expanded,*R);constant=s.ilcm(*[s.denom(z) for z in poly.coeffs()]) if not poly.is_zero else s.Integer(1)
        integer=s.Poly(s.expand(constant*expanded),*R)
        assert all(z.q==1 for z in integer.coeffs())
        assert integer.as_expr().subs(d,0)==0
        divided=integer.exquo(s.Poly(d,*R)) if not integer.is_zero else s.Poly(0,*R)
        content=s.igcd(*[abs(int(z)) for z in divided.coeffs()]) if not divided.is_zero else s.Integer(1)
        H=s.expand(divided.as_expr()/content);assert s.cancel(eq[idx].subs(sub,simultaneous=True)-content*d*H/(constant*D**m))==0
        assert H.subs(dict.fromkeys(R,0))==0
        numterms=encpoly(integer.as_expr(),R);hterms=encpoly(H,R);allterms+=len(numterms)+len(hterms)
        records.append({'original_equation_row':idx,'D_denominator_power':m,'constant_denominator':int(constant),
          'primitive_content_unit':int(content),'cleared_numerator':numterms,'divided_primitive_polynomial':hterms,
          'whole_source_identity':'Fi(sub)=content*d*H/(constant*D^power)','all_source_and_exact_d_division_identities_verified':True})
        quotients.append(H)
    rootsenc={str(q):{'numerator':encpoly(numer[q],R),'D_denominator_power':dpow[q]} for q in numer}
    allterms+=sum(len(r['numerator']) for r in rootsenc.values());assert allterms<=150000
    origin=dict.fromkeys(R,0);J=s.Matrix([[s.diff(q,z).subs(origin) for z in R] for q in quotients]);rank=J.rank()
    cols=list(J.rref()[1]);rows=list(J[:,cols].T.rref()[1]);minor=s.factor(J[rows,cols].det()) if rank else s.Integer(1)
    # Candidate pivots are only reported. Do not perform another elimination.
    pivotnames=[str(R[j]) for j in cols]
    nextdegree=[s.Poly(quotients[i],*[R[j] for j in cols]).total_degree() for i in rows] if rank else []
    proof={'original_source_file':str(pp.relative_to(ROOT)),'original14_variable_order':list(map(str,vars)),
      'remaining9_variable_order':list(map(str,R)),'unit_D':'1-d^2','unit_D_origin':'1',
      'five_explicit_rational_roots':rootsenc,'five_triangular_unit_identities':identities,
      'all5_original_selected_equations_substitute_to_zero':True,'new_unknown_denominator_inverses_used':False,
      'expansion_guard':guard,'all12_actual_remaining_and_divided_records':records,
      'whole_completed_E4_presentation':'QQ[[a,d,v2,...,v8]]/(d*H0,...,d*H11)',
      'localized_selected5_quotient_identity':'After inverting onlyD=1-d² and nonzeroQQconstants, selected5 ideal is EXACT graph of the five rational roots by triangular unit elimination.',
      'origin_Jacobian12_by9':matrix_strings(J),'origin_rank':rank,'candidate_unit_pivot_rows':rows,
      'candidate_unit_pivot_columns':cols,'candidate_unit_pivot_minor':str(minor),'candidate_pivot_names':pivotnames,
      'candidate_pivot_distribution':{'flag_a_d':sum(j<2 for j in cols),'CA_fiber':sum(j>=2 for j in cols)},
      'candidate_equation_degrees_in_joint_pivots':nextdegree,
      'no_next_elimination_GB_localstandardbasis_CM_or_D0_identification':True,
      'all_stored_expanded_terms':allterms}
    p=out/'actual12-flag-divided-polynomials.json';p.write_text(json.dumps(proof,separators=(',',':'))+'\n');assert p.stat().st_size<=32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pp,pr]
    rep={'id':STEM,'all_passed':True,'status':'exact_unit_rational5_elimination_and12_actual_d_divided_polynomials',
      'all5_triangular_and_substitution_identities_verified':True,'only_new_variable_denominator':'1-d^2',
      'conservative_expansion_terms_before_product':totalbound,'actual_expanded_stored_terms':allterms,
      'remaining_polynomials':12,'remaining_variables':9,'common_exact_flag_factor':'d',
      'divided_origin_rank':rank,'candidate_next_pivot_distribution':proof['candidate_pivot_distribution'],
      'candidate_next_pivots':pivotnames,'candidate_next_minor':str(minor),
      'no_next_elimination_GB_higherjets_D0_CM_or_depth_claim':True,'output_bytes':p.stat().st_size,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{p.name:sha(p),'expansion-guard.json':sha(out/'expansion-guard.json')},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
def matrix_strings(M):return [[str(c) for c in row] for row in M.tolist()]
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
