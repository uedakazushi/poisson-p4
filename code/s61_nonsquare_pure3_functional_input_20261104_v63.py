#!/usr/bin/env python3
"""Known-lambda divisions/source binding and functional-input counts only."""
import argparse
import hashlib
import json
import resource
import signal
import time
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_cubic_zw25_pair_unit_completion_20261104_v63 import Lazy

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_PURE3_FUNCTIONAL_INPUT_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text())
    ib=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63'
    cb=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_cubic_zw25_pair_unit_completion_20261104_v63.py',
           ib/'report.json',ib/'input.json',cb/'report.json',cb/'checkpoint.json',
           ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json',
           ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63_REVIEW/report.json',
           ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json']
    ir,data,cr,cp=[read(p)for p in paths[2:6]]
    assert sha(paths[3])==ir['input_json_sha256']and sha(paths[5])==cr['evidence_sha256']['checkpoint.json']
    assert all(read(p)['all_passed']for p in paths[6:8])
    rows=[10,13,16];relative=[0,3,7]
    assert [data['source25_labels'][i][1]for i in rows]==relative
    # Necessary actual source gate: prove zero BEFORE the genuine H*L pivot.
    assert all(not data['matrix25x32'][i][j]for i in rows for j in range(26,32))
    assert all(not cp['matrix25x32'][i][j][0]for i in rows for j in range(26,32))
    genuine=cp['all_completed_unit_steps'][-1]
    assert [genuine['pivot_row'],genuine['pivot_column']]==[11,26]
    assert not any(i in rows for i,_ in genuine['eliminations'])
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    C=100*a+18*e**2-9*e-2
    decp=lambda v:R.from_dict({tuple(m):s.QQ(c)for m,c in v})
    encp=lambda p:[[list(m),str(c)]for m,c in p.items()]
    enc=lambda p:[encp(p.numer),encp(p.denom)]
    lambdas=[decp(data['genuine16_source_clearing_units'][i])for i in relative]
    assert lambdas==[-64000*C]*3
    raw=[[decp(q)for q in data['matrix25x32'][i][:26]]for i in rows]
    normalized=[[K(q)/K(lam)for q in row]for row,lam in zip(raw,lambdas)]
    assert all(K(lam)*q==K(original)for row,orig,lam in zip(normalized,raw,lambdas)for q,original in zip(row,orig))
    # The divisions are by a known unit. Do not assume the results polynomial.
    denominator_exponents=[];cleared=[]
    for row in normalized:
        exponents=[]
        for q in row:
            den=q.denom;n=0
            while not den.rem(C):den=den.exquo(C);n+=1
            assert den.is_ground and den
            exponents.append(n)
        m=max(exponents);denominator_exponents.append(exponents)
        rr=[]
        for q in row:
            z=q*K(C)**m;assert z.denom.is_ground
            rr.append(z.numer.exquo(z.denom))
        cleared.append(rr)
    Q10=[[decp(q)for q in row[:26]]for row in data['matrix25x32'][:10]]
    assert sum(len(q)for row in Q10 for q in row)==317
    # Exact source cleanup recipe: all changes to these rows are Q10 multiples.
    lazy=Lazy();source_labels=list(range(10))+rows
    vectors={i:[int(i==j)for j in source_labels]for i in source_labels}
    coef_refs=[]
    for number,st in enumerate(cp['all_completed_unit_steps']):
        pivot=st['pivot_row']
        if pivot>=10:
            assert not any(i in source_labels for i,_ in st['eliminations']);continue
        inv=lazy.put(['inverse_old_named_pivot',number])
        vectors[pivot]=[lazy.mul(z,inv)for z in vectors[pivot]]
        for target,co in st['eliminations']:
            if target not in vectors:continue
            node=lazy.put(['old_elimination_coefficient',number,target]);coef_refs.append([number,target,co])
            vectors[target]=[lazy.sub(x,lazy.mul(node,y))for x,y in zip(vectors[target],vectors[pivot])]
    for i,row in enumerate(rows):
        assert vectors[row][10+i]==1
        assert all(vectors[row][10+j]==0 for j in range(3)if i!=j)
    F=500*a+90*e**2+3*e-24
    B=(50*a+18*e**2-9*e-2)*d+60*e-40
    J=5*a*d-6*e-4
    targets=[('1',R.one),('a',a),('d',d),('F',F),('B',B),('C',C),('J',J),('a*d',a*d)]
    shifts=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    def union_for(row,ss):
        return set().union(*[{tuple(x+y for x,y in zip(m,h))for m in p}for h in ss for p in row])
    equation_profiles=[]
    for i,row in enumerate(cleared):
        exponent=max(denominator_exponents[i]);unit=C**exponent
        targetpolys=[unit*p for _,p in targets]
        profiles=[]
        for ss,name in [(shifts[:1],'constant'),(shifts,'affine_1_a_e_d')]:
            qrows=[union_for(rr,ss)for rr in Q10]
            grows=union_for(row,ss)|set().union(*(set(p)for p in targetpolys))
            profile={'cofactor_class':name,'unknowns':26*len(ss),
                     'Q10_polynomial_equation_row_counts':[len(x)for x in qrows],
                     'genuine_polynomial_equation_row_count':len(grows),
                     'total_QQ_row_equations':sum(map(len,qrows))+len(grows),
                     'serialized_input_coefficient_terms':len(ss)*(317+sum(len(p)for p in row))+sum(len(p)for p in targetpolys)}
            profiles.append(profile)
        equation_profiles.append({'old_source_label':rows[i],'relative_genuine':relative[i],
                                  'normalized_equation':'Q10*w=0, Gi_normalized*w=u',
                                  'polynomial_clearing_C_exponent':exponent,
                                  'cleared_equation':'cleared_Gi*w=C^exponent*u',
                                  'profiles':profiles})
    proof={'old3_source_labels':rows,'relative16_genuine_indices':relative,
           'raw_genuine_source135_indices':[18,39,68],
           'all18_initial_heavy_entries_zero_BEFORE_G11':True,
           'all18_checkpoint_heavy_entries_zero_AFTER_G11':True,
           'actual_G11_source_subtraction_coefficients':['0','0','0'],
           'G11_pivot_remains_heavy_not_pure':[11,26],
           'normalized_known_lambdas':[encp(x)for x in lambdas],
           'all78_known_lambda_division_identities_verified':True,
           'normalized3x26_rationals':[[enc(q)for q in row]for row in normalized],
           'normalized_denominator_C_exponents':denominator_exponents,
           'cleared3x26_polynomials':[[encp(q)for q in row]for row in cleared],
           'Q10x26_polynomials':[[encp(q)for q in row]for row in Q10],
           'exact_source_cleanup_scalar_DAG':{'nodes':lazy.nodes,'source_basis_labels':source_labels,
               'old_checkpoint_reference':str(paths[5].relative_to(ROOT)),
               'coefficient_refs':[ [i,t,co]for i,t,co in coef_refs],
               'row3_source_cofactor_roots':[vectors[i]for i in rows]},
           'functional_descends_through_all_old_Q_operations':'Each old rowGi=rawGi+sum(beta_ij Qj); G11 coefficient is literally0. Thus Q10*w=0 implies oldGi*w=lambda_i*(Gi_normalized*w).',
           'known_targets':[[name,encp(p)]for name,p in targets],
           'cofactor_monomial_basis':shifts,'equation_profiles':equation_profiles,
           'success_would_prove':'A knownunit value on one normalizedGi and exactQ10annihilation gives a genuine source direction outside the entire quadratic image; pure/heavy tails remain in the actual source presentation.',
           'no_QQ_span_solve_rank_GB_point_or_newminor':True}
    (out/'certificate.json').write_text(json.dumps(proof,separators=(',',':'))+'\n')
    result={'id':ID,'actual_initial_and_post_G11_heavy_zero_gate_passed':True,
            'G11_cleanup_exactly_zero_not_omitted_by_assumption':True,
            'all78_known_lambda_divisions_passed':True,
            'all_source_changes_only_Q10_multiples':True,
            'raw3_total_terms':sum(len(q)for row in raw for q in row),
            'normalized_numerator_terms':[sum(len(q.numer)for q in row)for row in normalized],
            'normalized_denominator_max_C_power':[max(row)for row in denominator_exponents],
            'cleared_source_terms':[sum(len(p)for p in row)for row in cleared],
            'constant_and_affine_unknown_counts':[26,104],
            'equation_profiles':equation_profiles,
            'no_QQspan_solve_higher_degree_point_GB_or_207expansion':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),
            'seconds':time.monotonic()-began,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output.resolve())
