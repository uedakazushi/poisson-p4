#!/usr/bin/env python3
"""Independent source scaling/all35 determinants/fixed coprime pair audit.

No author native executable, solver, Gaussian, GB, point or prime is run.
Sparse Fraction products handle all fixed matrix identities; one exact QQ
multivariate gcd of the two source-bound last generators proves gcd22=1.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
import s61_nonsquare_genuine14_H5_review_20261104_v63 as P

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_MINORS_GCD_20261104_v63_REVIEW'
ZERO=(0,0,0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def determinant3(rows):
    out={}
    for j in range(3):
        k,l=(j+1)%3,(j+2)%3
        cofactor=P.add(P.mul(rows[1][k],rows[2][l]),P.scale(P.mul(rows[1][l],rows[2][k]),Q(-1)))
        out=P.add(out,P.mul(rows[0][j],cofactor))
    return out
def monic(p):
    if not p:return {}
    return P.scale(p,1/p[max(p)])
def key(p):return tuple(sorted(monic(p).items(),reverse=True))

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    costs=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63'
    minors=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63'
    gcds=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINOR_GCD_20261104_v63'
    hashes={}
    for base in (costs,minors,gcds):
        report=read(base/'report.json')
        for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
        for p,h in report['evidence_sha256'].items():assert sha(base/p)==h;hashes[str((base/p).relative_to(ROOT))]=h
        hashes[str((base/'report.json').relative_to(ROOT))]=sha(base/'report.json')
    gates=[
      'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_MINOR_OPEN_20261104_v63_REVIEW']
    for name in gates:
        p=ROOT/'reports'/name/'report.json';assert read(p)['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
    cost=read(costs/'certificate.json');cert=read(minors/'certificate.json');native=read(minors/'native-certificate.json')
    gproof=read(gcds/'native-certificate.json')
    rawpath=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63/certificate.json';raw=read(rawpath)
    rawdata=[[P.decode(p) for p in row] for row in raw['actual_block7x3']]
    A=[[P.decode(p) for p in row] for row in cost['prepared_integer7x3']]
    assert len(A)==len(rawdata)==7 and all(len(row)==3 for row in A+rawdata)
    assert sum(len(p) for row in rawdata for p in row)==908
    assert sum(len(p) for row in A for p in row)==345
    assert all(x.denominator==1 for row in A for p in row for x in p.values())
    assert cost['actual_original_block_rows14']==raw['all_other3_column_support_exact_rows']
    assert cost['actual_original_block_columns14']==raw['other3_original_columns']
    assert cost['variable_order']==['a','e','d']
    knownfile=ROOT/'reports/S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63/certificate.json'
    known=[P.decode(p) for p in read(knownfile)['known_factor_polynomials']]
    a={(1,0,0):Q(1)};e={(0,1,0):Q(1)};d={(0,0,1):Q(1)};one={ZERO:Q(1)}
    F=P.add(P.add(P.scale(a,500),P.scale(P.power(e,2),90)),P.add(P.scale(e,3),P.scale(one,-24)))
    C=P.add(P.add(P.scale(a,100),P.scale(P.power(e,2),18)),P.add(P.scale(e,-9),P.scale(one,-2)))
    B=P.add(P.mul(P.add(P.add(P.scale(a,50),P.scale(P.power(e,2),18)),P.add(P.scale(e,-9),P.scale(one,-2))),d),P.add(P.scale(e,60),P.scale(one,-40)))
    J=P.add(P.scale(P.mul(a,d),5),P.add(P.scale(e,-6),P.scale(one,-4)))
    assert known==[a,d,F,B,C,J]
    def unit(record):
        c=Q(record['constant']);assert c
        powers=record['a_d_F_B_C_J_powers'];assert len(powers)==6 and all(isinstance(n,int) and n>=0 for n in powers)
        p={ZERO:c}
        for f,n in zip(known,powers):p=P.mul(p,P.power(f,n))
        if 'polynomial' in record:assert p==P.decode(record['polynomial'])
        return p
    rowunits=list(map(unit,cost['known_row_unit_divisors']))
    colunits=list(map(unit,cost['known_column_unit_divisors']))
    multipliers=list(map(Q,cost['integer_column_multipliers_QQ_units']))
    assert len(rowunits)==7 and len(colunits)==len(multipliers)==3 and all(multipliers)
    for i in range(7):
        for j in range(3):
            assert P.mul(P.mul(rowunits[i],A[i][j]),colunits[j])==P.scale(rawdata[i][j],multipliers[j])
    common_column=one
    for u,c in zip(colunits,multipliers):common_column=P.mul(common_column,P.scale(u,1/c))
    records=native['records'];triples=list(itertools.combinations(range(7),3))
    assert len(records)==35 and [r['rows7'] for r in records]==list(map(list,triples))
    primitive=[];mapping=[];seen={};computed=[];total=0;zeros=0
    for index,(I,record) in enumerate(zip(triples,records)):
        det=determinant3([A[i] for i in I]);computed.append(det)
        assert det==P.decode(record['prepared_integer_minor'])
        source=P.decode(record['actual_source_minor'])
        factor=common_column
        for i in I:factor=P.mul(factor,rowunits[i])
        assert P.mul(factor,det)==source
        if not det:
            assert record['zero'];zeros+=1;mapping.append(None);continue
        assert not record['zero'];total+=len(det)
        stripped=P.decode(record['primitive_integer_polynomial'])
        assert all(c.denominator==1 for c in stripped.values())
        u=unit({'constant':record['named_unit_constant'],'a_d_F_B_C_J_powers':record['named_a_d_F_B_C_J_powers']})
        assert P.mul(u,stripped)==det
        k=key(stripped)
        if k not in seen:seen[k]=len(primitive);primitive.append(monic(stripped))
        mapping.append(seen[k])
        assert record['single_knownunit_minor']==(set(stripped)=={ZERO})
    assert zeros==7 and len(records)-zeros==28 and total==3525
    assert mapping==cert['minor35_to_distinct_generator_index']
    assert primitive==list(map(P.decode,cert['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']))
    assert len(primitive)==22 and sum(map(len,primitive))==2963
    assert max(sum(m) for p in primitive for m in p)==15
    assert mapping[33:]==[20,21]
    assert [len(primitive[i]) for i in (20,21)]==[88,199]
    assert [max(sum(m) for m in primitive[i]) for i in (20,21)]==[10,14]
    # TWO actual source determinants recomputed directly, rather than only
    # inferred from source entry scaling or from prepared determinant values.
    for index in (33,34):
        I=triples[index]
        actual=determinant3([rawdata[i] for i in I])
        assert actual==P.decode(records[index]['actual_source_minor'])
    selected=raw['selected_rows7'];selected_columns=raw['selected_cols3']
    assert len(set(selected))==3 and sorted(selected_columns)==[0,1,2]
    # The old selected Delta has NON-sorted row/column orders. Keep both
    # permutations rather than silently replacing them by lex order.
    selected_det=determinant3([[rawdata[i][j] for j in selected_columns] for i in selected])
    assert selected_det==P.decode(raw['nonzero_selected_minor_Delta'])
    row_sign=(-1)**sum(selected[i]>selected[j] for i in range(3) for j in range(i+1,3))
    column_sign=(-1)**sum(selected_columns[i]>selected_columns[j] for i in range(3) for j in range(i+1,3))
    selected_index=triples.index(tuple(sorted(selected)))
    assert selected_det==P.scale(P.decode(records[selected_index]['actual_source_minor']),Q(row_sign*column_sign))
    # One small fixed exact gcd of the SOURCE-BOUND primitive last pair.
    # No Groebner basis, factorization, modular point or native author job.
    ring_,ra,re,rd=ring('a,e,d',sp.QQ)
    rpoly=lambda p:ring_.from_dict({m:sp.QQ(c.numerator,c.denominator) for m,c in p.items()})
    pairgcd=rpoly(primitive[20]).gcd(rpoly(primitive[21]))
    assert pairgcd==ring_.one
    assert P.decode(gproof['gcd22'])==one and P.decode(gproof['remaining_h_coprime_U'])==one
    assert not gproof['localizer_removal_steps'] and P.decode(gproof['final_gcd_h_U'])==one
    assert len(gproof['all22_literal_quotients'])==22
    assert list(map(P.decode,gproof['all22_literal_quotients']))==primitive
    U=one
    for f in known:U=P.mul(U,f)
    assert U==P.decode(cert['FULLU_product'])==P.decode(native['FULLU_product'])
    assert len(U)==49 and max(map(sum,U))==11
    proof=dict(source_block21_entries_908terms_from_independent_rawH_gate=True,
      all21_unit_preconditioning_identities_and_integer345term_input_verified=True,
      all35_exact_compact_determinants_and_source_scaling_identities_verified=True,
      zero_minors=7,nonzero_minors=28,compact_minor_total_terms=3525,
      all_named_factor_remultiplications_and_QQ_only_dedup_verified=True,
      monic_generators=22,monic_generator_terms=2963,max_generator_degree=15,
      two_actual_source_determinants33_34_independently_recomputed=True,
      selected407Delta_exact_row_and_column_orders_preserved=True,
      coprime_pair_minor_indices=[33,34],coprime_pair_terms=[88,199],coprime_pair_degrees=[10,14],
      exact_QQ_multivariate_pair_gcd=1,all22_gcd_equals1=True,
      FULLU_localizer49terms_degree11_verified=True,
      B_rankbad_codimension_at_least2_or_empty_dimension_at_most1_only=True,
      no_unitideal_emptiness_or_global_FULLU_H7_claim=True,
      no_author_native_GB_point_prime_factorization_or_Gaussian=True,
      maximum_polynomial_intermediate_terms=P.MAXTERMS)
    (out/'independent-B35-gcd-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_nonsquare_genuine14_H5_review_20261104_v63.py'):
        hashes[str(p.relative_to(ROOT))]=sha(p)
    review=dict(id=ID,all_passed=True,input_sha256=hashes,
      evidence_sha256={'independent-B35-gcd-proof.json':sha(out/'independent-B35-gcd-proof.json')},
      all35_minors_source_units_and_gcd1_independently_accepted=True,
      classification='B-rankbad empty or codimension>=2/dimension<=1 on FULLU; no unit-ideal claim',
      seconds=time.monotonic()-start,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(review,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(review))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
