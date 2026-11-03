#!/usr/bin/env python3
"""Independent Fraction audit of all leading cases in the finite G2 cover."""
import argparse
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import resource
import shutil
import signal
import tarfile
import time
from s61_nonsquare_pair_empty_QQ_review_20261103_v63 import read, sha, decode, add, mul, scale, power, subtract

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_TERNARY_PENCIL_FINITE_20261103_v63_REVIEW'


def determinant(matrix, dim):
    answer = {}
    for perm in itertools.permutations(range(len(matrix))):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(len(perm)) for j in range(i+1, len(perm)))
        term = {(0,)*dim: Q(sign)}
        for i, j in enumerate(perm): term = mul(term, matrix[i][j])
        answer = add(answer, term)
    return answer


def main(out):
    started = time.monotonic(); resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000)); signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    parent = ROOT/'reports/S61_NONSQUARE_TERNARY_PENCIL_FINITE_20261103_v63'; author = read(parent/'report.json')
    assert author['frozen'] and author['all_actual_QQ_identities_passed']
    archive = ROOT/author['archive']; assert archive.stat().st_size == 207753 and sha(archive) == author['archive_sha256']
    snap = out/'immutable-input'; snap.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar:
            assert member.isfile(); dst = snap/member.name; assert dst.resolve().is_relative_to(snap.resolve())
            dst.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, dst.open('wb') as target: shutil.copyfileobj(src, target)
    for rel, h in author['input_sha256'].items(): assert sha(snap/rel if (snap/rel).exists() else ROOT/rel) == h, rel
    old_review = ROOT/'reports/S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW/report.json'
    accepted = read(old_review); assert accepted['all_passed']
    for rel, h in accepted['input_sha256'].items(): assert sha(ROOT/rel) == h, rel
    source = read(snap/'work/nonsquare-v62/pair-finite-proof-v76/certificate.json')
    profile = read(snap/'work/nonsquare-v62/quad78-small-gates-v83/pencil-profile.json')
    base = snap/'work/nonsquare-v62/quad78-pencil-pseudo-v84'
    reorder = lambda p: {(ex[0], ex[2], ex[1]): c for ex, c in p.items()}
    rows = [[reorder(decode(p, 3)) for p in row] for row in source['natural_rows']]  # t,H,v
    def symmetric(row):
        return [[row[0], scale(row[1], Q(1,2)), scale(row[3], Q(1,2))],
            [scale(row[1], Q(1,2)), row[2], scale(row[4], Q(1,2))],
            [scale(row[3], Q(1,2)), scale(row[4], Q(1,2)), row[5]]]
    m0, m3 = list(map(symmetric, rows)); pencil = [{} for _ in range(4)]
    for perm in itertools.permutations(range(3)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
        co = [{(0,0,0): Q(sign)}]
        for i, j in enumerate(perm):
            new = [{} for _ in range(len(co)+1)]
            for k, p in enumerate(co):
                new[k] = add(new[k], mul(p, m0[i][j])); new[k+1] = add(new[k+1], mul(p, m3[i][j]))
            co = new
        pencil = [add(p, q) for p, q in zip(pencil, co)]
    Ci = [reorder(decode(p, 3)) for p in profile['normalized_pencil_coefficients']]
    T = {(1,0,0): Q(1)}; P = {(0,1,0): Q(1), (1,0,0): Q(-8)}; J = {(0,1,0): Q(5), (0,0,1): Q(-2), (0,0,0): Q(-4)}
    one = {(0,0,0): Q(1)}
    for actual, ci, norm in zip(pencil, Ci, profile['normalization_records']):
        powers = norm['known_unit_powers_t_P_J']; assert len(powers) == 3 and min(powers) >= 0
        factor = {(0,0,0): Q(norm['QQcontent'])}; assert factor[(0,0,0)]
        for q, k in zip([T,P,J], powers): factor = mul(factor, power(q,k,3))
        assert actual == mul(ci, factor)
    assert list(map(len, Ci)) == [28,56,91,131]
    assert [max(sum(ex) for ex in p) for p in Ci] == [4,6,8,9]
    leading = read(base/'leading-coefficients.json'); LC, MID, CONST = [decode(leading[k], 2) for k in ['LC','MID','CONST']]
    lift = lambda p: {(0,ex[0],ex[1]): c for ex,c in p.items()}
    coefft = lambda p,k: {(ex[1],ex[2]):c for ex,c in p.items() if ex[0] == k}
    assert [coefft(Ci[0],k) for k in [2,1,0]] == [LC,MID,CONST] and max(ex[0] for ex in Ci[0]) == 2
    ab = []
    for i in [1,2,3]:
        item = read(base/f'remainder-{i}.json'); quo = decode(item['pseudo_quotient'],3); raw = decode(item['pseudo_remainder'],3)
        assert mul(power(lift(LC),item['LC_power'],3),Ci[i]) == add(mul(quo,Ci[0]),raw)
        normalized = decode(item['normalized_remainder'],3); pp = item['removed_units_t_P_J_LC']; assert len(pp)==4 and min(pp)>=0
        factor = {(0,0,0):Q(item['QQcontent'])}; assert factor[(0,0,0)]
        for q,k in zip([T,P,J,lift(LC)],pp): factor=mul(factor,power(q,k,3))
        assert raw == mul(factor,normalized) and all(ex[0] <= 1 for ex in normalized)
        ab.append((coefft(normalized,1),coefft(normalized,0)))
    N = []
    J2 = {(ex[1],ex[2]):c for ex,c in J.items()}
    for i,j in [(0,1),(0,2)]:
        item=read(base/f'eliminate-t-{i}-{j}.json'); q=decode(item['polynomial'],2); pj,pl=item['removed_J_LC']; assert min(pj,pl)>=0
        factor=scale(mul(power(J2,pj,2),power(LC,pl,2)),Q(item['QQcontent']))
        assert subtract(mul(ab[i][0],ab[j][1]),mul(ab[j][0],ab[i][1])) == mul(q,factor); N.append(q)
    cert=read(base/'finite-certificate.json'); checks=[]
    def evalv(p,value):
        result={}
        for (ih,iv),c in p.items(): result=add(result,{(ih,):c*value**iv})
        return result
    def coprime(data,F,G):
        assert decode(data['F'],2)==F and decode(data['G'],2)==G
        coefficients=[]
        for idx,degree in data['coefficient_labels']:
            assert idx in [0,1]; coefficients.append({(ex[1],):c for ex,c in [F,G][idx].items() if ex[0]==degree})
        weights=[decode(p,1) for p in data['coefficient_content_Bezout']]; assert len(weights)==len(coefficients)
        identity={}
        for p,q in zip(weights,coefficients): identity=add(identity,mul(p,q))
        assert identity == {(0,):Q(1)}
        value=Q(data['special_v']); f,g=evalv(F,value),evalv(G,value)
        assert max(ex[0] for ex in F)==max(ex[0] for ex in f)==data['F_H_degree']
        assert max(ex[0] for ex in G)==max(ex[0] for ex in g)==data['G_H_degree']
        u,w=[decode(p,1) for p in data['specialized_QQ_Bezout']]
        assert add(mul(u,f),mul(w,g)) == {(0,):Q(1)}
        checks.append(dict(H_degrees=[data['F_H_degree'],data['G_H_degree']], specialization=str(value),content_and_specialized_sum1=True))
    coprime(cert['generic_LC_nonzero_projection'],*N)
    linear=read(base/'LC-zero-MID-nonzero.json'); Rlin=decode(linear['linear_pseudo_remainder'],2); q=decode(linear['quotient'],3)
    pj,pm=linear['removed_J_MID']; assert min(pj,pm)>=0
    raw=scale(mul(mul(power(J,pj,3),power(lift(MID),pm,3)),lift(Rlin)),Q(linear['QQcontent']))
    divisor=add(mul(lift(MID),T),lift(CONST)); degree=max(ex[0] for ex in Ci[1])
    assert mul(power(lift(MID),degree,3),Ci[1]) == add(mul(q,divisor),raw)
    coprime(cert['LC_zero_MID_nonzero_projection'],LC,Rlin)
    def resultant(f,g):
        m,n=max(ex[0] for ex in f),max(ex[0] for ex in g)
        aa=[{(ex[1],):c for ex,c in f.items() if ex[0]==k} for k in range(m,-1,-1)]
        bb=[{(ex[1],):c for ex,c in g.items() if ex[0]==k} for k in range(n,-1,-1)]
        mat=[[{}]*i+aa+[{}]*(n-i-1) for i in range(n)]+[[{}]*i+bb+[{}]*(m-i-1) for i in range(m)]
        return determinant(mat,1)
    last=cert['LC_MID_zero_case']; R1,R2=[decode(last[k],1) for k in ['resultant_R1','resultant_R2']]
    assert R1==resultant(LC,MID) and R2==resultant(LC,CONST)
    u,w=[decode(p,1) for p in last['resultant_gcd_Bezout']]
    fthird={(2,):Q(1),(1,):Q(13,3),(0,):Q(4)}
    assert add(mul(u,R1),mul(w,R2))==decode(last['resultant_gcd'],1)==fthird
    assert len(last['special_branches'])==2; branches=[]
    def evalhv(p,hh,vv):
        ans={}
        for (it,ih,iv),c in p.items(): ans=add(ans,{(it,):c*hh**ih*vv**iv})
        return ans
    for branch,(vv,hh) in zip(last['special_branches'],[(Q(-3),Q(8)),(Q(-4,3),Q(-14,3))]):
        assert Q(branch['v'])==vv and Q(branch['forced_H'])==hh
        polys=[evalv(p,vv) for p in [LC,MID,CONST]]; weights=[decode(p,1) for p in branch['leading_Bezout']]
        identity={}
        for p,q in zip(weights,polys): identity=add(identity,mul(p,q))
        expected=decode(branch['leading_H_polynomial'],1); assert identity==expected
        assert max(ex[0] for ex in expected)==1 and sum(c*hh**ex[0] for ex,c in expected.items())==0
        higher=[evalhv(p,hh,vv) for p in Ci[1:]]; assert higher==[decode(p,1) for p in branch['higher_t_polynomials']]
        weights=[decode(p,1) for p in branch['higher_t_Bezout']]; identity={}
        for p,q in zip(weights,higher): identity=add(identity,mul(p,q))
        assert identity=={(0,):Q(1)}; branches.append({'v':str(vv),'H':str(hh),'leading_and_higher_literal_QQ_Bezout_verified':True})
    evidence=dict(actual_ternary_pencil_all4_coefficients_rebuilt_from_accepted_actual_pair=True,
        known_normalizations_only_t_P_J_QQ=True, four_profile_term_counts=[28,56,91,131],
        all3_pseudo_quotient_remainder_identities_and_2_linear_eliminants_checked=True,
        leading_LC_nonzero_projection_finite=True, leading_LC_zero_MID_unit_projection_finite=True,
        exact_content_and_degree_preserving_good_specialization_certificates=checks,
        independent_small_Sylvester_determinants_and_univariate_sum_f_over3=True,
        leading_allzero_branches=branches, no_global_inverse_of_LC_or_MID=True,
        common_quadratic_factor_case_excluded_by_prior_whole_pair_rank2=True,
        singular_pencil_necessary_locus_and_G2_failure_dimension_at_most0=True,
        no_finite_Zp_proper_model_empty_locus_actual_rank78_or_fullquartic_claim=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n'); shutil.rmtree(snap)
    result=dict(id=STEM,all_passed=True,all_actual_pencil_and_saved_identities_checked_by_independent_Fraction_arithmetic=True,
        all_leading_coefficient_cases_retained=True,QQ_singular_pencil_necessary_locus_dimension_at_most0=True,
        QQ_actual_ternary_commonfactor_locus_dimension_at_most0=True,open='a*P*J!=0',
        finite_not_empty_and_no_arithmetic_Nakayama_promotion=True,
        no_conditional_syzygy_gap_reinstatement_or_quad78_claim=True,
        no_GB_new_field_prime_or_heavy_source_replay=True,seconds=time.monotonic()-started)
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',parent/'report.json',archive,old_review,
        ROOT/'research/notes/S61_NONSQUARE_PAIR_OMIT_BSQUARE_CONDITIONAL_SYZYGY_GAP_20261103_v63.md']
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p) for p in paths}
    result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items() if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
