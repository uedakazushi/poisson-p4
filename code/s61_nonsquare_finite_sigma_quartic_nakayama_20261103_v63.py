#!/usr/bin/env python3
"""Exact split provenance for the observed finite-Sigma quartic rank test.

No new rank computation is run: the one frozen GF rank receipt is bound.
All16 QQ Q*L identities and all6*126 K*M coefficient identities are checked.
"""
import argparse
from fractions import Fraction as Q
import gzip
import itertools
import json
from pathlib import Path
import resource
import signal
import time
from s61_nonsquare_structural_reductions_review_20261103_v63 import read,sha,add,mul

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_FINITE_SIGMA_QUARTIC_NAKAYAMA_20261103_v63'
PCHAR=32003


def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    probe=ROOT/'reports/S61_NONSQUARE_FINITE_SIGMA_QUARTIC_PROBE_20261103_v63';pr=read(probe/'report.json')
    assert pr['native_completed']and pr['single_nmod_mat_rank_call']and pr['rank']==2184 and pr['one_rank_call_only']
    assert pr['regular_representation_shape']==[2352,3528]
    for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    published=read(probe/'publish-receipt.json')
    for rel,h in published['copied_immutable_execution_files'].items():assert sha(probe/rel)==h,rel
    assert published['generated_block_sha256']==pr['actual_block_u32le_sha256']
    finite=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW/report.json'
    fr=read(finite);assert fr['all_passed']and fr['finite_V_algebra_without_flatness_reviewed']and fr['FULLU_unit_on_ENTIRE_special_fiber_and_finite_V_model_reviewed']
    for rel,h in fr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    split_review=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v63_REVIEW/report.json'
    sr=read(split_review);assert sr['all_passed']and sr['independent_full_Fraction_right_inverse']
    for rel,h in sr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';data=read(source)
    def mons(n,degree):return[(degree,)]if n==1 else[(i,)+v for i in range(degree+1)for v in mons(n-1,degree-i)]
    m2,m3,m4=mons(6,2),mons(6,3),mons(6,4);assert[len(m2),len(m3),len(m4)]==[21,56,126]
    quadrics=[];count=0
    for row,shift in zip(data['cleared_qfour_sparse'],[0,2,2,1]):
        q={}
        for mon,coeffs in row:
            assert len(mon)==7 and mon[6]==0 and sum(mon)==2
            p={}
            for ex,c in coeffs:
                c=Q(c);assert c.denominator%PCHAR;new=(ex[0]-ex[2]+shift,ex[1],ex[2]);assert min(new)>=0
                p=add(p,{new:c});count+=1
            q[tuple(mon[:6])]=p
        quadrics.append(q)
    assert count==804
    Qmatrix=[[q.get(mon,{})for mon in m2]for q in quadrics]
    # Verify the source used by the native module, not a Q2prime row variant.
    raw=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW/own-raw-modp.txt'
    tokens=iter(raw.read_text().split());assert[int(next(tokens))for _ in range(3)]==[24,56,3]
    literal=[]
    for i in range(24):
        row=[]
        for j in range(56):
            p={}
            for _ in range(int(next(tokens))):
                c=int(next(tokens));ex=tuple(int(next(tokens))for _ in range(3))
                if c:p[ex]=c
            row.append(p)
        literal.append(row)
    assert next(tokens,None)is None
    for i in range(4):
        for k in range(6):
            for col,mon in enumerate(m3):
                lower=tuple(x-int(j==k)for j,x in enumerate(mon))
                q=quadrics[i].get(lower,{})if min(lower)>=0 else{}
                expected={ex:c.numerator*pow(c.denominator,-1,PCHAR)%PCHAR for ex,c in q.items()}
                expected={ex:c for ex,c in expected.items()if c}
                assert literal[6*i+k][col]==expected
    split=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/right_inverse.txt.gz'
    assert sha(split)=='eaa64c439ca14d7ce5d274e5141e02f41f3c27ae25e2f2a24ace305f1d9e3e43'
    tokens=iter(gzip.decompress(split.read_bytes()).decode().split());assert[int(next(tokens))for _ in range(3)]==[21,4,0]
    L=[];ncofs=0
    for i in range(21):
        row=[]
        for j in range(4):
            p={}
            for _ in range(int(next(tokens))):
                c=Q(next(tokens));ex=tuple(int(next(tokens))for _ in range(3));assert c and min(ex)>=0 and ex not in p
                assert c.denominator%PCHAR;p[ex]=c;ncofs+=1
            row.append(p)
        L.append(row)
    assert next(tokens,None)is None and ncofs==2461
    g={(2,0,3):Q(-5),(2,1,2):Q(6),(2,0,2):Q(4)}
    for i in range(4):
        for j in range(4):
            actual={}
            for k in range(21):actual=add(actual,mul(Qmatrix[i][k],L[k][j]))
            assert actual==(g if i==j else{})
    # K's exact six left syzygies have the SAME quadratic-cofactor ordering as M.
    pairs=list(itertools.combinations(range(4),2));km_entries=0
    for i,j in pairs:
        for gamma in m4:
            actual={}
            for alpha in m2:
                beta=tuple(x-y for x,y in zip(gamma,alpha))
                if min(beta)<0:continue
                actual=add(actual,mul(quadrics[j].get(alpha,{}),quadrics[i].get(beta,{})))
                actual=add(actual,{ex:-c for ex,c in mul(quadrics[i].get(alpha,{}),quadrics[j].get(beta,{})).items()})
            assert not actual;km_entries+=1
    assert km_entries==6*126
    # Split K without a 78x78 determinant or Artin socle inference:
    # column(a,b) has L[:,b]/g in block a, all other blocks zero.
    split_signs=[]
    for i,j in pairs:
        row=[]
        for a,b in pairs:
            value=int(j==b)if a==i else(-int(i==b)if a==j else 0)
            row.append(value)
        split_signs.append(row)
    assert split_signs==[[int(i==j)for j in range(6)]for i in range(6)]
    assert 84*28-6*28==pr['rank']==2184
    evidence=dict(actual804_quad_terms_and_all24_native_source_rows_matched=True,
        original_quartic_rows='quadric_i times m2[shift], row=i*21+shift; m2 lexicographic exponent order',
        original_quartic_columns='m4 lexicographic exponent order',
        exact_all16_QQ_Q_times_L_equals_gI4=True,right_inverse_terms=2461,
        all_quad_and_L_coefficients_32003_integral=True,g='-a^2*P^2*(5P-6e-4)',
        g_unit_on_finite_V_model_from_existing_FULLU_unit_proof=True,
        actual_K6by84_M84by126_all756_QQ_entries_zero=True,
        K_right_inverse='column(a,b) puts L[:,b]/g in block a, for a<b',K_right_inverse_product_I6=split_signs,
        regular_representation_orientation='M is S^126 -> S^84; coefficient operators use output-row/input-column orientation',
        kerK_special_dimension=2184,observed_native_rank=2184,
        imageMbar_equals_kerKbar_by_exact_dimension=True,
        split_kernel_commutes_with_reduction_modp_without_S_flatness=True,
        finite_V_cokernel_zero_by_Nakayama=True,
        raw_quartic_rank78_on_entire_QQ_finite_Sigma=True,
        original_light146_and_weighted152_on_Sigma_via_accepted_degree3_colon_and_ZW19=True,
        no_second_rank_large_inverse_new_field_or_GB=True,
        no_global_quad78_or_C122_rank64_or_fullquartic_claim=True)
    (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    result=dict(id=STEM,all_exact_QQ_split_and_source_checks_passed=True,
        completed_frozen_GF_rank2184_bound=True,proper_finite_V_plus_split_kernel_Nakayama_bridge_accepted=True,
        actual_raw_quad78_proved_on_ENTIRE_finite_Sigma=True,
        original_light146_weighted152_on_Sigma_corollary=True,
        other_finite_quad_exceptions_global_quad78_C122_rank64_not_proved=True,
        no_second_rank_or_automatic_heavy_launch=True,independent_review_pending=True,seconds=time.monotonic()-start)
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',source,split,finite,split_review,raw,probe/'report.json',probe/'native.json',probe/'native.log',probe/'publish-receipt.json']
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths}
    result['evidence_sha256']={'certificate.json':sha(out/'certificate.json')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
