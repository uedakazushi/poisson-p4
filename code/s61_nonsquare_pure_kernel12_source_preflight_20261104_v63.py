#!/usr/bin/env python3
"""Read-only intrinsic-heavy-chart12x19 and prePOWER10x17 source comparison."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return {tuple(m):Q(v)for m,v in p}
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    cb=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63';cp=read(cb/'checkpoint.json');assert sha(cb/'checkpoint.json')==read(cb/'report.json')['evidence_sha256']['checkpoint.json']
    hb=ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63';hc=read(hb/'certificate.json')
    gb=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';gc=read(gb/'certificate.json')
    qb=ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63';qc=read(qb/'certificate.json')
    pb=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63';pc=read(pb/'certificate.json')
    rows,columns=cp['remaining_rows'],cp['remaining_columns'];qr=[i for i in rows if i<10];gr=[i for i in rows if i>=10];pure=[j for j in columns if j<26];heavy=[j for j in columns if j>=26]
    assert qr==[4,5,9]and len(gr)==14 and len(pure)==19 and heavy==[27,28,29,30,31]
    assert gr==hc['final14_source_labels']and heavy==hc['final5_heavy_target_labels']
    assert all(not cp['matrix25x32'][i][j][0]for i in qr for j in heavy)
    assert hc['old_checkpoint70_rationals']==[[cp['matrix25x32'][i][j]for j in heavy]for i in gr]
    pureterms=sum(len(cp['matrix25x32'][i][j][0])for i in rows for j in pure);qterms=sum(len(cp['matrix25x32'][i][j][0])for i in qr for j in pure)
    N=hc['raw_compact14x5'];nterms=sum(len(p)for row in N for p in row);assert pureterms==70399 and qterms==692 and nterms==2407
    nz=[i for i,row in enumerate(N)if any(row)];assert len(nz)==11
    potential=[]
    for I in itertools.combinations(nz,5):
        if any(all(N[I[k]][p[k]]for k in range(5))for p in itertools.permutations(range(5))):potential.append(list(I))
    assert len(potential)==371
    rawmap={r:i for i,(original,relative,r)in enumerate(hc['actual_genuine_raw_checkpoint_index_map'])}
    pair=[rawmap[i]for i in gc['sparse_pair_rows']];brows=[rawmap[i]for i in gc['all_other3_column_support_exact_rows']]
    assert pair==[1,3]and brows==[6,8,9,10,11,12,13]
    assert all(poly(N[r][1+j])==poly(gc['actual_block7x3'][i][j])for i,r in enumerate(brows)for j in range(3))
    fixed35=[sorted(pair+list(I))for I in itertools.combinations(brows,3)];assert len(set(map(tuple,fixed35)))==35
    heavydata=dict(source17_rows=rows,Q3_source_original25_labels=qr,G14_source_original25_labels=gr,
        pure19_target_original32_labels=pure,heavy5_target_original32_labels=heavy,
        source_shape=[17,24],Q3_heavy_entries_all_literal_zero=True,pure_source_terms=pureterms,Q3_pure_terms=qterms,
        compact_heavy_terms=nterms,active_parameter_terms_before_new_products=pureterms+nterms,
        compact_heavy_max_entry_terms=max(len(p)for row in N for p in row),
        compact_heavy_max_degrees=[max(m[j]for row in N for p in row for m,v in p)for j in range(3)],
        normalization='H_actual=diag(rho14)*N14x5*diag(1,1,1,1,gamma), allrho/gamma FULLUunits by accepted70 exactidentities. ScaleG14sources byrho^-1 and5Htargets bygamma^-1; pureG14 become rho^-1*P, Q3 unchanged.',
        no_POWER_or_old26g_normalizer_appears=True,
        all_H5_support_potential_charts=potential,nonzero_heavy_rows=nz,
        support_chart_count=371,these_determinants_NOT_computed_or_claimed_nonzero=True,
        fixed_sparse_pair_rows14=pair,fixed_B7_rows14=brows,Dp_35_H5_row_charts=fixed35,
        Dp_coverage='OnD(p),pairdet=-(6/5)*B*p^2 is unit. Thus rankN=2+rankactualB7x3;35BminorscoverH5good there. Up to explicit row/col permutation sign,theH5minor ispairdet*B3minor.',
        p0_warning='No implicationBgood=>H7 globally. Onp=0,H7maypersist and requiresotheractualH5minors orseparateTheta branch data. DoNOTdropthissurface.',
        H7_source_condition='UseactualrawH7 plus acceptedCdelta/HLunit-splits toinferNrank5; doNOTinferitfromBgood alone.')
    recipe=dict(choose='Anyactual5rowminorDelta=det(N_I) nonzero; no genericGaussian normalization.',
        source='3unchangedQrows +9genuine kernelrows. lambda_u=N_u*adj(N_I).',
        pure12x19='Q3 unchanged; foru outsideI useDelta*rho_u^-1*P_u-sum_i(lambda_ui*rho_i^-1*P_i).',
        source_basis_minor='Delta^9 with3unit diagonalentries; allsource/targetchangesparameteronly.',
        exact_rank_formula='rank(actual17x24)=5+rank(this12x19) onD(Delta). Fullrowrank12 iff17x24rank17.',
        full57_only_sufficient='rank12=>acceptedaug25rank25=>fullcubic57. SmallerprojectedrankdoesNOTrefutefull57.',
        inverse_factors='ONLYknownFULLUunits andchosenactualH5minor; no4oldnormalizerremainders/noa29cofactor.',
        proposed_next='Ifbudgeted, finiteDAG constructionwithsmall5x5adjugateuniversalidentities andexternal70399purecoeffs. No suchconstructionornewcoefficient/minorhere.',
        no_NEW_field_evaluation_minor_rank_prime_point_GB=True)
    # A smaller optional source uses two already accepted GLOBAL pure pivots.
    p15=qc['new15_source_original25_labels'];p22=qc['new22_target_original32_columns'];pure17=[j for j in p22 if j<26]
    assert p15==[9]+gr and len(pure17)==17 and [j for j in p22 if j>=26]==heavy
    hs=pc['actual_h9'];labels=pc['actual_h9_original_target_labels'];assert len(hs)==len(labels)==9 and all(j in pure17 for j in labels)
    h9refs=[qc['remaining15x22_scalar_roots'][0][p22.index(j)]for j in labels]
    hdata=dict(source_h9=hs,original_target_labels=labels,source_prePOWER_DAG_roots=h9refs,
        source_prePOWER_certificate=str((qb/'certificate.json').relative_to(ROOT)),
        source_h9_binding_certificate=str((pb/'certificate.json').relative_to(ROOT)),
        source_h9_terms=sum(map(len,hs)),source_h9_degrees=[max(m[j]for p in hs for m,v in p)for j in range(3)],
        source_h9_max_QQ_bits=max(max(abs(Q(v).numerator).bit_length(),Q(v).denominator.bit_length())for p in hs for m,v in p),
        pure_tail8_labels=[j for j in pure17 if j not in labels],all8pure_tails_must_be_retained=True,
        covering_proof='Acceptedliteral sumHi*hi=a29 andainvertible=>nineD(hi)coverFULLU. Reuseproofhash; doNOTexpand/readnewdegree207cofactors.',
        no_POWER_cofactor_evaluation_performed=True)
    alternative=dict(source_shape=[15,22],source15_labels=p15,target22_labels=p22,
        why_legal='TwoacceptedGLOBALpureQcompletions from17x24, source/targetinverses onlyaandQQunits; finalPOWERGL10notused.',
        on_H5_chart='KeeponeQ9row and9genuineheavykernelrows ->10x17 purematrix.',
        on_one_of_nine_Q9_charts='Withhi nonzero, cleartheQ9purepivot retainingall8otherpuretails; get9x16. InvertONLYactualH5minorandactualhi.',
        exact_rank='rank(M17)=2+rank(M15)=2+5+rank(K10)=8+rank(N9). Same projected detector, NOTfullsourcecokernel.',
        chart_cover='ActualH5minorscoverH7, ninehi coverFULLU. No old4g denominatorartifact needed; actualrank testsstillOPEN.',
        no_hypersurface_equation_or_globalFULLU57_claim=True)
    for name,data in [('H5-source-and-chart-profile.json',heavydata),('kernel12-recipe.json',recipe),('prePOWER-h9-source.json',hdata),('smaller10x17-alternative.json',alternative)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n')
    deps=['S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63_REVIEW']
    assert all(read(ROOT/'reports'/x/'report.json')['all_passed']for x in deps)
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',cb/'report.json',cb/'checkpoint.json',hb/'report.json',hb/'certificate.json',gb/'report.json',gb/'certificate.json',qb/'report.json',qb/'certificate.json',pb/'report.json',pb/'certificate.json']+[ROOT/'reports'/x/'report.json'for x in deps]
    report=dict(id=ID,all_passed=True,source_split_Q3_G14_pure19_H5_verified=True,active_parameter_terms=72806,
        support_potential_H5_charts=371,fixedDp_H5_charts=35,optional_h9_terms=373,
        no_new_parameter_coefficients_minor_rank_expansion_point_prime_or_GB=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
