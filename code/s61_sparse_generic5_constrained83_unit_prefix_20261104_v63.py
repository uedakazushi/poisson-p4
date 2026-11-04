#!/usr/bin/env python3
"""Bounded deterministic constant/parameter-monomial prefix: orchestration.

The native worker forecasts every update's support, visible coefficient
scratch, rational height and complete serialized output before a pivot.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT_PREFIX_20261104_v63'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_RAW_BLOCK_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def prepare(out,work):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    sr=ROOT/'reports'/SOURCE/'report.json';rp=read(sr)
    assert rp['all_passed'] and sha(sr)=='3b9d32cf070644c7651adf382155f2379702faa967d4e12047a0d730befdc303'
    mp=sr.parent/'actual-236x83-matrix.json';assert sha(mp)==rp['evidence_sha256'][mp.name]
    m=read(mp);assert m['all_passed'] and m['shape']==[236,83] and m['total_matrix_terms']==49039
    ip=work/'input.txt'
    with ip.open('w') as f:
        f.write(f"5 236 83 {len(m['sparse_entries'])}\n")
        for row,col,p in m['sparse_entries']:
            f.write(f"{row} {col} {len(p['terms'])}\n")
            for t in p['terms']:f.write(str(Q(p['content'])*int(t[0]))+' '+' '.join(map(str,t[1:]))+'\n')
    assert ip.stat().st_size<32*1024*1024
    plan={
        'id':ID,'phase':'frozen_pre_arithmetic_policy',
        'source_shape':[236,83],'source_total_coefficients':49039,
        'source_column_order':'c,alpha20,beta20,n42; no column is dropped except a proved unit pair',
        'five_parameter_order':['H','L','s0','q0','t0'],
        'allowed_pivots':'QQ constants first whenever present; otherwise one Laurent parameter monomial with nonzero QQ coefficient. No nonmonomial pivot.',
        'selection':'Minimum exact structural fill=sum over target rows of #pivot-row active columns absent in that row; tie by term-pair majorant, coefficient height, original row, original column.',
        'forecast_before_any_pivot':'For each proposed update, multiply support sets of the old column factor and normalized pivot entry; union with old entry. This proves k product terms and b result terms before arithmetic. The whole prospective pivot is forecast before any normalization/update.',
        'visible_live_bound':'Let L be current/predicted owned coefficient slots, k product support bound and b union bound. L+2k+3b+16 bounds product, two aligned add/sub scratch polynomials, result and literal back-add identity check. Support-set scratch is counted separately against150000. After each update L is advanced by b-old_length; cleared column and retired pivot row are subtracted.',
        'normalization_bound':'Monomial scalar/shift scaling in place changes no support; no polynomial inverse or GCD is used.',
        'serialized_bound':'Before a pivot, normalized row/factor receipts plus the entire predicted residual are bounded using exact supports, rational numerator/denominator bit bounds and signed exponent digit bounds, with1MiB metadata reserve; must be<=32MiB.',
        'time_policy':'Stop before a new pivot at47 native seconds to reserve time for final residual serialization; parent process hard60seconds/2GB.',
        'stop_policy':'No allowed pivot, next deterministic pivot forecast over any cap, or time budget: serialize the exact current residual and all completed operations; no larger retry or nonmonomial pivot.',
        'receipt':'Stream normalized pivot row and all target factors with original source row/column identifiers; retain per-update forecast records and literal back-add checks. Forward row operations and reverse triangular substitutions define a factored inverse instruction, not a dense236-square inverse.',
        'all_original_parsed_source_arrays_released':'Native reads only the new49039-term matrix. Python JSON is not present in the native address space.',
        'no_T_P_selector_Gram_or_dense_inverse':True,
        'hard_live_terms':150000,'hard_artifact_bytes':32*1024*1024,'hard_seconds':60,'hard_memory_bytes':2_000_000_000,
        'input_sha256':{str(sr.relative_to(ROOT)):sha(sr),str(mp.relative_to(ROOT)):sha(mp)},
        'work_input_sha256':sha(ip),'work_input_bytes':ip.stat().st_size,
        'prepare_seconds':time.monotonic()-started}
    (out/'preflight.json').write_text(json.dumps(plan,indent=2)+'\n')
    print(json.dumps({'preflight_passed':True,'coefficients':49039,'input_bytes':ip.stat().st_size,'preflight_sha256':sha(out/'preflight.json')}))
def finish(out,work):
    signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));assert not(out/'report.json').exists()
    p=read(out/'preflight.json');assert sha(work/'input.txt')==p['work_input_sha256']
    for f,h in p['input_sha256'].items():assert sha(ROOT/f)==h
    n=read(out/'native-summary.json');assert n['all_passed']
    residual=read(out/'residual-matrix.json');assert residual['all_passed']
    assert len(residual['row_original_ids'])==236-n['pivot_count']
    assert len(residual['column_original_ids'])==83-n['pivot_count']
    receipts=[];phase_counts={'QQ_constant':0,'parameter_monomial':0}
    for k in range(n['pivot_count']):
        f=out/f'pivot-{k:02d}.json';j=read(f);assert j['completed'] and j['step']==k
        assert j['forecast']['peak_live_bound']<=150000 and j['forecast']['output_bytes_upper']<=32*1024*1024
        assert j['all_row_update_back_add_identities_checked']
        phase_counts[j['phase']]+=1
        receipts.append({'step':k,'row':j['pivot_original_row'],'column':j['pivot_original_column'],'phase':j['phase'],'sha256':sha(f)})
    assert phase_counts==n['phase_counts']
    assert n['peak_live_terms']<=150000 and n['max_forecast_terms']<=150000
    assert sum(f.stat().st_size for f in out.iterdir()if f.is_file())<32*1024*1024
    proof={
        'coefficient_field':'QQ(H,L,s0,q0,t0)',
        'exact_equivalence':'Each normalized pivot and streamed subtraction is an invertible row operation over the coefficient field. Unit pair removal gives an identity block plus the exact residual; all alpha/beta constraints are retained through these row operations.',
        'left_inverse_instructions':'Apply the same streamed forward operations to an arbitrary236-row RHS; solve the residual if a later justified method succeeds; reverse the pivot list with x_col=normalized_RHS−sum(normalized_pivot_row[col_j]*x_j). This is a finite factored instruction, not an expanded dense left inverse.',
        'first42_target':'For the raw top134 RHS and zero bottom102 constraints, reverse substitution in the actual n42 source columns gives the first42 inverse directions once the residual is solved. No right-hand-side/source variable is inverted.',
        'current_scope':'Only this exact prefix and residual; no nonmonomial pivot, residual inverse, remaining84-direction inverse, ideal membership, or exhaustion.',
        'receipts':receipts,'stop_reason':n['stop_reason'],
        'retained_parameter_names':['H','L','s0','q0','t0'],
        'all_inverse_scalars':'Nonzero QQ constants and explicitly recorded parameter monomials ONLY',
        'unit_variables_inverted':n['unit_variables_inverted'],
        'residual_shape':residual['shape'],
        'residual_rank_from_existing_rank83_only':83-n['pivot_count'],
        'new_rank_search':False}
    (out/'exact-prefix-and-factored-inverse-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[str(Path(__file__).relative_to(ROOT)),f'code/{Path(__file__).stem}.cpp',f'research/notes/{ID}.md',*p['input_sha256']]
    r={'id':ID,'all_passed':True,'exact_unit_pivots':n['pivot_count'],'phase_counts':phase_counts,'stop_reason':n['stop_reason'],'exact_residual_shape':residual['shape'],'exact_residual_terms':residual['total_terms'],'peak_live_terms':n['peak_live_terms'],'max_forecast_terms':n['max_forecast_terms'],'parameter_monomial_units_only':True,'no_nonmonomial_pivot_Gram_dense_inverse_T_selector_or_exhaustion':True,'input_sha256':{f:sha(ROOT/f)for f in paths},'evidence_sha256':{f.name:sha(f)for f in out.iterdir()if f.is_file()},'work_input_sha256':sha(work/'input.txt'),'work_executable_sha256':sha(work/'native'),'native_seconds':n['native_seconds']}
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',required=True,type=Path);a.add_argument('--work',required=True,type=Path);a.add_argument('--finish',action='store_true');v=a.parse_args();(finish if v.finish else prepare)(v.output.resolve(),v.work.resolve())
