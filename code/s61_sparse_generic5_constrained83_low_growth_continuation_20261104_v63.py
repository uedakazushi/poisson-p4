#!/usr/bin/env python3
"""New bounded continuation from the immutable24-step residual, not a retry."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_LOW_GROWTH_CONTINUATION_20261104_v63'
PARENT='S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT_PREFIX_20261104_v63'
TIGHT='S61_SPARSE_GENERIC5_CONSTRAINED83_TIGHT_FORECAST_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def prepare(out,work):
    signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports'/PARENT;pr=read(parent/'report.json')
    assert pr['all_passed'] and sha(parent/'report.json')=='31ba1bf17d5f78015125e070147b2250353c3054504ef16f1e34e12f12c14e92'
    tp=ROOT/'reports'/TIGHT/'report.json';tr=read(tp);assert tr['all_passed'] and sha(tp)=='fd7deb2016bd20e84c08aaab6993d5025225c86b580cd73ac2a583676c90d0f7'
    mp=parent/'residual-matrix.json';assert sha(mp)==pr['evidence_sha256'][mp.name]
    m=read(mp);assert m['shape']==[212,59] and m['total_terms']==124840 and m['nonzero_entries']==6657
    first=next(p for i,j,p in m['sparse_entries']if(i,j)==(192,22))
    assert Q(first['content'])*int(first['terms'][0][0])==-1 and first['terms']==[['1',0,0,0,0,1]] and first['shift']==[0]*5
    inp=work/'input.txt'
    with inp.open('w')as f:
        f.write('5 236 83 212 59 6657\n');f.write(' '.join(map(str,m['row_original_ids']))+'\n');f.write(' '.join(map(str,m['column_original_ids']))+'\n')
        for row,col,p in m['sparse_entries']:
            f.write(' '.join(map(str,[row,col,len(p['terms']),*p['shift']]))+'\n')
            for t in p['terms']:f.write(str(Q(p['content'])*int(t[0]))+' '+' '.join(map(str,t[1:]))+'\n')
    assert inp.stat().st_size<32*1024*1024
    pre={
      'id':ID,'phase':'frozen_before_new_arithmetic','initial_shape':[212,59],'initial_coefficients':124840,
      'old24_steps_replayed':False,'original_row_column_ids_retained':True,
      'forced_first_pivot':{'row':192,'column':22,'value':'-t0','prior_forecast':127259,'new_numeric_scratch_reserve':64,'new_whole_live_bound':127307},
      'policy_after_first':'Consider all literal QQ-constant OR Laurent parameter-monomial entries. Choose minimum completed, safe whole coefficient/scratch forecast; tie by original source row,column. Candidates enumerated by cheap growth only to find an incumbent. Monotone partial forecast may prune candidates proved unable to beat it. An incomplete time-limited selection executes NO pivot.',
      'row_update_preforecast':'For each entire prospective pivot, exact support sums/unions bound k product terms,b output terms and whole-pivot evolving live L. L+2k+3b+64 includes aligned add/sub scratch, product,next,back-add check and ephemeral scalar LCM/height integers. Support-set scratch is separately counted. Guard150000 applies before any normalization/update.',
      'RAW_bytes':'Exact content denominator LCM, integer-product height bound from primitive coefficient bit bounds, digits(g)+sum digits(H_i/g)<=n digits(M)+1, exact signed/shifted exponent lengths. Full receipts+selection records+exact residual+1MiB reserve<=32MiB. No gzip.',
      'negative_exponents':'A value is x^shift times an exact fmpq_mpoly. Division by q*x^e changes scalar content by1/q and shift by-e. Every add/sub aligns monomials exactly. Exponent/intermediate representation bounds are checked before arithmetic, not evaluated at a numerical parameter.',
      'whole_worker_induction':'Native retains only the current residual124840-term matrix, no original49039 terms or old24 operation values. Before EACH next pivot enforce all support/RAW bounds; stop at the first state with no safe eligible pivot or before47s. Original-source references remain external and hash-bound.',
      'first_changed_method_RAW_bound':'Previous exact-M bound7.06MB; replacing coefficient maxima by2^bits-1 multiplies M by<4, adding at most one digit per term.8MiB is conservative for the first new pivot, before exact native reforecast.',
      'hard_live_terms':150000,'hard_RAW_bytes':32*1024*1024,'hard_seconds':60,'hard_memory_bytes':2_000_000_000,
      'no_nonmonomial_pivot_T_P_selector_dense_inverse_Gram_series_or_source_point':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [parent/'report.json',mp,tp,parent/'stopped-candidate.json']},
      'algorithm_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),ROOT/f'code/{Path(__file__).stem}.cpp']},
      'work_input_sha256':sha(inp),'work_input_bytes':inp.stat().st_size}
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');print(json.dumps({'preflight_passed':True,'preflight_sha256':sha(out/'preflight.json'),'input_bytes':inp.stat().st_size,'first_pivot':'-t0','first_live_bound':127307}))
def finish(out,work):
    signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));assert not(out/'report.json').exists()
    pre=read(out/'preflight.json');assert sha(work/'input.txt')==pre['work_input_sha256']
    for f,h in pre['input_sha256'].items():assert sha(ROOT/f)==h
    for f,h in pre['algorithm_source_sha256'].items():assert sha(ROOT/f)==h
    n=read(out/'native-summary.json');m=read(out/'residual-matrix.json');assert n['all_passed'] and m['all_passed']
    k=n['pivot_count'];assert m['shape']==[212-k,59-k] and n['completed_original_prefix_steps_not_replayed']==24
    assert n['peak_live_terms']<=150000 and n['max_forecast_terms']<=150000 and n['numeric_bound_scratch_reserve']==64
    for step in range(k):
        p=read(out/f'pivot-{step:02d}.json');s=read(out/f'selection-{step:02d}.json')
        assert p['completed'] and p['all_row_update_back_add_identities_checked'] and p['forecast']['peak_live_bound']<=150000 and p['forecast']['output_bytes_upper']<=32*1024*1024
        assert s['selected']==[p['pivot_original_row'],p['pivot_original_column']]
        if step==0:
            assert s['selected']==[192,22] and p['pivot_QQ_coefficient']=='-1' and p['pivot_parameter_exponent']==[0,0,0,0,1]
        else:
            assert s['all_eligible_considered_or_dominated_before_arithmetic']
            value=(p['forecast']['peak_live_bound'],*s['selected'])
            for row,col,status,peak,b in s['records']:
                if status in ('safe','dominated'):assert value<=(peak,row,col)
    proof={
      'exact_scope':'New continuation from immutable24-step residual; both the old cap stop and all source/constraint identities are preserved.',
      'parameter_field_units':'Every pivot is a nonzero QQ scalar times a Laurent monomial in H,L,s0,q0,t0. The first is -t0. No unknown normal coordinate or numerical parameter is inverted.',
      'negative_exponent_semantics':'Each saved polynomial is content*x^shift*integerPolynomial. Shifted supports and coefficient bounds include negatives exactly; no series/truncation or polynomial GCD is used.',
      'changed_policy':'After explicitly authorized first pivot, all eligible constants/monomials are either fully forecast or proved dominated/unsafe. The chosen safe whole forecast is minimum, with source ID ties.',
      'row_operation_exactness':'All streamed updates satisfy literal new+factor*normalizedPivot=old. Unit normalization and pair removal are invertible operations over the same K5; alpha/beta constraints are transformed, not discarded.',
      'factored_inverse':'Apply both the immutable24 receipts and these new receipts to a RHS. Reverse all pivot rows after a later justified residual solve; this is a finite factored inverse instruction, not a computed full left inverse.',
      'no_new_rank_confirmation':'Residual rank=59-newPivotCount is only the corollary of accepted original rank83 and exact unit pair removals.',
      'stop':n['stop_reason'],'continuation_unit_pivots':k,'total_unit_pivots_including_old':24+k,'remaining_shape':m['shape'],
      'no_nonmonomial_pivot_Gram_series_source_point_or_exhaustion':True}
    (out/'exact-continuation-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[str(Path(__file__).relative_to(ROOT)),f'code/{Path(__file__).stem}.cpp',f'research/notes/{ID}.md',*pre['input_sha256']]
    result={'id':ID,'all_passed':True,'new_unit_pivots':k,'total_including_old24':24+k,'phase_counts':n['phase_counts'],'stop_reason':n['stop_reason'],'exact_residual_shape':m['shape'],'exact_residual_terms':m['total_terms'],'peak_live_terms':n['peak_live_terms'],'max_selected_forecast':n['max_forecast_terms'],'all_back_add_identity_checks':n['literal_row_update_identity_checks'],'unit_variables_inverted':n['unit_variables_inverted'],'old24_pivots_and_stop_unchanged':True,'no_nonmonomial_pivot_Gram_series_or_exhaustion':True,'input_sha256':{f:sha(ROOT/f)for f in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'work_input_sha256':sha(work/'input.txt'),'work_executable_sha256':sha(work/'native'),'native_seconds':n['native_seconds']}
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path);p.add_argument('--finish',action='store_true');a=p.parse_args();(finish if a.finish else prepare)(a.output.resolve(),a.work.resolve())
