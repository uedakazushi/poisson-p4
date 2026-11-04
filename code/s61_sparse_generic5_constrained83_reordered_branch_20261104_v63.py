#!/usr/bin/env python3
"""ONE reordered exact branch, original49039 source, hard47seconds total."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_REORDERED_BRANCH_20261104_v63'
BASE='S61_SPARSE_GENERIC5_CONSTRAINED83_RAW_BLOCK_20261104_v63'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_REORDERED_START_PREFLIGHT_20261104_v63'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT_PREFIX_20261104_v63'
LOW='S61_SPARSE_GENERIC5_CONSTRAINED83_LOW_GROWTH_CONTINUATION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def run(out,work):
    started=time.monotonic();signal.alarm(47);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    gp=ROOT/'reports'/GATE/'report.json';gate=read(gp);assert gate['all_passed'] and sha(gp)=='b486f3d7893fc5ed70d5fca7814eecccabe2408f2beb883bc4d9092e47d400c6'
    rp=ROOT/'reports'/BASE/'report.json';r=read(rp);assert r['all_passed'] and sha(rp)=='3b9d32cf070644c7651adf382155f2379702faa967d4e12047a0d730befdc303'
    mp=rp.parent/'actual-236x83-matrix.json';assert sha(mp)==r['evidence_sha256'][mp.name];m=read(mp)
    assert m['shape']==[236,83] and m['total_matrix_terms']==49039 and m['nonzero_entries']==9431
    cpp=ROOT/f'code/{Path(__file__).stem}.cpp';inp=work/'input.txt'
    with inp.open('w')as f:
        f.write('5 236 83 236 83 9431\n');f.write(' '.join(map(str,range(236)))+'\n');f.write(' '.join(map(str,range(83)))+'\n')
        for row,col,p in m['sparse_entries']:
            f.write(f'{row} {col} '+str(len(p['terms']))+' 0 0 0 0 0\n')
            for t in p['terms']:f.write(str(Q(p['content'])*int(t[0]))+' '+' '.join(map(str,t[1:]))+'\n')
    del m
    oldpaths=[ROOT/'reports'/name/'report.json'for name in (OLD,LOW)]
    for p in oldpaths:assert read(p)['all_passed']
    pre={'id':ID,'source_original_236_by83_terms':49039,'source_original_all_constraints':True,
         'force_first':[192,22,'-t0'],'subsequent_policy':'minimum safe whole coefficient/scratch forecast across every constant or parameter monomial, then source row,column; full or monotone-dominated candidate certificates before arithmetic',
         'no_prior24_or29_replay':True,'numeric_scratch_reserve':64,
         'precise_coefficient_guard':'Before the whole pivot, per update L+2k+3b+64 bounds owned state, product, aligned sum copies,next,back-add and scalar LCM heights; support-set scratch separately bounded.',
         'RAW_guard':'LCM/content/digit bound including signed exponents. Include current selection/pivot files, exact residual,1MiB metadata PLUS128*remainingRows*remainingCols bytes reserved for the NEXT incomplete selection frontier.',
         'Laurent_semantics':'content*x^shift*primitive integer polynomial; unit isQQ*parameter monomial only; negative exponents exact, exponent representation checked before arithmetic.',
         'hard_total_seconds':47,'native_soft_seconds':39,'hard_bytes':32*1024*1024,'hard_live_terms':150000,'hard_memory_bytes':2_000_000_000,
         'soft_time_tightened_from_gate_to_allow_compile_and_final_save':True,
         'no_nonmonomial_pivot_point_prime_gcd_Gram_series_or_T_P':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [gp,rp,mp,*oldpaths]},
         'algorithm_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),cpp]},'work_input_sha256':sha(inp)}
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
    binary=work/'native'
    status={}
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,46-(time.monotonic()-started)))
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(inp),str(out)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,46-(time.monotonic()-started)))
        summary=read(out/'native-summary.json');assert summary['all_passed']
        status={'one_worker_passed':True,'compile_plus_worker_total_seconds':time.monotonic()-started}
    except Exception as e:
        status={'one_worker_passed':False,'failure':repr(e),'compile_plus_worker_total_seconds':time.monotonic()-started,'no_retry':True}
        raise
    finally:
        (out/'execution-status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps({**status,'pivot_count':summary['pivot_count'],'stop_reason':summary['stop_reason']}))
def finish(out,work):
    signal.alarm(47);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));assert not(out/'report.json').exists()
    pre=read(out/'preflight.json');ex=read(out/'execution-status.json');assert ex['one_worker_passed'] and ex['compile_plus_worker_total_seconds']<47
    assert sha(work/'input.txt')==pre['work_input_sha256']
    for section in ('input_sha256','algorithm_source_sha256'):
        for f,h in pre[section].items():assert sha(ROOT/f)==h
    n=read(out/'native-summary.json');m=read(out/'residual-matrix.json');assert n['all_passed'] and m['all_passed']
    k=n['pivot_count'];assert m['shape']==[236-k,83-k] and n['source_initial_terms']==49039 and n['prior_24_or29_steps_replayed']==0
    assert n['peak_live_terms']<=150000 and n['max_forecast_terms']<=150000 and n['future_frontier_RAW_reserve_used']
    for step in range(k):
        p=read(out/f'pivot-{step:02d}.json');s=read(out/f'selection-{step:02d}.json')
        assert p['completed'] and p['all_row_update_back_add_identities_checked']
        assert p['forecast']['peak_live_bound']<=150000 and p['forecast']['output_bytes_upper']<=32*1024*1024
        assert s['selected']==[p['pivot_original_row'],p['pivot_original_column']]
        if step==0:assert s['selected']==[192,22] and p['pivot_QQ_coefficient']=='-1' and p['pivot_parameter_exponent']==[0,0,0,0,1]
        else:
            assert s['all_eligible_considered_or_dominated_before_arithmetic']
            for row,col,status,peak,b in s['records']:
                if status in('safe','dominated'):assert(p['forecast']['peak_live_bound'],*s['selected'])<=(peak,row,col)
    proof={'source':'ORIGINAL49039-term constrained236x83; separate elimination path from immutable old24/29',
           'all_parameters_constraints_preserved':'SameK=QQ(H,L,s0,q0,t0), exact invertible unit row operations and unit pair removals; no constraint dropped.',
           'inverse_instructions':'Apply these source-ID row operations to a236-componentRHS; if a later justified method solves the residual, reverse the unit pivots. No full dense inverse computed.',
           'remaining_rank_from_accepted_rank83':'83-new_branch_pivot_count only; no new rank test',
           'row_update_arithmetic':'All new+factor*normalizedPivot=old identities checked in exact Laurent arithmetic; unit denominators are recorded coordinate monomials only.',
           'separate_branch_pivots':k,'residual_shape':m['shape'],'stop_reason':n['stop_reason'],
           'no_old_prefix_replay_gcd_nonmonomial_or_series':True}
    (out/'exact-reordered-branch-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[str(Path(__file__).relative_to(ROOT)),f'code/{Path(__file__).stem}.cpp',f'research/notes/{ID}.md',*pre['input_sha256']]
    result={'id':ID,'all_passed':True,'separate_original_matrix_branch':True,'unit_pivots':k,'phase_counts':n['phase_counts'],'stop_reason':n['stop_reason'],'exact_residual_shape':m['shape'],'exact_residual_terms':m['total_terms'],'peak_live_terms':n['peak_live_terms'],'max_selected_forecast':n['max_forecast_terms'],'literal_back_add_checks':n['literal_row_update_identity_checks'],'unit_parameters_inverted':n['unit_variables_inverted'],'prior24_and29_unchanged':True,'no_nonmonomial_pivot_gcd_Gram_series_point_prime_or_exhaustion':True,'input_sha256':{f:sha(ROOT/f)for f in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'work_input_sha256':sha(work/'input.txt'),'work_executable_sha256':sha(work/'native'),'native_seconds':n['native_seconds'],'compile_plus_worker_seconds':ex['compile_plus_worker_total_seconds']}
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path);p.add_argument('--finish',action='store_true');a=p.parse_args();(finish if a.finish else run)(a.output.resolve(),a.work.resolve())
