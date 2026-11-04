#!/usr/bin/env python3
"""ONE conditional continuation from36, bound to the first FULL-safe frontier."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FRONTIER_SAFE_CONTINUATION_20261104_v63'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_SHORTLIST_CONTINUATION_20261104_v63'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_REMAINING_UNIT_FRONTIER_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def build_input(out,work):
    sd=ROOT/'reports'/SOURCE;gd=ROOT/'reports'/GATE;sr=sd/'report.json';gr=gd/'report.json'
    s=read(sr);g=read(gr)
    assert s['all_passed'] and sha(sr)=='836bbe72ed6c9dcae7c67e748f6400dcca20e61a93f02421c31b2cbe3a1ac6fb'
    assert g['all_passed'] and g['FULL_safe_found'] and sha(gr)=='a53b79a89e7a203cc6aa3f206941881d721a6bd44448586b4a6db54c9f67d464'
    cp=gd/'chosen-candidate.txt';fp=gd/'forecast-029.json';mp=sd/'residual-matrix.json'
    for p in(cp,fp):assert sha(p)==g['evidence_sha256'][p.name]
    assert sha(mp)==s['evidence_sha256'][mp.name]
    assert cp.read_text().split()==['55','78','16/405','0','0','0','0','1']
    f=read(fp);assert f['FULL_safe'] and f['peak_forecast']==148418 and f['complete_RAW_forecast']==13477441
    m=read(mp);assert m['shape']==[200,47] and m['total_terms']==143477 and m['nonzero_entries']==5483
    first=next(p for i,j,p in m['sparse_entries']if(i,j)==(55,78))
    assert Q(first['content'])*int(first['terms'][0][0])==Q(16,405)
    assert [first['shift'][v]+first['terms'][0][v+1]for v in range(5)]==[0,0,0,0,1]
    inp=work/'input.txt'
    with inp.open('w')as o:
        o.write('5 236 83 200 47 5483\n');o.write(' '.join(map(str,m['row_original_ids']))+'\n');o.write(' '.join(map(str,m['column_original_ids']))+'\n')
        for i,j,p in m['sparse_entries']:
            o.write(' '.join(map(str,[i,j,len(p['terms']),*p['shift']]))+'\n')
            for t in p['terms']:o.write(str(Q(p['content'])*int(t[0]))+' '+' '.join(map(str,t[1:]))+'\n')
    return inp,[sr,gr,cp,fp,mp]
def run(out,work):
    started=time.monotonic();signal.alarm(47);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);inp,paths=build_input(out,work);gc.collect()
    cpp=ROOT/f'code/{Path(__file__).stem}.cpp'
    pre={'id':ID,'conditional_FULL_safe_gate':True,'initial_terms':143477,'prior36_not_replayed':True,
         'first_source_candidate':[55,78,'16*t0/405'],'full_source_forecast_terms':148418,'full_source_RAW_bound':13477441,
         'fresh_native_reforecast_before_arithmetic':True,'later_policy':'at most16 deterministic cheap-fill constant/monomial candidates; bestSAFE among examined, NOT global minimum',
         'guards':'L+2k+3b+64 whole-pivot coefficient/scratch; exact signed Laurent; RAW input included, numericLCM/content/exponent byte bound, full future-frontier reserve',
         'hard_total_seconds':47,'native_soft_seconds':39,'hard_live_terms':150000,'hard_RAW_bytes':32*1024*1024,'hard_memory_bytes':2_000_000_000,
         'source_arrays_released_before_native':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
         'algorithm_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),cpp]},'work_input_sha256':sha(inp)}
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');binary=work/'native';status={}
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,46-(time.monotonic()-started)))
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(inp),str(out)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,46-(time.monotonic()-started)))
        n=read(out/'native-summary.json');assert n['all_passed']
        status={'one_authorized_arithmetic_worker_passed':True,'input_compile_worker_total_seconds':time.monotonic()-started}
    except Exception as e:
        status={'one_authorized_arithmetic_worker_passed':False,'failure':repr(e),'no_retry':True};raise
    finally:(out/'execution-status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps({**status,'new_pivots':n['pivot_count'],'stop_reason':n['stop_reason']}))
def finish(out,work):
    assert not(out/'report.json').exists();pre=read(out/'preflight.json');ex=read(out/'execution-status.json')
    assert ex['one_authorized_arithmetic_worker_passed'] and ex['input_compile_worker_total_seconds']<47
    for sec in('input_sha256','algorithm_source_sha256'):
        for f,h in pre[sec].items():assert sha(ROOT/f)==h
    assert sha(work/'input.txt')==pre['work_input_sha256']
    n=read(out/'native-summary.json');m=read(out/'residual-matrix.json');k=n['pivot_count']
    assert n['all_passed'] and m['all_passed'] and m['shape']==[200-k,47-k] and n['prior36_steps_replayed']==0
    assert n['source_initial_terms']==143477 and n['peak_live_terms']<=150000 and n['max_forecast_terms']<=150000
    for step in range(k):
        p=read(out/f'pivot-{step:02d}.json');s=read(out/f'selection-{step:02d}.json')
        assert p['completed'] and p['all_row_update_back_add_identities_checked'] and p['forecast']['peak_live_bound']<=150000 and p['forecast']['output_bytes_upper']<=32*1024*1024
        if step==0:assert s['selected']==[55,78] and p['pivot_QQ_coefficient']=='16/405' and p['pivot_parameter_exponent']==[0,0,0,0,1]
        else:
            assert len(s['records'])<=16 and s['not_global_minimum_claim'] and s['all_shortlisted_considered_or_dominated_before_arithmetic']
            for i,j,status,peak,b in s['records']:
                if status in('safe','dominated'):assert(p['forecast']['peak_live_bound'],*s['selected'])<=(peak,i,j)
    proof={'actual_source':'Frozen exact36-state200x47; hash-bound first FULL-safe source-only frontier, all source parameters/constraints unchanged.',
           'composite_path':'Same reordered30 + shortlist6 + this continuation; prior36 never replayed.',
           'exact_unit_operations':'Each new unit isQQ times a parameter Laurent monomial. Exact back-add; no unknown normal or numeric parameter inverted.',
           'selection_scope':'After first gated candidate only≤16 cheap-fill candidates; no global minimum or all-unit unsafe claim.',
           'inverse_instructions':'Apply the frozen36 receipts then this list to RHS; reverse after a later justified residual solve. No full inverse/exhaustion.',
           'new_pivots':k,'total_on_reordered_path':36+k,'residual_shape':m['shape'],'stop':n['stop_reason']}
    (out/'exact-frontier-safe-continuation-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[str(Path(__file__).relative_to(ROOT)),f'code/{Path(__file__).stem}.cpp',f'research/notes/{ID}.md',*pre['input_sha256']]
    r={'id':ID,'all_passed':True,'source_gate_bound_before_arithmetic':True,'new_unit_pivots':k,'total_on_reordered_path':36+k,
       'phase_counts':n['phase_counts'],'stop_reason':n['stop_reason'],'exact_residual_shape':m['shape'],'exact_residual_terms':m['total_terms'],
       'peak_live_terms':n['peak_live_terms'],'max_selected_forecast':n['max_forecast_terms'],'literal_back_add_checks':n['literal_row_update_identity_checks'],
       'unit_parameters_inverted':n['unit_variables_inverted'],'all_prior_results_immutable':True,
       'no_nonmonomial_pivot_gcd_series_Gram_point_prime_or_exhaustion':True,
       'input_sha256':{f:sha(ROOT/f)for f in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
       'work_input_sha256':sha(work/'input.txt'),'work_executable_sha256':sha(work/'native'),'native_seconds':n['native_seconds'],
       'input_compile_worker_seconds':ex['input_compile_worker_total_seconds']}
    assert (work/'input.txt').stat().st_size+sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path);p.add_argument('--finish',action='store_true')
    a=p.parse_args();(finish if a.finish else run)(a.output.resolve(),a.work.resolve())
