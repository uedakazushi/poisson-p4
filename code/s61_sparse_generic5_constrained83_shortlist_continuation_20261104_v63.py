#!/usr/bin/env python3
"""Source-only gate, then ONE conditionally authorized bounded continuation."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,gc,hashlib,json,resource,signal,subprocess,time
from s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63 import stat,supp,common_denominator_primitive_bytes,encoded_bytes
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_SHORTLIST_CONTINUATION_20261104_v63'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_REORDERED_BRANCH_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def build_gate(out,work):
    folder=ROOT/'reports'/SOURCE;rp=folder/'report.json';r=read(rp)
    assert r['all_passed'] and sha(rp)=='218e79735f178bf6e1cabeb869e7cb1167872896a3be3e9a8e22f91f880a4e8f'
    mp=folder/'residual-matrix.json';assert sha(mp)==r['evidence_sha256'][mp.name]
    fp=folder/'selection-frontier.json';assert sha(fp)==r['evidence_sha256'][fp.name];front=read(fp)
    saved_safe=sorted([x for x in front['records']if x[2]=='safe'],key=lambda x:(x[3],x[0],x[1]));assert len(saved_safe)==4
    row,col,_,savedpeak,_=saved_safe[0];assert(row,col)==(12,67) and savedpeak==90349
    m=read(mp);assert m['shape']==[206,53] and m['total_terms']==85459 and m['nonzero_entries']==6078
    A={(i,j):p for i,j,p in m['sparse_entries']};rr=m['row_original_ids'];cc=m['column_original_ids'];live=85459
    first=A[row,col];pe=supp(first)[0];unit=Q(first['content'])*int(first['terms'][0][0]);assert unit==Q(-32,405) and pe==(0,0,0,0,1)
    pcols=[j for j in cc if j!=col and(row,j)in A];targets=[i for i in rr if i!=row and(i,col)in A]
    R=sum(len(A[row,j]['terms'])for j in cc if(row,j)in A);logical=live;peak=live
    receipt=4096;residual=8192+32*(len(rr)-1)*(len(cc)-1);updates=[]
    for(i,j),p in A.items():
        if i!=row and j!=col:residual+=encoded_bytes(p)
    for j in pcols:
        p=A[row,j];sh=[p['shift'][v]-pe[v]for v in range(5)];receipt+=encoded_bytes(p)-len(p['content'])+len(str(Q(p['content'])/unit))+encoded_bytes(sh)-encoded_bytes(p['shift'])+128
    for i in targets:
        factor=A[i,col];nf,cf,mf,_,_=stat(factor);mf=(1<<mf.bit_length())-1;fs=supp(factor);receipt+=encoded_bytes(factor)+512
        for j in pcols:
            old=A.get((i,j));a,ca,ma,_,_=stat(old);ma=(1<<ma.bit_length())-1
            np,cp,mpm,_,_=stat(A[row,j]);mpm=(1<<mpm.bit_length())-1;cp/=unit
            ps=[tuple(e[v]-pe[v]for v in range(5))for e in supp(A[row,j])]
            prod={tuple(f[v]+p[v]for v in range(5))for f in fs for p in ps};union=set(supp(old))|prod;k=len(prod);b=len(union)
            scratch=len(fs)+len(ps)+a+k+b;bound=max(logical+2*k+3*b+64,live+scratch+64);peak=max(peak,bound)
            cprod=cf*cp;D=ca.denominator*cprod.denominator//gcd(ca.denominator,cprod.denominator)
            M=abs(ca.numerator)*(D//ca.denominator)*ma+abs(cprod.numerator)*(D//cprod.denominator)*min(nf,np)*mf*mpm
            os=[0]*5 if old is None else old['shift'];base=[min(os[v],factor['shift'][v]+A[row,j]['shift'][v]-pe[v])for v in range(5)]
            # FLINT sizeinbase may overestimate decimal digits by one: add a
            # full digit per numerator term and denominator to the exact bound.
            newbytes=common_denominator_primitive_bytes(union,base,M,D)+b+1
            residual+=newbytes-(encoded_bytes(old)if old is not None else 0)
            meta=[i,j,a,k,b,logical,bound,M.bit_length(),D.bit_length()];updates.append(meta);receipt+=encoded_bytes(meta)+1;logical+=b-a
        logical-=nf
    logical-=R
    rawbound=receipt+residual+4096+2*1024*1024+128*(len(rr)-1)*(len(cc)-1)
    gate={'source_shape':[206,53],'initial_terms':85459,'saved_FULL_safe_records':saved_safe,
          'chosen_best_among_saved_examined_not_global':[row,col],'pivot_QQ':str(unit),'pivot_parameter_exponent':list(pe),
          'fresh_FULL_whole_coefficient_and_scratch_forecast':peak,'fresh_final_terms_union_upper':logical,
          'fresh_complete_RAW_bound':rawbound,'all_row_update_forecasts':updates,
          'no_actual_coefficient_update_in_gate':True,
          'source_only_gate_passed':peak<=150000 and rawbound<=32*1024*1024}
    assert peak==savedpeak
    (out/'source-only-initial-gate.json').write_text(json.dumps(gate,indent=2)+'\n')
    if not gate['source_only_gate_passed']:return gate,None
    inp=work/'input.txt'
    with inp.open('w')as f:
        f.write('5 236 83 206 53 6078\n');f.write(' '.join(map(str,rr))+'\n');f.write(' '.join(map(str,cc))+'\n')
        for i,j,p in m['sparse_entries']:
            f.write(' '.join(map(str,[i,j,len(p['terms']),*p['shift']]))+'\n')
            for t in p['terms']:f.write(str(Q(p['content'])*int(t[0]))+' '+' '.join(map(str,t[1:]))+'\n')
    gate['RAW_input_bytes_included']=inp.stat().st_size
    gate['fresh_complete_RAW_bound_including_input']=rawbound+inp.stat().st_size
    gate['source_only_gate_passed']=peak<=150000 and gate['fresh_complete_RAW_bound_including_input']<=32*1024*1024
    (out/'source-only-initial-gate.json').write_text(json.dumps(gate,indent=2)+'\n')
    if not gate['source_only_gate_passed']:return gate,None
    return gate,inp
def run(out,work):
    started=time.monotonic();signal.alarm(47);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    gate,inp=build_gate(out,work)
    if inp is None:print(json.dumps({'gate_failed_no_arithmetic':True}));return
    gc.collect() # All parsed source coefficient arrays were local to build_gate.
    cpp=ROOT/f'code/{Path(__file__).stem}.cpp';folder=ROOT/'reports'/SOURCE
    helper=ROOT/'code/s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63.py'
    pre={'id':ID,'gate_passed':True,'initial_source_terms':85459,'prior30_not_replayed':True,
         'first_pivot_from_saved_FULL_safe_examined':[12,67],
         'subsequent_policy':'Enumerate literal constant/parameter-monomial entries by structural_fill,term_pair_majorant,source row,column. FULL-forecast at most first16; choose best SAFE completed among examined. NO global minimum claim. Stop if none of this shortlist is safe, not all-unit unsafe.',
         'unchanged_exact_guards':'Whole-row support/scratch L+2k+3b+64, literal back-add, exact signed Laurent shifts/exponent guards, tight numericLCM/content/exponentRAW bound and full future-frontier reserve. RAW input file is additionally included in32MiB.',
         'hard_total_seconds':47,'native_soft_seconds':39,'hard_live_terms':150000,'hard_RAW_bytes':32*1024*1024,'hard_memory_bytes':2_000_000_000,
         'source_parsed_arrays_released_before_native':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [folder/'report.json',folder/'residual-matrix.json',folder/'selection-frontier.json',helper]},
         'algorithm_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),cpp]},'work_input_sha256':sha(inp)}
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
    print(json.dumps({'source_gate_passed':True,'peak':gate['fresh_FULL_whole_coefficient_and_scratch_forecast'],'RAW_bound_with_input':gate['fresh_complete_RAW_bound_including_input']}),flush=True)
    status={};binary=work/'native'
    try:
        with(out/'compile.log').open('wb')as log:subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,46-(time.monotonic()-started)))
        with(out/'native.log').open('wb')as log:subprocess.run([str(binary),str(inp),str(out)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,46-(time.monotonic()-started)))
        n=read(out/'native-summary.json');assert n['all_passed'];status={'one_conditionally_authorized_worker_passed':True,'gate_compile_worker_total_seconds':time.monotonic()-started}
    except Exception as e:
        status={'one_conditionally_authorized_worker_passed':False,'failure':repr(e),'no_retry':True};raise
    finally:(out/'execution-status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps({**status,'pivot_count':n['pivot_count'],'stop_reason':n['stop_reason']}))
def finish(out,work):
    signal.alarm(47);assert not(out/'report.json').exists();p=read(out/'preflight.json');e=read(out/'execution-status.json');assert e['one_conditionally_authorized_worker_passed'] and e['gate_compile_worker_total_seconds']<47
    for section in('input_sha256','algorithm_source_sha256'):
        for f,h in p[section].items():assert sha(ROOT/f)==h
    assert sha(work/'input.txt')==p['work_input_sha256']
    n=read(out/'native-summary.json');m=read(out/'residual-matrix.json');k=n['pivot_count'];assert n['all_passed'] and m['all_passed'] and m['shape']==[206-k,53-k]
    assert n['source_initial_terms']==85459 and n['prior30_steps_replayed']==0 and n['shortlist_limit']==16 and n['not_global_minimum_claim']
    assert n['peak_live_terms']<=150000 and n['max_forecast_terms']<=150000
    for step in range(k):
        q=read(out/f'pivot-{step:02d}.json');s=read(out/f'selection-{step:02d}.json');assert q['completed'] and q['all_row_update_back_add_identities_checked'] and s['not_global_minimum_claim']
        assert q['forecast']['peak_live_bound']<=150000 and q['forecast']['output_bytes_upper']<=32*1024*1024
        if step==0:assert s['selected']==[12,67] and q['pivot_QQ_coefficient']=='-32/405' and q['pivot_parameter_exponent']==[0,0,0,0,1]
        else:
            assert len(s['records'])<=16 and s['all_shortlisted_considered_or_dominated_before_arithmetic']
            for row,col,status,peak,b in s['records']:
                if status in('safe','dominated'):assert(q['forecast']['peak_live_bound'],*s['selected'])<=(peak,row,col)
    proof={'actual_source':'Exact30-step reordered85459-term residual, all206rows53cols; no prior row operation replay.',
           'selection_scope':'First best FULL-safe saved examined record; later at most16 cheap-fill candidates only. No global minimum or all-unit unsafe conclusion.',
           'exactness':'All row operations invertible over the same five-parameter coefficient field; streamed literal back-add and exact negative Laurent exponents. Constraints transformed, never discarded.',
           'factored_inverse':'Apply source30 receipts then these new receipts to RHS; reverse both lists only after a later justified residual solve. No full inverse or exhaustion.',
           'completed_unit_pairs_this_continuation':k,'total_on_reordered_path':30+k,'residual_shape':m['shape'],'stop':n['stop_reason']}
    (out/'exact-shortlist-continuation-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    files=[str(Path(__file__).relative_to(ROOT)),f'code/{Path(__file__).stem}.cpp',f'research/notes/{ID}.md',*p['input_sha256']]
    r={'id':ID,'all_passed':True,'source_only_gate_passed_before_arithmetic':True,'new_unit_pivots':k,'total_on_reordered_path':30+k,'phase_counts':n['phase_counts'],'stop_reason':n['stop_reason'],'exact_residual_shape':m['shape'],'exact_residual_terms':m['total_terms'],'peak_live_terms':n['peak_live_terms'],'max_selected_forecast':n['max_forecast_terms'],'literal_back_add_checks':n['literal_row_update_identity_checks'],'unit_parameters_inverted':n['unit_variables_inverted'],'no_global_minimum_or_all_unit_unsafe_claim':True,'prior_results_immutable':True,'no_nonmonomial_pivot_gcd_series_Gram_point_prime_or_exhaustion':True,'input_sha256':{f:sha(ROOT/f)for f in files},'evidence_sha256':{x.name:sha(x)for x in out.iterdir()if x.is_file()},'work_input_sha256':sha(work/'input.txt'),'work_executable_sha256':sha(work/'native'),'native_seconds':n['native_seconds'],'gate_compile_worker_seconds':e['gate_compile_worker_total_seconds']}
    assert (work/'input.txt').stat().st_size+sum(x.stat().st_size for x in out.iterdir()if x.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path);p.add_argument('--finish',action='store_true');a=p.parse_args();(finish if a.finish else run)(a.output.resolve(),a.work.resolve())
