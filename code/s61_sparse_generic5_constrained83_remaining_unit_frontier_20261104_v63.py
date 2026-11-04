#!/usr/bin/env python3
"""ONE40s source-only remaining312 scan, no arithmetic row operation."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_REMAINING_UNIT_FRONTIER_20261104_v63'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_SHORTLIST_CONTINUATION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def prepare_input(out,work):
    folder=ROOT/'reports'/SOURCE;rp=folder/'report.json';r=read(rp)
    assert r['all_passed'] and sha(rp)=='836bbe72ed6c9dcae7c67e748f6400dcca20e61a93f02421c31b2cbe3a1ac6fb'
    mp=folder/'residual-matrix.json';fp=folder/'selection-frontier.json'
    for p in(mp,fp):assert sha(p)==r['evidence_sha256'][p.name]
    m=read(mp);f=read(fp)
    assert m['shape']==[200,47] and m['total_terms']==143477 and f['eligible_candidates']==328 and len(f['records'])==16 and all(x[2]=='terms'for x in f['records'])
    inp=work/'input.txt'
    with inp.open('w')as o:
        o.write('5 236 83 200 47 5483\n');o.write(' '.join(map(str,m['row_original_ids']))+'\n');o.write(' '.join(map(str,m['column_original_ids']))+'\n')
        for row,col,p in m['sparse_entries']:
            o.write(' '.join(map(str,[row,col,len(p['terms']),*p['shift']]))+'\n')
            for t in p['terms']:o.write(str(Q(p['content'])*int(t[0]))+' '+' '.join(map(str,t[1:]))+'\n')
        o.write('16\n')
        for row,col,*_ in f['records']:o.write(f'{row} {col}\n')
    return inp,[rp,mp,fp]
def run(out,work):
    started=time.monotonic();signal.alarm(40);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    inp,paths=prepare_input(out,work);gc.collect();cpp=ROOT/f'code/{Path(__file__).stem}.cpp'
    engine=ROOT/'code/s61_sparse_generic5_constrained83_shortlist_continuation_20261104_v63.cpp'
    pre={'id':ID,'source_only_no_row_operations':True,'source_terms':143477,'exclude_exact_previous16_positions':True,
         'remaining312_once_only':True,'stop_first_FULL_safe':True,'save_every_completed_forecast_and_cursor':True,
         'native_soft_seconds':32,'hard_total_seconds':40,'owned_polynomial_terms_limit':150000,
         'RAW_input_plus_outputs_limit':32*1024*1024,'memory_limit':2_000_000_000,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
         'algorithm_source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),cpp,engine]},
         'work_input_sha256':sha(inp)}
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');binary=work/'source-scan'
    with(out/'compile.log').open('wb')as log:
        subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint) -lgmp','compile',str(cpp),str(binary)],
                       check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,39-(time.monotonic()-started)))
    with(out/'native.log').open('wb')as log:
        subprocess.run([str(binary),str(inp),str(out)],check=True,stdout=log,stderr=subprocess.STDOUT,
                       timeout=max(1,39-(time.monotonic()-started)))
    n=read(out/'source-frontier-summary.json')
    assert n['all_passed'] and n['source_owned_terms_unchanged']==n['peak_owned_polynomial_terms']==143477
    ex={'source_worker_passed':True,'input_compile_scan_total_seconds':time.monotonic()-started}
    (out/'execution-status.json').write_text(json.dumps(ex,indent=2)+'\n');print(json.dumps({**ex,**n}))
def finish(out,work):
    assert not(out/'report.json').exists();pre=read(out/'preflight.json');n=read(out/'source-frontier-summary.json');ex=read(out/'execution-status.json')
    assert ex['source_worker_passed'] and ex['input_compile_scan_total_seconds']<40
    for sec in('input_sha256','algorithm_source_sha256'):
        for f,h in pre[sec].items():assert sha(ROOT/f)==h
    assert sha(work/'input.txt')==pre['work_input_sha256']
    for k in range(n['completed_remaining']):
        f=read(out/f'forecast-{k:03d}.json');c=read(out/f'cursor-{k:03d}.json')
        assert f['completed_source_forecast'] and f['no_coefficient_update_executed'] and c['next_cursor']==k+1
    files=[str(Path(__file__).relative_to(ROOT)),f'code/{Path(__file__).stem}.cpp',f'research/notes/{ID}.md',
           *pre['input_sha256'],*pre['algorithm_source_sha256']]
    r={'id':ID,'all_passed':True,'source_only':True,'exact36_result_immutable':True,
       'completed_remaining312':n['completed_remaining'],'FULL_safe_found':n['FULL_safe_found'],
       'all328_unsafe_under_frozen_forecast':n['all328_unsafe_under_frozen_forecast'],'stop_reason':n['stop_reason'],
       'peak_owned_terms':143477,'no_pivot_gcd_series_or_new_source_coefficients':True,
       'input_sha256':{f:sha(ROOT/f)for f in files},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
       'work_input_sha256':sha(work/'input.txt'),'work_executable_sha256':sha(work/'source-scan'),
       'total_seconds':ex['input_compile_scan_total_seconds']}
    assert (work/'input.txt').stat().st_size+sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path)
    p.add_argument('--finish',action='store_true');a=p.parse_args()
    (finish if a.finish else run)(a.output.resolve(),a.work.resolve())
