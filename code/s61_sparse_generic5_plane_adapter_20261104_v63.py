#!/usr/bin/env python3
"""ONE capped symbolic18 kernel and exact plane/normal42 affine adapter."""
import argparse,gzip,json,resource,shlex,subprocess,time
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_PLANE_ADAPTER_20261104_v63'
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    pref=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';gate=read(pref/'report.json');assert gate['all_passed']
    for n,h in gate['input_sha256'].items():assert sha(ROOT/n)==h,n
    data=read(pref/'known-five-input.json.gz');sel=data['saved_selector'];inp=work/'input.txt'
    with inp.open('w')as f:
        f.write('5 51 20 45\n')
        polys=[p for row in data['L51x20']for p in row]+[p for row in data['N45x20']for p in row]+[p for key in['q_center','Q52_minus_center','Q51','Q52']for p in data[key]]
        for p in polys:
            f.write(str(len(p))+'\n')
            for e,c in p:f.write(c+' '+' '.join(map(str,e))+'\n')
        for key in['L_rows','L_pivot_columns','L_free_columns']:f.write(' '.join(map(str,sel[key]))+'\n')
    with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
    native=ROOT/'code/s61_sparse_generic5_plane_adapter_20261104_v63.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'plane-adapter-native';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    result=None;error=None
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        result=read(out/'native-certificate.json');assert result['all_passed']and result['rank_L_exact']==18 and result['kernel51x2_verified']and result['plane_rank3']
        assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',pref/'report.json',pref/'certificate.json',pref/'known-five-input.json.gz']
    checkpoint=read(out/'kernel20x2.json')if(out/'kernel20x2.json').exists()else None
    rank18=bool(checkpoint and checkpoint['all_passed']and checkpoint['all51x2_kernel_equations_zero']and checkpoint['T_free_equals_DeltaI2'])
    report=dict(id=ID,all_passed=error is None,error=error,status='exact_rank18_plane3_and_normal42_adapter'if result else'bounded_frontier',
        completed_exact_rank18=rank18,completed_all51_kernel=rank18,plane_and_normal42_complete=bool(result),native_log_tail=(out/'native.log').read_text().splitlines()[-3:]if(out/'native.log').exists()else[],
        completed_blocks=[name for name in['old18.json','kernel20x2.json','plane-normal42-adapter.json','native-certificate.json']if(out/name).exists()],
        no_point_prime_GB_J0_or134_normal_substitution=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],error=error,completed=report['completed_blocks'],seconds=report['seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
