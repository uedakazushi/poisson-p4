#!/usr/bin/env python3
"""ONE QQ-unit-priority prefix, then at most one field pivot and small plane."""
import argparse,gzip,json,resource,shlex,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTANT_UNIT_PREFIX_20261104_v63'
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    pref=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';old=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_ADAPTER_20261104_v63';k=read(old/'kernel20x2.json');assert k['all_passed']and k['all51x2_kernel_equations_zero'];data=read(pref/'known-five-input.json.gz')
    for path in[pref/'report.json',old/'report.json']:
        for n,h in read(path)['input_sha256'].items():assert sha(ROOT/n)==h,n
    inp=work/'input.txt'
    def write(f,p):
        f.write(str(len(p))+'\n')
        for e,c in p:f.write(c+' '+' '.join(map(str,e))+'\n')
    def primitive(p):return [[t[1:],str(Q(p['content'])*int(t[0]))]for t in p['terms']]
    with inp.open('w')as f:
        f.write('5 51 20 45\n')
        for p in[p for row in data['L51x20']for p in row]+[p for row in data['N45x20']for p in row]+data['q_center']+data['Q52_minus_center']:write(f,p)
        write(f,primitive(k['Delta']))
        for row in k['T_numerators20x2']:
            for p in row:write(f,primitive(p))
        f.write(' '.join(map(str,k['free_columns']))+'\n')
    with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
    native=ROOT/'code/s61_sparse_generic5_constant_unit_prefix_20261104_v63.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'constant-unit-prefix';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    result=None;error=None
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        result=read(out/'native-certificate.json');assert result['all_passed'];assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',pref/'report.json',pref/'known-five-input.json.gz',old/'report.json',old/'kernel20x2.json']
    prefix=read(out/'constant-prefix.json')if(out/'constant-prefix.json').exists()else None
    tail=None
    if prefix:
        k=prefix['constant_QQ_pivots'];ps=[p for row in prefix['reduced51x20'][k:]for p in row[k:]]
        tail=dict(shape=[51-k,20-k],nonzero=sum(bool(p['terms'])for p in ps),terms=sum(len(p['terms'])for p in ps),maxterms=max(len(p['terms'])for p in ps),degree=max((sum(t[1:])for p in ps for t in p['terms']),default=0))
    report=dict(id=ID,all_passed=error is None,error=error,status=result['status']if result else'bounded_frontier',constant_QQ_pivots=result['constant_QQ_pivots']if result else None,
        residual_tail_statistics=tail,
        kernel_complete=bool(result and result['kernel_complete']),plane_normal42_complete=bool(result and result['plane_complete']),
        generic_coefficient_field_only_not_FULLU=True,no_oldCramer_retry_gcd_point_prime_GB_or134_normal_substitution=True,
        completed_blocks=[n for n in['constant-prefix.json','small-kernel.json','plane-normal42-adapter.json','native-certificate.json']if(out/n).exists()],
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],pivots=report['constant_QQ_pivots'],kernel=report['kernel_complete'],plane=report['plane_normal42_complete'],seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
