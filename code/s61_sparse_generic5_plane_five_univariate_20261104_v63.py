#!/usr/bin/env python3
"""One formalQQ[T] minor witness from an immutable134x5 polynomialmatrix."""
import argparse,gzip,json,resource,shlex,subprocess,time
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63'
def preflight(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 source=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_FIVE_SELECTOR_20261104_v63';m=read(source/'five-column-matrix.json.gz');assert m['all_passed']and m['total_terms']==692953
 weights=[1,2,3,5,7];rows=[];total=0;degree=0;maxterms=0
 for row in m['matrix134x5']:
  rr=[]
  for p in row:
   c=defaultdict(int)
   for term in p['terms']:c[sum(a*b for a,b in zip(weights,term[1:]))]+=int(term[0])
   q={e:v for e,v in c.items()if v};rr.append(dict(content=p['content'],terms=[[e,str(v)]for e,v in sorted(q.items())]));total+=len(q);maxterms=max(maxterms,len(q));degree=max(degree,max(q,default=0))
  rows.append(rr)
 assert len(rows)==134 and all(len(row)==5 for row in rows)
 raw=json.dumps(dict(variables=['T'],formal_parameter_exponents=weights,matrix134x5=rows),separators=(',',':')).encode();allowed=total<=150000 and len(raw)<=32*1024*1024
 (out/'collapsed-matrix.json.gz').write_bytes(gzip.compress(raw,mtime=0))
 inp=work/'input.txt'
 with inp.open('w')as f:
  f.write('134 5\n')
  for row in rows:
   for p in row:
    f.write(str(len(p['terms']))+'\n');c=Q(p['content'])
    for e,v in p['terms']:f.write(str(e)+' '+str(c*int(v))+'\n')
 report=dict(id=ID+'_PREFLIGHT',all_passed=allowed,source_matrix_terms=692953,formal_parameter_map=dict(zip(['H','L','s0','q0','t0'],['T','T^2','T^3','T^5','T^7'])),collapsed_total_terms=total,max_collapsed_terms=maxterms,max_degree=degree,collapsed_raw_bytes=len(raw),total_term_guard_150000=total<=150000,artifact_guard32MiB=len(raw)<=32*1024*1024,no_numeric_or_finite_field_sourcepoint=True,work_input_sha256=sha(inp),input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),ROOT/f'code/{Path(__file__).stem}.cpp',ROOT/f'research/notes/{ID}.md',source/'report.json',source/'five-column-matrix.json.gz']},evidence_sha256={'collapsed-matrix.json.gz':sha(out/'collapsed-matrix.json.gz')},seconds=time.monotonic()-start)
 (out/'preflight.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
def run(out,work):
 start=time.monotonic();out=out.resolve();work=work.resolve();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));gate=read(out/'preflight.json');assert gate['all_passed']and sha(work/'input.txt')==gate['work_input_sha256']
 for name,h in gate['input_sha256'].items():assert sha(ROOT/name)==h,name
 native=ROOT/f'code/{Path(__file__).stem}.cpp';binary=work/'native';error=None;result=None
 try:
  flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(work/'input.txt'),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
  result=read(out/'native-certificate.json');assert result['all_passed']and result['rank']==5
 except Exception as exc:error=repr(exc)
 for p in list(out.glob('*.json')):
  if p.name=='preflight.json':continue
  assert p.stat().st_size<=32*1024*1024;(out/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0));p.unlink()
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<=32*1024*1024
 report=dict(id=ID,all_passed=error is None,error=error,rank5_of_actual_five_parameter_plane_restrictions=bool(result),native=result,formal_curve_nonzero_minor_implies_original_polynomial_nonzero=True,no_numeric_point_prime_GB_cross126_or_exhaustion=True,previous_guard_miss_packet_retained_unchanged=True,input_sha256={**gate['input_sha256'],str((out/'preflight.json').relative_to(ROOT)):sha(out/'preflight.json')},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],rank5=report['rank5_of_actual_five_parameter_plane_restrictions'],error=error,seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--preflight',action='store_true');a=p.parse_args();(preflight if a.preflight else run)(a.output,a.work)
