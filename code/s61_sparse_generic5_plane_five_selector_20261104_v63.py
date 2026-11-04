#!/usr/bin/env python3
"""ONE point-free five-functional plane selector, field pivots at mostfive."""
import argparse,gzip,json,resource,shlex,subprocess,time
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_PLANE_FIVE_SELECTOR_20261104_v63'
def dec(p):return{tuple(t[1:]):Q(t[0])for t in p}
def np(p):return{tuple(t[1:]):Q(p['content'])*int(t[0])for t in p['terms']}
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 a=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63';b=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';c=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63';full=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63';assert read(a/'report.json')['all_passed']and read(full/'report.json')['all_passed']
 original=read(a/'plane-W-original.json.gz');sh=read(c/'shared-denominator-kernel-plane.json.gz');data=read(b/'pulled134-input.json.gz');five=read(b/'known-five-input.json.gz');P=dec(sh['common_denominator_P']);C=list(map(np,original['center45']));W=[[dec(row[j])for j in range(2)]for row in original['W45x3']];d=[{tuple(e):Q(v)for e,v in p}for p in five['Q52_minus_center']]
 equations=[]
 for encoded in data['equations']:
  row={}
  for m,e,v in encoded:row.setdefault(tuple(m),{})[tuple(e)]=Q(v)
  equations.append(row)
 assert len(equations)==134
 # A centered, noncollinear triangle can be sent to(1,0),(0,1),(-1,-1).
 # Coefficient order1,a,b,sigma,a2,ab,aSigma,b2,bSigma,sigma2.
 ideal=[[1,0,0,1,-1,1,0,-1,0,-2],[0,1,0,0,-1,2,0,0,0,0],[0,0,1,0,0,2,0,-1,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,1,0]];points=[(Q(1),Q(0),Q(0)),(Q(0),Q(1),Q(0)),(Q(-1),Q(-1),Q(0)),(Q(0),Q(0),Q(-1,2)),(Q(0),Q(0),Q(1))];select=[0,1,2,6,8]
 for row in ideal:
  for x,y,z in points:assert sum(v*t for v,t in zip(row,[1,x,y,z,x*x,x*y,x*z,y*y,y*z,z*z]))==0
 assert [[row[j]for j in select]for row in ideal]==[[int(i==j)for j in range(5)]for i in range(5)]
 # Conservative scalar multiplication count, not a coefficient support/rank assertion.
 pairs=0
 for row in equations:
  for m,p in row.items():
   if len(m)==1:pairs+=len(p)*len(C[m[0]])
   elif len(m)==2:
    i,j=m;pairs+=len(p)*len(C[i])*len(C[j])
    for h in range(2):pairs+=len(p)*(len(W[i][h])*len(d[j])+len(W[j][h])*len(d[i]))
 pre=dict(coefficient_field='QQ(H,L,s0,q0,t0)',functional_columns=['constant','linear_a','linear_b','a_sigma','b_sigma'],five_functionals_I5_on_universal_ideal_basis=True,abstract_ideal_basis5x10=ideal,coefficient_indices=select,points_for_formal_identity_not_source_samples=[[str(v)for v in p]for p in points],centered_triangle_trace_dependency='TORIC_E2_SCALAR_INVERSE_20261101',actual_source_genericity='sparse6Sp4xGm dominance plus invariant noncollinearity/discriminant; no old fixed row minor inference',raw_scalar_multiplication_pair_majorant_partial=pairs,scope_of_pair_count='center constant and polar columns before collecting; J(C)W linear columns counted by actual native construction',hard_seconds=60,hard_memory_bytes=2_000_000_000,hard_artifact_bytes=32*1024*1024,hard_terms_per_poly_or_field_num_den=150000,field_pivots_at_most=5,no_points_primes_GB_cross126=True)
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 inp=work/'input.txt'
 with inp.open('w')as f:
  f.write('5 134 45\n')
  def put(p):
   f.write(str(len(p))+'\n')
   for e,v in sorted(p.items()):f.write(str(v)+' '+' '.join(map(str,e))+'\n')
  put(P)
  for p in C+d:put(p)
  for row in W:
   for p in row:put(p)
  for row in equations:
   f.write(str(len(row))+'\n')
   for m,p in sorted(row.items()):f.write(str(len(m))+' '+(' '.join(map(str,m))+' 'if m else'')+'\n');put(p)
 (out/'input.txt.gz').write_bytes(gzip.compress(inp.read_bytes(),mtime=0));native=ROOT/f'code/{Path(__file__).stem}.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'native';error=None;result=None
 try:
  flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
  result=read(out/'selector-certificate.json');assert result['all_passed']and result['rank']==5
 except Exception as exc:error=repr(exc)
 # Preserve complete bounded prefix only. No rerun or different candidate family.
 oversized=[]
 for p in list(out.glob('*.json')):
  if p.name=='preflight.json':continue
  n=p.stat().st_size
  if n>32*1024*1024:oversized.append(dict(name=p.name,bytes=n,sha256=sha(p)));p.unlink()
  else:(out/(p.name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0));p.unlink()
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
 paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',a/'report.json',a/'plane-W-original.json.gz',b/'report.json',b/'pulled134-input.json.gz',b/'known-five-input.json.gz',c/'shared-denominator-kernel-plane.json.gz',full/'report.json',ROOT/'reports/TORIC_E2_SCALAR_INVERSE_20261101.json',ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json',ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json',ROOT/'reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63/report.json']
 report=dict(id=ID,all_passed=error is None,error=error,rank5_proved=bool(result),selector=result,complete_matrix_saved=(out/'five-column-matrix.json.gz').exists(),abstract_five_functional_injectivity_proved=True,no_inference_from_length5_alone=True,original134_equations_preserved=True,no_full_fiber_exhaustion_or_cross126_rank=True,oversized_artifacts_removed=oversized,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],rank5=report['rank5_proved'],error=error,seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
