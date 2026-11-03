#!/usr/bin/env python3
"""Fresh exact bounded unit lift; no all-point rank/coker claim."""
from pathlib import Path
from fractions import Fraction
import argparse,hashlib,json,os,re,shlex,subprocess,tarfile,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_BINARY_NONSQUARE_BOUNDED_UNIT_LIFT_20261103_v62'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def poly(line):
 z=line.split();n=int(z[0]);assert len(z)==1+4*n
 return {tuple(map(int,z[2+4*t:5+4*t])):Fraction(z[1+4*t])for t in range(n)}
def enc(p):
 return str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in sorted(p.items(),reverse=True))+'\n'
def add(A,B):
 C=A.copy()
 for m,c in B.items():
  C[m]=C.get(m,0)+c
  if not C[m]:del C[m]
 return C
def mul(A,B):
 C={}
 for m,c in A.items():
  for n,d in B.items():
   k=tuple(a+b for a,b in zip(m,n));C[k]=C.get(k,0)+c*d
 return {m:c for m,c in C.items()if c}
def parse_matrix(path,metadata=False):
 with Path(path).open()as f:
  n,m,k=map(int,f.readline().split());meta=[]
  if metadata:meta=[list(map(int,f.readline().split())),list(map(int,f.readline().split()))];assert len(meta[0])==k and len(meta[1])==n
  A=[[poly(f.readline())for j in range(m)]for i in range(n)];assert not f.read().strip()
 return n,m,k,meta,A
def run(cmd,out,expected=0,timeout=120):
 r=subprocess.run([str(x)for x in cmd],cwd=ROOT,capture_output=True,text=True,timeout=timeout)
 out.write_text(r.stdout+r.stderr);assert r.returncode==expected,(cmd,r.returncode,r.stdout[-300:],r.stderr[-300:]);return out.read_text()
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 parent=ROOT/'reports/S61_BINARY_NONSQUARE_PRODUCT_ALT_20261103_v62';r=json.loads((parent/'report.json').read_text());assert r['all_passed'];member='work/nonsquare-product-alt-v62/mod-ell/chart-3-ad-residual.txt';unitmember='work/nonsquare-product-alt-v62/mod-ell/chart-3-ad-units.txt'
 with tarfile.open(parent/'evidence.tar.gz')as t:
  for name,dest in [(member,out/'source14.txt'),(unitmember,out/'units.txt')]:
   b=t.extractfile(name).read();assert hashlib.sha256(b).hexdigest()==r['archive_member_manifest'][name]['sha256'];dest.write_bytes(b)
 n,m,k,meta,A=parse_matrix(out/'source14.txt',True);assert(n,m,k)==(14,46,10)
 keep=[j for j in range(m)if A[1][j]];assert len(keep)==19
 negtarget={(2,0,1):Fraction(-5),(1,1,0):Fraction(6),(1,0,0):Fraction(4)}
 scalar=out/'scalar-source.txt';scalar.write_text('20 1 20\n'+''.join(enc(A[1][j])for j in keep)+enc(negtarget));bounds=out/'scalar.bounds';bounds.write_text(' '.join(['3']*19+['0'])+'\n')
 cpp=ROOT/'code/polynomial_bounded_unit_lift_20261103_v62.cpp';onecpp=ROOT/'code/polynomial_one_witness_pivot_20261103_v62.cpp';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));binary=out/'bounded-unit-lift';one=out/'one-pivot'
 for src,dest in [(cpp,binary),(onecpp,one)]:run(['c++','-O2','-std=c++17',src,'-o',dest,*flags],out/(dest.name+'.compile.log'))
 cert=out/'scalar.kernel';log=run([binary,scalar,3,cert,bounds],out/'scalar.log');assert'EXACT_FULL_QQ_UNIT_LIFT_PASS 1 DEGREE 3'in log and'EXACT_QQ_DIXON_SOLVE_DONE'in log
 kn,km,kk,_,KA=parse_matrix(cert,True);assert(kn,km,kk)==(1,20,0)and KA[0][-1]=={(0,0,0):1};co=KA[0][:-1];assert all(sum(ex)<=3 for p in co for ex in p)
 total={}
 for p,c in zip([A[1][j]for j in keep],co):total=add(total,mul(p,c))
 assert add(total,negtarget)=={},'pure Fraction full original scalar identity'
 witness=[{}for j in range(m)]
 for j,c in zip(keep,co):witness[j]=c
 assert sum(map(len,witness))==154;wit=out/'scalar.witness';wit.write_text('46\n'+''.join(enc(p)for p in witness))
 residual=out/'residual13.txt';plog=run([one,out/'source14.txt',residual,out/'units.txt',1,wit,1],out/'pivot.log');assert'EXACT_WITNESS_REPLACE 22 TARGET 1'in plog and'PASS 1'in plog
 rn,rm,rk,_,RA=parse_matrix(residual,True);terms=sum(len(p)for row in RA for p in row);assert(rn,rm,rk,terms)==(13,45,1,225656)
 # Fresh exact fixtures: small reconstruction, forced Dixon, unlucky-prime refusal.
 fixtures={}
 p61=2**61-1
 for tag,number,expect in [('small',1,0),('dixon',123456789012345678901234567890123456789,0),('unlucky_prime',p61,4)]:
  f=out/(tag+'.input');f.write_text('2 1 2\n1 '+str(number)+' 1 0 0\n1 -1 1 0 0\n');bf=out/(tag+'.bounds');bf.write_text('0 0\n');dest=out/(tag+'.kernel');flog=run([binary,f,0,dest,bf],out/(tag+'.log'),expected=expect)
  if expect==0:
   _,_,_,_,S=parse_matrix(dest,True);assert S[0][0]=={(0,0,0):Fraction(1,number)}and S[0][1]=={(0,0,0):1}
   if tag=='dixon':assert'EXACT_QQ_DIXON_SOLVE_DONE'in flog
  else:assert'NO_CERTIFICATE_AT_THIS_PRIME_AND_BOUND'in flog
  fixtures[tag]=dict(expected_exit=expect,passed=True)
 # Stronger isolated target: no polynomial solution within degree <=5.
 direct=out/'isolated-source.txt'
 direct.write_text('47 14 47\n'+''.join(enc(A[i][j])for j in range(m)for i in range(n))+''.join(enc(negtarget if i==1 else{})for i in range(n)))
 db=out/'isolated.bounds';db.write_text(' '.join(['5']*46+['0'])+'\n');dlog=run([binary,direct,5,out/'isolated.kernel',db],out/'isolated.log',expected=4,timeout=120)
 assert'TARGET_FREE 0'in dlog and'SPARSE_COMPLETE_RANK 2577 NULLITY 0'in dlog
 note=ROOT/'research/notes/S61_BINARY_NONSQUARE_BOUNDED_UNIT_LIFT_20261103_v62.md'
 report=dict(id=ID,all_passed=True,meaning='Exact scalar identity, one known-unit replacement pivot, and finite degree bound only; no all-core rank theorem.',known_open='a*d*F*B*C*J!=0; J=5ad-6e-4',source_member=member,source_sha256=sha(out/'source14.txt'),source_shape=[14,46],scalar_target='a*(5*a*d-6*e-4)',target_zero_based_row=1,coefficient_degree_bound=3,witness_terms=154,finite_coefficient_matrix_shape=[423,381],selected_exact_QQ_square_shape=[258,258],native_full_QQ_identity=True,independent_python_Fraction_identity=True,one_invertible_column_replacement=22,additional_unit_pivots=1,total_mod_ell_unit_pivots=11,residual_shape=[13,45],residual_terms=terms,residual_sha256=sha(residual),residual_bytes=residual.stat().st_size,large_residual_omitted=True,isolated_target='(aJ)*e_1',isolated_multiplier_degree_bound=5,isolated_finite_matrix_columns=2577,isolated_modular_rank=2577,isolated_QQ_nonexistence_within_bound=True,fixtures=fixtures,all_point_product_rank41=False,all_point_cubic_rank57=False,quartic_or_ambient_exhaustion=False,flint_version=subprocess.check_output(['pkg-config','--modversion','flint'],text=True).strip(),input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__).resolve(),cpp,onecpp,note,parent/'report.json',parent/'evidence.tar.gz']},seconds=time.monotonic()-started)
 files=[p for p in out.iterdir()if p.is_file()and p not in [binary,one,residual]];report['evidence']={p.name:dict(bytes=p.stat().st_size,sha256=sha(p))for p in files}
 with tarfile.open(out/'evidence.tar.gz','w:gz')as t:
  for p in files:t.add(p,arcname=p.name)
 report['archive_sha256']=sha(out/'evidence.tar.gz');(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 for p in list(out.iterdir()):
  if p.name not in {'evidence.tar.gz','report.json'}:p.unlink()
 print('PASS',ID,report['seconds'])
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();main(args.output.resolve())
