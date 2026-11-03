#!/usr/bin/env python3
"""Independent sparse-Fraction replay of the first degree-one J0 unit.

This does not repeat the search.  It proves the full polynomial row identity,
replays the saved one-pivot operation, and independently verifies all rational
diagonal scalings coefficient by coefficient.  Large regenerable matrices are
not copied into the public certificate.
"""
import argparse,gzip,hashlib,json,resource,shlex,shutil,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def poly(line):
 s=line.split();n=int(s[0]);assert len(s)==1+4*n
 p={}
 for k in range(n):
  ex=tuple(map(int,s[2+4*k:5+4*k]));assert ex[2]==0 and ex not in p
  p[ex]=F(s[1+4*k]);assert p[ex]
 return p
def encoded(p):return [[list(e),str(c)]for e,c in sorted(p.items())]
def decoded(p):return {tuple(e):F(c)for e,c in p}
def identity(rows,cofs):
 out={}
 for p,q in zip(rows,cofs):
  for e,c in p.items():
   for d,b in q.items():
    k=tuple(x+y for x,y in zip(e,d));out[k]=out.get(k,F(0))+c*b
 return {e:c for e,c in out.items()if c}
def replay_payload(p):
 rows=list(map(decoded,p['row_polynomials']));cofs=list(map(decoded,p['cofactor_polynomials']))
 target={(0,2,0):F(9),(0,1,0):F(12),(0,0,0):F(4)}
 assert identity(rows,cofs)==target
 return {'all_passed':True,'columns':len(rows),'active_cofactors':sum(bool(q)for q in cofs),'cofactor_terms':sum(map(len,cofs)),'target':'(3e+2)^2'}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path);ap.add_argument('--replay',type=Path);args=ap.parse_args()
 if args.replay:
  with gzip.open(args.replay,'rt')as f:p=json.load(f)
  print(json.dumps(replay_payload(p),indent=2));return
 assert args.output;out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 source=ROOT/'work/nonsquare-J0-pid-v63/primitive-chain-v1/stage-1/residual.txt';w=ROOT/'work/nonsquare-J0-pid-v63/degree1-row4-v1';norm=w/'scalar-primitive';previous=ROOT/'reports/S61_NONSQUARE_J0_SECOND_UNIT_20261103_v63/report.json'
 assert sha(source)==json.loads(previous.read_text())['new_residual_sha256']
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[21,55,1];f.readline();f.readline();selected=[]
  for i in range(21):
   for j in range(55):
    s=f.readline()
    if i==4:selected.append(s)
 assert (w/'row.txt').read_text()=='1 55 0\n'+''.join(selected)
 rows=list(map(poly,selected))
 with(w/'witness.txt').open()as f:
  assert int(f.readline())==55;cofs=[poly(f.readline())for _ in range(55)];assert not f.read().strip()
 assert max(sum(ex)for q in cofs for ex in q)==1
 payload={'source_shape':[21,55],'target_row':4,'base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]','row_polynomials':list(map(encoded,rows)),'cofactor_polynomials':list(map(encoded,cofs)),'target':'(3e+2)^2'}
 checked=replay_payload(payload);assert checked['active_cofactors']==39 and checked['cofactor_terms']==45
 with gzip.open(out/'identity.json.gz','wt',compresslevel=9)as f:json.dump(payload,f,separators=(',',':'))
 # A new saved-witness replay checks the actual cokernel reduction, not merely
 # the scalar row identity.  No coefficient-field or Smith search is repeated.
 native=ROOT/'code/polynomial_one_witness_pivot_20261103_v62.cpp';units=ROOT/'work/nonsquare-J0-pid-v63/units.txt'
 flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
 exe=out/'replay-pivot';subprocess.run(['c++','-std=c++17','-O3',str(native),*flags,'-lgmp','-o',str(exe)],check=True)
 replay=out/'replayed-residual.txt'
 def limits():resource.setrlimit(resource.RLIMIT_AS,(1800000000,1800000000));resource.setrlimit(resource.RLIMIT_CPU,(40,41))
 t=time.monotonic()
 with(out/'fresh-pivot.log').open('w')as f:
  p=subprocess.run([str(exe),str(source),str(replay),str(units),'1',str(w/'witness.txt'),'4'],stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits,timeout=45)
 assert p.returncode==0 and 'EXACT_WITNESS_REPLACE'in(out/'fresh-pivot.log').read_text()and'PASS 1'in(out/'fresh-pivot.log').read_text()
 assert sha(replay)==sha(w/'residual.txt')=='80179a8145d79c1046764da2aabc0c51dfed774e67cdd8d0b969b697732c2c00'
 pivotseconds=time.monotonic()-t;replay.unlink();exe.unlink()
 # Streaming full QQ diagonal identity N_ij=r_i M_ij c_j.  This is independent
 # Fraction arithmetic, and does not assume numerical factors are p-units.
 with(norm/'factors.txt').open()as f:
  assert list(map(int,f.readline().split()))==[20,54];rf=[];cf=[]
  for i in range(20):s=f.readline().split();assert s[:2]==['ROW',str(i)];rf.append(F(s[2]))
  for j in range(54):s=f.readline().split();assert s[:2]==['COLUMN',str(j)];cf.append(F(s[2]))
  assert not f.read().strip()
 assert all(rf+cf);terms=0;t=time.monotonic()
 with(w/'residual.txt').open()as f,(norm/'residual.txt').open()as g:
  for _ in range(3):assert f.readline()==g.readline()
  for i in range(20):
   for j in range(54):
    p=poly(f.readline());q=poly(g.readline());scale=rf[i]*cf[j]
    assert q=={e:c*scale for e,c in p.items()};assert all(c.denominator==1 for c in q.values());terms+=len(p)
  assert not f.read().strip()and not g.read().strip()
 assert terms==351404
 for src,name in [(w/'witness.txt','witness.txt'),(w/'solver.log','search.log'),(w/'pivot.log','producer-pivot.log'),(norm/'factors.txt','QQ-factors.txt'),(norm/'solver.log','producer-QQ-scaling.log'),(units,'units.txt')]:shutil.copyfile(src,out/name)
 inputs=[Path(__file__),source,previous,w/'witness.txt',w/'row.txt',native,ROOT/'code/s61_nonsquare_J0_degree1_unit_span_20261103_v63.cpp',ROOT/'code/polynomial_QQ_scalar_primitive_20261103_v63.cpp',norm/'factors.txt']
 report={'id':'S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_V63','all_passed':True,**checked,'base_ring':payload['base_ring'],'J_not_inverted':True,'unknown_core_factor_inverted':False,'source_shape':[21,55],'target_row':4,'new_residual_shape':[20,54],'total_safe_pivots_from57basis':37,'independent_sparse_QQ_full_row_identity':True,'fresh_saved_witness_pivot_replay_seconds':pivotseconds,'fresh_pivot_output_byte_identical':True,'QQ_diagonal_identity_independently_verified_all_terms':terms,'QQ_diagonal_check_seconds':time.monotonic()-t,'numeric_scalings_are_only_QQ_units_no_mixed_characteristic_claim':True,'raw_residual_sha256':sha(w/'residual.txt'),'raw_residual_bytes':(w/'residual.txt').stat().st_size,'normalized_residual_sha256':sha(norm/'residual.txt'),'normalized_residual_bytes':(norm/'residual.txt').stat().st_size,'large_outputs_omitted_as_fast_reproducible':True,'full_coker_unit_proved':False,'J0_defect_support_finite_proved':False,'loader_failure_before_activated_run_excluded':True,'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','new_residual_shape','active_cofactors','cofactor_terms','QQ_diagonal_identity_independently_verified_all_terms','seconds']},indent=2))
if __name__=='__main__':main()
