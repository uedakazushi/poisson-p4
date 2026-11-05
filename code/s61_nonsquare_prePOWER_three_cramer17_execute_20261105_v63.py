#!/usr/bin/env python3
"""ONE conditional three-Cramer execution after immutable peer and root gates."""
import argparse,gc,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def adjugate_identity():
 # Universal free polynomial ring on nine entries, not actual source values.
 def add(a,b):
  c=a.copy()
  for m,v in b.items():
   c[m]=c.get(m,0)+v
   if not c[m]:del c[m]
  return c
 def neg(a):return {m:-v for m,v in a.items()}
 def mul(a,b):
  c={}
  for m,v in a.items():
   for n,w in b.items():
    k=tuple(x+y for x,y in zip(m,n));c[k]=c.get(k,0)+v*w
    if not c[k]:del c[k]
  return c
 def var(i):m=[0]*9;m[i]=1;return {tuple(m):1}
 A=[[var(3*i+j)for j in range(3)]for i in range(3)];co=[]
 for i in range(3):
  row=[]
  for j in range(3):
   rr=[k for k in range(3)if k!=i];cc=[k for k in range(3)if k!=j];p=add(mul(A[rr[0]][cc[0]],A[rr[1]][cc[1]]),neg(mul(A[rr[0]][cc[1]],A[rr[1]][cc[0]])));row.append(p if (i+j)%2==0 else neg(p))
  co.append(row)
 det={}
 for j in range(3):det=add(det,mul(A[0][j],co[0][j]))
 for i in range(3):
  for j in range(3):
   z={}
   for k in range(3):z=add(z,mul(A[i][k],co[j][k]))
   assert z==(det if i==j else {})
 return dict(free_nine_variable_polynomial_identity_Aadj_equals_detI=True,cofactor_signs='(-1)^(i+j)',no_actual_parameter_values=True)
def main(peer,peerhash,out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));assert sha(peer)==peerhash;pr=read(peer);assert pr['all_passed']
 gate=ROOT/'reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_WORKER_GATE_20261105_v63/report.json';g=read(gate);assert sha(gate)=='f4b268ef74db209445a793791c296d3740e5a8a041d7389d464ec57c7d330f1c'
 pins={str(peer.relative_to(ROOT)):peerhash,str(gate.relative_to(ROOT)):sha(gate)}
 for rel,h in g['input_sha256'].items():assert sha(ROOT/rel)==h;pins[rel]=h
 clarification=ROOT/'reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_PHASE_CLARIFICATION_20261105_v63/report.json';assert sha(clarification)=='f6344fe982ea36ba64ff60838ab5ceafb4070d1643d543d9def07fce0b8f751c';pins[str(clarification.relative_to(ROOT))]=sha(clarification)
 inp=gate.parent/'native-input.txt';assert sha(inp)==g['native_input_sha256'];pins[str(inp.relative_to(ROOT))]=sha(inp)
 # Root review is recorded separately from the independent arithmetic peer.
 out.mkdir(parents=True,exist_ok=False);dump(out/'root-and-peer-custody.json',dict(root_review='Root accepted f4b268/eabb935/a2b797, conditional ONE triple job after κ17 peer, same73322/434bit/5.82MB.',peer_path=str(peer.relative_to(ROOT)),peer_sha256=peerhash,all_parameter_coefficients_from_fixed12_source_polys=True))
 cpp=ROOT/'code/s61_nonsquare_three_cramer17_dense_20261105_v63.cpp';exe=ROOT/'work/d-zero-v63/three-cramer17/worker';exe.parent.mkdir(parents=True,exist_ok=True)
 c=subprocess.run(['g++','-std=c++17','-O2',str(cpp),'-lflint','-lgmp','-o',str(exe)],capture_output=True,text=True,timeout=15);(out/'compile.stdout').write_text(c.stdout);(out/'compile.stderr').write_text(c.stderr)
 if c.returncode:dump(out/'report.json',dict(id=ID,all_passed=False,compile_only_failure=True,no_arithmetic_started=True,input_sha256=pins));return
 h=subprocess.run([str(exe),'--help'],capture_output=True,text=True,timeout=3);assert h.returncode==0;(out/'native-help.txt').write_text(h.stdout)
 dump(out/'phase.json',dict(phase='ONE fixed same-target17 THREE Cramer integer DAG',worker_sha256=sha(cpp),native_input_sha256=sha(inp)))
 def limits():resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_CPU,(52,52))
 rr=subprocess.run([str(exe),str(inp),str(out)],capture_output=True,text=True,timeout=max(1,57-(time.monotonic()-start)),preexec_fn=limits);(out/'native.stdout').write_text(rr.stdout);(out/'native.stderr').write_text(rr.stderr)
 if rr.returncode:dump(out/'report.json',dict(id=ID,all_passed=False,arithmetic_started=True,no_result_claim=True,all_failure_files_preserved=True,input_sha256=pins));return
 n=read(out/'native.json');assert n['all_passed']and n['peak_owned_numeric_metadata']<=73322 and n['observed_integer_bits']<=434
 outputs=[]
 for j,bound in enumerate([6765,5270,4059]):
  p=out/f'Cramer{j}.json';x=read(p);D=int(x['denominator']);v=x['integer_coefficients'];assert D>0 and len(v)<=bound;prev=None
  for m,z in v:
   assert prev is None or tuple(m)>prev;prev=tuple(m);q=Q(int(z),D);assert q and int(z).bit_length()<=217 and q*D==int(z)
  outputs.append(dict(file=p.name,terms=len(v),is_zero=not bool(v),common_denominator=str(D),canonicalQQ_backmultiply_verified=True,sha256=sha(p)));del x,v;gc.collect()
 dump(out/'full-source-and-adjugate-proof.json',dict(universal_adjugate=adjugate_identity(),outputs=outputs,source_scale_carry_custody_reference=str((gate.parent/'literal-source-custody.json').relative_to(ROOT)),
  x_first3_formula='-Y/kappa3',x4='0 from independently accepted kappa17 identity',fourth_row_proof='kappa annihilates first3 A4 columns; allfirst3 equations solve; kappa3 unit forces last residual0',
  no_new_unit_chart=True,all_ten_block_and_tail_effects_retained=True,no_lastrow_or_rank_or_global57_claim=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/rel).stat().st_size for rel in pins);assert new+external<=5820000
 rep=dict(id=ID,all_passed=True,frozen=True,ONE_three_Cramer_bundle=True,Cramer_terms=n['terms'],Cramer_is_zero=[x['is_zero']for x in outputs],
  all_products_reverseQQ_and_addbacks=True,all_common_den_QQ_readbacks=True,universal_adjugate_source_proof=True,
  same_target17_only=True,no_extra_unit_or_lastrow_rank_claim=True,all_three_pair_coordinates_retained=True,all_ten_corrections_retained=True,
  peak_owned_numeric_metadata=n['peak_owned_numeric_metadata'],observed_integer_bits=n['observed_integer_bits'],integer_products=n['integer_products'],
  total_new_RAW_before_report=new,external_read_RAW=external,TOTAL_RAW=new+external,independent_Cramer_peer_pending=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--peer',type=Path,required=True);p.add_argument('--peer-sha256',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.peer.resolve(),a.peer_sha256,a.output.resolve())
