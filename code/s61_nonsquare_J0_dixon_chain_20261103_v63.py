#!/usr/bin/env python3
"""Modular-guided exact-Dixon degree1 pivots, checkpointing EVERY unit."""
import argparse,hashlib,json,math,shlex,shutil,subprocess,time
from pathlib import Path
from s61_nonsquare_J0_primitive_chain_20261103_v63 import job,profile
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True,type=Path);ap.add_argument('--expected-sha256',required=True);ap.add_argument('--output',required=True,type=Path)
 ap.add_argument('--budget',type=int,default=300);ap.add_argument('--per-row-cap',type=int,default=30);args=ap.parse_args();assert 0<args.budget<=300
 args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic();source=args.source.resolve();assert sha(source)==args.expected_sha256
 units=ROOT/'work/nonsquare-J0-pid-v63/units.txt';searchsrc=ROOT/'code/s61_nonsquare_J0_degree1_dixon_20261103_v63.cpp';pivotsrc=ROOT/'code/polynomial_one_witness_pivot_20261103_v62.cpp';normsrc=ROOT/'code/polynomial_QQ_scalar_primitive_20261103_v63.cpp'
 flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
 exes={}
 for name,cpp in [('search',searchsrc),('pivot',pivotsrc),('normalize',normsrc)]:
  exe=args.output/name;subprocess.run(['c++','-std=c++17','-O3',str(cpp),*flags,'-lgmp','-o',str(exe)],check=True);exes[name]=exe
 stages=[];state=source;initial,_=profile(state);total=0;stop='budget_or_no_witness'
 while args.budget-(time.monotonic()-started)>8:
  if shutil.disk_usage(ROOT).free<1500000000:stop='workspace_free_space_below1p5GB';break
  if state.stat().st_size>1200000000:stop='matrix_bytes_above1p2GB';break
  shape,rows=profile(state)
  if not shape[0]:stop='full_coker_unit';break
  rows.sort(key=lambda r:(r['max_degree'],r['terms'],r['row']));folder=args.output/f'stage-{len(stages)+1}';folder.mkdir()
  stage={'source':str(state.relative_to(ROOT)),'source_shape':shape[:2],'source_sha256':sha(state),'profiles':rows,'attempts':[]};found=None
  # Retain only candidate row text, never reparsing the full matrix natively.
  for row in rows:
   remaining=args.budget-(time.monotonic()-started)
   if remaining<10:break
   extracted=folder/f'row{row["row"]}.txt'
   with state.open()as f:
    f.readline();f.readline();f.readline();selected=[]
    for i in range(shape[0]):
     for j in range(shape[1]):
      s=f.readline()
      if i==row['row']:selected.append(s)
   extracted.write_text(f'1 {shape[1]} 0\n'+''.join(selected));witness=folder/f'row{row["row"]}.witness';log=folder/f'row{row["row"]}.log'
   cpu=max(1,min(args.per_row_cap,math.floor(remaining)-5));attempt=job([str(exes['search']),str(extracted),str(units),'0','0',str(witness),'8'],log,cpu,min(remaining-4,cpu+2))
   text=log.read_text();attempt.update(row_profile=row,original_polynomial_unit_verified='EXACT_KNOWN_UNIT_DENOMINATOR_LIFT_PASS'in text,extracted_sha256=sha(extracted));stage['attempts'].append(attempt)
   print('D1_STAGE',len(stages)+1,'ROW',row['row'],attempt['status'],'UNIT',attempt['original_polynomial_unit_verified'],'SECONDS',time.monotonic()-started,flush=True)
   if attempt['status']=='completed'and attempt['original_polynomial_unit_verified']and witness.exists():found=(row,witness);break
  remaining=args.budget-(time.monotonic()-started)
  if found and remaining>4:
   row,witness=found;raw=folder/'raw-residual.txt';cpu=max(1,math.floor(remaining)-1)
   applied=job([str(exes['pivot']),str(state),str(raw),str(units),'1',str(witness),str(row['row'])],folder/'pivot.log',cpu,remaining);stage['pivot']=applied
   if applied['status']=='completed':
    newshape,_=profile(raw);assert newshape[0]==shape[0]-1 and newshape[1]in[shape[1]-1,shape[1]]and newshape[2]==1
    log=(folder/'pivot.log').read_text();assert'PASS 1'in log;total+=1;state=raw
    stage.update(raw_shape=newshape[:2],raw_sha256=sha(raw),raw_bytes=raw.stat().st_size,witness=str(witness.relative_to(ROOT)),witness_sha256=sha(witness),original_column_replaced='EXACT_WITNESS_REPLACE'in log,redundant_column_appended='EXACT_WITNESS_APPEND'in log)
    latest={'additional_exact_unit_pivots':total,'shape':newshape[:2],'state':str(state.relative_to(ROOT)),'sha256':sha(state),'source_only_known_units':True}
    (args.output/'latest-confirmed.json').write_text(json.dumps(latest,indent=2)+'\n');print('D1_EXACT_UNIT_EXPORTED',total,'SHAPE',newshape[:2],'SECONDS',time.monotonic()-started,flush=True)
    # Numeric QQ diagonal scaling only; failure leaves the raw exact state safe.
    remaining=args.budget-(time.monotonic()-started)
    if remaining>8:
     normalized=folder/'residual.txt';cpu=max(1,min(60,math.floor(remaining)-1));norm=job([str(exes['normalize']),str(raw),str(normalized),str(folder/'QQ-factors.txt'),'1'],folder/'normalize.log',cpu,min(remaining,cpu+2));stage['QQ_normalization']=norm
     if norm['status']=='completed':
      assert'EXACT_QQ_DIAGONAL_SCALING_PASS'in(folder/'normalize.log').read_text();assert profile(normalized)[0]==newshape;state=normalized
      stage.update(normalized_sha256=sha(state),normalized_bytes=state.stat().st_size)
      latest.update(state=str(state.relative_to(ROOT)),sha256=sha(state),QQ_diagonal_inverse_checks_passed=True);(args.output/'latest-confirmed.json').write_text(json.dumps(latest,indent=2)+'\n')
  stages.append(stage);(folder/'receipt.json').write_text(json.dumps(stage,indent=2)+'\n')
  if'raw_shape'not in stage:stop='no_further_witness_in_budget';break
 final=profile(state)[0][:2]
 report={'id':'S61_NONSQUARE_J0_DIXON_CHAIN_20261103_V63','base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]','J_not_inverted':True,'unknown_core_factor_inverted':False,
  'primitive_cofactor_basis':['1','a','e'],'unit_target_factor_bound':8,'export_after_every_unit':True,'QQ_scalar_diagonal_optimization_only':True,
  'initial_shape':initial[:2],'additional_exact_unit_pivots':total,'last_confirmed_shape':final,'last_confirmed_state':str(state.relative_to(ROOT)),'last_confirmed_sha256':sha(state),
  'stages':stages,'stop_reason':stop,'full_coker_unit_proved':final[0]==0,'J0_rank_defect_finite_proved':final[0]==0,
  'budget_seconds':args.budget,'seconds':time.monotonic()-started,'per_row_cap_seconds':args.per_row_cap,'AS_limit_bytes':1800000000,
  'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in [Path(__file__),source,units,searchsrc,pivotsrc,normsrc,ROOT/'code/s61_nonsquare_J0_primitive_chain_20261103_v63.py']}}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['additional_exact_unit_pivots','last_confirmed_shape','seconds']},indent=2))
if __name__=='__main__':main()
