#!/usr/bin/env python3
"""One cached-full9 actual honest T2 minimal-generator computation."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM=H.STEM
def parse(v):return H.parse(v)
def checked(v):return H.guard(v)
def main(out,gate,sourcepeer=None):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized_cached_value_worker_started'})
 gp=ROOT/gate;g=json.loads(gp.read_text());assert g['all_passed'] and not g['actual741_Y_or_rank_started'];assert g['code_sha256']==H.sha(Path(H.__file__))
 ep=gp.parent/'actual-T2-mu-input-source-and-cache-gate.json';assert H.sha(ep)==g['evidence_sha256'][ep.name];d=json.loads(ep.read_text())
 pins={gate:H.sha(gp),str(ep.relative_to(ROOT)):H.sha(ep)}
 for p,s in g['input_sha256'].items():assert H.sha(ROOT/p)==s;pins[p]=s
 if sourcepeer:
  pp=ROOT/sourcepeer;peer=json.loads(pp.read_text());assert peer['all_passed'];pins[sourcepeer]=H.sha(pp)
 # Existing cached ds plus exact original84 restoration arithmetic, no old8 re-eval.
 memo={b['original_node']:parse(b['origin_and_ds']) for b in d['all_exact_cached_leaf_bindings_origin_and_ds']}
 arithmetic=dict(d['old_exact84_original_arithmetic_nodes'])
 for i,n in sorted(arithmetic.items()):
  if n[0]=='c':v={0:F(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=add(memo[n[1]],memo[n[2]])
  elif n[0]=='n':v=scale(memo[n[1]],-1)
  else:assert n[0]=='m';v=H.pairmul(memo[n[1]],memo[n[2]])
  memo[i]=checked(v)
 atom={i:memo[i] for i in d['old79_atom_ids']};typed={}
 for i,n in d['old640_typed_closure_nodes']:
  if n[0]=='c':v={0:F(n[1],n[2])} if n[1] else {}
  elif n[0]=='actual_C0_scalar_function':v=atom[n[1]]
  elif n[0]=='sum':
   v={}
   for j,a,b in n[1]:v=add(v,typed[j],F(a,b))
  else:
   assert n[0]=='product';v={0:F(n[1],n[2])} if n[1] else {}
   for j in n[3]:v=H.pairmul(v,typed[j])
  typed[i]=checked(v)
 P=H.selected(H.PATHS['quartic'],'all640_exact_D0_firstjet81_columns');assert len(P)==640
 X0=H.selected(H.PATHS['origin'],'full_X0_56_by768');newell=[];compiler_counts=[]
 for j,raw in enumerate(d['old640_source_columns']):
  v={};orig={}
  for row,i in raw:
   p=typed[i];c=p.get(1,F(0))/5
   if c:v[row]=c
   c=p.get(0,F(0))
   if c:orig[row]=c
  for row in range(8):
   c=F(P[j][9*(row+1)+8])/5
   if c:v[48+row]=c
  assert all(orig.get(row,F(0))==F(X0[row][j]) for row in range(48))
  assert all(F(X0[row][j])==0 for row in range(48,56));newell.append(checked(v))
 compiler_live=sum(len(p) for p in memo.values())+sum(len(p) for p in typed.values())+sum(map(len,newell))+sum(len(c) for c in P)+sum(map(len,atom.values()))+10000
 assert compiler_live<150000;del P,memo,typed,atom,arithmetic,orig
 ell_evidence={'coordinate':'ell atfixedz7,rA, ds/5','actual_first640_ell_columns':[srec(v) for v in newell],
  'all640_oldB_origins_match_X0':True,'all_P8_ell_from_original_full9_raw640_firstjet81':True,
  'old8_derivatives_reused_without_recalculation':True,'source_compiler_live_slot_bound':compiler_live}
 p1=out/'phase1-actual-first640-ell-from-existing-full9-cache.json';H.dump(p1,ell_evidence)
 old640=H.selected(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9');old128=H.selected(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9')
 old=[]
 for col in old640+old128:
  q={}
  for row,p in col:
   for k,c in parse(p).items():
    if k:q[row*9+k-1]=c
  old.append(q)
 assert sum(map(len,old))<=4271;del old640,old128
 L=list(map(parse,d['L0_Res4_Ntranspose_columns']));S0=list(map(parse,d['native_S0_columns']));S1=list(map(parse,d['native_S1ell_columns']))
 # These128 intermediates are deliberately not claimed individualactual X1ell.
 effective=[checked(mv(L,v)) for v in S1];Xcols=[{i:F(X0[i][j]) for i in range(56) if F(X0[i][j])} for j in range(768)];del X0
 M=list(map(parse,d['Mnew0_columns']));G2=list(map(parse,d['full_G2_linear234_by102_columns']));norms={p:parse(v) for p,v in d['active42_normalized_basis']};reps={p:parse(v) for p,v in d['active42_original_source_representations']}
 ranknorm={};rankrepr={};originalbasis=[];selected=[];factors=[];receipts=[];equation_hash=__import__('hashlib').sha256();peak=compiler_live
 fixed_slots=sum(map(len,old))+sum(map(len,newell))+sum(map(len,effective))+sum(map(len,M))+sum(map(len,G2))+sum(map(len,norms.values()))+sum(map(len,reps.values()))+sum(map(len,L))+sum(map(len,S0))+sum(map(len,S1))+759+97+10000
 H.dump(out/'phase.json',{'phase':'all640ell_complete_before741_unique_G2_solves','compiler_slots':compiler_live,'fixed_rank_inputs_bound':fixed_slots})
 for k,record in enumerate(d['all741_zero_origin_source_combinations']):
  combo=parse(record['original768']);assert not mv(Xcols,combo)
  cn={j-640:c for j,c in combo.items() if j>=640};assert not mv(S0,cn)
  x1=mv(old,combo);ell={}
  for j,c in combo.items():ell=add(ell,newell[j] if j<640 else effective[j-640],c)
  for row,c in ell.items():x1=add(x1,{row*9+8:c})
  rhs={}
  for i,c in x1.items():
   row,var=divmod(i,9)
   for target,a in M[row].items():rhs=add(rhs,{target*9+var:a*c})
  checked(rhs);y=H.solve_G2(rhs,G2,norms,reps)
  equation_hash.update(json.dumps([k,srec(rhs),srec(y)],separators=(',',':')).encode())
  # Original-column rank certificate, not a second retained denseY matrix.
  v=dict(y);representation={};steps=[]
  while v and min(v) in ranknorm:
   p=min(v);c=v[p];v=checked(add(v,ranknorm[p],-c));representation=checked(add(representation,rankrepr[p],c));steps.append([p,str(c)])
  if v:
   p=min(v);a=v[p];j=len(originalbasis);originalbasis.append(y);selected.append(k)
   ranknorm[p]=checked(scale(v,1/a));rankrepr[p]=checked(scale(add({j:F(1)},representation,-1),1/a))
   assert mv(originalbasis,rankrepr[p])==ranknorm[p]
   receipts.append({'kernel_column':k,'original_source_column':record['source_column'],'basis_column':j,'pivot_row':p,'divide_by':str(a),'subtract':steps})
   factor={j:F(1)}
  else:factor=representation
  assert mv(originalbasis,factor)==y;factors.append(factor)
  live=fixed_slots+sum(map(len,ranknorm.values()))+sum(map(len,rankrepr.values()))+sum(map(len,originalbasis))+sum(map(len,factors))+sum(len(r['subtract']) for r in receipts)+2*len(x1)+2*len(rhs)+4*len(y)+1024
  peak=max(peak,live);assert live<=150000,('actual_live_slot_guard',k,live)
  if (k+1)%64==0:H.dump(out/'phase.json',{'phase':'741_unique_G2_and_rank_stream','completed_columns':k+1,'rank':len(originalbasis),'peak_slots':peak})
 rank=len(originalbasis);mu=105-rank;assert 94<=rank<=102 and 3<=mu<=11
 prows=[r['pivot_row'] for r in receipts];rowidx={p:i for i,p in enumerate(sorted(prows))}
 minor=[{rowidx[i]:c for i,c in v.items() if i in rowidx} for v in originalbasis]
 sign=(-1)**sum(prows[i]>prows[j] for i in range(rank) for j in range(i+1,rank))
 evidence={'source_column_order':'741 fixed triangular zero-origin combinations, same index as frozeninputgate',
  'full_G2_equations_verified':741*234,'all741_intrinsic_nativeS0_zero_rechecked':True,
  'all741_unique_Y0_source_equation_fingerprint':equation_hash.hexdigest(),
  'actual_selected_Y0_source_columns':selected,'actual_selected_Y0_columns102':[srec(v) for v in originalbasis],
  'fixed_rank_receipts':receipts,'normalized_basis_rows':[[p,srec(v)] for p,v in sorted(ranknorm.items())],
  'normalized_basis_original_selected_source_representations':[[p,srec(v)] for p,v in sorted(rankrepr.items())],
  'all741_factorizations_in_original_selected_Y0_basis':[srec(v) for v in factors],
  'fixed_original_minor_rows':sorted(prows),'fixed_original_minor_columns':selected,
  'fixed_original_minor_sparse_columns':[srec(v) for v in minor],
  'fixed_original_minor_determinant_nonzero_product':{'sign':sign,'QQ_factors':[r['divide_by'] for r in receipts]},
  'actual_Y0_rank':rank,'actual_mu_B_T2':mu,'formula':'mu_B(T2)=132-27-rank(Y0)=105-rank(Y0)',
  'actual_T2_length_or_K_exponents_computed':False,'new_implicit_physical_or_N1_coefficients':False,
  'peak_owned_QQ_slot_bound':peak,'limits':d['limits']}
 p2=out/'actual741-Y0-fixed-rank-and-T2-minimal-generators.json';H.dump(p2,evidence)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'actual_Y0_rank':rank,'actual_mu_B_T2':mu,'all741_unique_G2linear_equations_checked':True,
  'all741_original_source_column_factorizations_checked':True,'new_implicit_physical_or_N1_coefficients':False,
  'actual_T2_length_or_K_exponents_computed':False,'input_sha256':pins,
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{p.name:H.sha(p) for p in [p1,p2]},'peak_owned_QQ_slot_bound':peak,
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 H.dump(out/'phase.json',{'phase':'FINAL_741_actual_origin_values_and_exact_minimal_generators','rank':rank,'mu':mu});H.dump(out/'report.json',report)
 print(json.dumps({k:report[k] for k in ['all_passed','actual_Y0_rank','actual_mu_B_T2','peak_owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer');a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)
