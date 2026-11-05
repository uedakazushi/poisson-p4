#!/usr/bin/env python3
"""Independent formal source fold / symbolic rell AD. No new values or solves."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time,resource,signal
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63';STEM=A+'_REVIEW'
O='AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63';G='AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_INPUT_GATE_20261105_v63'
P={'fold_report':f'reports/{A}/report.json','fold':f'reports/{A}/actual-two-rell-source-folds-and-retained-native-first-operators.json',
'ad_report':f'reports/{O}/report.json','ad':f'reports/{O}/actual-two-by619-old-formal-rell-second-sensitivities.json',
'gate_report':f'reports/{G}/report.json','gate':f'reports/{G}/actual-old-rell-AD-sparse-dependencies-and-budget.json',
'adjoint_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/report.json',
'adjoint':'reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/independent-two-full918-adjoints-and-source-rell-support.json',
'known':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',
'known_peer':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json',
'mini':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',
'mu_gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
'mu_peer':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json',
'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
'old640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
'ell640':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',
'G8':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
'S':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def dot(a,b):return sum((c*b.get(i,Q(0)) for i,c in a.items()),Q(0))
def kr(p):return [[list(k),str(c)] for k,c in sorted(p.items())]
def jetmul(a,b):
 # Sparse square-free two-variable algebra with literal second-unknowns.
 # Pairing two firsts is represented directly, not via producer gradients.
 out={}
 for (r,e),v in a.items():
  for (s,f),w in b.items():
   if r+s<=1 and e+f<=1:out=add(out,{(r+s,e+f):v*w})
 return out

def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'SOURCE_FOLD_AND_FORMAL_AD_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k not in {'C0','G8','S'}}
 for k in ['fold_report','ad_report','gate_report','adjoint_peer','known_peer','mu_peer']:assert D[k]['all_passed']
 for k,stem in [('fold_report',A),('ad_report',O),('gate_report',G)]:
  for p,h in D[k]['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in D[k]['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/p)==h
 F=D['fold'];E=D['ad'];Gamma=list(map(sp,D['adjoint']['actual_Gamma_on_source56_rell']));assert list(map(enc,Gamma))==F['two_Gamma_on56_actual_X2']
 L=list(map(sp,D['mu_gate']['L0_Res4_Ntranspose_columns']));GL=[{j:dot(g,col) for j,col in enumerate(L) if dot(g,col)} for g in Gamma];assert list(map(enc,GL))==F['two_Gamma_L0_on48_intrinsic_native_S2']
 S=fields(ROOT/P['S'],['native_S12_by32'])['native_S12_by32'];G8=fields(ROOT/P['G8'],['G8'])['G8']
 sr=[{} for _ in range(128)]
 for row,col,k,sign in S['koszul28_in_first8_rows']:
  for i in range(4):
   for mask,c in enumerate(map(Q,G8[k][32:36])):
    if c and not mask&i:sr[4*col+i]=add(sr[4*col+i],{4*row+(mask|i):sign*c})
 s0=list(map(sp,D['mu_gate']['native_S0_columns']));se=list(map(sp,D['mu_gate']['native_S1ell_columns']));piv=D['mini']['original27_X_source_columns']
 s0p=[s0[i-640] if i>=640 else {} for i in piv];srp=[sr[i-640] if i>=640 else {} for i in piv];sep=[se[i-640] if i>=640 else {} for i in piv]
 ops=D['chain']['operator_scalar_nodes'];B=[dict(v) for v in D['chain']['B48_by640_sparse_columns']];Pcols=fields(ROOT/P['C0'],['honest_D0_relation_columns'])['honest_D0_relation_columns']
 oldr=[{i:sp(v).get(8,Q(0)) for i,v in col if sp(v).get(8,Q(0))} for col in D['old640']['source_columns640_rows56_sparse_dual9']];olde=list(map(sp,D['ell640']['actual_first640_ell_columns']))
 needed=set();atoms=set();omega=0
 for record,stored in zip(D['known']['all_zero_bottom_and_Y0_source_combinations'],F['all619_source_folds']):
  c=sp(record['original768_combination']);a1=sp(record['actual_a1_27_times9']);ar={i:a1[i*9+7] for i in range(27) if i*9+7 in a1};ae={i:a1[i*9+8] for i in range(27) if i*9+8 in a1};cn={j-640:v for j,v in c.items() if j>=640}
  assert not mv(s0,cn);omr=add(mv(sr,cn),scale(mv(s0p,ar),-1));ome=add(mv(se,cn),scale(mv(s0p,ae),-1));omega+=bool(omr)
  assert enc(ar)==stored['a1_r'] and enc(ae)==stored['a1_ell'] and enc(cn)==stored['native_intrinsic_S2_source_combination'] and enc(omr)==stored['native_Omega_r'] and enc(ome)==stored['native_Omega_ell']
  correction=add(mv(srp,ae),mv(sep,ar))
  for channel,g in enumerate(Gamma):
   formal={}
   for j,v in c.items():
    if j>=640:continue
    for row,w in g.items():
     if row>=48:key=(1,Pcols[j][row-47])
     else:
      if row not in B[j]:continue
      root=B[j][row];n=ops[root]
      if n[0]=='c':continue
      key=(1,n[1]) if n[0]=='actual_C0_scalar_function' else (0,root)
      if not key[0]:needed.add(root)
     formal=add(formal,{key:v*w})
   offset=-dot(GL[channel],correction)
   for i,j in enumerate(piv):
    if j<640:offset-=dot(g,oldr[j])*ae.get(i,Q(0))+dot(g,olde[j])*ar.get(i,Q(0))
   assert kr(formal)==stored['two_old_typedB_or_originalP_rell_root_functionals'][channel]
   assert offset==Q(stored['two_known_oldpivot_and_L0S1_corrections'][channel])
 assert omega==12
 # Exact transitive source closure independently follows source IDs.
 def visit(i):
  if i in closure:return
  n=ops[i]
  children=[] if n[0] in ['c','actual_C0_scalar_function'] else [a[0] for a in n[1]] if n[0]=='sum' else n[3]
  for j in children:visit(j)
  closure.add(i)
 closure=set()
 for i in needed:visit(i)
 assert [[i,ops[i]] for i in sorted(closure)]==F['needed_typed_closure_nodes']
 atoms={ops[i][1] for i in closure if ops[i][0]=='actual_C0_scalar_function'};assert sorted(atoms)==F['needed_original_oldB_atom_labels'] and len(closure)==629 and len(atoms)==59
 # Reconstruct firsts and original84 restoration arithmetic (not new AD).
 first={i:sp(v) for i,v in D['atomfirst']['all79_atom_values']};memo={b['original_node']:sp(b['origin_and_ds']) for b in D['mu_gate']['all_exact_cached_leaf_bindings_origin_and_ds']}
 for i,n in D['mu_gate']['old_exact84_original_arithmetic_nodes']:
  if n[0]=='c':v={0:Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=add(memo[n[1]],memo[n[2]])
  elif n[0]=='n':v=scale(memo[n[1]],-1)
  else:
   a,b=memo[n[1]],memo[n[2]];v={0:a.get(0,Q(0))*b.get(0,Q(0)),1:a.get(0,Q(0))*b.get(1,Q(0))+a.get(1,Q(0))*b.get(0,Q(0))};v={i:c for i,c in v.items() if c}
  memo[i]=v
 jets={};sens={};claimedj=dict(E['typed_node_known4_coefficients']);claimeds={i:sp(v) for i,v in E['typed_node_original_second_gradients']}
 for i in sorted(closure):
  n=ops[i]
  if n[0]=='c':v={(0,0):Q(n[1],n[2])} if n[1] else {};g={}
  elif n[0]=='actual_C0_scalar_function':
   a=n[1];v={m:c for m,c in [((0,0),first[a].get(0,Q(0))),((1,0),first[a].get(8,Q(0))),((0,1),memo[a].get(1,Q(0))/5)] if c};g={a:Q(1)}
  elif n[0]=='sum':v={};g={};
  else:v={(0,0):Q(n[1],n[2])};g={}
  if n[0]=='sum':
   for j,num,den in n[1]:v=add(v,scale(jets[j],Q(num,den)));g=add(g,scale(sens[j],Q(num,den)))
  elif n[0]=='product':
   for j in n[3]:g=add(scale(g,jets[j].get((0,0),Q(0))),scale(sens[j],v.get((0,0),Q(0))));v=jetmul(v,jets[j])
  assert [v.get(m,Q(0)) for m in [(0,0),(1,0),(0,1),(1,1)]]==list(map(Q,claimedj[i])) and g==claimeds[i]
  jets[i]=v;sens[i]=g
 unknown=set();count=0
 for f,r in zip(F['all619_source_folds'],E['all619_results']):
  for ch,func in enumerate(f['two_old_typedB_or_originalP_rell_root_functionals']):
   known=Q(f['two_known_oldpivot_and_L0S1_corrections'][ch]);g={}
   for (kind,i),c in func:
    if kind:g=add(g,{i:Q(c)})
    else:known+=Q(c)*jets[i].get((1,1),Q(0));g=add(g,scale(sens[i],Q(c)))
   assert g==sp(r['two_formal_old_functionals'][ch]['literal_original_C0_second_rell_sensitivity']) and known==Q(r['two_formal_old_functionals'][ch]['known_lower_first_response_only']);unknown.update(g);count+=len(g)
 assert count==1478 and len(unknown)==104
 ep=o/'independent619-source-fold-and-old-formal-rell-AD.json';dump(o,ep.name,{'inputs':inputs,'all619_original768_source_folds_and_both_old_native_corrections_verified':True,'all_nativeS0cN_zero':True,'native_Omega_r_nonzero_count':omega,'native_L1ell_Omega_r_retained':True,'all629_typed_nodes_and59_literal_atom_firsts_verified':True,'all1238_formal_old_second_functionals_verified':True,'unknown_original_C0_second_labels':sorted(unknown),'sensitivity_nonzeros':count,'native_expansion_exact':'L0*Srell*cN+Lr*(Sell*cN-S0p*aell)+Lell*(Sr*cN-S0p*ar)-L0*(Spr*aell+Spell*ar)','ordinary_rell_product':'a0*b_rell+b0*a_rell+a_r*b_ell+a_ell*b_r; nofactor2','no_actual_second_values_or_HF_or_N1_or_gauge_computed':True,'next_cache_bridge_not_accepted_here':'Each originalC0 scalar label must be independently bound to coherent full9 ancestor second response. Univariate rr cannot substitute for rell.'})
 report={'id':STEM,'all_passed':True,'both619_source_folds_and_all1238_formal_old_AD_accepted':True,'actual_second_or_T2_HF1_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
