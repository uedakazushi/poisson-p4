#!/usr/bin/env python3
"""SOURCE8 root/M1 fold and sparse-AD budget. No second values or AD yet."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
MU='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
MINI='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json'
def count(v):
 if isinstance(v,F) or isinstance(v,int):return 1
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 return 0
def bounds(x):return x.denominator,abs(x.numerator)
def bsum(a,b):
 d=math.lcm(a[0],b[0]);return d,a[1]*(d//a[0])+b[1]*(d//b[0])
def bmul(a,b):return a[0]*b[0],a[1]*b[1]
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'SOURCE8_pre_AD_knownarrays_no_second_values'});hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 fact=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63/actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json')
 assert read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63_REVIEW/report.json')['all_passed']
 adj=read('reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63/actual-two-full918-source-adjoint-certificates-or-obstructions.json')
 lam=[]
 for r in adj['results']:
  q=H.parse(r['lambda52']);assert all(i%2 for i in q);lam.append({i//2:c for i,c in q.items()})
 M=list(map(H.parse,field(GATE,'Mnew0_columns')));G2=list(map(H.parse,field(GATE,'full_G2_linear234_by102_columns')))
 Gamma=[{j:sum((a*q.get(i,F(0)) for i,a in v.items()),F(0)) for j,v in enumerate(M) if sum((a*q.get(i,F(0)) for i,a in v.items()),F(0))} for q in lam]
 ops=dict(field(GATE,'old640_typed_closure_nodes'));B=[dict(v) for v in field(GATE,'old640_source_columns')]
 P=field(H.PATHS['C0'],'honest_D0_relation_columns')
 combs=[H.parse(r['original768']) for r in field(GATE,'all741_zero_origin_source_combinations')]
 selected=field(MU,'actual_selected_Y0_source_columns');Y0=list(map(H.parse,field(MU,'actual_selected_Y0_columns102')))
 selectedcombs=[combs[j] for j in selected];del combs
 L=list(map(H.parse,field(GATE,'L0_Res4_Ntranspose_columns')));S0=list(map(H.parse,field(GATE,'native_S0_columns')));Se=list(map(H.parse,field(GATE,'native_S1ell_columns')));effectiveell=[mv(L,v) for v in Se]
 oldr=[];X0=[]
 for path,key in [(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9'),(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9')]:
  cols=field(path,key)
  for col in cols:
   v0={};r={}
   for row,p in col:
    v=H.parse(p)
    if v.get(0):v0[row]=v[0]
    if v.get(8):r[row]=v[8]
   oldr.append(r);X0.append(v0)
  del cols
 oldell=list(map(H.parse,field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json','actual_first640_ell_columns')))
 piv=field(MINI,'original27_X_source_columns');C=field(MINI,'free30_V_indices');inv=list(map(H.parse,field(MINI,'fixed_bottom30_inverse_columns')));Xp=[X0[j] for j in piv]
 Xeff=[];a1=[]
 for n,c in enumerate(selectedcombs):
  cn={j-640:x for j,x in c.items() if j>=640};assert not mv(S0,cn) and not mv(X0,c)
  r=H.guard(mv(oldr,c));e={}
  for j,x in c.items():e=H.guard(add(e,oldell[j] if j<640 else effectiveell[j-640],x))
  ar=mv(inv,{i:r[row] for i,row in enumerate(C) if r.get(row)});ae=mv(inv,{i:e[row] for i,row in enumerate(C) if e.get(row)})
  ar={i:x for i,x in ar.items() if i<27};ae={i:x for i,x in ae.items() if i<27}
  er=H.guard(add(r,mv(Xp,ar),-1));ee=H.guard(add(e,mv(Xp,ae),-1));Xeff.append([er,ee]);a1.append([ar,ae])
  for var,v in [(7,er),(8,ee)]:
   rhs={}
   for k,x in Y0[n].items():rhs=add(rhs,{i//9:y*x for i,y in G2[k].items() if i%9==var})
   assert mv(M,v)==rhs
 # No dense X0 or source matrices are retained past this point.
 del oldr,oldell,effectiveell,X0,Xp,inv
 allweights=[]
 for r in fact['all_four_pairs_of_formal_psi_contractions']:
  for name in ['A','B']:allweights.append((r['source_extra_index'],name,H.parse(r[f'formal_{name}_weights_on99psi2'])))
 folds=[];needed=set();owned_fixed=count([ops,B,P,selectedcombs,Y0,L,S0,Se,M,G2,Xeff,a1,Gamma,lam,allweights])+6000
 assert owned_fixed<60000,('pre_fold_owned_guard',owned_fixed)
 for extra,name,w in allweights:
  source=[{},{}];vecr=[{},{}];vece=[{},{}];apr=[{},{}];ape=[{},{}];y=[{},{}]
  for k,x in w.items():
   j,ch=divmod(k,2);source[ch]=H.guard(add(source[ch],selectedcombs[j],x));vecr[ch]=H.guard(add(vecr[ch],Xeff[j][0],x));vece[ch]=H.guard(add(vece[ch],Xeff[j][1],x))
   apr[ch]=H.guard(add(apr[ch],a1[j][0],x));ape[ch]=H.guard(add(ape[ch],a1[j][1],x));y[ch]=H.guard(add(y[ch],Y0[j],x))
  old={};native=[]
  for ch,c in enumerate(source):
   cn={j-640:x for j,x in c.items() if j>=640};assert not mv(S0,cn);native.append(cn)
   for j,x in c.items():
    if j>=640:continue
    for row,a in Gamma[ch].items():
     if row<48:
      root=B[j].get(row)
      if root is None:continue
      node=ops[root]
      if node[0]=='c':continue
      key=(1,node[1]) if node[0]=='actual_C0_scalar_function' else (0,root)
     else:key=(1,P[j][row-48+1])
     v=old.get(key,F(0))+a*x
     if v:old[key]=v
     else:old.pop(key,None)
  for kind,root in old:
   if not kind:needed.add(root)
  folds.append({'extra_active_source':extra,'contraction':name,'source_original768_two_channels':[srec(c) for c in source],
   'old_typed_or_literal_C0_root_functional':[[list(k),str(v)] for k,v in sorted(old.items())],
   'native_S2rell_two_zero_origin_combinations':[srec(v) for v in native],
   'M1_known_Xeff_r_two_channels':[srec(v) for v in vecr],'M1_known_Xeff_ell_two_channels':[srec(v) for v in vece],
   'known_Xp1_correction_a_r_two_channels':[srec(v) for v in apr],'known_Xp1_correction_a_ell_two_channels':[srec(v) for v in ape],
   'pure_known_G2quadratic_Y0_two_channels':[srec(v) for v in y],
   'native_S0_zero_both_channels':True})
  assert owned_fixed+count(folds)<60000,('source_fold_live_guard',owned_fixed+count(folds))
 closure=set()
 def visit(i):
  if i in closure:return
  closure.add(i);n=ops[i]
  if n[0]=='sum':
   for j,a,b in n[1]:visit(j)
  elif n[0]=='product':
   for j in n[3]:visit(j)
 for i in needed:visit(i)
 nodes=[[i,ops[i]] for i in sorted(closure)];atoms=sorted({n[1] for i,n in nodes if n[0]=='actual_C0_scalar_function'})
 source_phase_peak=owned_fixed+count(folds)+count(nodes);assert source_phase_peak<60000,('source_before_lastuse_guard',source_phase_peak)
 # Release broad earlier arrays before source first bindings/positive profile.
 del ops,B,P,selectedcombs,Y0,L,S0,Se,Xeff,a1,M,G2
 # Actual old first values are source inputs, not second coefficients.
 old={i:H.parse(v) for i,v in field('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json','all79_atom_values')}
 leaf={v['original_node']:H.parse(v['origin_and_ds']) for v in field(GATE,'all_exact_cached_leaf_bindings_origin_and_ds')}
 for i,n in field(GATE,'old_exact84_original_arithmetic_nodes'):
  if n[0]=='c':v={0:F(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=add(leaf[n[1]],leaf[n[2]])
  elif n[0]=='n':v=scale(leaf[n[1]],-1)
  else:v=H.pairmul(leaf[n[1]],leaf[n[2]])
  leaf[i]=H.guard(v)
 atomvalues={i:[old[i].get(0,F(0)),old[i].get(8,F(0)),leaf[i].get(1,F(0))/5,F(0)] for i in atoms}
 source_phase_peak=max(source_phase_peak,count([folds,nodes,lam,Gamma,old,leaf,atomvalues])+10000);assert source_phase_peak<60000
 del old,leaf
 centers={};deps={};positive={};gbound={};highest=0
 def check(b):
  nonlocal highest
  highest=max(highest,b[0].bit_length(),b[1].bit_length());return b
 def plus4(a,b):return [check(bsum(x,y)) for x,y in zip(a,b)]
 def times4(a,b):return [check(bmul(a[0],b[0])),check(bsum(bmul(a[0],b[1]),bmul(a[1],b[0]))),check(bsum(bmul(a[0],b[2]),bmul(a[2],b[0]))),check(bsum(bsum(bmul(a[0],b[3]),bmul(a[3],b[0])),bsum(bmul(a[1],b[2]),bmul(a[2],b[1]))))]
 for i,n in nodes:
  if n[0]=='c':v=F(n[1],n[2]);ss=set();p=[bounds(v),(1,0),(1,0),(1,0)];g=(1,0)
  elif n[0]=='actual_C0_scalar_function':v=atomvalues[n[1]][0];ss={n[1]};p=list(map(bounds,atomvalues[n[1]]));g=(1,1)
  elif n[0]=='sum':
   v=F(0);ss=set();p=[(1,0)]*4;g=(1,0)
   for j,a,b in n[1]:
    c=F(a,b);v+=c*centers[j];ss.update(deps[j]);p=plus4(p,[check(bmul(bounds(c),x)) for x in positive[j]]);g=check(bsum(g,bmul(bounds(c),gbound[j])))
  else:
   v=F(n[1],n[2]);ss=set();p=[bounds(v),(1,0),(1,0),(1,0)];g=(1,0)
   for j in n[3]:
    g=check(bsum(bmul(g,bounds(centers[j])),bmul(bounds(v),gbound[j])));p=times4(p,positive[j]);v*=centers[j]
   for k,j in enumerate(n[3]):
    c=F(n[1],n[2])
    for ell,t in enumerate(n[3]):
     if ell!=k:c*=centers[t]
    if c:ss.update(deps[j])
  centers[i]=v;deps[i]=ss;positive[i]=p;gbound[i]=g
  live=count([folds,nodes,lam,Gamma,atomvalues,centers,deps,positive,gbound])+10000
  source_phase_peak=max(source_phase_peak,live);assert live<60000,('positive_profile_owned',i,live)
 rowbound=0;boundrows=[]
 for r in folds:
  labels=set();p=(1,0);g=(1,0)
  for (kind,i),cs in r['old_typed_or_literal_C0_root_functional']:
   c=F(cs)
   if kind:labels.add(i);g=check(bsum(g,bounds(c)))
   else:labels.update(deps[i]);p=check(bsum(p,bmul(bounds(c),positive[i][3])));g=check(bsum(g,bmul(bounds(c),gbound[i])))
  rowbound+=len(labels);boundrows.append({'old_known_cross_bound':list(map(str,p)),'old_gradient_L1_bound':list(map(str,g))})
 # The input is deliberately pruned for a later AD job; large earlier arrays
 # are not reimported by that job.
 slots={'folded_source_and_M1_factored_records':count(folds),'typed_nodes_and_dependency_indices':count(nodes)+count(deps),
  'typed_known4_and_gradients':4*len(nodes)+sum(map(len,deps.values())),
  'eight_formal_output_sensitivities':rowbound,'source_atom_firsts':count(atomvalues),
  'retained_lambda_Gamma_and_source_gauge_recipes':count(lam+Gamma)+4000,'scratch_metadata_and_bounds':10000}
 whole=sum(slots.values());assert whole<60000,('future_AD_whole',whole)
 assert highest<2048,('positive_AD_height_bound',highest)
 raw=whole*(2*len(str(2**highest))+24)+2_000_000;assert raw<32*1024*1024
 ep=out/'actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json';H.dump(ep,{
  'eight_source_folds':folds,'lambda26_two_channels':[srec(v) for v in lam],'Gamma56_two_channels':[srec(v) for v in Gamma],
  'needed_typed_nodes':nodes,'actual_atom_origin_r_ell_firsts':[[i,list(map(str,v))] for i,v in sorted(atomvalues.items())],
  'typed_constants':[[i,str(v)] for i,v in sorted(centers.items())],'typed_structural_second_dependency_labels':[[i,sorted(v)] for i,v in sorted(deps.items())],
  'positive_output_bound_rows':boundrows,'positive_AD_height_bound_bits':highest,'next_AD_whole_QQ_bound':whole,'next_AD_slot_breakdown':slots,'next_AD_raw_byte_bound':raw,
  'source_phase_owned_upper':source_phase_peak,
  'AD_started':False,'actual_second_values_queried':False,
  'M1_source_formula':'sum_psi lambda_psi*(M1_r*Xeff_ell +M1_ell*Xeff_r), withthetwofactoredsourcevectorsforeach8',
  'known_pivot_formula':'minus sum_psi Gamma_psi*(Xp1_r*a_ell +Xp1_ell*a_r); typedfirst-onlyresponse, notevaluatedhere',
  'native_formula':'sum_psi Gamma_psi[L0*S2rell*cN+Lr*Sell*cN+Lell*Sr*cN]; L2dropsONLYexactS0cN0',
  'known_G2_formula':'minus sum_psi lambda_psi*G2quadratic_rell*weightedY0_psi; noqueryhere',
  'gauge_formula':'M1_normalized=M1_canonical+G2linear*A0; keepcanonicalmu2/typednormalizingA0. No claimA0=0 orremovablebeyondprovedgauge.',
  'formal_AD_scope':'opaqueC0secondrell symbolsindependent; noactualphysical/copy/honestsecondvalues; new8sensitivitymustberecomputed',
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'QQ_slots':60000,'raw_bytes':32*1024*1024,'bits':2048}})
 report={'id':STEM,'all_passed':True,'AD_started':False,'actual_second_values_queried':False,'source_contractions':8,
  'typed_closure_nodes':len(nodes),'atom_count':len(atoms),'next_AD_whole_QQ_bound':whole,'positive_AD_height_bound_bits':highest,'next_AD_raw_byte_bound':raw,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','typed_closure_nodes','atom_count','next_AD_whole_QQ_bound','positive_AD_height_bound_bits','next_AD_raw_byte_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
