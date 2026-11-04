#!/usr/bin/env python3
"""Independent fixed619 high-order source folds, no scalar values."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63_REVIEW'
A='AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63'
PATHS={
 'report':f'reports/{A}/report.json','fold':f'reports/{A}/all619-high-order-root-adjoints.json',
 'profile':f'reports/{A}/preflight.json','typed':f'reports/{A}/whole619-typed-source-and-order-profile.json',
 'actual619':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/actual619-quartic-source-order-and-unknown-gates.json',
 'bridge_report':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/report.json',
 'ownGamma':'reports/AMBIENT_B18_T_CORRECTED168_QUARTIC_VALUES_20261104_v63_REVIEW/independent-c4-and-all168-source-postverification.json',
 'ownGamma_report':'reports/AMBIENT_B18_T_CORRECTED168_QUARTIC_VALUES_20261104_v63_REVIEW/report.json',
 'S':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'gates':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(x):return {int(i):Q(c) for i,c in x if Q(c)}
def enc(x):return [[list(i) if isinstance(i,tuple) else i,str(c)] for i,c in sorted(x.items()) if c]
def put(a,k,c):
 z=a.get(k,Q(0))+c
 if z:a[k]=z
 else:a.pop(k,None)
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'independent_high_order_folds_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['report','bridge_report','ownGamma_report'])
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 G=D['ownGamma'];gamma=rec(G['Gamma0_independent'])
 for i,c in rec(G['Gamma_b3_independent']).items():put(gamma,i,Q(-42215,3456)*c)
 for i,c in rec(G['Gamma_c4_independent']).items():put(gamma,i,Q(4565,68832)*c)
 assert enc(gamma)==D['fold']['Gamma4_56_times120']
 mons2=list(itertools.combinations_with_replacement(range(8),2));i2={m:i for i,m in enumerate(mons2)}
 mons3=list(itertools.combinations_with_replacement(range(8),3))
 assert D['fold']['monomials2_36']==list(map(list,mons2)) and D['fold']['monomials3_120']==list(map(list,mons3))
 xp=list(map(rec,D['S']['first27_source_combinations']));pcols=D['C0']['honest_D0_relation_columns'];nodes=D['C0']['DAG_scalar_nodes']
 bcols=list(map(dict,D['chain']['B48_by640_sparse_columns']));ops=D['chain']['operator_scalar_nodes']
 def label(j,row,order,mon):
  if j>=640:return (2,j-640,row,order,mon) if row<48 else None
  if row>=48:
   node=pcols[j][row-48+1];return None if nodes[node][0]=='c' else (1,node,order,mon)
  node=bcols[j].get(row)
  if node is None or ops[node][0]=='c':return None
  return (1,ops[node][1],order,mon) if ops[node][0]=='actual_C0_scalar_function' else (0,node,order,mon)
 roots3=set();roots2=set();p4roots=set();xp3roots=set();native3=set();native2=set();support=0;check619=[]
 assert len(D['fold']['all619'])==619
 for source,row in zip(D['actual619']['all619'],D['fold']['all619']):
  j=source['constant619_index'];assert row['kernel619_index']==j
  c=rec(source['actual_original768_constant_combination_c0']);a1=rec(source['actual_Xp27_linear_a1_sparse27_times8'])
  assert row['actual_original768_constant_combination']==enc(c) and row['actual_Xp27_linear_a1']==enc(a1)
  p4={};xp3={};f3={};f2={}
  for k,a in c.items():
   if k<640 and nodes[pcols[k][8]][0]!='c':put(p4,pcols[k][8],a)
  for p in range(27):
   a=a1.get(8*p+7,Q(0))
   if a:
    for root,b in D['gates']['normalized_Xp27_source55_rr'][p]['rawP8_root_coefficients']:put(xp3,root,-a*Q(b))
  # Forward multiplication by each exact homogeneous source coefficient.
  # The Xp2 terms are actual coefficients, preserving every native L2*S0.
  for index,g in gamma.items():
   src,m=divmod(index,120);mon=mons3[m]
   for k,a in c.items():
    lab=label(k,src,3,m)
    if lab is not None:put(f3,lab,-g*a)
   for v in set(mon):
    rest=list(mon);rest.remove(v);m2=i2[tuple(rest)]
    for p,cp in enumerate(xp):
     a=a1.get(8*p+v,Q(0))
     if a:
      for k,b in cp.items():
       lab=label(k,src,2,m2)
       if lab is not None:put(f2,lab,g*a*b)
  assert row['rawP55_pure_r4_roots']==enc(p4) and row['Xp55_pure_r3_times_a1r_roots']==enc(xp3)
  assert row['actual_X3_high_order_source_functional']==enc(f3)
  assert row['actual_Xp2_times_a1_source_functional']==enc(f2)
  assert set(row['lower_source_operators_not_dropped'])=={'b1','l1','a2','linear_oldP3_corrections','Gamma4_Xstar2_b1','Gamma4_XL2_l1','Gamma4_Xp1_a2','Xp55_2_a2rr'}
  roots3.update(f3);roots2.update(f2);p4roots.update(p4);xp3roots.update(xp3)
  native3.update(k[1] for k in f3 if k[0]==2);native2.update(k[1] for k in f2 if k[0]==2)
  support+=len(f3)+len(f2)+len(p4)+len(xp3)+len(c)+len(a1);assert support<150000
  check619.append([j,len(f3),len(f2),len(p4),len(xp3)])
 orders3=sorted({mons3[k[-1]] for k in roots3});orders2=sorted({mons2[k[-1]] for k in roots2})
 assert orders3==[(0,7,7),(1,7,7),(4,7,7),(5,7,7),(7,7,7)]
 profile=D['profile']
 assert len(native3)==profile['native_X3_source_columns']==110
 assert len(native2)==profile['native_ACTUAL_X2_source_columns']==21
 assert len(p4roots)==profile['rawP4_root_count']==90 and len(xp3roots)==profile['Xp3_root_count']==4
 assert profile['high_order_X3_monomials']==list(map(list,orders3))
 assert profile['lower_Xp2_a1_monomials']==list(map(list,orders2))
 assert not D['typed']['lower_operator_computations_executed']
 ev={'source_manifest':inputs,'own_Gamma4':enc(gamma),'all619_source_folds_and_signs_checked':True,
  'all619_fold_support_counts':check619,'actual_X3_monomials':orders3,'actual_Xp2_a1_monomials':orders2,
  'native_ACTUAL_X3_columns':len(native3),'native_ACTUAL_X2_columns':len(native2),
  'rawP4_roots':sorted(p4roots),'Xp3_roots':sorted(xp3roots),'lower_operators_retained_and_not_evaluated':True,
  'no_effective_S2_replacement_for_individual_native_coefficients':True,
  'native_L3_drop_and_remaining_L2_gate':'S0c0=0, but S1c0+S0c1 not assumedzero; complete c1 includes the unevaluated b1/l1 corrections',
  'no_complete_coefficient_closure_budget_or_source6_sufficiency_claimed':True,
  'no_new_source_coefficients_or_actual619_quartic_values_or_HF4_computed':True,'review_QQ_support':support}
 ep=out/'independent619-high-order-source-folds.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'all619_high_order_adjoint_and_source_labels_accepted':True,
  'lower_operators_complete_source_evaluation_and_new_HF4_not_claimed':True,
  'input_hashes':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
