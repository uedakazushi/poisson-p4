#!/usr/bin/env python3
"""Independent source-Jacobian0/r and coherent-copy high-adjoint replay."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,json,resource,signal,time
from ambient_B18_T_source4_adjoint9_review_20261104_v63 import PATHS as P9
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63';STEM=A+'_REVIEW'
PATHS={k:P9[k] for k in ['C0','cache2','cache3','FJ','copy','atom','high']}
PATHS.update({'primary':f'reports/{A}/report.json','jacC':f'reports/{A}/phase2-C0-source-Jacobians0-r.json',
 'jacB':f'reports/{A}/phase3-typed-high-old-source-Jacobians0-r.json','folds':f'reports/{A}/phase4-all619-high-old-source-sensitivities.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'old':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
 'source4_peer':'reports/AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63_REVIEW/report.json'})
ZERO=(Q(0),Q(0));ONE=(Q(1),Q(0))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def scale(a,c):return (c*a[0],c*a[1])
def times(a,b):return (a[0]*b[0],a[0]*b[1]+a[1]*b[0])
def gad(a,b,c=ONE):
 d=dict(a)
 for k,v in b.items():
  x=add(d.get(k,ZERO),times(c,v))
  if x!=ZERO:d[k]=x
  else:d.pop(k,None)
 return d
def product(a,b):return times(a[0],b[0]),gad(gad({},a[1],b[0]),b[1],a[0])
def summation(a,b,c=Q(1)):return add(a[0],scale(b[0],c)),gad(a[1],b[1],(c,Q(0)))
def encode(g):return [[list(k),list(map(str,p))] for k,p in sorted(g.items())]
def push(out,k,v):
 a=out.get(k,Q(0))+v
 if a:out[k]=a
 else:out.pop(k,None)
def dump(out,name,x):
 (out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()};dump(out,'phase.json',{'phase':'SOURCE_ONLY_AD_INPUTS_SAVED','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()};assert D['primary']['all_passed'] and D['source4_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 C=D['C0'];nodes=C['DAG_scalar_nodes'];FJ=D['FJ'];HI=D['high'];ops=D['chain']['operator_scalar_nodes']
 mons=[tuple(m) for m in D['cache2']['base20_monomials_original8']];ri=mons.index((7,));bound={}
 def val2(p):
  v={int(i):Q(c) for i,c in p};return (v.get(0,Q(0)),v.get(ri,Q(0)))
 for i,p in D['cache2']['all105_actual_source2_values']:bound[i]=(val2(p),{('physical',i):ONE})
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in FJ[key].items():bound[int(i)]=((Q(p[0]),Q(p[8])+Q(p[9])),{('honest',int(i)):ONE})
 for row in D['atom']['all35_ordered_CA_prefix_firsts']:
  v={int(i):Q(c) for i,c in row['first9_bar']};i=row['node'];bound[i]=((v.get(0,Q(0)),v.get(8,Q(0))),{('honest',i):ONE})
 for i,p in D['atom']['all8_g_quotient_overrides']:
  v={int(j):Q(c) for j,c in p};bound[i]=((v.get(0,Q(0)),v.get(8,Q(0))),{('honest',i):ONE})
 for i in C['actual_C0_basis8_scalar_nodes']:
  if i not in bound:bound[i]=(ZERO,{})
 cb=FJ['both57_source2']['CB'];ys=cb['F9_nodes'][7:9];mi={tuple(m):i for i,m in enumerate(FJ['source2_monomials'])}
 own={i:list(map(Q,p)) for i,p in cb['all57_source2jets']};copy={j:i for i,j in D['copy']['explicit57_source_copy_bindings']}
 cv={i:val2(p) for i,p in D['cache2']['all57_actual_copy_source2_values']}
 croots=set();troots=set()
 for row in HI['all619']:
  croots.update(i for i,c in row['rawP55_pure_r4_roots']+row['Xp55_pure_r3_times_a1r_roots'])
  for key,c in row['actual_X3_high_order_source_functional']:
   if key[0]==0:troots.add(key[1])
   elif key[0]==1:croots.add(key[1])
 tc=set();todo=list(troots)
 while todo:
  i=todo.pop()
  if i in tc:continue
  tc.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  elif n[0]=='actual_C0_scalar_function':croots.add(n[1])
  else:assert n[0]=='c'
 def children(i):
  if i in bound:return []
  if i in copy:return ys
  n=nodes[i];assert n[0] in ['c','a','m','n','i'],('UNBOUND_SOURCE',i,n)
  return n[1:] if n[0]!='c' else []
 needed=set();todo=list(croots)
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);todo.extend(children(i))
 uses=collections.Counter(j for i in needed for j in children(i));uses.update(croots)
 values={};peak=0;live=0;copychecks=[]
 for i in sorted(needed):
  n=nodes[i];o=n[0]
  if i in bound:v=bound[i]
  elif i in copy:
   p=own[copy[i]];gradient={};hessian=[]
   for a,y in enumerate(ys):
    h0=p[mi[(7+a,)]];h1=Q(0)
    for b,x in enumerate(ys):
     pair=tuple(sorted((7+a,7+b)));h=p[mi[pair]]*(2 if a==b else 1)
     h1+=h*values[x][0][1]
    gradient=gad(gradient,values[y][1],(h0,h1));hessian.append([str(h0),str(h1)])
   assert cv[i][0]==p[0] and cv[i][1]==sum(p[mi[(7+a,)]]*values[y][0][1] for a,y in enumerate(ys))
   v=(cv[i],gradient);copychecks.append({'copy_node':i,'own_source_node':copy[i],'gradient_YB0_r':hessian})
  elif o=='c':v=((Q(n[1],n[2]),Q(0)),{})
  elif o=='a':v=summation(values[n[1]],values[n[2]])
  elif o=='n':v=(scale(values[n[1]][0],Q(-1)),{k:scale(a,Q(-1)) for k,a in values[n[1]][1].items()})
  elif o=='m':v=product(values[n[1]],values[n[2]])
  elif o=='i':
   x,g=values[n[1]];assert x[0];iv=(1/x[0],-x[1]/x[0]**2);assert times(iv,x)==ONE
   v=(iv,{k:scale(times(times(iv,iv),a),Q(-1)) for k,a in g.items()})
  else:raise AssertionError(n)
  values[i]=v;live+=2+2*len(v[1]);peak=max(peak,live);assert peak<125000
  for j in children(i):
   uses[j]-=1
   if not uses[j]:live-=2+2*len(values[j][1]);del values[j]
 reconstructed=[[i,list(map(str,values[i][0])),encode(values[i][1])] for i in sorted(croots)]
 assert reconstructed==D['jacC']['all_C0_roots']
 dump(out,'phase.json',{'phase':'COHERENT_COPY_PSI2_AND_ALL_C0_JACOBIANS_PASS','source_roots':len(croots),'peak_QQ':peak})
 tu=collections.Counter()
 for i in tc:
  n=ops[i]
  if n[0]=='sum':tu.update(j for j,a,b in n[1])
  elif n[0]=='product':tu.update(n[3])
 tu.update(troots);tv={};tlive=0;tpeak=0
 for i in sorted(tc):
  n=ops[i];o=n[0]
  if o=='c':v=((Q(n[1],n[2]),Q(0)),{})
  elif o=='actual_C0_scalar_function':v=values[n[1]]
  elif o=='sum':
   v=(ZERO,{})
   for j,a,b in n[1]:v=summation(v,tv[j],Q(a,b))
  elif o=='product':
   v=((Q(n[1],n[2]),Q(0)),{})
   for j in n[3]:v=product(v,tv[j])
  else:raise AssertionError(n)
  tv[i]=v;tlive+=2+2*len(v[1]);tpeak=max(tpeak,tlive);assert peak+tpeak<140000
  ch=[j for j,a,b in n[1]] if o=='sum' else n[3] if o=='product' else []
  for j in ch:
   tu[j]-=1
   if not tu[j]:tlive-=2+2*len(tv[j][1]);del tv[j]
 assert [[i,list(map(str,tv[i][0])),encode(tv[i][1])] for i in sorted(troots)]==D['jacB']['typed_roots']
 old={i:(Q(v),{int(k):Q(c) for k,c in p}) for i,v,gg,p,q in D['old']['node_data']}
 assert all(tv[i][0]==(old[i][0],old[i][1].get(7,Q(0))) for i in troots if i in old)
 all3=set();all4=set();folds=[]
 for row,claimed in zip(HI['all619'],D['folds']['all619']):
  t3={};t4={}
  for root,c in row['rawP55_pure_r4_roots']:
   for key,p in values[root][1].items():push(t4,key,Q(c)*p[0]);push(t3,key+(119,),Q(c)*p[1])
  for root,c in row['Xp55_pure_r3_times_a1r_roots']:
   for key,p in values[root][1].items():push(t3,key+(119,),Q(c)*p[0])
  for key,c in row['actual_X3_high_order_source_functional']:
   if key[0]==2:continue
   gradient=(tv if key[0]==0 else values)[key[1]][1]
   for tag,p in gradient.items():push(t3,tag+(key[-1],),Q(c)*p[0])
  x3=[[list(k),str(v)] for k,v in sorted(t3.items())];x4=[[list(k),str(v)] for k,v in sorted(t4.items())]
  assert x3==claimed['old_high_source3_sensitivity'] and x4==claimed['old_high_source4_sensitivity']
  assert row['kernel619_index']==claimed['kernel619_index'];all3.update(t3);all4.update(t4)
  folds.append({'kernel619_index':row['kernel619_index'],'source3':x3,'source4':x4})
 assert len(all3)==71 and len(all4)==9 and all(k[0]=='physical' for k in all4)
 honest=sorted({k[1] for k in all3 if k[0]=='honest'});assert honest==[330164,330289,330524,330594,566400,566401,566402,566403]
 ep=out/'independent-source-Jacobian-and-all619-high-folds.json'
 dump(out,ep.name,{'source_manifest':inputs,'all_C0_and_typed_Jacobians0_r_independently_match':True,
  'coherent_copy_psi1_psi2_records':copychecks,'all619_old_high_folds':folds,
  'active_source3_pairs':len(all3),'physical_source4_leaves':len(all4),'only_honest_higher_nodes':honest,
  'source_AD_closure':len(needed),'typed_AD_closure':len(tc),'peak_live_QQ':peak+tpeak,
  'third_callback_monomials':'index in combinations_with_replacement(range(8),3), ordinary coefficients',
  'native_G3_not_combined_or_evaluated':True,'new_source_or_gauge_coefficients':0,
  'no_necessity_HF4_or_actual619_quartic_claim':True})
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'all619_high_old_source_adjoints_and_coherent_copy_rule_independently_accepted':True,
  'source3_active_pairs':71,'source4_physical_targets':9,'new_coefficients_computed':0,
  'native_third_combined_or_actual619_quartics':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-st,
  'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
