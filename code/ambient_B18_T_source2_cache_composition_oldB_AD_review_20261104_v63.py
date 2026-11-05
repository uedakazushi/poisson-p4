#!/usr/bin/env python3
"""Independent nonlinear FULL9 transport and oldB honest-second cancellation."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time
from ambient_B18_T_gdiff_second_univariate5_review_20261104_v63 import Symbolic,closure,rec,enc,add,put
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW'
A='AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD'
PATHS={
 'report':f'reports/{A}/report.json','oldvalues':f'reports/{A}/phase3-oldB-common-source2-values-and-symbols.json',
 'parts':f'reports/{A}/phase4-all619-lower-knownparts-and-symbols.json',
 'source_cache':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'inverse':'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'copy2':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
 'composite':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
 'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'lower':'reports/AMBIENT_B18_T_COMMON_LOWER_ADJOINT_20261104_v63/phase3-all619-common-D-contractions.json',
 'D_peer':'reports/AMBIENT_B18_T_COMMON_LOWER_ADJOINT_20261104_v63_REVIEW/report.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,name,v):(out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'independent_source2_transport_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert D['report']['all_passed'] and D['D_peer']['all_passed']
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 active=[0,1,4,5,7];mon=[()]+[(i,) for i in active]+[m for m in itertools.combinations_with_replacement(active,2) if m!=(1,1)]
 def ex(m):return tuple(m.count(i) for i in range(8))
 ring=Symbolic([ex(m) for m in mon]);idx={m:i for i,m in enumerate(mon)};quad=[i for i,m in enumerate(mon) if len(m)==2]
 assert D['source_cache']['base20_monomials_original8']==list(map(list,mon))
 ownmons=[tuple(m) for m in D['FJ']['source2_monomials']]
 assert ownmons==[tuple(m) for m in D['inverse']['source2_monomials']]
 variables=[({idx[(i,)]:Q(1)}, {}) if i in active else ring.c(0) for i in range(8)];YA=variables[:7]+[variables[7],variables[7]]
 def comp(p,args):
  z=ring.c(0)
  for m,c in zip(ownmons,p):
   if not Q(c):continue
   t=ring.c(Q(c))
   for i in m:t=ring.mul(t,args[i])
   z=ring.add(z,t)
  return z
 qa,qb=[D['inverse']['families'][f] for f in ['CA','CB']]
 XA=[comp(p,YA) for p in qa['physical_R0_readback_jets_in_own_F9']]
 # Independent nonlinear inverse formula from XB1 inverse and its raw
 # quadratic part; do not reuse the saved inverseCB polynomial table.
 Bi=[[Q(c) for c in row] for row in qb['fixed_inverse9']]
 Bl=[[Q(c) for c in row] for row in qb['R0_by_F9_unit_Jacobian']]
 assert all(sum(Bi[i][k]*Bl[k][j] for k in range(9))==Q(i==j) for i in range(9) for j in range(9))
 def linear(p):return ({i:c for i,c in p[0].items() if len(mon[i])==1},{})
 Ylin=[]
 for row in Bi:
  z=ring.c(0)
  for a,p in zip(row,XA):z=ring.add(z,ring.scale(linear(p),a))
  Ylin.append(z)
 Xq=[]
 for p in qb['physical_R0_readback_jets_in_own_F9']:
  pp=[Q(c) if len(m)==2 else Q(0) for m,c in zip(ownmons,p)];Xq.append(comp(pp,Ylin))
 YB=[]
 for i,row in enumerate(Bi):
  z=Ylin[i]
  for a,p,q in zip(row,XA,Xq):z=ring.add(z,ring.scale(ring.add(ring.add(p,linear(p),-1),q,-1),a))
  YB.append(z)
 assert [comp(p,YB) for p in qb['physical_R0_readback_jets_in_own_F9']]==XA
 assert [comp(p,XA) for p in qa['inverse_F9_in_common_R0_through2']]==YA
 ca={i:comp(p,YA) for i,p in D['FJ']['both57_source2']['CA']['all57_source2jets']}
 cb={i:comp(p,YB) for i,p in D['FJ']['both57_source2']['CB']['all57_source2jets']}
 assert len(set(ca)&set(cb))==9 and all(ca[i]==cb[i] for i in set(ca)&set(cb))
 common=ca|cb;assert len(common)==105
 # Copied regular-function section is deliberately the accepted CB zero7
 # section, distinct from the full common-CB48 physical base-change above.
 copyargs=[ring.c(0)]*7+YB[7:9];labels=dict(D['copy']['explicit57_source_copy_bindings'])
 copied={labels[i]:comp(p,copyargs) for i,p in D['FJ']['both57_source2']['CB']['all57_source2jets']}
 cache=D['source_cache'];assert [[i,enc(p[0])] for i,p in sorted(common.items())]==cache['all105_actual_source2_values']
 assert [[i,enc(p[0])] for i,p in sorted(copied.items())]==cache['all57_actual_copy_source2_values']
 assert list(map(lambda p:enc(p[0]),YB))==cache['actual_YB9_nonlinear_through2']
 bmons=[tuple(i for i,n in enumerate(m) for _ in range(n)) for m in D['B11']['B11_basis_monomials8']];bi={m:i for i,m in enumerate(bmons)}
 def project(p):return {bi[mon[i]]:c for i,c in p[0].items() if mon[i] in bi}
 assert {i:project(p) for i,p in common.items()}=={i:rec(p) for i,p in D['B11']['all105_actual_B11_root_values']}
 assert {i:project(p) for i,p in copied.items()}=={i:rec(p) for i,p in D['copy2']['all57_copied_unknown_values_through2']}
 nodes=D['C0']['DAG_scalar_nodes'];assert nodes[:len(D['source']['DAG_scalar_nodes'])]==D['source']['DAG_scalar_nodes']
 def evaluate(roots,bind):
  order=closure(nodes,roots,set(bind));uses=collections.Counter()
  for i in order:
   if i in bind:continue
   n=nodes[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
  uses.update(roots);values={};peak=0
  for i in order:
   n=nodes[i];o=n[0]
   if i in bind:p=bind[i]
   elif o=='c':p=ring.c(Q(n[1],n[2]))
   elif o=='a':p=ring.add(values[n[1]],values[n[2]])
   elif o=='n':p=ring.scale(values[n[1]],-1)
   elif o=='m':p=ring.mul(values[n[1]],values[n[2]])
   elif o=='i':p=ring.inv(values[n[1]])
   else:raise AssertionError((i,n))
   values[i]=p;peak=max(peak,sum(len(a)+len(b) for a,b in values.values()));assert peak<150000
   if i not in bind:
    for j in n[1:3] if o in ['a','m'] else n[1:] if o in ['n','i'] else []:
     uses[j]-=1
     if not uses[j]:values.pop(j,None)
  return {i:values[i] for i in roots},peak
 eq105=D['joint']['all105_equation_original_node_order'];eqcopy=D['copy']['all57_copied_equation_indices']
 eqvalues,peak=evaluate(eq105+eqcopy,common|copied)
 assert all(eqvalues[i]==(YA[j-48] if 48<=j<57 else ring.c(0)) for j,i in enumerate(eq105))
 assert all(eqvalues[i]==ring.c(0) for i in eqcopy)
 dump(out,'phase1-independent-nonlinear-source2-transport.json',{'all105_and57_equations_checked':True,'all105_and57_B11_projections_exact':True,
  'actual_YB9_nonlinear':list(map(lambda p:enc(p[0]),YB)),'all105_source_values_match':True,'copy_zero7_section_not_common_CB48':True})
 # Second coefficients of genuine source quotient functions remain
 # independent for each quadratic monomial, rather than one untagged sum.
 opaque={}
 def from_first(p):
  p=list(map(Q,p));v={0:p[0]} if p[0] else {}
  for k in active:
   c=p[8]+p[9] if k==7 else p[1+k]
   if c:v[idx[(k,)]]=c
  return v
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in D['FJ'][key].items():opaque[int(i)]=from_first(p)
 for row in D['atomfirst']['all35_ordered_CA_prefix_firsts']:
  q=rec(row['first9_bar']);opaque[row['node']]={0:q[0]} if q.get(0) else {}
  for k in active:
   if q.get(1+k):opaque[row['node']][idx[(k,)]]=q[1+k]
 for i,p in D['atomfirst']['all8_g_quotient_overrides']:
  q=rec(p);opaque[i]={0:q[0]} if q.get(0) else {}
  for k in active:
   if q.get(1+k):opaque[i][idx[(k,)]]=q[1+k]
 bind=common|copied
 for i,v in opaque.items():bind[i]=(v,{(i,mon[q],q):Q(1) for q in quad})
 for j,i in enumerate(D['C0']['actual_C0_basis8_scalar_nodes']):bind[i]=variables[j] if j<7 else ring.c(0)
 oldlabels={tuple(k) for row in D['lower']['all619_source_functionals'] for k,c in row['X2_source_functional'] if k[0]<2}
 troot={k[1] for k in oldlabels if k[0]==0};croot={k[1] for k in oldlabels if k[0]==1};ops=D['chain']['operator_scalar_nodes'];tc=set();todo=list(troot)
 while todo:
  i=todo.pop()
  if i in tc:continue
  tc.add(i);n=ops[i];o=n[0]
  if o=='actual_C0_scalar_function':croot.add(n[1])
  elif o=='sum':todo.extend(j for j,a,b in n[1])
  elif o=='product':todo.extend(n[3])
  elif o in ['a','m']:todo.extend(n[1:3])
  elif o in ['n','i']:todo.append(n[1])
  else:assert o=='c'
 oldcv={i:(rec(v),rec(g)) for i,v,g in D['composite']['all108_composite_values_and_opaque_second_gradients']};croot.update(oldcv)
 cv,p2=evaluate(sorted(croot),bind);tv={}
 for i in sorted(tc):
  n=ops[i];o=n[0]
  if o=='c':p=ring.c(Q(n[1],n[2]))
  elif o=='actual_C0_scalar_function':p=cv[n[1]]
  elif o=='sum':
   p=ring.c(0)
   for j,a,b in n[1]:p=ring.add(p,tv[j],Q(a,b))
  elif o=='product':
   p=ring.c(Q(n[1],n[2]))
   for j in n[3]:p=ring.mul(p,tv[j])
  elif o=='a':p=ring.add(tv[n[1]],tv[n[2]])
  elif o=='n':p=ring.scale(tv[n[1]],-1)
  elif o=='m':p=ring.mul(tv[n[1]],tv[n[2]])
  elif o=='i':p=ring.inv(tv[n[1]])
  else:raise AssertionError((i,n))
  tv[i]=p
 def check_saved(i,p,saved):
  assert enc(p[0])==saved[0]
  g=rec(saved[1])
  for q in quad:assert {(node):c for (node,tag,out),c in p[1].items() if out==q}==g
 for i,v,g in D['oldvalues']['C0_composite_values_and_opaque_second_gradients']:check_saved(i,cv[i],[v,g])
 for i,v,g in D['oldvalues']['typedB_root_values_and_opaque_second_gradients']:check_saved(i,tv[i],[v,g])
 for i,(v,g) in oldcv.items():assert project(cv[i])==v
 pairs=list(itertools.combinations_with_replacement(range(8),2));known=[];sens=[];native=[]
 for row in D['lower']['all619_source_functionals']:
  z=Q(row['known_Xp1_a1']);g={};nn={}
  for k,c in row['X2_source_functional']:
   k=tuple(k);c=Q(c)
   if k[0]==2:put(nn,k,c);continue
   p=tv[k[1]] if k[0]==0 else cv[k[1]];q=idx[pairs[k[-1]]];z+=c*p[0].get(q,Q(0))
   for (node,tag,out),a in p[1].items():
    if out==q:put(g,(node,k[-1]),c*a)
  assert not g;known.append(z);sens.append(g);native.append(nn)
 assert list(map(str,known))==D['parts']['all619_oldB_known_parts']
 assert list(map(enc,native))==D['parts']['all619_uncomputed_native_second_source_functionals']
 assert sum(bool(v) for v in native)==70 and sum(map(len,native))==28302
 ev={'source_manifest':inputs,'independent_method':'nonlinear inverse from XB1/XB2, sparse monomial symbolic algebra with separate honest-second symbol for every leaf and quadratic pair',
  'all105_shared_source_and57_zero7_copy_equations_and_B11_checked':True,
  'all_oldB_numeric_values_and_gradients_match':True,'all619_oldB_known_parts':list(map(str,known)),
  'all619_honest_second_sensitivities_exact_zero':True,'remaining_native_sources':70,'remaining_native_QQ_coefficients':28302,
  'B20_meaning':'m³ truncation in z0,z1,z4,z5,r with z1²=0; distinct from the conductor H20 algebra',
  'no_native_extra_second_gauge_values_or_rawP4_or_X3_or_newsource_solve':True,
  'not_complete619_lower_or_quartic_values':True,'QQ_source_peak':max(peak,p2)}
 ep=out/'independent-source2-transport-and619-oldB-sensitivity.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'nonlinear_FULL9_transport_and105_57_B11_checks_accepted':True,
  'all619_oldB_honest_second_sensitivities_zero_accepted':True,'native70_sources_still_symbolic':True,
  'complete619_quartics_or_HF4_not_claimed':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep),'phase1-independent-nonlinear-source2-transport.json':sha(out/'phase1-independent-nonlinear-source2-transport.json')},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
