#!/usr/bin/env python3
"""Independent118 source reconstruction, symbolic ordinary P3, forward B11 AD.

No author worker/import, source solve, new coefficient, point or rank search.
"""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time
from ambient_B18_T_gdiff_second_univariate5_review_20261104_v63 import (
 Symbolic,closure,numeric_eval,series_mul,put,add,rec,enc)

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63_REVIEW'
A='AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63'
PATHS={
 'report':f'reports/{A}/report.json','author_phase1':f'reports/{A}/phase1-Gamma-and-b3-source.json',
 'author_phase2':f'reports/{A}/phase2-all118-source-functionals.json',
 'author_values':f'reports/{A}/actual-b3-and-all118-pure-quartic-values.json',
 'theorem':'reports/AMBIENT_B18_T_RESTRICTED118_QUARTIC_PREFLIGHT_20261104_v63/restricted118-quartic-source-operator-preflight.json',
 'theorem_report':'reports/AMBIENT_B18_T_RESTRICTED118_QUARTIC_PREFLIGHT_20261104_v63/report.json',
 'fifth_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
 'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'cache':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
 'cache_peer':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'copy1':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'copy2':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
 'fourth':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
 'native_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'w0':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
 'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,name,v):(out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'independent_source_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['report','theorem_report','fifth_peer','cache_peer','native_peer'])
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 T,M=D['theorem'],D['m1'];ops=D['chain']['operator_scalar_nodes'];nodes=D['C0']['DAG_scalar_nodes']
 bcols=[dict(c) for c in D['chain']['B48_by640_sparse_columns']];pcols=D['C0']['honest_D0_relation_columns']
 m2=list(itertools.combinations_with_replacement(range(8),2));i2={m:i for i,m in enumerate(m2)}
 m3=list(itertools.combinations_with_replacement(range(8),3));i3={m:i for i,m in enumerate(m3)}
 qnorm=list(map(rec,D['low']['normalized84_unitcycle_Q1_columns_sparse']))
 wb=list(map(rec,D['w0']['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):j for j,v in enumerate(wb)}
 ys=[Q(r['known_rr']) for r in T['known_Xp_and_Ystar_source55_second']['normalized_Ystar84_source55_rr']]
 J=[];H=[]
 for w in range(84):
  rem={w:Q(1)};coeff={}
  while rem:
   k=min(rem);p=wl[k];a=rem[k];put(coeff,p,a);rem=add(rem,wb[p],-a)
  y={};h=Q(0)
  for p,a in coeff.items():y=add(y,qnorm[p],a);h+=a*ys[p]
  J.append(y);H.append(h)
 m0=[list(map(Q,r)) for r in D['V']['all24_by56_m0']]
 U=M['free_kernel_rows32'];pi=[{U[j]:c for j,c in rec(r).items()} for r in M['pi_free32_rows']]
 # Forward derivative/Euler scattering on each virtual input coefficient.
 # This independently reconstructs Gamma, rather than reading the author's
 # contraction rows or differentiating the author's root formula.
 def gamma(degree,theta,include_H):
  mon=m2 if degree==2 else m3;idx=i2 if degree==2 else i3;n=len(mon);virt={};outadj={}
  for w,(sector,pair) in enumerate(D['w0']['literal_wedge_basis84']):
   i,j=pair;terms=[(8*sector+7,i,1,1)] if j==7 else [(8*sector+j,i,1,2),(8*sector+i,j,-1,2)]
   for row,var,sign,shift in terms:
    for q,mm in enumerate(mon):
     if var not in mm:continue
     small=list(mm);small.remove(var)
     weight=H[w] if include_H and small==[7]*(degree-1) else Q(0)
     for tv,c in J[w].items():
      target,v=divmod(tv,8);weight-=c*theta.get(n*target+idx[tuple(sorted(small+[v]))],Q(0))
     if weight:put(virt,n*row+q,Q(sign*mm.count(var),sum(v<7 for v in small)+shift)*weight)
  for index,c in theta.items():
   target,q=divmod(index,n)
   for src,a in pi[target].items():put(outadj,n*src+q,c*a)
  for index,c in virt.items():
   row,q=divmod(index,n)
   for src,a in enumerate(m0[row]):
    if a:put(outadj,n*src+q,c*a)
  return outadj
 theta0=rec(T['Theta3_on_cubic600']['constant_part']);thetaB=rec(T['Theta3_on_cubic600']['canonical_qB_pure_cubic_coefficient_part'])
 g0=gamma(3,theta0,True);gB=gamma(3,thetaB,False)
 assert enc(g0)==D['author_phase1']['Gamma0_56_times120'] and enc(gB)==D['author_phase1']['Gammab_56_times120']
 theta2={71:Q(-55495,20736),107:Q(27785,1536)};g2adj=gamma(2,theta2,True)
 assert g0=={120*(i//36)+i3[m2[i%36]+(7,)]:c for i,c in g2adj.items()}
 assert gB=={120*29+119:Q(1,4),120*30+119:Q(1)}
 xp1=list(map(rec,M['actual_pivot_X1_by27_sparse56_times8']))
 xproots=[rec(r['rawP8_root_coefficients']) for r in T['known_Xp_and_Ystar_source55_second']['normalized_Xp27_source55_rr']]
 src37=list(map(rec,T['actual_selected37_original768_constant_combinations']))
 def label(j,row,mon):
  if j>=640:return (2,j-640,row,mon) if row<48 else None
  if row>=48:
   root=pcols[j][row-48+1];return None if nodes[root][0]=='c' else (1,root,mon)
  root=bcols[j].get(row)
  if root is None or ops[root][0]=='c':return None
  return (1,ops[root][1],mon) if ops[root][0]=='actual_C0_scalar_function' else (0,root,mon)
 def P3(c):
  z={}
  for j,a in c.items():
   if j<640 and nodes[pcols[j][8]][0]!='c':put(z,pcols[j][8],a)
  return z
 def XP(a,stride):
  z={}
  for p,roots in enumerate(xproots):
   ar=a.get(stride*p+stride-1,Q(0))
   for root,b in roots.items():
    if nodes[root][0]!='c':put(z,root,-ar*b)
  return z
 def fold(c1,a2,adj):
  f={};known=Q(0)
  for index,c in adj.items():
   row,q=divmod(index,120);mon=m3[q]
   for v in set(mon):
    pair=list(mon);pair.remove(v);q2=i2[tuple(pair)]
    for j,a in c1[v].items():
     lab=label(j,row,q2)
     if lab is not None:put(f,lab,-c*a)
    for p,first in enumerate(xp1):known+=c*first.get(8*row+v,Q(0))*a2.get(36*p+q2,Q(0))
  return f,known
 # Reconstruct b3 in the exact accepted quadratic-source normalization.
 c=rec(M['620_actual_source_combinations_in_original768'][563]);a1=rec(M['620_actual_linear_pivot_coefficients_a1_by27_times8'][563])
 bf={};bk=Q(0)
 for index,cg in g2adj.items():
  row,q=divmod(index,36);u,v=m2[q]
  for p,first in enumerate(xp1):
   bk+=cg*a1.get(8*p+u,Q(0))*first.get(8*row+v,Q(0))
   if u!=v:bk+=cg*a1.get(8*p+v,Q(0))*first.get(8*row+u,Q(0))
  for j,a in c.items():
   lab=label(j,row,q)
   if lab is not None:put(bf,lab,-cg*a)
 bformula=(P3(c),XP(a1,8),bf,bk)
 assert {'selected620':563,'divide':'-250/3','rawP3':enc(bformula[0]),'XP2':enc(bformula[1]),'X2':enc(bf),'known':str(bk)}==D['author_phase1']['b3_source']
 formulas=[];usedlabels=set(bf);rawroots=set(bformula[0])|set(bformula[1])
 for ki,combo in enumerate(T['actual118_constant_combinations_of119']):
  coeffs=[{} for _ in range(8)];a2={}
  for syz,scale in rec(combo).items():
   for index,cg in rec(T['all119_linear_source_syzygies296'][syz]['coefficients296_sparse']).items():
    j,v=divmod(index,8);coeffs[v]=add(coeffs[v],src37[j],scale*cg)
   a2=add(a2,rec(T['all119_actual_a2_sparse27_times36'][syz]),scale)
  f0,k0=fold(coeffs,a2,g0);fb,kb=fold(coeffs,a2,gB);p3=P3(coeffs[7]);xp=XP(a2,36)
  r=D['author_phase2']['all118'][ki]
  assert r=={'kernel118_index':ki,'combination_in119':combo,'rawP3':enc(p3),'XP2':enc(xp),'X2_constant':enc(f0),'X2_times_b3':enc(fb),'known_constant':str(k0),'known_times_b3':str(kb)}
  formulas.append((p3,xp,f0,fb,k0,kb));usedlabels.update(f0);usedlabels.update(fb);rawroots.update(p3);rawroots.update(xp)
 dump(out,'phase1-independent-Gamma-b3-and118-source.json',{'Gamma0':enc(g0),'Gamma_b3':enc(gB),'Gamma2':enc(g2adj),'all118_source_combinations_and_corrections_matched':True,'b3_source_matched':True})
 dump(out,'phase.json',{'phase':'all118_source_operators_independently_reconstructed','seconds':time.monotonic()-start})
 # Independent ordinary-source substitution, with ALL absent regular quotient
 # second/third coefficients symbolic until the final functionals contract.
 source3=D['joint'];common={}
 for i,p in source3['all105_unknown_jets']:
  q=[Q(0)]*4
  for m,cg in zip(source3['monomials'],p):q[sum(m)]+=Q(cg)
  common[i]=tuple(q)
 yroots=[i for i,p in D['copy2']['actual_YB_last2_from_saved105']];eq105=source3['all105_equation_original_node_order']
 oldv,units=numeric_eval(nodes,eq105+yroots,common,4)
 assert all(oldv[i]==(0,1 if j in [55,56] else 0,0,0) for j,i in enumerate(eq105))
 copybind=dict(D['copy1']['explicit57_source_copy_bindings']);copyvals={}
 for i,p in D['fourth']['families']['CB']['all57_source4_jets']:
  q=(Q(0),)*4
  for m,cg in zip(D['fourth']['monomials_bivariate_through4'],p):
   if sum(m)>3 or not Q(cg):continue
   t=(Q(cg),0,0,0)
   for k,n in enumerate(m):
    for _ in range(n):t=series_mul(t,oldv[yroots[k]],4)
   q=tuple(a+b for a,b in zip(q,t))
  copyvals[copybind[i]]=q
 ring=Symbolic([(i,) for i in range(4)]);bind={i:({k:c for k,c in enumerate(p) if c},{}) for i,p in (common|copyvals).items()};opaque={}
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in D['FJ'][key].items():opaque[int(i)]=(Q(p[0]),Q(p[8])+Q(p[9]))
 for r in D['first640']['all35_ordered_CA_prefix_firsts']:
  q=rec(r['first9_bar']);opaque[r['node']]=(q.get(0,Q(0)),q.get(8,Q(0)))
 for i,p in D['first640']['all8_g_quotient_overrides']:
  q=rec(p);opaque[i]=(q.get(0,Q(0)),q.get(8,Q(0)))
 for i,(c,a) in opaque.items():bind[i]=({k:v for k,v in enumerate([c,a]) if v},{(i,(7,7),2):Q(1),(i,(7,7,7),3):Q(1)})
 ga=[Q(0)]*4
 for m,cg in zip(source3['monomials'],source3['g_CA_all10_coefficients']):ga[sum(m)]+=Q(cg)
 for j,i in enumerate(D['C0']['actual_C0_basis8_scalar_nodes']):bind[i]=({k:c for k,c in enumerate(ga) if c},{}) if j==7 else ring.c(0)
 copyeq=D['copy1']['all57_copied_equation_indices'];roots=sorted(rawroots|set(copyeq));order=closure(nodes,roots,set(bind));uses=collections.Counter()
 for i in order:
  if i in bind:continue
  n=nodes[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
 uses.update(roots);vv={};peak=0
 for i in order:
  n=nodes[i];op=n[0]
  if i in bind:p=bind[i]
  elif op=='c':p=ring.c(Q(n[1],n[2]))
  elif op=='a':p=ring.add(vv[n[1]],vv[n[2]])
  elif op=='n':p=ring.scale(vv[n[1]],-1)
  elif op=='m':p=ring.mul(vv[n[1]],vv[n[2]])
  elif op=='i':p=ring.inv(vv[n[1]]);units+=1
  else:raise AssertionError((i,n))
  vv[i]=p;peak=max(peak,sum(len(x)+len(y) for x,y in vv.values()));assert peak<150000
  if i not in bind:
   for j in n[1:3] if op in ['a','m'] else n[1:] if op in ['n','i'] else []:
    uses[j]-=1
    if not uses[j]:vv.pop(j,None)
 assert all(vv[i]==ring.c(0) for i in copyeq)
 raw={i:vv[i] for i in rawroots}
 # Forward typed B11 reader with symbolic absent atom2 coefficients.
 def ex(m):return tuple(m.count(i) for i in range(8))
 bm=[(0,)*8,ex((7,)),ex((7,7)),ex((0,)),ex((1,)),ex((4,)),ex((5,)),ex((0,7)),ex((1,7)),ex((4,7)),ex((5,7))];br=Symbolic(bm)
 cv={i:rec(p) for i,p,g in D['cache']['all108_composite_values_and_opaque_second_gradients']};assert all(not g for i,p,g in D['cache']['all108_composite_values_and_opaque_second_gradients'])
 atom={i:rec(p) for i,p in D['first640']['all79_atom_values']};lin={0:0,8:1,1:3,2:4,5:5,6:6};todo=[k[1] for k in usedlabels if k[0]==0];needed=set()
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
 bv={}
 for i in sorted(needed):
  n=ops[i]
  if n[0]=='c':p=br.c(Q(n[1],n[2]))
  elif n[0]=='actual_C0_scalar_function':
   node=n[1]
   if node in cv:p=(cv[node],{})
   else:p=({b:atom[node][j] for j,b in lin.items() if atom[node].get(j)},{(node,tuple(k for k,cg in enumerate(bm[b]) for _ in range(cg)),b):Q(1) for b in [2,7,8,9,10]})
  elif n[0]=='sum':
   p=br.c(0)
   for j,a,b in n[1]:p=br.add(p,bv[j],Q(a,b))
  elif n[0]=='product':
   p=br.c(Q(n[1],n[2]))
   for j in n[3]:p=br.mul(p,bv[j])
  else:raise AssertionError(n)
  bv[i]=p
 native=[{r:rec(p) for r,p in c} for c in D['native']['source_columns128_rows48_sparse_B11_effective_second']]
 pairB={(0,7):7,(1,7):8,(4,7):9,(5,7):10,(7,7):2}
 assert all(m2[k[-1]] in pairB for k in usedlabels)
 actualg2={i:Q(c) for i,c in D['fifth_peer']['actual_eight_regular_g2_verified']}
 def coefficient(p,k):return p[0].get(k,Q(0)),{(i,tag):c for (i,tag,out),c in p[1].items() if out==k}
 def value(p3,xp,f,k0):
  z=k0;g={}
  for roots,degree in [(p3,3),(xp,2)]:
   for i,c in roots.items():a,b=coefficient(raw[i],degree);z+=c*a;g=add(g,b,c)
  for key,c in f.items():
   b=pairB[m2[key[-1]]]
   if key[0]==2:a=native[key[1]].get(key[2],{}).get(b,Q(0));gg={}
   elif key[0]==1:a=cv[key[1]].get(b,Q(0));gg={}
   else:a,gg=coefficient(bv[key[1]],b)
   z+=c*a;g=add(g,gg,c)
  assert all(tag==(7,7) and node in actualg2 for node,tag in g),('unresolved_real_coefficient',g)
  return z+sum(c*actualg2[node] for (node,tag),c in g.items()),z,g
 braw,bknown,bgrad=value(*bformula);b3=braw/Q(-250,3);assert b3==Q(-42215,3456)
 answers=[];known=[];sens=[]
 for p3,xp,f0,fb,k0,kb in formulas:
  a,ak,ag=value(p3,xp,f0,k0);b,bk,bg=value({}, {},fb,kb);answers.append(a+b3*b);known.append(ak+b3*bk);sens.append(add(ag,bg,b3))
 assert all(a==0 for a in answers)
 assert D['author_values']['canonical_qB_pure_cubic_b3']==str(b3)
 assert D['author_values']['all118_actual_corrected_pure_quartic_values']==list(map(str,answers))
 ev={'source_manifest':inputs,'independent_method':'forward derivative/Euler source scattering, own118 original768 composition, symbolic ordinary source3 replay and forward B11 typed reader',
  'Gamma0_independent':enc(g0),'Gamma_b3_independent':enc(gB),'b3_corrected_F620_563':str(braw),'canonical_b3':str(b3),
  'b3_independent_known_part':str(bknown),'b3_independent_missing_g2_sensitivity':enc(bgrad),
  'all118_source_functionals_and_correction_signs_independently_matched':True,
  'all105_original_and57_copied_source_equations_through3_checked':True,
  'all_ordinary_quotient_second_and_third_sensitivities_retained_then_canceled_except_accepted_g2':True,
  'all118_independent_known_parts':list(map(str,known)),'all118_independent_g2_sensitivities':list(map(enc,sens)),
  'all118_actual_corrected_values':list(map(str,answers)),'all118_zero':True,
  'source_unit_inverse_checks':units,'symbolic_live_QQ_peak':peak,
  'negative_scope_only_restricted118_family':True,'actual_HF4_or_exact_T_length_computed':False,
  'no_new_source_coefficient_solve_jet_point_prime_or_rank':True}
 ep=out/'independent-b3-and-all118-source-postverification.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'canonical_b3_independently_verified':str(b3),
  'all118_source_and_all118_zero_values_independently_verified':True,
  'actual_HF4_or_exact_T_length_claimed':False,'scope':'negative118 family only',
  'input_hashes':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep),'phase1-independent-Gamma-b3-and118-source.json':sha(out/'phase1-independent-Gamma-b3-and118-source.json')},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)

if __name__=='__main__':main()
