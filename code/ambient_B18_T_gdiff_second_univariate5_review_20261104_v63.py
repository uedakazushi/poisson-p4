#!/usr/bin/env python3
"""Read-only provided114 source5 verification and ALL619 symbolic cubic test.

No unknown solve, author import, new source jet, or numerical matrix search.
The scalar DAG is evaluated with all missing regular quotient second/third
coefficients symbolic until the final619 contractions cancel them.
"""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW'
A='AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63'
PATHS={
 'author_report':f'reports/{A}/report.json',
 'CA':f'reports/{A}/CA-source-and-regular-quotient-second.json',
 'CB':f'reports/{A}/CB-source-and-regular-quotient-second.json',
 'result':f'reports/{A}/actual-gdiff-second-and-complete619-pure-cubic-values.json',
 'input':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63/exact-two57-univariate5-input.json',
 'input_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63_REVIEW/report.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'complete':'reports/AMBIENT_B18_T_PURE_CUBIC619_COMPLETENESS_PREFLIGHT_20261104_v63/independent-pure619-Schur-Theta-and-completeness-preflight.json',
 'complete_report':'reports/AMBIENT_B18_T_PURE_CUBIC619_COMPLETENESS_PREFLIGHT_20261104_v63/report.json',
 'adjoint':'reports/AMBIENT_B18_T_CONSTANT_QUADRATIC_KERNEL_R3_ADJOINT_20261104_v63/phase1-theta-and-pure-r3-second-adjoint.json',
 'functionals':'reports/AMBIENT_B18_T_CONSTANT_QUADRATIC_KERNEL_R3_ADJOINT_20261104_v63/phase2-all619-pure-cubic-source-functionals.json',
 'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
 'economy':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase2-known-parts-and-missing-independent-projection.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'S':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'w0':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
 'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
 'fourth':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'copy1':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'copy2':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'atom1':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'cache':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
 'cache_peer':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json',
 'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
 'native_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(x):return {int(i):Q(c) for i,c in x if Q(c)}
def enc(x):return [[list(k) if isinstance(k,tuple) else k,str(c)] for k,c in sorted(x.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for k,v in b.items():
  z=a.get(k,Q(0))+c*v
  if z:a[k]=z
  else:a.pop(k,None)
 return a
def put(a,k,c):
 z=a.get(k,Q(0))+c
 if z:a[k]=z
 else:a.pop(k,None)
def series_mul(a,b,n):return tuple(sum(a[i]*b[k-i] for i in range(k+1)) for k in range(n))
def series_inv(a,n):
 assert a[0];q=[1/a[0]]
 for k in range(1,n):q.append(-sum(a[i]*q[k-i] for i in range(1,k+1))/a[0])
 q=tuple(q);assert series_mul(a,q,n)==(1,)+(0,)*(n-1);return q
def closure(nodes,roots,stops):
 seen=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in stops:continue
  p=nodes[i]
  if p[0] in ['a','m']:todo.extend(p[1:3])
  elif p[0] in ['n','i']:todo.append(p[1])
  else:assert p[0]=='c',('uncached_leaf',i,p)
 return sorted(seen)
def numeric_eval(nodes,roots,bind,n,order=None):
 order=order or closure(nodes,roots,set(bind));uses=collections.Counter()
 for i in order:
  if i in bind:continue
  p=nodes[i];uses.update(p[1:3] if p[0] in ['a','m'] else p[1:] if p[0] in ['n','i'] else [])
 uses.update(roots);vals={};inverse_checks=0
 for i in order:
  p=nodes[i];op=p[0]
  if i in bind:v=bind[i]
  elif op=='c':v=(Q(p[1],p[2]),)+(Q(0),)*(n-1)
  elif op=='a':v=tuple(a+b for a,b in zip(vals[p[1]],vals[p[2]]))
  elif op=='n':v=tuple(-a for a in vals[p[1]])
  elif op=='m':v=series_mul(vals[p[1]],vals[p[2]],n)
  elif op=='i':v=series_inv(vals[p[1]],n);inverse_checks+=1
  else:raise AssertionError(p)
  vals[i]=v
  if i not in bind:
   for j in p[1:3] if op in ['a','m'] else p[1:] if op in ['n','i'] else []:
    uses[j]-=1
    if not uses[j]:vals.pop(j,None)
 return {i:vals[i] for i in roots},inverse_checks

class Symbolic:
 """Finite monomial algebra plus formal order>=2 coefficient perturbations."""
 def __init__(self,mons):
  self.mons=mons;self.index={m:i for i,m in enumerate(mons)};self.table={(i,j):self.index.get(tuple(a+b for a,b in zip(m,n))) for i,m in enumerate(mons) for j,n in enumerate(mons)}
 def c(self,x):return ({0:Q(x)} if x else {},{})
 def add(self,a,b,c=Q(1)):return add(a[0],b[0],c),add(a[1],b[1],c)
 def scale(self,a,c):return ({i:c*v for i,v in a[0].items() if c*v},{i:c*v for i,v in a[1].items() if c*v})
 def mul(self,a,b):
  p={};g={}
  for i,c in a[0].items():
   for j,d in b[0].items():
    k=self.table[i,j]
    if k is not None:put(p,k,c*d)
  for x,y in [(a,b),(b,a)]:
   for (node,tag,i),c in x[1].items():
    for j,d in y[0].items():
     k=self.table[i,j]
     if k is not None:put(g,(node,tag,k),c*d)
  # Unknown*unknown has order>=4 and is zero in both algebras here.
  return p,g
 def inv(self,a):
  c=a[0].get(0,Q(0));assert c
  h=self.scale(self.add(a,self.c(-c)),-1/c);v=self.c(1);t=self.c(1)
  for _ in range(3):t=self.mul(t,h);v=self.add(v,t)
  v=self.scale(v,1/c);assert self.mul(a,v)==self.c(1);return v

def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'independent_postverification_inputs_saved_no_unknown_solve','inputs':inputs},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['author_report','input_peer','complete_report','joint_peer','cache_peer','native_peer'])
 for p,h in D['author_report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['author_report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 # FIRST PHASE: provided114 fifth coefficients substituted, not recomputed.
 src=D['source']['DAG_scalar_nodes'];g2={};source_checks=[]
 for f in D['input']['families']:
  name=f['family'];actual=D[name];u=f['original_unknown57_order'];E=f['original_equations57'];g=f['actual_g_original_node'];bn=f['actual_bbar4_original_nodes'];candidate={int(i):tuple(map(Q,p)) for i,p in actual['all57_source_jets_through5']}
  assert list(candidate)==u
  old={int(i):tuple(map(Q,p)) for i,p in f['cached_univariate_source_through4']}
  assert all(candidate[i][:5]==old[i] for i in u)
  values,units=numeric_eval(src,E+[g]+bn,candidate,6,f['source_needed_original_DAG_nodes']);slopes=list(map(Q,f['last2_target_slopes']))
  for j,i in enumerate(E):assert values[i]==(0,slopes[j-55] if j>=55 else 0,0,0,0,0)
  assert [[i,list(map(str,values[i]))] for i in E]==actual['all57_source_equations_through5']
  gp=values[g];assert gp[:3]==(0,0,0) and gp[3] and not gp[4]
  vals=[]
  for b,h1 in zip(bn,map(Q,f['cached_four_quotient_first_line'])):
   p=values[b];assert p[:4]==(0,0,0,0) and p[4]==gp[3]*h1
   h2=p[5]/gp[3];assert series_mul(gp,(Q(0),h1,h2,Q(0),Q(0),Q(0)),6)==p;vals.append(h2)
  assert list(map(str,vals))==actual['new_four_regular_quotient_second_line']
  ids=[566400,566401,566402,566403] if name=='CA' else [330164,330289,330524,330594]
  g2.update(zip(ids,vals));source_checks.append({'family':name,'provided57_fifth_coefficients':[[i,str(candidate[i][5])] for i in u],
   'all57_equations_all6_coefficients_checked':True,'actual_g':list(map(str,gp)),'actual_bbar4':[list(map(str,values[i])) for i in bn],
   'actual_four_h2':list(map(str,vals)),'unit_inverse_checks':units,'source5_uniqueness':'same saved lower4 coefficients and accepted invertibleJ0 imply the fifth coefficients are unique'})
 (out/'phase1-independent114-source5-and-g2.json').write_text(json.dumps({'source_checks':source_checks,'actual_eight_g2':[[i,str(v)] for i,v in sorted(g2.items())],'no_unknown_solve':True},separators=(',',':'))+'\n')
 (out/'phase.json').write_text(json.dumps({'phase':'all114_provided_values_source_equations_verified','seconds':time.monotonic()-start})+'\n')
 # Remaining phases reconstruct all619 source functionals and retain every
 # unknown ordinary regular coefficient until final contraction.
 verify619(D,g2,out,inputs,start)

def verify619(D,g2,out,inputs,start):
 C,S,M=D['complete'],D['S'],D['m1'];theta=rec(C['Theta_extension180_sparse']);ys=list(map(Q,C['actual_Ystar2_last_rr84']))
 wb=list(map(rec,D['w0']['complete_rank_certificate']['basis_normalized_columns']));wl={min(p):j for j,p in enumerate(wb)};qnorm=list(map(rec,D['low']['normalized84_unitcycle_Q1_columns_sparse']))
 # Reconstruct the quadratic source adjoint FORWARDS from Theta, U0, Y1,Y2.
 J=[];H=[]
 for w in range(84):
  v={w:Q(1)};coeff={}
  while v:
   k=min(v);p=wl[k];a=v[k];put(coeff,p,a);v=add(v,wb[p],-a)
  y={};h=Q(0)
  for p,a in coeff.items():y=add(y,qnorm[p],a);h+=a*ys[p]
  J.append(y);H.append(h)
 ms2=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(ms2)}
 va={}
 for w,(sector,pair) in enumerate(D['w0']['literal_wedge_basis84']):
  i,j=pair;terms=[(8*sector+7,i,1,1)] if j==7 else [(8*sector+j,i,1,2),(8*sector+i,j,-1,2)]
  weights={}
  for v in range(8):
   z=H[w] if v==7 else Q(0)
   for tv,c in J[w].items():
    t,k=divmod(tv,8);z-=c*theta.get(36*t+mi[tuple(sorted((v,k)))],Q(0))
   if z:weights[v]=z
  for row,var,sign,shift in terms:
   for n,mon in enumerate(ms2):
    if var not in mon:continue
    small=list(mon);small.remove(var);v=small[0]
    if v in weights:put(va,36*row+n,Q(sign*mon.count(var),int(v<7)+shift)*weights[v])
 m0=[list(map(Q,r)) for r in D['V']['all24_by56_m0']];U=D['X0']['free_kernel_rows32'] if 'X0' in D else D['m1']['free_kernel_rows32']
 pi=[{U[j]:c for j,c in rec(r).items()} for r in S['constant_free32_quotient_projection5_rows']];adj={}
 for index,c in theta.items():
  t,n=divmod(index,36)
  for src,a in pi[t].items():put(adj,36*src+n,c*a)
 for index,c in va.items():
  row,n=divmod(index,36)
  for src,a in enumerate(m0[row]):
   if a:put(adj,36*src+n,c*a)
 assert enc(adj)==D['adjoint']['theta_plus_Ystar2_beta_quadratic_source_adjoint56_times36']
 # Independently reconstruct every rawP3 / Xp2 / second-source contraction.
 pcols=D['C0']['honest_D0_relation_columns'];nodes=D['C0']['DAG_scalar_nodes'];ops=D['chain']['operator_scalar_nodes'];bcols=[dict(x) for x in D['chain']['B48_by640_sparse_columns']]
 source620=list(map(rec,M['620_actual_source_combinations_in_original768']));a1=list(map(rec,M['620_actual_linear_pivot_coefficients_a1_by27_times8']));xp1=list(map(rec,M['actual_pivot_X1_by27_sparse56_times8']))
 xp_roots=[]
 for c in S['first27_source_combinations']:
  r={}
  for j,a in rec(c).items():
   if j<640:put(r,pcols[j][8],a)
  xp_roots.append(r)
 all620=[]
 for c,a in zip(source620,a1):
  p3={};xp2={};x2={};known=Q(0)
  for j,b in c.items():
   if j<640:put(p3,pcols[j][8],b)
  for p in range(27):
   ar=a.get(8*p+7,Q(0))
   for root,b in xp_roots[p].items():put(xp2,root,-ar*b)
  for index,t in adj.items():
   row,n=divmod(index,36);u,v=ms2[n]
   for p,first in enumerate(xp1):
    known+=t*a.get(8*p+u,Q(0))*first.get(8*row+v,Q(0))
    if u!=v:known+=t*a.get(8*p+v,Q(0))*first.get(8*row+u,Q(0))
   for j,b in c.items():
    key=None
    if j>=640:
     if row<48:key=(2,j-640,row,n)
    elif row>=48:
     root=pcols[j][row-48+1]
     if nodes[root][0]!='c':key=(1,root,n)
    else:
     root=bcols[j].get(row)
     if root is not None and ops[root][0]!='c':key=(1,ops[root][1],n) if ops[root][0]=='actual_C0_scalar_function' else (0,root,n)
    if key is not None:put(x2,key,-t*b)
  all620.append(({i:c for i,c in p3.items() if nodes[i][0]!='c'},{i:c for i,c in xp2.items() if nodes[i][0]!='c'},x2,known))
 final=[]
 for index,c in enumerate(C['actual619_constant_quadratic_kernel_combinations620']):
  p3={};xp2={};x2={};known=Q(0)
  for j,a in rec(c).items():
   p,q,r,s=all620[j];p3=add(p3,p,a);xp2=add(xp2,q,a);x2=add(x2,r,a);known+=a*s
  row=D['functionals']['all619'][index]
  assert enc(p3)==row['rawP8_last_cubic_roots'] and enc(xp2)==row['Xp_second_r2_roots'] and enc(x2)==row['quadratic_source_root_correction'] and str(known)==row['known_Xp1_a1_contraction']
  final.append((p3,xp2,x2,known))
 (out/'phase.json').write_text(json.dumps({'phase':'all619_source_operators_independently_reconstructed','seconds':time.monotonic()-start})+'\n')
 # Actual existing common105 and copiedCB57 source3, evaluated without solve.
 source3=D['joint'];common={}
 for i,p in source3['all105_unknown_jets']:
  q=[Q(0)]*4
  for m,c in zip(source3['monomials'],p):q[sum(m)]+=Q(c)
  common[i]=tuple(q)
 yroots=[i for i,p in D['copy2']['actual_YB_last2_from_saved105']];eq105=source3['all105_equation_original_node_order']
 oldv,unit105=numeric_eval(nodes,eq105+yroots,common,4)
 for j,i in enumerate(eq105):assert oldv[i]==(0,1 if j in [55,56] else 0,0,0)
 cb=D['fourth']['families']['CB'];copybind=dict(D['copy1']['explicit57_source_copy_bindings']);copyvalues={};ownmons=D['fourth']['monomials_bivariate_through4']
 for i,p in cb['all57_source4_jets']:
  q=(Q(0),)*4
  for m,c in zip(ownmons,p):
   if sum(m)>3 or not Q(c):continue
   t=(Q(c),0,0,0)
   for k,n in enumerate(m):
    for _ in range(n):t=series_mul(t,oldv[yroots[k]],4)
   q=tuple(a+b for a,b in zip(q,t))
  copyvalues[copybind[i]]=q
 ring=Symbolic([(i,) for i in range(4)])
 knownbind={i:({k:c for k,c in enumerate(p) if c},{}) for i,p in (common|copyvalues).items()}
 opaque={}
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in D['FJ'][key].items():opaque[int(i)]=(Q(p[0]),Q(p[8])+Q(p[9]))
 for r in D['atom1']['all35_ordered_CA_prefix_firsts']:
  q=rec(r['first9_bar']);opaque[r['node']]=(q.get(0,Q(0)),q.get(8,Q(0)))
 for i,p in D['atom1']['all8_g_quotient_overrides']:
  q=rec(p);opaque[i]=(q.get(0,Q(0)),q.get(8,Q(0)))
 for i,(c,a) in opaque.items():knownbind[i]=({k:v for k,v in enumerate([c,a]) if v},{(i,(7,7),2):Q(1),(i,(7,7,7),3):Q(1)})
 ga=[Q(0)]*4
 assert source3['g_CA_all10_coefficients']==D['C0']['cached_actual_CA_gA_all_degree3_coefficients']
 for m,c in zip(source3['monomials'],source3['g_CA_all10_coefficients']):ga[sum(m)]+=Q(c)
 for j,i in enumerate(D['C0']['actual_C0_basis8_scalar_nodes']):knownbind[i]=({k:c for k,c in enumerate(ga) if c},{}) if j==7 else ring.c(0)
 proots={i for p,q,r,k in final for i in set(p)|set(q)};copyeq=D['copy1']['all57_copied_equation_indices'];roots=sorted(proots|set(copyeq))
 order=closure(nodes,roots,set(knownbind));uses=collections.Counter()
 for i in order:
  if i in knownbind:continue
  n=nodes[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
 uses.update(roots);v={};unit3=0;peak=0
 for i in order:
  n=nodes[i];op=n[0]
  if i in knownbind:p=knownbind[i]
  elif op=='c':p=ring.c(Q(n[1],n[2]))
  elif op=='a':p=ring.add(v[n[1]],v[n[2]])
  elif op=='n':p=ring.scale(v[n[1]],-1)
  elif op=='m':p=ring.mul(v[n[1]],v[n[2]])
  elif op=='i':p=ring.inv(v[n[1]]);unit3+=1
  else:raise AssertionError((i,n))
  v[i]=p;peak=max(peak,sum(len(a)+len(b) for a,b in v.values()));assert peak<150000
  if i not in knownbind:
   for j in n[1:3] if op in ['a','m'] else n[1:] if op in ['n','i'] else []:
    uses[j]-=1
    if not uses[j]:v.pop(j,None)
 assert all(v[i]==ring.c(0) for i in copyeq)
 raw={i:v[i] for i in proots}
 # Forward typed-root B11 evaluation with the SAME symbolic rr coefficient
 # labels as the common source3 calculation. All other mixed seconds remain.
 def ex(m):return tuple(m.count(i) for i in range(8))
 bmons=[(0,)*8,ex((7,)),ex((7,7)),ex((0,)),ex((1,)),ex((4,)),ex((5,)),ex((0,7)),ex((1,7)),ex((4,7)),ex((5,7))];br=Symbolic(bmons)
 cv={i:rec(p) for i,p,g in D['cache']['all108_composite_values_and_opaque_second_gradients']};assert all(not g for i,p,g in D['cache']['all108_composite_values_and_opaque_second_gradients'])
 atoms={i:rec(p) for i,p in D['atom1']['all79_atom_values']};lin={0:0,8:1,1:3,2:4,5:5,6:6};needed=set();todo=[k[1] for p,q,r,z in final for k in r if k[0]==0]
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
   else:p=({b:atoms[node][j] for j,b in lin.items() if atoms[node].get(j,Q(0))},{(node,tuple(k for k,c in enumerate(bmons[b]) for _ in range(c)),b):Q(1) for b in [2,7,8,9,10]})
  elif n[0]=='sum':
   p=br.c(0)
   for j,a,b in n[1]:p=br.add(p,bv[j],Q(a,b))
  elif n[0]=='product':
   p=br.c(Q(n[1],n[2]))
   for j in n[3]:p=br.mul(p,bv[j])
  else:raise AssertionError(n)
  bv[i]=p
 native=[{i:rec(p) for i,p in col} for col in D['native']['source_columns128_rows48_sparse_B11_effective_second']]
 pairB={(0,7):7,(1,7):8,(4,7):9,(5,7):10,(7,7):2}
 def coefficient(jet,k):return jet[0].get(k,Q(0)),{(node,tag):c for (node,tag,output),c in jet[1].items() if output==k}
 known=[];gradients=[];answers=[]
 for j,(p3,xp2,x2,k0) in enumerate(final):
  z=k0;g={}
  for roots,degree in [(p3,3),(xp2,2)]:
   for i,c in roots.items():a,b=coefficient(raw[i],degree);z+=c*a;g=add(g,b,c)
  for key,c in x2.items():
   b=pairB[ms2[key[-1]]]
   if key[0]==2:a=native[key[1]].get(key[2],{}).get(b,Q(0));gg={}
   elif key[0]==1:a=cv[key[1]].get(b,Q(0));gg={}
   else:a,gg=coefficient(bv[key[1]],b)
   z+=c*a;g=add(g,gg,c)
  assert all(tag==(7,7) and node in g2 for node,tag in g),'uncancelled ordinary missing derivative'
  primitive={node:c for (node,tag),c in g.items()}
  assert primitive==rec(D['missing']['all619'][j]['opaque_second_r2_sensitivity'])
  assert str(z)==D['economy']['all619_known_part_missing_coefficients_split'][j]
  actual=z+sum(c*g2[node] for node,c in primitive.items());assert actual==0
  known.append(z);gradients.append(primitive);answers.append(actual)
 assert list(map(str,answers))==D['result']['all619_actual_pure_cubic_relation_values']
 assert [[i,str(c)] for i,c in sorted(g2.items())]==D['result']['actual_eight_regular_g_second_values']
 # Optional economy rows are checked as true symbolic identities too.
 for r in D['economy']['all_missing_annihilator_combinations']:
  c=rec(r['combination619']);gg={};z=Q(0)
  for j,a in c.items():gg=add(gg,gradients[j],a);z+=a*known[j]
  assert not gg and z==0
 assert len(D['economy']['all_missing_annihilator_combinations'])==615
 evidence={'actual_source114_substitution_phase':'phase1-independent114-source5-and-g2.json',
  'all619_source_functionals_independently_reconstructed':True,'Theta_source_adjoint31_independent':enc(adj),
  'all105_original_source_equations_through3_verified':True,'all57_copy_source_equations_through3_verified':True,
  'actual_YB_last2_through3':[[i,list(map(str,oldv[i]))] for i in yroots],
  'all615_economy_annihilator_rows_independently_zero':True,
  'all619_independent_known_parts':list(map(str,known)),'all619_independent_missing_second_sensitivities':list(map(enc,gradients)),
  'all_ordinary_Hadamard_second_third_and_g_third_sensitivities_cancel':True,
  'all619_actual_values':list(map(str,answers)),
  'all619_actual_values_zero':True,'actual_four_old_minus_new_g2':[str(g2[a]-g2[b]) for a,b in zip([330164,330289,330524,330594],[566400,566401,566402,566403])],
  'source_unit_inverse_checks':unit105+unit3,'symbolic_source_peak_QQ_support':peak,
  'actual_HF_through_degree3':[5,3,2,1],'actual_T_length_interval':[11,17],
  'actual_T_total_length_exact':None,'full_nil_cube_not_inferred':True,
  'no_new_unknown_source_solve_jet_or_rank_search':True,'source_manifest':inputs}
 ep=out/'independent-complete619-symbolic-source-postverification.json';blob=json.dumps(evidence,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 p1=out/'phase1-independent114-source5-and-g2.json'
 report={'id':STEM,'all_passed':True,'provided114_source5_coefficients_verified_without_solving':True,
  'all114_source_equations_all6_coefficients_and_uniqueness_J0_accepted':True,
  'actual_eight_regular_g2_verified':[[i,str(c)] for i,c in sorted(g2.items())],
  'all619_source_functionals_and_all619_zero_values_independently_accepted':True,
  'all615_economy_symbolic_rows_independently_zero':True,
  'actual_HF_through_degree3_accepted':[5,3,2,1],'actual_T_length_interval_accepted':[11,17],
  'exact_total_length_or_full_nil_cube_claimed':False,'input_hashes':inputs,
  'evidence_sha256':{p1.name:sha(p1),ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
