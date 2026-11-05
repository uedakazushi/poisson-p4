#!/usr/bin/env python3
"""Prepared independent remaining-quadratic source projector.
Only fixed polynomial stencil operations; no producer homotopy imports.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
from ambient_B18_honest_nil_T_gauge120_review_20261105_v63 import Base
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
PAIRS=list(itertools.combinations(range(8),2));MON2=list(itertools.combinations_with_replacement(range(8),2))
def homotopy_first_stencil(sector,row,mon):
 """Linear output of the ORIGINAL ordered-Hadamard/Euler formula on v2.
 The inadmissible pure-r last component is excluded by the source-cycle proof.
 """
 if row==7 and mon==(7,7):return{}
 w=[{}for _ in range(7)];cross={}
 if row<7:w[row][mon]=F(1)
 else:
  zi=max(i for i in mon if i<7);other=list(mon);other.remove(zi);cross[zi]=other[0]
  if sector<2:w[zi][tuple(sorted((7,other[0])))]=F(1)
 out={}
 for k,(i,j)in enumerate(PAIRS):
  if j==7:
   if i in cross:out[28*sector+k,cross[i]]=F(1)
   continue
  for inp,var,sign in[(w[j],i,F(1)),(w[i],j,F(-1))]:
   for p,c in inp.items():
    count=p.count(var)
    if count:
     other=list(p);other.remove(var);remaining=other[0];Q.put(out,(28*sector+k,remaining),sign*c*count/(2+int(remaining<7)))
 return out
def rebuild_second_adjoints(D):
 """Derivative of full fixed Schur/operator formula with actual C1 carry.
 Caller must first verify all actual C1 data and source frame constraints.
 """
 pi=[Q.pairs(v)for v in D['quad']['exact_quotient_duals_rows_sparse']]
 M=D['m0'];m0=[Q.sparse(v)for v in M['m0_24_by44']];free=[i for i in range(44)if i not in M['selected_source_columns24']]
 I=[free[i]for i in D['minimal']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in D['minimal']['fixed_minor_data']['inverse']]
 xp=D['X1']['X0_pivot_columns'];raw1=[{int(i):Q.pairs(v)for i,v in col}for col in D['X1']['all780_source_columns44_sparse_dual9']]
 mu1=list(map(F,D['first2']['true_actual_m1_row7_r_coefficients44']));hi=[[F(c)for c in row]for row in D['first3']['all54_actual_pivot_h0_readback']]
 AF=[[F(c)for c in row]for row in D['first2']['all48_actual_adapter_row20_first9']]
 L=[{(h,v):F(c)for h,v,c in row}for row in D['carry']['all6_actual_effective_H1_callbacks_row84_var8_sparse']]
 H={(s,r,q):homotopy_first_stencil(s,r,m)for s in range(3)for r in range(8)for q,m in enumerate(MON2)}
 AD=[[{}for _ in MON2]for _ in range(6)];AB=[[{}for _ in MON2]for _ in range(6)]
 for f in range(6):
  for q,mon in enumerate(MON2):
   for j in range(44):
    alpha={k:row.get(I.index(j),F(0))for k,row in enumerate(Ji)}if j in I else{}
    v=F(0)
    for target,rr in enumerate([42,43]):
     v+=pi[f].get(36*target+q,F(0))*(F(j==rr)-sum((raw1[p].get(rr,{}).get(0,F(0))*alpha.get(k,F(0))for k,p in enumerate(xp)),F(0)))
    for n,zi in enumerate([0,4,5]):
     c=pi[f].get(36*(n+2)+35,F(0))
     if mon==(zi,7):v+=c*mu1[j]
     if mon==(7,7):v-=c*sum((hi[n][k]*alpha.get(k,F(0))for k in range(18)),F(0))
    for s in range(3):
     for r in range(8):
      a=m0[8*s+r].get(j,F(0))
      if a:v-=a*sum((c*H[s,r,q].get(k,F(0))for k,c in L[f].items()),F(0))
    if v:AD[f][q][j]=v
   for row in range(48):
    v=F(0)
    for n,zi in enumerate([0,4,5]):
     c=F(5,2)*pi[f].get(36*(n+2)+35,F(0))
     if mon==(7,7):v+=c*AF[row][zi+1]
     if mon==(zi,7):v+=c*AF[row][8]
    if v:AB[f][q][row]=v
 return AD,AB
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 primary='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63';rep=load(primary);assert pins[f'reports/{primary}/report.json']=='6d9ffb9faa0b88fd7cea7db859cdd840933274f9f009f93cbc16fbec91819e5d'
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 assert load('AMBIENT_B18_HONEST_NIL_T_GAUGE120_20261105_v63_REVIEW')['all_passed']
 assert load('AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63_REVIEW')['all_passed']
 names={
 'physics':('AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63','phase1-common105-and-physical-source-base25.json'),
 'proper':('AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63','phase3-actual-proper-G8-base25.json'),
 'honest':('AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63','phase3-full661-honest-second-sensitivities.json'),
 'honestgrad':('AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63','phase2-allB-and-P-true-honest-second-gradients.json'),
 'source2':('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63','full9-firstjet-certificate.json'),
 'adapter2':('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63','actual-special-fiber20-coefficients.json'),
 'C0':('AMBIENT_B18_T0_HONEST_C0_20261104_v63','honest-C0-D0-D4-DAG.json'),
 'atomfirst':('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','phase1-actual-atom-firsts.json'),
 'copy':('AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63','phase1-exact-leaf-label-transport.json'),
 'oldcache':('AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD','phase3-oldB-common-source2-values-and-symbols.json'),
 'chain':('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63','actual-UP-QB-640-typed-operator.json'),
 'radical':('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json'),
 'GL':('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63','actual-E4-nil-two-generators-and-relations.json'),
 'X1':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json'),
 'm0':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json'),
 'minimal':('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json'),
 'first2':('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase2-row20-surrogates-and-true-m1row7.json'),
 'first3':('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json'),
 'carry':('AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63','actual81-pivot-first-carries-and-six-H1-source-adjoints.json'),
 'higher':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63','actual661-first-corrections-m2-replacement-and-native-highest-sensitivities.json'),
 'free43':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_FREE43_CACHED_20261105_v63','actual-one-by661-free43-quadratic-values-and-source-sensitivities.json'),
 'operator':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63','actual-complete-six-by661-quadratic-source-operator.json'),
 'quad':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63','known-quadratic-closure-and-complete814-source-projection.json'),
 'Koszul':('AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63','actual-reduced-Cbar-three-sector-Koszul.json')}
 D={k:load(*p)for k,p in names.items()};p1=load(primary,'phase1-all-source-root-second-representatives.json');p2=load(primary,'phase2-actual-native-and-m1-source-inputs.json');p3=load(primary,'phase3-all661-actual-six-quadratic-values-and-rank.json')
 save('phase.json',dict(phase='all immutable actual source/frame peers pinned; fixed-rank source replay starts'))
 B=Base(D['physics']['base25_monomials']);MI=B.idx;MON=B.mon
 # Independent scalar source2 composition for copied CB; common105 itself
 # was already independently checked in the gauge120 peer.
 YA=[{MI[(i,)]:F(1)}if(i,)in MI else{}for i in range(7)]+[{MI[(7,)]:F(1)}]*2;YB=[]
 for lin,quad in zip(D['adapter2']['YB_in_YA_linear'],D['adapter2']['YB_in_YA_quadratic']):
  z={}
  for i,c in enumerate(lin):z=Q.plus(z,Q.scaled(YA[i],F(c)))
  for e,c in quad:
   q={0:F(c)}if F(c)else{}
   for i,n in enumerate(e):
    for _ in range(n):q=B.mul(q,YA[i])
   z=Q.plus(z,q)
  YB.append(z)
 def compose(v,args):
  z={}
  for mon,c in zip(D['source2']['source2_monomials'],v):
   if not F(c):continue
   q={0:F(c)}
   for i in mon:q=B.mul(q,args[i])
   z=Q.plus(z,q)
  return z
 bind={i:(Q.sparse(v),{})for i,v in D['physics']['all105_actual_source2_values']};cbmap=dict(D['copy']['explicit57_source_copy_bindings'])
 for i,v in D['source2']['both57_source2']['CB']['all57_source2jets']:bind[cbmap[i]]=(compose(v,[{}]*7+YB[7:]),{})
 def firstpoly(v,full9=False):
  v=list(map(F,v));z={0:v[0]}if v[0]else{}
  for i in range(8):
   if(i,)in MI:
    a=v[8]+v[9]if full9 and i==7 else v[i+1]
    if a:z[MI[(i,)]]=a
  return z
 for key in['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,v in D['source2'][key].items():bind[int(i)]=(firstpoly(v,True),{int(i):F(1)})
 for row in D['atomfirst']['all35_ordered_CA_prefix_firsts']:
  v=Q.pairs(row['first9_bar']);bind[row['node']]=(firstpoly([v.get(i,F(0))for i in range(9)]),{row['node']:F(1)})
 for i,raw in D['atomfirst']['all8_g_quotient_overrides']:
  v=Q.pairs(raw);bind[i]=(firstpoly([v.get(j,F(0))for j in range(9)]),{i:F(1)})
 # These are actual coordinate functions and the source-proved gA in m³;
 # not arbitrary unprovided second coefficients declared zero.
 for i,node in enumerate(D['C0']['actual_C0_basis8_scalar_nodes']):bind[node]=(YA[i]if i<7 else{}, {})
 R=D['radical'];rn=R['new_Cscalar_and_R0component_nodes'];Pcols=R['actual_P8_by640_columns'];tn=D['chain']['operator_scalar_nodes'];Bcols=D['chain']['B48_by640_sparse_columns'];broots={i for col in Bcols for row,i in col};proots={rn[i][1]for col in Pcols for i in col if i}
 tc=set();cr=set(proots);todo=list(broots)
 while todo:
  i=todo.pop()
  if i in tc:continue
  tc.add(i);n=tn[i];o=n[0]
  if o=='actual_C0_scalar_function':cr.add(n[1])
  elif o=='sum':todo.extend(j for j,a,b in n[1])
  elif o=='product':todo.extend(n[3])
  elif o in['a','m','n','i']:todo.extend(n[1:])
  else:assert o=='c'
 cn=D['C0']['DAG_scalar_nodes'];extra=D['copy']['all57_copied_equation_indices'];cc=set();todo=list(cr)+extra
 while todo:
  i=todo.pop()
  if i in cc:continue
  cc.add(i)
  if i in bind:continue
  n=cn[i]
  if n[0]in['a','m','n','i']:todo.extend(n[1:])
  else:assert n[0]=='c',('UNCACHED_SOURCE_LEAF',i,n)
 def add(x,y):return Q.plus(x[0],y[0]),Q.plus(x[1],y[1])
 def scale(x,c):return Q.scaled(x[0],c),Q.scaled(x[1],c)
 def mult(x,y):return B.mul(x[0],y[0]),Q.plus(Q.scaled(x[1],y[0].get(0,F(0))),Q.scaled(y[1],x[0].get(0,F(0))))
 def inverse(x):return B.inv(x[0]),Q.scaled(x[1],-1/x[0][0]**2)
 uses=collections.Counter(list(cr)+extra)
 for i in cc:
  if i not in bind and cn[i][0]in['a','m','n','i']:uses.update(cn[i][1:])
 cv={};peak=0
 for i in sorted(cc):
  n=cn[i];o=n[0]
  if i in bind:v=bind[i]
  elif o=='c':v=({0:F(n[1],n[2])}if n[1]else{}, {})
  elif o=='a':v=add(cv[n[1]],cv[n[2]])
  elif o=='m':v=mult(cv[n[1]],cv[n[2]])
  elif o=='n':v=scale(cv[n[1]],-1)
  else:assert o=='i';v=inverse(cv[n[1]])
  cv[i]=v;peak=max(peak,sum(len(p)+len(g)for p,g in cv.values()));assert peak<100000
  if i not in bind and o in['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del cv[j]
 assert all(cv[i]==({}, {})for i in extra)
 tv={}
 for i in sorted(tc):
  n=tn[i];o=n[0]
  if o=='c':v=({0:F(n[1],n[2])}if n[1]else{}, {})
  elif o=='actual_C0_scalar_function':v=cv[n[1]]
  elif o=='sum':
   v=({}, {})
   for j,a,b in n[1]:v=add(v,scale(tv[j],F(a,b)))
  elif o=='product':
   v=({0:F(n[1],n[2])}, {})
   for j in n[3]:v=mult(v,tv[j])
  elif o=='a':v=add(tv[n[1]],tv[n[2]])
  elif o=='m':v=mult(tv[n[1]],tv[n[2]])
  elif o=='n':v=scale(tv[n[1]],-1)
  else:assert o=='i';v=inverse(tv[n[1]])
  tv[i]=v
 assert [[i,[[k,str(c)]for k,c in sorted(tv[i][0].items())]]for i in sorted(broots)]==p1['all1239_Broot_second_representatives']
 assert [[i,[[k,str(c)]for k,c in sorted(cv[i][0].items())]]for i in sorted(cr)]==p1['all768_C0root_second_representatives']
 for i,z,g in D['honestgrad']['all1239_B_root_origin_and_gradients']:assert tv[i][0].get(0,F(0))==F(z)and tv[i][1]==Q.pairs(g)
 for i,old,z,g in D['honestgrad']['all723_P_root_origin_and_gradients']:assert rn[i][1]==old and cv[old][0].get(0,F(0))==F(z)and cv[old][1]==Q.pairs(g)
 for key,values in[('C0_composite_values_and_opaque_second_gradients',cv),('typedB_root_values_and_opaque_second_gradients',tv)]:
  for i,v,g in D['oldcache'][key]:
   if i in values:assert {j:c for j,c in values[i][0].items()if j<20}==Q.pairs(v)
 save('phase1-independent-full-source-representatives-and-gradients.json',dict(all1239_B_and768_C0_representatives=True,all_original_true_opaque_gradients=True,all57_copied_equations=True,all_old_B20_overlaps=True,peak_sparse_QQ=peak,no_individual_unknown_seconds_claimed=True))
 # Direct finite-algebra source polynomials; no canonical gauge recomputation.
 exp=[tuple(m.count(i)for i in range(8))for m in MON];nu=[[F(c)for c in v]for v in D['physics']['nu_A_B_base25']];C=Q.FiniteC(exp,nu);P=D['proper']
 core=[Q.sparse(P['kappa9_C100'][i])for i in[0,1,2,6]];vs=[Q.sparse(v)for v in P['actualG8_C100'][:5]]+[C.polynomial(D['GL']['forward_GL4'][s],core)for s in['h','c','r','l']]
 native=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  polys={row:C.polynomial(p,vs)for row,p in col}
  for mask in range(4):
   z={}
   for row,p in polys.items():
    q=C.mul(p,{mask:F(1)})
    for k in range(4):
     v=C.component(q,k)
     if v:z[4*row+k]=v
   native.append(z)
 qv=[];d=C.polynomial(R['d_polynomial'],vs)
 for p in R['R9_generators_polynomials9']:
  v=C.mul(d,C.polynomial(p,vs));qv += [C.mul(v,{mask:F(1)})for mask in range(4)]
 for i in range(8):qv.append({4*MI[(i,)]+1:F(-1)}if i<7 and(i,)in MI else{})
 vg=[[C.component(v,1),C.component(v,3),Q.plus(C.component(v,2),Q.scaled(C.component(v,1),F(-11,12)))]for v in qv]
 m0=[Q.sparse(v)for v in D['m0']['m0_24_by44']];mu1=list(map(F,D['first2']['true_actual_m1_row7_r_coefficients44']))
 for j,v in enumerate(vg):
  for s in range(3):
   for i in range(7):
    if(i,)in MI:assert v[s].get(MI[(i,)],F(0))==m0[8*s+i].get(j,F(0))
   assert v[s].get(MI[(7,)]if s<2 else MI[(7,7)],F(0))==m0[8*s+7].get(j,F(0))
 assert [v[0].get(MI[(7,7)],F(0))for v in vg]==mu1
 triples=[(i,7,7)for i in[0,4,5]];mon28=MON+triples;exp28=exp+[tuple(m.count(i)for i in range(8))for m in triples];C3=Q.FiniteC(exp28,[v+[F(0)]*3 for v in nu]);n3=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  p=next((p for row,p in col if row==5),[]);v=C3.polynomial(p,vs)
  for mask in range(4):
   q=C3.component(C3.mul(v,{mask:F(1)}),0);n3.append([q.get(25+i,F(0))for i in range(3)])
 assert n3==[[F(c)for c in row]for row in p2['all140_native_row20_third_lower_surrogates']]
 # Source-highest responses are checked independently at C4+epsilon.
 responses=[];epsmons=[[0]*8,[1]+[0]*7]
 for which in range(2):
  E=Q.FiniteC(epsmons,[[0,int(k==which)]for k in range(2)]);v0=[{i:c for i,c in v.items()if i<4}for v in vs];z=[F(0)]*640
  for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
   p=next((p for row,p in col if row==5),[]);v=E.polynomial(p,v0)
   z += [E.mul(v,{mask:F(1)}).get(4,F(0))for mask in range(4)]
  responses.append(z)
 assert responses==[[F(c)for c in row]for row in D['higher']['native_highest_monic_responses_raw2_by780']]
 # u5 even-highest sensitivities; other core sensitivities have even part0
 # since their differentiated quadratic coefficients remain in (tauA,tauB).
 E0=Q.FiniteC(epsmons,[[0,0],[0,0]]);vs0=[{i:c for i,c in v.items()if i<4}for v in vs];u5resp=[]
 for k in range(5):
  vv=[p.copy()for p in vs0];vv[k]=Q.plus(vv[k],{4:F(1)});z=[F(0)]*640
  for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
   p=next((p for row,p in col if row==5),[]);v=E0.polynomial(p,vv)
   z += [E0.mul(v,{mask:F(1)}).get(4,F(0))for mask in range(4)]
  u5resp.append(z)
 for k in range(5,9):
  for slot in range(4):
   vv=[p.copy()for p in vs0];vv[k]=Q.plus(vv[k],{4+slot:F(1)})
   for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
    v=E0.polynomial(next((p for row,p in col if row==5),[]),vv)
    assert all(not E0.mul(v,{mask:F(1)}).get(4,F(0))for mask in range(4))
 # First-unit and homotopy coefficients, all on literal accepted sources.
 raw1=[{int(i):Q.pairs(v)for i,v in col}for col in D['X1']['all780_source_columns44_sparse_dual9']];xp=D['X1']['X0_pivot_columns'];free=[i for i in range(44)if i not in D['m0']['selected_source_columns24']];I=[free[i]for i in D['minimal']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in D['minimal']['fixed_minor_data']['inverse']]
 L=[{(h,v):F(c)for h,v,c in row}for row in D['carry']['all6_actual_effective_H1_callbacks_row84_var8_sparse']];queries=sorted(set().union(*(set(v)for v in L)));Pi=[Q.pairs(v)for v in D['quad']['exact_quotient_duals_rows_sparse']]
 def m1(row,j,var):
  s,i=divmod(row,8);assert i!=7,('LAST_M1_MUST_REMAIN_A_CONTRACTION',row)
  if var!=7 and var>i:return F(0)
  mon=tuple(sorted((i,var)));assert mon in MI,('UNAVAILABLE_CANONICAL_M1_PAIR',row,j,var,mon)
  return vg[j][s].get(MI[mon],F(0))
 def first(p,row,v):return raw1[p].get(row,{}).get(v+1,F(0))
 def origin(p,row):return raw1[p].get(row,{}).get(0,F(0))
 def mxp1(p,row,var):
  z=sum((a*first(p,j,var)for j,a in m0[row].items()),F(0))
  if row%8==7:return F(0)if var==7 else z
  return z+sum((m1(row,j,var)*origin(p,j)for j in range(44)),F(0))
 h0={}
 for h,var in queries:
  s,pid=divmod(h,28);i,j=PAIRS[pid];h0[h]=[mxp1(p,8*s+7,i)if j==7 else(mxp1(p,8*s+j,i)-mxp1(p,8*s+i,j))/2 for p in xp]
 assert [[h,[[i,str(v)]for i,v in enumerate(vals)if v]]for h,vals in []]==[] # no hidden homotopy evaluation
 assert [[h,list(map(str,v))]for h,v in sorted(h0.items())]==p2['all_needed_hXp0']
 hi=[[F(c)for c in row]for row in D['first3']['all54_actual_pivot_h0_readback']]
 last={r:[-sum((a*first(p,j,7)for j,a in m0[r].items()),F(0))for p in xp]for r in[7,15,23]}
 assert [[r,list(map(str,v))]for r,v in last.items()]==p2['all54_last_m1_Xp0_contractions']
 AF=[[F(c)for c in row]for row in D['carry']['all48_actual_ResAT_row20_B11']];AF1=[[F(c)for c in row]for row in D['first2']['all48_actual_adapter_row20_first9']]
 def pivot2(p,mon):
  if p>=640:return native[p-640].get(20,{}).get(MI[mon],F(0))
  z=F(0)
  for row,node in Bcols[p]:
   q=tv[node][0]
   if mon==(7,7):z+=AF[row][2]*q.get(0,F(0))+AF1[row][8]*q.get(MI[(7,)],F(0))
   else:
    i=mon[0];b={0:7,4:9,5:10}[i];z+=AF[row][b]*q.get(0,F(0))+AF1[row][i+1]*q.get(MI[(7,)],F(0))+AF1[row][8]*q.get(MI[(i,)],F(0))
  return z
 hp1=[[F(5,2)*pivot2(p,(i,7))+sum((mu1[j]*first(p,j,i)for j in range(44)),F(0))for p in xp]for i in[0,4,5]]
 mu2=[-sum((mu1[j]*first(p,j,7)for j in range(44)),F(0))-F(5,2)*pivot2(p,(7,7))for p in xp]
 assert hp1==[[F(c)for c in row]for row in p2['all54_hXp1_r']]
 assert mu2==list(map(F,p2['all18_mu2_Xp0_replacement_values']))
 HS={q:{}for q in queries}
 for s in range(3):
  for row in range(8):
   for mi,mon in enumerate(MON):
    if len(mon)!=2:continue
    h=homotopy_first_stencil(s,row,mon)
    for q in queries:
     if h.get(q,F(0)):HS[q][8*s+row,mi]=h[q]
 assert [[h,var,[[row,mi,str(c)]for(row,mi),c in sp.items()]]for(h,var),sp in HS.items()]==p2['all_constant_H2_stencils']
 AD,AB=rebuild_second_adjoints(D)
 assert AD==[[Q.pairs(v)for v in row]for row in D['honest']['adjoint_deltaX2_6_by36_by44_sparse']]
 assert AB==[[Q.pairs(v)for v in row]for row in D['honest']['adjoint_deltaBold2_6_by36_by48_sparse']]
 save('phase2-independent-canonical-m1-native-highest-and-pivot-source-proof.json',dict(direct_native140_C100_and_lower3=True,source_highest_two_nu_and_five_u5_responses_bound=True,unknown_core_highest_even_responses_zero_C4=True,all_actual_m0_and_mu1_readbacks=True,all_needed_hXp0_and54_hXp1=True,all18_mu2_and54_last_m1_contractions=True,all_H2_stencils_and_second_adjoints=True,no_unavailable_prefix_pair_or_last_m1_assigned_zero=True))
 # Whole-source fold BEFORE transport: no individual unknown second is used.
 AA=[]
 for row in D['X1']['new_A12_by9_C36']:
  rr=[]
  for v in row:
   z={k:F(c)for k,c in enumerate(v[:4])if F(c)}
   for var in range(8):
    if(var,)in MI:
     for k in range(4):
      c=F(v[4*(var+1)+k])
      if c:z[4*MI[(var,)]+k]=c
   rr.append(z)
  AA.append(rr)
 sources=D['free43']['all661_constant_original780_source_combinations'];original=D['operator']['all661_source_columns_in_original849'];d3=D['Koszul']['all8_differentials_sparse_columns'][2];a1saved=D['higher']['all661_a1_18_by8_sparse'];allfinal=[];alldirect=[];allcarry=[];allrow3=[];checks=0
 for sc,src in enumerate(sources):
  co=Q.pairs(src);bp=[{}for _ in range(48)];bg=[{}for _ in range(48)];pp=[{}for _ in range(8)];pg=[{}for _ in range(8)];nc={};x1=[{}for _ in range(44)]
  for j,a in co.items():
   for row,p in raw1[j].items():x1[row]=Q.plus(x1[row],{v-1:a*c for v,c in p.items()if v})
   if j<640:
    for row,node in Bcols[j]:bp[row]=Q.plus(bp[row],Q.scaled(tv[node][0],a));bg[row]=Q.plus(bg[row],Q.scaled(tv[node][1],a))
    for row,node in enumerate(Pcols[j]):
     if node:pp[row]=Q.plus(pp[row],Q.scaled(cv[rn[node][1]][0],a));pg[row]=Q.plus(pg[row],Q.scaled(cv[rn[node][1]][1],a))
   else:nc[j-640]=a
  assert all(not p for p in bg+pg),('UNCANCELED_TRUE_OPAQUE_SECOND',sc)
  assert all(not p.get(0,F(0))for p in bp),('B0_MUST_CANCEL_BEFORE_A2_DROPPED',sc)
  assert all(sum((a*v[j]for j,a in co.items()),F(0))==0 for v in responses+u5resp),('UNKNOWN_HIGHEST_COEFFICIENT_NOT_CANCELED',sc)
  bc=[{4*mi+mask:c for mask in range(4)for mi,c in bp[4*i+mask].items()}for i in range(12)];x2=[{}for _ in range(44)]
  for gen in range(9):
   v=Q.plus(*(C.mul(AA[i][gen],bc[i])for i in range(12)if bc[i]))
   for mask in range(4):x2[4*gen+mask]=C.component(v,mask)
  for i in range(8):x2[36+i]=pp[i]
  for j,a in nc.items():
   for row,p in native[j].items():x2[row]=Q.plus(x2[row],Q.scaled(p,a))
  for row in range(44):
   for var in[0,1,4,5,6,7]:assert x2[row].get(MI[(var,)],F(0))==x1[row].get(var,F(0))
  a1=Q.mprod(Ji,[x1[i]for i in I]);assert [[8*k+v,str(a)]for k,p in enumerate(a1)for v,a in sorted(p.items())]==a1saved[sc]
  for row in range(44):assert x1[row]==Q.plus(*(Q.scaled(a1[k],origin(p,row))for k,p in enumerate(xp)))
  def prod(v,w,mon):
   i,j=mon;return v.get(i,F(0))*w.get(j,F(0))+(v.get(j,F(0))*w.get(i,F(0))if i!=j else F(0))
  a2cache={}
  def a2(mon):
   if mon not in a2cache:
    z={k:x2[row].get(MI[mon],F(0))-sum((prod({v:first(p,row,v)for v in range(8)},a1[j],mon)for j,p in enumerate(xp)),F(0))for k,row in enumerate(I)}
    a2cache[mon]=Q.mvec(Ji,z)
   return a2cache[mon]
  r3=[]
  for n,i in enumerate([0,4,5]):
   v=sum((AF1[r][i+1]*bp[r].get(MI[(7,7)],F(0))+AF1[r][8]*bp[r].get(MI[(i,7)],F(0))+AF[r][2]*bp[r].get(MI[(i,)],F(0))+AF[r][[7,9,10][n]]*bp[r].get(MI[(7,)],F(0))for r in range(48)),F(0))
   r3.append(v+sum((a*n3[j][n]for j,a in nc.items()),F(0)))
  direct=[F(0)]*6
  for f in range(6):
   for index,c in Pi[f].items():
    target,m=divmod(index,36);mon=MON2[m]
    if target<2:
     row=[42,43][target];aa=a2(mon);v=x2[row].get(MI[mon],F(0))-sum((origin(p,row)*aa.get(k,F(0))+prod({i:first(p,row,i)for i in range(8)},a1[k],mon)for k,p in enumerate(xp)),F(0))
    else:
     assert mon==(7,7);n=target-2;i=[0,4,5][n];aa=a2(mon)
     v=F(5,2)*r3[n]+sum((mu1[j]*x2[j].get(MI[(i,7)],F(0))for j in range(44)),F(0))+sum((mu2[k]*a1[k].get(i,F(0))for k in range(18)),F(0))-sum((hi[n][k]*aa.get(k,F(0))+hp1[n][k]*a1[k].get(7,F(0))for k in range(18)),F(0))
    direct[f]+=c*v
  v2cache={}
  def v2(row,mi):
   key=row,mi
   if key in v2cache:return v2cache[key]
   mon=MON[mi];v=sum((a*x2[j].get(mi,F(0))for j,a in m0[row].items()),F(0))
   if row%8==7:
    if 7 in mon:v+=sum((last[row][k]*a1[k].get(mon[0],F(0))for k in range(18)),F(0))
   else:
    i,j=mon;v+=sum((m1(row,k,i)*x1[k].get(j,F(0))+(m1(row,k,j)*x1[k].get(i,F(0))if i!=j else F(0))for k in range(44)),F(0))
   v2cache[key]=v;return v
  b1={}
  for h,var in queries:
   v=sum((a*v2(row,mi)for(row,mi),a in HS[h,var].items()),F(0))-sum((h0[h][k]*a1[k].get(var,F(0))for k in range(18)),F(0))
   for j,a in original[sc]:
    if j>=681:
     for rr,e,c in d3[j-681]:
      if rr==h and sum(e)==1 and e[var]==1:v+=F(a)*F(c)
   b1[h,var]=v
  for row in[7,15,23]:assert v2(row,MI[(7,7)])==0;checks+=1
  carry=[-sum((a*b1[k]for k,a in L[f].items()),F(0))for f in range(6)];final=[a+b for a,b in zip(direct,carry)]
  assert direct==list(map(F,p3['all661_direct_six_values'][sc]))and carry==list(map(F,p3['all661_actual_carry_six_values'][sc]))
  assert final==list(map(F,p3['all661_actual_quadratic_six_values'][sc]))and r3==list(map(F,p3['all661_row20_third_functionals3'][sc]))
  assert final[2]==0 and final[3]==F(D['free43']['all661_actual_quadratic_dual_values'][sc])
  allfinal.append(final);alldirect.append(direct);allcarry.append(carry);allrow3.append(r3)
  assert sum(map(len,tv.values()))<150000 # tuples contain the two sparse components, not dense coefficients
 assert checks==1983
 cp=p3['fixed_rank_columns'];rr=p3['fixed_rank_rows'];assert cp==[280,570,586]and rr==[0,1,3]
 A=[{j:allfinal[c][r]for j,c in enumerate(cp)if allfinal[c][r]}for r in rr];inv=[Q.sparse(v)for v in p3['fixed_minor_inverse']]
 assert Q.mprod(A,inv)==Q.mprod(inv,A)==Q.eye(3)and Q.det_fixed(A)==F(p3['fixed_minor_determinant'])==F(4504862780165234375,1935581132388096)
 for col,w in zip(allfinal,p3['all661_rank_factorizations']):
  c=Q.mvec(inv,{i:col[r]for i,r in enumerate(rr)if col[r]});assert c==Q.sparse(w)
  assert [sum((allfinal[cp[k]][i]*a for k,a in c.items()),F(0))for i in range(6)]==col
 save('independent-all661-direct-carry-values-and-fixed-rank-proof.json',dict(all661_raw_Bold0_and_true_opaque_second_gradients_cancel=True,
  all_highest_monic_u5_core_sensitivities_checked_before_surrogates=True,all_actual_first_source_and18_a1_readbacks=True,
  original_ordered_H2_and15_mixed_d3_sources_retained=True,all661_direct_carry_six_values=allfinal,
  all1983_last_component_pure_r2_cycles=True,literal_zero_row2_and_independent_free43_row3=True,
  fixed_rows=rr,fixed_columns=cp,fixed_minor_determinant=F(4504862780165234375,1935581132388096),all661_fixed_rank3_factorizations=True,
  actual_rank=3,actual_HF_Tnil=[5,5,3],no_new_source_jet_solve_X3_or_nilcube=True))
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,actual_complete6x661_rank=3,actual_HF_Tnil_through2=[5,5,3],all_source_direct_carry_and_fixed_rank_gates_verified=True,no_new_source_gauge_jet_X3_or_nilcube=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
