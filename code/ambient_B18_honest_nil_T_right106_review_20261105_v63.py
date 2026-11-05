#!/usr/bin/env python3
"""Independent fixed-source RIGHT106 replay; no producer imports or solves."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
from ambient_B18_honest_nil_T_gauge120_review_20261105_v63 import Base
from ambient_B18_honest_nil_T_quadratic_remaining_review_20261105_v63 import homotopy_first_stencil
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_HONEST_NIL_T_RIGHT106_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
M2=list(itertools.combinations_with_replacement(range(8),2));M3=list(itertools.combinations_with_replacement(range(8),3));P8=list(itertools.combinations(range(8),2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 sid='AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63';vid='AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63'
 sr=load(sid);vr=load(vid);assert sr['all_passed']and vr['all_passed']
 for stem,rep in[(sid,sr),(vid,vr)]:
  for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
  for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/n)==h
 PF=load(sid,'right106-full-first-span-and-cache-footprint.json');v1=load(vid,'phase1-all35-nine-convolution-projections.json');v2=load(vid,'phase2-actual106-cubic-values-and-rank.json')
 names={
 'cub':('AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63','complete-cubic917-typed-source-and-three-NF-adjoints.json'),
 'first':('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json'),
 'first2':('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase2-row20-surrogates-and-true-m1row7.json'),
 'X':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json'),
 'W':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase3-actual-W84-rank-and-mu.json'),
 'm0':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json'),
 'carry':('AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63','actual81-pivot-first-carries-and-six-H1-source-adjoints.json'),
 'grad':('AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63','phase2-allB-and-P-true-honest-second-gradients.json'),
 'roots':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63','phase1-all-source-root-second-representatives.json'),
 'lower':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63','phase2-actual-native-and-m1-source-inputs.json'),
 'higher':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63','actual661-first-corrections-m2-replacement-and-native-highest-sensitivities.json'),
 'rad':('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json'),
 'chain':('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63','actual-UP-QB-640-typed-operator.json'),
 'C0':('AMBIENT_B18_T0_HONEST_C0_20261104_v63','honest-C0-D0-D4-DAG.json'),
 'proper':('AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63','phase3-actual-proper-G8-base25.json'),
 'physics':('AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63','phase1-common105-and-physical-source-base25.json'),
 'GL':('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63','actual-E4-nil-two-generators-and-relations.json'),
 'minimal':('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json')}
 D={k:load(*p)for k,p in names.items()}
 for p in['AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_GAUGE120_20261105_v63_REVIEW','AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63_REVIEW']:assert load(p)['all_passed']
 save('phase.json',dict(phase='fixed actual source and accepted independent cache/filtration proofs pinned'))
 X,W=D['X'],D['W'];raw=[{int(i):Q.pairs(v)for i,v in col}for col in X['all780_source_columns44_sparse_dual9']];xp=X['X0_pivot_columns'];m0=[Q.sparse(v)for v in D['m0']['m0_24_by44']]
 coeff=[list(map(F,v))for v in X['all780_X0_factorization_coefficients']];O=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[Q.sparse(v)for v in W['W_fixed_minor_inverse']];final=[i for i in range(762)if i not in wp];R35=D['first']['rank_minor_columns'];assert R35==PF['selected35_actual_first_relations']and len(R35)==35
 weights=[];x1=[]
 for rel in R35:
  k=final[rel];w=Q.pairs(W['all_W84_columns_sparse'][k]);a=Q.mvec(Wi,{i:w[r]for i,r in enumerate(wr)if r in w});src={O[k]:F(1)}
  for p,c in zip(xp,coeff[O[k]]):Q.put(src,p,-c)
  for i,c in a.items():
   Q.put(src,O[wp[i]],-c)
   for p,b in zip(xp,coeff[O[wp[i]]]):Q.put(src,p,c*b)
  assert src==Q.pairs(PF['all35_constant_original780_source_combinations'][len(weights)])
  v=[{}for _ in range(44)]
  for j,c in src.items():
   for row,z in raw[j].items():v[row]=Q.plus(v[row],Q.scaled(z,c))
  assert all(not z.get(0,F(0))for z in v)
  assert all(not sum((c*v[j].get(var+1,F(0))for j,c in row.items()),F(0))for row in m0 for var in range(8))
  weights.append(src);x1.append([[z.get(var+1,F(0))for var in range(8)]for z in v])
 cp=PF['selected20_contraction_basis_columns'];rr=PF['fixed20_minor_rows'];assert cp==list(range(18))+[53,54]
 ncol=D['m0']['selected_source_columns24'];ni=[Q.sparse(v)for v in D['m0']['selected_minor_inverse']];k42={42:F(1)}
 aa=Q.mvec(ni,{j:row.get(42,F(0))for j,row in enumerate(m0)})
 for j,c in aa.items():Q.put(k42,ncol[j],-c)
 assert k42==Q.sparse(PF['extra_kernel_k42_vector44'])
 cols=[{i:raw[p].get(i,{}).get(0,F(0))for i in range(44)if raw[p].get(i,{}).get(0,F(0))}for p in xp]+[{i:F(i==43)for i in range(44)if i==43},k42]
 A=[{j:col.get(i,F(0))for j,col in enumerate(cols)if col.get(i,F(0))}for i in rr];iv=[Q.sparse(v)for v in PF['fixed20_minor_inverse']];assert Q.det_fixed(A)==F(PF['fixed20_minor_determinant'])
 for i in range(20):
  assert Q.mvec(A,{j:iv[j].get(i,F(0))for j in range(20)if iv[j].get(i,F(0))})=={i:F(1)}
  assert Q.mvec(iv,{j:A[j].get(i,F(0))for j in range(20)if A[j].get(i,F(0))})=={i:F(1)}
 assert all(not Q.mvec(m0,col)for col in cols)
 factors=[]
 for z in x1:
  fs=[]
  for var in range(8):
   a=Q.mvec(iv,{j:z[i][var]for j,i in enumerate(rr)if z[i][var]});re={}
   for j,c in a.items():re=Q.plus(re,Q.scaled(cols[j],c))
   assert re=={i:z[i][var]for i in range(44)if z[i][var]};fs.append([a.get(j,F(0))for j in range(20)])
  factors.append(fs)
 assert factors==[[list(map(F,v))for v in row]for row in PF['all35_by8_X1_factorizations_into20_basis']]
 assert all(not f[7][19]for f in factors) # all R35 pure-r vectors lie in the known rank19 span
 bc=D['chain']['B48_by640_sparse_columns'];pc=D['rad']['actual_P8_by640_columns'];bg={i:Q.pairs(g)for i,c,g in D['grad']['all1239_B_root_origin_and_gradients']};pg={i:Q.pairs(g)for i,root,c,g in D['grad']['all723_P_root_origin_and_gradients']};tv={i:Q.pairs(v)for i,v in D['roots']['all1239_Broot_second_representatives']};cv={i:Q.pairs(v)for i,v in D['roots']['all768_C0root_second_representatives']}
 for src in weights:
  origin={};grad={}
  for j,c in src.items():
   if j>=640:continue
   for row,node in bc[j]:
    Q.put(origin,row,c*tv[node].get(0,F(0)))
    for leaf,b in bg[node].items():Q.put(grad,(row,leaf),c*b)
   for row,node in enumerate(pc[j]):
    if node:
     for leaf,b in pg[node].items():Q.put(grad,(48+row,leaf),c*b)
  assert not origin and not grad
 # Independently convolve the cubic duals, and differentiate the full81
 # target Schur carry before the original ordered-Hadamard/Euler lift.
 delta=[Q.pairs(v)for v in D['cub']['cubic_NF_dual3_by600_sparse']];mi3={m:i for i,m in enumerate(M3)};duals=[];eff={};C1=D['carry']['all81_actual_visible_C1_top2_by8_and_wedge3_r']
 for f,dd in enumerate(delta):
  for var in range(8):
   pi={36*t+i:dd[120*t+mi3[tuple(sorted(mon+(var,)))]]for t in range(5)for i,mon in enumerate(M2)if dd.get(120*t+mi3[tuple(sorted(mon+(var,)))],F(0))}
   if not pi:continue
   duals.append((f,var,pi));adj={}
   for index,c in pi.items():
    target,q=divmod(index,36);a,b=M2[q]
    for k in range(81):
     if target<2:
      Q.put(adj,(k,b),c*F(C1[k][target][a]))
      if a!=b:Q.put(adj,(k,a),c*F(C1[k][target][b]))
     else:
      assert(a,b)==(7,7);Q.put(adj,(k,7),c*F(C1[k][2][target-2]))
   z={}
   for i,h in enumerate(wr):
    for v in range(8):Q.put(z,(h,v),sum((adj.get((k,v),F(0))*Wi[k].get(i,F(0))for k in range(81)),F(0)))
   eff[f,var]=z
 assert [[f,v,[[i,str(c)]for i,c in p.items()]]for f,v,p in duals]==PF['all_nonzero_convolution_duals']
 assert [[f,v,[[h,w,str(c)]for(h,w),c in p.items()]]for(f,v),p in eff.items()]==PF['all_actual_effective_C1_H2_callbacks']
 HS={};rawpairs=set();m1pairs=set()
 for z in eff.values():
  for h,v in z:HS[h,v]={}
 for s in range(3):
  for row in range(8):
   for q,mon in enumerate(M2):
    stencil=homotopy_first_stencil(s,row,mon)
    for h,v in HS:
     if stencil.get((h,v),F(0)):
      HS[h,v][8*s+row,mon]=stencil[h,v];rawpairs.add(mon)
      if row<7:
       for a,b in[mon]+([mon[::-1]]if mon[0]!=mon[1]else[]):
        if a==7 or a<=row:
         if any(z[j][b]for z in x1 for j in range(44)):m1pairs.add(tuple(sorted((row,a))))
 for f,v,pi in duals:rawpairs.update(M2[i%36]for i in pi if i//36<2)
 assert sorted(rawpairs)==[tuple(x)for x in PF['required_raw_X2_monomials']] and sorted(m1pairs)==[tuple(x)for x in PF['required_canonical_m1_pairs']]
 MON=[tuple(m)for m in D['physics']['base25_monomials']];assert rawpairs|m1pairs<=set(MON)
 sensu=[];sensv=[]
 for syz in D['cub']['all106_first_linear_syzygies280_sparse_columns']:
  u=[F(0)]*3;v=[F(0)]*3
  for index,c in syz:
   p,var=divmod(index,8);c=F(c)
   for f in range(3):
    if var==7:
     for n,zi in enumerate([0,4,5]):u[f]+=c*factors[p][zi][19]*delta[f].get(120*(n+2)+mi3[(7,7,7)],F(0))
    for i in range(7):v[f]-=c*factors[p][i][19]*eff.get((f,var),{}).get((56+P8.index((i,7)),7),F(0))
  sensu.append(u);sensv.append(v)
 assert sensu==[list(map(F,v))for v in PF['all106_U_sensitivities3']] and sensv==[list(map(F,v))for v in PF['all106_V_sensitivities3']] and not any(any(v)for v in sensu+sensv)
 save('phase1-independent-full-source20-and-two-formal-sensitivities.json',dict(all35_actual_source_combinations=True,both20_inverse_products_and_all12320_factorizations=True,known19_plus_k42=True,all35_Bold0_opaque_second_gradients_zero=True,all9_convolution_duals_and_ordered_H2_callbacks=True,all106_by3_U_and_V_sensitivities_zero=True,all_required_source_pairs_in_accepted_base25=True))
 # Our sparse finite algebra, previously independently checked source25.
 # No canonical solve or new source polynomial is evaluated here.
 B=Base(MON);MI=B.idx;exps=[tuple(m.count(i)for i in range(8))for m in MON];nu=[list(map(F,v))for v in D['physics']['nu_A_B_base25']];C=Q.FiniteC(exps,nu);P=D['proper'];R=D['rad'];rn=R['new_Cscalar_and_R0component_nodes']
 core=[Q.sparse(P['kappa9_C100'][i])for i in[0,1,2,6]];vs=[Q.sparse(v)for v in P['actualG8_C100'][:5]]+[C.polynomial(D['GL']['forward_GL4'][k],core)for k in['h','c','r','l']]
 native=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  vv=[(row,C.polynomial(pol,vs))for row,pol in col]
  for mask in range(4):
   z={}
   for row,pol in vv:
    pol=C.mul(pol,{mask:F(1)})
    for k in range(4):
     p=C.component(pol,k)
     if p:z[4*row+k]=p
   native.append(z)
 qv=[];dd=C.polynomial(R['d_polynomial'],vs)
 for pol in R['R9_generators_polynomials9']:
  pp=C.mul(dd,C.polynomial(pol,vs));qv.extend(C.mul(pp,{mask:F(1)})for mask in range(4))
 for i in range(8):qv.append({4*MI[(i,)]+1:F(-1)}if i<7 and(i,)in MI else{})
 vg=[[C.component(v,1),C.component(v,3),Q.plus(C.component(v,2),Q.scaled(C.component(v,1),F(-11,12)))]for v in qv]
 mu1=list(map(F,D['first2']['true_actual_m1_row7_r_coefficients44']));last15=[v[1].get(MI[(7,7)],F(0))for v in vg]
 assert [v[0].get(MI[(7,7)],F(0))for v in vg]==mu1
 for j,v in enumerate(vg):
  for s in range(3):
   for i in range(7):
    if(i,)in MI:assert v[s].get(MI[(i,)],F(0))==m0[8*s+i].get(j,F(0))
   assert v[s].get(MI[(7,)]if s<2 else MI[(7,7)],F(0))==m0[8*s+7].get(j,F(0))
 # Unknown highest native coefficients remain independent directions;
 # their response is tested on all35 original source combinations.
 eps=[[0]*8,[1]+[0]*7];vs0=[{i:c for i,c in v.items()if i<4}for v in vs];responses=[]
 for typ in range(7):
  E=Q.FiniteC(eps,[[0,int(k==typ)]for k in range(2)]if typ<2 else[[0,0],[0,0]]);ev=[v.copy()for v in vs0]
  if typ>=2:ev[typ-2]=Q.plus(ev[typ-2],{4:F(1)})
  vals=[F(0)]*640
  for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
   v=E.polynomial(next((p for row,p in col if row==5),[]),ev);vals.extend(E.mul(v,{mask:F(1)}).get(4,F(0))for mask in range(4))
  responses.append(vals)
  assert all(not sum((c*vals[j]for j,c in src.items()),F(0))for src in weights)
 assert responses[:2]==[list(map(F,v))for v in D['higher']['native_highest_monic_responses_raw2_by780']]
 AA=[]
 for row in X['new_A12_by9_C36']:
  rr1=[]
  for v in row:
   p={i:F(v[i])for i in range(4)if F(v[i])}
   for var in range(8):
    if(var,)in MI:
     for k in range(4):Q.put(p,4*MI[(var,)]+k,F(v[4*(var+1)+k]))
   rr1.append(p)
  AA.append(rr1)
 AF=[list(map(F,v))for v in D['carry']['all48_actual_ResAT_row20_B11']];AF1=[list(map(F,v))for v in D['first2']['all48_actual_adapter_row20_first9']]
 def origin(p,row):return raw[p].get(row,{}).get(0,F(0))
 def first(p,row,var):return raw[p].get(row,{}).get(var+1,F(0))
 def pivot2(p,mon):
  if p>=640:return native[p-640].get(20,{}).get(MI[mon],F(0))
  z=F(0)
  for row,node in bc[p]:
   qq=tv[node]
   if mon==(7,7):z+=AF[row][2]*qq.get(0,F(0))+AF1[row][8]*qq.get(MI[(7,)],F(0))
   else:
    i=mon[0];z+=AF[row][{0:7,4:9,5:10}[i]]*qq.get(0,F(0))+AF1[row][i+1]*qq.get(MI[(7,)],F(0))+AF1[row][8]*qq.get(MI[(i,)],F(0))
  return z
 mutwo=[-sum((mu1[j]*first(p,j,7)for j in range(44)),F(0))-F(5,2)*pivot2(p,(7,7))for p in xp]
 hxp1=[[F(5,2)*pivot2(p,(i,7))+sum((mu1[j]*first(p,j,i)for j in range(44)),F(0))for p in xp]for i in[0,4,5]]
 assert mutwo==list(map(F,D['lower']['all18_mu2_Xp0_replacement_values'])) and hxp1==[list(map(F,v))for v in D['lower']['all54_hXp1_r']]
 g3=sum(map(F,D['C0']['cached_actual_CA_gA_all_degree3_coefficients']));assert g3==F(125,24)
 mu20=mutwo+[-g3,F(0)];l23=[-sum((m0[23].get(j,F(0))*first(p,j,7)for j in range(44)),F(0))for p in xp]+[F(11,12)*g3,F(0)]
 assert l23[:18]==list(map(F,dict(D['lower']['all54_last_m1_Xp0_contractions'])[23]))
 def dott(a,p,var):return sum((c*d for c,d in zip(a,factors[p][var])),F(0))
 def m1(row,j,var):
  s,i=divmod(row,8);assert i<7
  if var!=7 and var>i:return F(0)
  mon=tuple(sorted((i,var)));assert mon in MI,('UNCACHED_PAIR',row,var)
  return vg[j][s].get(MI[mon],F(0))
 def mxp1(p,row,var):
  z=sum((a*first(p,j,var)for j,a in m0[row].items()),F(0))
  if row%8==7:return F(0)if var==7 else z
  return z+sum((m1(row,j,var)*origin(p,j)for j in range(44)),F(0))
 h0={}
 for h,v in HS:
  s,pair=divmod(h,28);i,j=P8[pair];h0[h]=[mxp1(p,8*s+7,i)if j==7 else(mxp1(p,8*s+j,i)-mxp1(p,8*s+i,j))/2 for p in xp]
 hi=[list(map(F,v))for v in D['first']['all54_actual_pivot_h0_readback']];free=[i for i in range(44)if i not in ncol];I=[free[i]for i in D['minimal']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in D['minimal']['fixed_minor_data']['inverse']];n3=[list(map(F,v))for v in D['lower']['all140_native_row20_third_lower_surrogates']]
 projected=[];savedinputs=[]
 def convolution(a,b,mon):
  i,j=mon;return a[i]*b[j]+(a[j]*b[i]if i!=j else F(0))
 for p,src in enumerate(weights):
  bp=[{}for _ in range(48)];pp=[{}for _ in range(8)];nc=[]
  for j,c in src.items():
   if j<640:
    for row,node in bc[j]:bp[row]=Q.plus(bp[row],Q.scaled(tv[node],c))
    for row,node in enumerate(pc[j]):
     if node:pp[row]=Q.plus(pp[row],Q.scaled(cv[rn[node][1]],c))
   else:nc.append((j-640,c))
  assert all(not z.get(0,F(0))for z in bp)
  bC=[{4*mi+mask:c for mask in range(4)for mi,c in bp[4*j+mask].items()}for j in range(12)]
  x2=[{}for _ in range(44)]
  for g in range(9):
   z={}
   for j in range(12):z=Q.plus(z,C.mul(AA[j][g],bC[j]))
   for mask in range(4):x2[4*g+mask]=C.component(z,mask)
  for row in range(8):x2[36+row]=pp[row]
  for j,c in nc:
   for row,v in native[j].items():x2[row]=Q.plus(x2[row],Q.scaled(v,c))
  for row in range(44):
   for var in range(8):
    if(var,)in MI:assert x2[row].get(MI[(var,)],F(0))==x1[p][row][var]
  a1=[{v:c for v in range(8)if(c:=sum((a*x1[p][I[k]][v]for k,a in row.items()),F(0)))}for row in Ji]
  a2={}
  for mon in rawpairs|{(7,7)}:
   rhs={j:x2[row].get(MI[mon],F(0))-sum((convolution([first(z,row,v)for v in range(8)],[a1[k].get(v,F(0))for v in range(8)],mon)for k,z in enumerate(xp)),F(0))for j,row in enumerate(I)}
   a2[mon]=Q.mvec(Ji,rhs)
  x3=[]
  for n,i in enumerate([0,4,5]):
   z=sum((AF1[row][i+1]*bp[row].get(MI[(7,7)],F(0))+AF1[row][8]*bp[row].get(MI[(i,7)],F(0))+AF[row][2]*bp[row].get(MI[(i,)],F(0))+AF[row][[7,9,10][n]]*bp[row].get(MI[(7,)],F(0))for row in range(48)),F(0));x3.append(z+sum((c*n3[j][n]for j,c in nc),F(0)))
  vv={}
  for stens in HS.values():
   for row,mon in stens:
    if(row,mon)in vv:continue
    z=sum((a*x2[j].get(MI[mon],F(0))for j,a in m0[row].items()),F(0));s,i=divmod(row,8)
    if i==7:
     if 7 in mon:
      other=mon[0]
      z+=sum((mu1[j]*x1[p][j][other]for j in range(44)),F(0))if row==7 else sum((last15[j]*x1[p][j][other]for j in range(44)),F(0))if row==15 else dott(l23,p,other)
    else:
     i,j=mon;z+=sum((m1(row,k,i)*x1[p][k][j]+(m1(row,k,j)*x1[p][k][i]if i!=j else F(0))for k in range(44)),F(0))
    vv[row,mon]=z
  bh={(h,v):sum((c*vv[key]for key,c in sten.items()),F(0))-sum((h0[h][k]*a1[k].get(v,F(0))for k in range(18)),F(0))for(h,v),sten in HS.items()}
  values=[]
  for f,var,pi in duals:
   val=F(0)
   for index,c in pi.items():
    target,q=divmod(index,36);mon=M2[q]
    if target<2:
     row=[42,43][target];v=x2[row].get(MI[mon],F(0))-sum((origin(z,row)*a2[mon].get(k,F(0))+convolution([first(z,row,w)for w in range(8)],[a1[k].get(w,F(0))for w in range(8)],mon)for k,z in enumerate(xp)),F(0))
    else:
     assert mon==(7,7);n=target-2;i=[0,4,5][n];v=F(5,2)*x3[n]+sum((mu1[j]*x2[j].get(MI[(i,7)],F(0))for j in range(44)),F(0))+dott(mu20,p,i)-sum((hi[n][k]*a2[(7,7)].get(k,F(0))+hxp1[n][k]*a1[k].get(7,F(0))for k in range(18)),F(0))
    val+=c*v
   val-=sum((c*bh[key]for key,c in eff[f,var].items()),F(0));values.append(val)
  projected.append(values);savedinputs.append(dict(source_relation=R35[p],actual_scalar_X3_zir2=x3,mu2_X1_surrogate8=[dott(mu20,p,v)for v in range(8)],last23_X1_surrogate8=[dott(l23,p,v)for v in range(8)]))
 assert projected==[list(map(F,v))for v in v1['all35_projected_M2_lower_representatives']]
 assert enc(savedinputs)==v1['all35_source_input']
 labels={(f,var):j for j,(f,var,pi)in enumerate(duals)};actual=[]
 for syz in D['cub']['all106_first_linear_syzygies280_sparse_columns']:
  z=[F(0)]*3
  for index,c in syz:
   p,var=divmod(index,8)
   for f in range(3):
    if(f,var)in labels:z[f]+=F(c)*projected[p][labels[f,var]]
  actual.append(z)
 assert actual==[list(map(F,v))for v in v2['all106_actual_cubic_columns3']] and all(not any(v)for v in actual) and v2['actual_right106_rank']==0
 save('phase2-independent-cache35-contractions-and-all106-values.json',dict(all35_nine_lower_projections=projected,all35_source_input=savedinputs,all106_actual_cubic_columns3=actual,all7_highest_native_response_source_directions_cancel=True,all18_mu2_and54_hXp1_recomputed=True,actual_gA3=g3,two_formal_unknowns_cancel_before_quotient_values=True,actual_right106_rank=0,allzero_is_only_partial_negative=True,no_actual_full_tail_or_nilcube_decided=True))
 for p in[Path(__file__),Path(Q.__file__),ROOT/'code/ambient_B18_honest_nil_T_gauge120_review_20261105_v63.py',ROOT/'code/ambient_B18_honest_nil_T_quadratic_remaining_review_20261105_v63.py',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,actual_right106_rank=0,all106_actual_values_zero=True,full_source20_two_unknown_and_cache_value_proof=True,allzero_partial_negative_only=True,no_new_source_jet_or_higher_callback=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20;save('report.json',rep);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
