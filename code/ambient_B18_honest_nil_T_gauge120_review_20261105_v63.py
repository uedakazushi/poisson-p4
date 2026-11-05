#!/usr/bin/env python3
"""Independent PROVIDED120 canonical-gauge verifier (prepared before FINAL).
No implicit source solve, canonical coefficient solve, X3 or source3.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
from ambient_B18_T_native_S2_B11_review_20261104_v63 import B44,V
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_HONEST_NIL_T_GAUGE120_20261105_v63_REVIEW'
OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
class Base:
 def __init__(self,mon):self.mon=[tuple(m)for m in mon];self.idx={m:i for i,m in enumerate(self.mon)}
 def mul(self,a,b):
  z={}
  for i,c in a.items():
   for j,d in b.items():
    k=self.idx.get(tuple(sorted(self.mon[i]+self.mon[j])))
    if k is not None:Q.put(z,k,c*d)
  return z
 def inv(self,p):
  c=p[0];v=Q.plus(Q.scaled(p,-1/c),{0:F(1)})
  ans=Q.scaled(Q.plus({0:F(1)},v,self.mul(v,v)),1/c)
  assert self.mul(p,ans)=={0:F(1)}
  return ans
class Cover(B44):
 def __init__(self,base,nu):super().__init__(nu);self.b=base
 def mul(self,a,b):
  z={}
  for(i,u),c in a.items():
   for(j,v),d in b.items():
    p=self.b.mul({i:c},{j:d})
    if u&v&1:p=self.b.mul(p,self.nus[0])
    if u&v&2:p=self.b.mul(p,self.nus[1])
    z=self.add(z,{(k,u^v):x for k,x in p.items()})
  return z
def source_inputs(D):
 """Recompose all105 leaves and physical coefficients in the25-base ring."""
 mon=D['preflight']['base25_original8_monomials'];B=Base(mon)
 assert len(B.mon)==25 and all(tuple(m[:j]+m[j+1:])in B.idx for m in B.mon for j in range(len(m)))
 assert B.mon[:20]==[tuple(m)for m in D['source20']['base20_monomials_original8']]
 J=D['source2'];own={f:dict(J['both57_source2'][f]['all57_source2jets'])for f in['CA','CB']}
 def fold(v):
  z={}
  for m,c in zip(J['source2_monomials'],v):
   key=tuple(sorted(7 if i==8 else i for i in m));k=B.idx.get(key)
   if k is not None and F(c):Q.put(z,k,F(c))
  return z
 Y=[]
 for lin,quad in zip(D['adapter2']['YB_in_YA_linear'],D['adapter2']['YB_in_YA_quadratic']):
  z={}
  for i,c in enumerate(lin):
   k=B.idx.get((7 if i==8 else i,))
   if k is not None and F(c):Q.put(z,k,F(c))
  for e,c in quad:
   m=tuple(sorted(7 if i==8 else i for i,n in enumerate(e)for _ in range(n)));k=B.idx.get(m)
   if k is not None and F(c):Q.put(z,k,F(c))
  Y.append(z)
 def cb(v):
  z={}
  for mon,c in zip(J['source2_monomials'],v):
   if not F(c):continue
   p={0:F(c)}
   for i in mon:p=B.mul(p,Y[i])
   z=Q.plus(z,p)
  return z
 leaves={i:fold(v)for i,v in own['CA'].items()}
 for i,v in own['CB'].items():
  z=cb(v)
  if i in leaves:assert z==leaves[i]
  else:leaves[i]=z
 assert len(leaves)==105
 old={i:Q.pairs(v)for i,v in D['source20']['all105_actual_source2_values']}
 assert all({k:c for k,c in v.items()if k<20}==old[i]for i,v in leaves.items())
 src=D['source'];nodes=src['DAG_scalar_nodes'];roots=[i for f in src['families'].values()for pair in f['complete_normalized_physical_w_coefficients'][0]for i in pair]+[f['actual_tracefree_nu']for f in src['families'].values()]+sum((J['both57_source2'][f]['equations57']for f in['CA','CB']),[])
 need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i)
  if i in leaves:continue
  op=nodes[i][0]
  if op in['a','m','n','i']:todo.extend(nodes[i][1:])
  else:assert op=='c'
 uses=collections.Counter(roots)
 for i in need:
  if i not in leaves and nodes[i][0]in['a','m','n','i']:uses.update(nodes[i][1:])
 vals={};peak=0
 for i in sorted(need):
  n=nodes[i];op=n[0]
  if i in leaves:v=leaves[i]
  elif op=='c':v={0:F(n[1],n[2])}if n[1]else{}
  elif op=='a':v=Q.plus(vals[n[1]],vals[n[2]])
  elif op=='m':v=B.mul(vals[n[1]],vals[n[2]])
  elif op=='n':v=Q.scaled(vals[n[1]],-1)
  elif op=='i':v=B.inv(vals[n[1]])
  else:raise AssertionError(op)
  vals[i]=v;peak=max(peak,sum(map(len,vals.values())));assert peak+86820<150000
  if i not in leaves and op in['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j]
 pairs={f:[[vals[a],vals[b]]for a,b in d['complete_normalized_physical_w_coefficients'][0]]for f,d in src['families'].items()}
 nu=[vals[d['actual_tracefree_nu']]for d in src['families'].values()]
 YA=[{B.idx[(i,)]:F(1)}if(i,)in B.idx else{}for i in range(7)]+[{B.idx[(7,)]:F(1)}]*2
 for f,args in[('CA',YA),('CB',Y)]:
  eq=[vals[i]for i in J['both57_source2'][f]['equations57']];assert eq[:48]==[{}]*48 and eq[48:]==args
 return B,leaves,pairs,nu,peak
def canonical_context(D,B,leaves,pairs,nu):
 ring=Cover(B,nu);ring.install();I=D['interface'];P=D['proper'];GL=D['old_preflight']['GL18_source_parameter_positions'];src=D['source']
 def frame(family):
  data=src['families'][family];t=ring.tau(0 if family=='CA'else 1);p=ring.add(t,ring.scale(ring.base(leaves[data['p_quadratic_relation']['s']]),F(1,2)));M=V.eye(5)
  for name,(a,b)in data['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);M[i][j]=ring.add(M[i][j],ring.base(leaves[a]),ring.mul(ring.base(leaves[b]),p))
  C=V.eye(5)
  if family=='CA':C[0][3]=ring.scale(p,-2);C[1][0]=p;C[1][3]=ring.neg(ring.mul(p,p))
  else:C[0][3]=ring.scale(p,F(1,2));C[1][0]=p;C[1][3]=ring.scale(ring.mul(p,p),F(1,4))
  return V.matmul(M,C)
 def physical(f):return[ring.add(ring.base(a),ring.mul(ring.base(b),ring.tau(0 if f=='CA'else 1)))for a,b in pairs[f]]
 HA,HB=frame('CA'),frame('CB');PA,PB=physical('CA'),physical('CB');center=list(map(F,P['actual_CA_source19_center']));rows=P['zero_coordinate_rows36']
 def own(moved):return[r[0]for r in V.matmul(V.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in P['CA_constant_left_inverse_rows19']])]
 def gauge(ev,od):
  G=V.eye(5)
  for(i,j),a,b in zip(I['actual_group15_positions'],ev,od):G[i][j]=ring.add(G[i][j],ring.base(a),ring.mul(ring.base(b),ring.tau(0)))
  return G
 def constraints(G):
  Fm=V.normalized_flag_frame(V.matmul(G,HA));co=own(V.tensor_push(PA,V.matmul(Fm,G)))
  return[Fm[1][0],Fm[1][3],Fm[1][4],Fm[0][3],Fm[0][4]]+[ring.add(co[j],ring.c(-center[j]))for j in P['CA_fiber10_fixed_coefficient_indices']]
 return ring,HA,HB,PA,PB,own,gauge,constraints
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 primary='AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63';r=load(primary);assert pins[f'reports/{primary}/report.json']=='b4cff89d643c79c4ca4bd2f8c094c2590d66d4ccba551b74ceedea91bd3b64d7'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 names={
 'preflight':('AMBIENT_B18_HONEST_NIL_T_GAUGE120_SOURCE_PREFLIGHT_20261105_v63','four-pair-base25-and-last-m1-contraction-preflight.json'),
 'source':('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json'),
 'source2':('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63','full9-firstjet-certificate.json'),
 'adapter2':('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63','actual-special-fiber20-coefficients.json'),
 'source20':('AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63','phase1-actual105-source2-and-copy.json'),
 'interface':('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63','actual-source-physical-firstjet-interface.json'),
 'proper':('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63','literal-proper-source14-preflight.json'),
 'second5':('AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63','actual4variable7core-polynomials.json'),
 'preparation':('AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63','preparation.json'),
 'old_preflight':('AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63','actual-native128-first-input-and-mixed-index-profile.json'),
 'gauge1':('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase3-gauge30-firsts.json'),
 'gauge20':('AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63','phase2-canonical270-zz-gauge.json'),
 'proper20':('AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63','phase3-actual-proper-G8-base20.json'),
 'J30':('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase2-fixed30-Jacobian.json'),
 'b11':('AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63','actual-B11-source-and-gauge-preflight.json'),
 'oldgauge':('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase2-gauge30-B11-correction.json'),
 'oldproper':('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase3-actual-proper-G8-B11.json')}
 D={k:load(*p)for k,p in names.items()}
 assert load('AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63_REVIEW')['all_passed']
 p1=load(primary,'phase1-common105-and-physical-source-base25.json');p2=load(primary,'phase2-canonical120-four-pair-gauge.json');p3=load(primary,'phase3-actual-proper-G8-base25.json')
 save('phase.json',dict(phase='immutable PROVIDED120 and source/gauge baseline peers pinned; no solve'))
 B,leaves,pairs,nu,peak=source_inputs(D)
 assert [list(m)for m in B.mon]==p1['base25_monomials']
 assert [[i,[v.get(k,F(0))for k in range(25)]]for i,v in sorted(leaves.items())]==[[i,list(map(F,v))]for i,v in p1['all105_actual_source2_values']]
 assert {f:[[[p.get(i,F(0))for i in range(25)]for p in pair]for pair in v]for f,v in pairs.items()}=={f:[[[F(c)for c in p]for p in pair]for pair in v]for f,v in p1['physical150_tracepairs_base25'].items()}
 assert [[p.get(i,F(0))for i in range(25)]for p in nu]==[[F(c)for c in p]for p in p1['nu_A_B_base25']]
 save('phase1-independent105-and-physical25-proof.json',dict(all105_source2_composed_in_physical_frame=True,all114_source_equations=True,all150_tracepairs_and_two_nu25_checked=True,source_peak_sparse_QQ=peak))
 ring,HA,HB,PA,PB,own,gauge,constraints=canonical_context(D,B,leaves,pairs,nu)
 ev0=[Q.sparse(v)for v in D['gauge20']['actual_even15_base20']];od0=[Q.sparse(v)for v in D['gauge20']['actual_odd15_base20']]
 for i in range(15):
  for dest,src in[(ev0,D['gauge1']['even15_first9']),(od0,D['gauge1']['odd15_first9'])]:
   c=F(src[i][7])
   if c:dest[i][20]=c
 residual=constraints(gauge(ev0,od0));newq=[21,22,23,24]
 assert all(k in newq and mask<2 for p in residual for k,mask in p)
 rhs=[[p.get((i,mask),F(0))for i in newq]for mask in range(2)for p in residual];assert rhs==[[F(c)for c in row]for row in p2['RHS30_by4']]
 ev=[Q.sparse(v)for v in p2['actual_even15_base25']];od=[Q.sparse(v)for v in p2['actual_odd15_base25']];corr=[[F(c)for c in row]for row in p2['correction30_by4']]
 for j in range(15):
  for k in range(25):
   assert ev[j].get(k,F(0))-ev0[j].get(k,F(0))==(corr[j][newq.index(k)]if k in newq else F(0))
   assert od[j].get(k,F(0))-od0[j].get(k,F(0))==(corr[15+j][newq.index(k)]if k in newq else F(0))
 J0=[[F(c)for c in row]for row in D['J30']['A0']];J1=[[F(c)for c in row]for row in D['J30']['A1_eta_delta']]
 J=[[F(0)]*30 for _ in range(30)]
 for i in range(15):
  for j in range(15):J[i][j]=J0[i][j];J[15+i][j]=J1[i][j];J[15+i][15+j]=J0[i][j]
 Ji=[Q.sparse(v)for v in D['J30']['J30_inverse']];Js=[Q.sparse(v)for v in J];assert Q.mprod(Js,Ji)==Q.mprod(Ji,Js)==Q.eye(30)
 assert all(sum((J[i][j]*corr[j][k]for j in range(30)),F(0))==-rhs[i][k]for i in range(30)for k in range(4))
 G=gauge(ev,od);assert not any(constraints(G));assert V.normalized_flag_frame(V.matmul(G,HA))==V.eye(5)
 FB=V.normalized_flag_frame(V.matmul(G,HB));co=own(V.tensor_push(PA,G));P=D['proper'];center=list(map(F,P['actual_CA_source19_center']))
 assert all(co[j]==ring.c(center[j])for j in P['CA_fiber10_fixed_coefficient_indices'])
 actual=[ring.add(co[j],ring.c(-center[j]))for j in P['CA_free_fiber9_coefficient_indices']]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j]for j in D['interface']['retained_source9_indices']]
 F12=[ring.add(PA[j],ring.neg(PB[j]))for j in[18,19,28,29,60,67,72,63,64,68,73,74]]
 du=ring.add(ring.ONE,ring.neg(ring.mul(kappa[1],kappa[1])));pu=ring.add(ring.c(5),ring.mul(kappa[1],kappa[1]));core=[kappa[j]for j in[0,1,2,6]];K=D['second5'];roots=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];roots.append(ring.tm(V.polynomial(p['numerator'],core),ring.power(du,-p['D_power']),ring.power(pu,-p['P_power'])))
 axes=[K['source9_variable_order'].index(n)for n in K['selected5_pivot_order']]
 gs=[ring.add(kappa[j],ring.neg(p))for j,p in zip(axes,roots)]+[V.polynomial(p,core)for p in D['preparation']['generators3']]
 oldidx=[B.idx[tuple(i for i,n in enumerate(e)for _ in range(n))]for e in D['b11']['B11_basis_monomials8']]
 for name,values in[('actual14',actual),('kappa9',kappa),('physicalF12',F12),('actualG8',gs)]:
  assert [[p.get((i,t),F(0))for i in range(25)for t in range(4)]for p in values]==[[F(c)for c in row]for row in p3[name+'_C100']]
  assert [[p.get((i,t),F(0))for i in range(20)for t in range(4)]for p in values]==[[F(c)for c in row]for row in D['proper20'][name+'_C80']]
  assert [[p.get((i,t),F(0))for i in oldidx for t in range(4)]for p in values]==[[F(c)for c in row]for row in D['oldproper'][name+'_C_B44']]
 for values,key in[(ev,'even15_B11'),(od,'odd15_B11')]:assert [[p.get(i,F(0))for i in oldidx]for p in values]==[[F(c)for c in row]for row in D['oldgauge'][key]]
 proof=dict(independent_own57_full_source2_and_CB_F9_transport=True,all114_source_equations_all105_values=True,all150_physical_tracepairs_two_nu=True,
  provided120_correction_equations_checked_without_solving=True,fixedJ30_both_inverses=True,all30_full_nonlinear_canonical_constraints_zero=True,
  actual_auxiliary_etaA_trace_frame_and_nu_carry=True,all_actual14_kappa_G8_F12_C100_checked=True,all_base20_prior270_B11_firstz6_overlaps=True,
  uniqueness='The new degree2 ideal annihilates the base maximal ideal; the trace-origin J30 is invertible, so the provided four corrections are unique.',
  no_new_source_gauge_solve_X3_or_source3=True,remaining4x661_NOT_reviewed_in_this_packet=True)
 save('independent-provided120-and-actual-C100-proof.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,provided120_independently_verified_without_solving=True,all30_constraints_114_source_equations_and_prior_overlaps=True,remaining4x661_values_NOT_yet_reviewed=True,no_new_source_jet_or_X3=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
