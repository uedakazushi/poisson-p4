#!/usr/bin/env python3
"""Independent selected135 Jacobian/RHS/native/adjoint postverification."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time
from ambient_B18_T_native_S2_B11_review_20261104_v63 import B44,V
from ambient_B18_T_gdiff_second_univariate5_review_20261104_v63 import Symbolic,closure,rec,put
from ambient_B18_T_source4_adjoint9_review_20261104_v63 import forward
from ambient_B18_T_native_zz_gauge_cached_review_20261104_v63 import PATHS as P20
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63';STEM=A+'_REVIEW'
PATHS={k:P20[k] for k in ['source','I','P','K','prep','PF','cache','gauge','proper','Op']}
PATHS.update({
 'primary':f'reports/{A}/report.json','J':f'reports/{A}/phase1-actual-J135-blocks-and-target-Jacobian.json',
 'rhs':f'reports/{A}/phase2-known-third-RHS-and-proper-response.json','factors':f'reports/{A}/phase3-native-requested-functional-factors.json',
 'values':f'reports/{A}/phase4-selected-joint135-third-adjoint-values.json','ring':f'reports/{A}/phase2-whole-live-gate.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'N0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
 'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
 'oldhigh':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/phase4-all619-high-old-source-sensitivities.json',
 'oldhigh_peer':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63_REVIEW/report.json',
 'method':'reports/AMBIENT_B18_JOINT135_SELECTED_FUNCTIONAL_METHOD_PREFLIGHT_20261104_v63/report.json',
 'MS_gate':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63/report.json'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,n,x):
 (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
def sp(p):return {int(i):Q(c) for i,c in p if Q(c)}
def enc(p):return [[i,str(c)] for i,c in sorted(p.items()) if c]
def adds(a,b,c=Q(1)):
 z=dict(a)
 for i,v in b.items():put(z,i,c*v)
 return z
def rowmul(v,rows):
 out={}
 for i,c in v.items():out=adds(out,rows[i],c)
 return out
class Cover(B44):
 def __init__(self,mons,nus):
  super().__init__(nus);self.mons=mons;idx={m:i for i,m in enumerate(mons)}
  self.btable={(i,j):idx.get(tuple(sorted(a+b))) for i,a in enumerate(mons) for j,b in enumerate(mons)}
 def base_mul(self,a,b):
  out={}
  for i,c in a.items():
   for j,d in b.items():
    k=self.btable[i,j]
    if k is not None:put(out,k,c*d)
  return out
 def mul(self,a,b):
  out={}
  for (i,u),c in a.items():
   for (j,v),d in b.items():
    z=self.base_mul({i:c},{j:d});carry=u&v;mask=u^v
    if carry&1:z=self.base_mul(z,self.nus[0])
    if carry&2:z=self.base_mul(z,self.nus[1])
    for k,q in z.items():put(out,(k,mask),q)
  return out
def proper(ring,D,phi,pairs,ev,od):
 """Literal independent physical tensor/flags; no author Proper evaluator."""
 ring.install();S,I,P,K,prep=[D[k] for k in ['source','I','P','K','prep']]
 GL=D['PF']['GL18_source_parameter_positions'];rows=P['zero_coordinate_rows36'];center=list(map(Q,P['actual_CA_source19_center']))
 def frame(family):
  f=S['families'][family];p=ring.add(ring.tau(0 if family=='CA' else 1),ring.scale(ring.base(phi[f['p_quadratic_relation']['s']]),Q(1,2)));M=V.eye(5)
  for name,(a,b) in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma') and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);M[i][j]=ring.add(M[i][j],ring.base(phi[a]),ring.mul(ring.base(phi[b]),p))
  C=V.eye(5)
  if family=='CA':C[0][3]=ring.scale(p,-2);C[1][0]=p;C[1][3]=ring.neg(ring.mul(p,p))
  else:C[0][3]=ring.scale(p,Q(1,2));C[1][0]=p;C[1][3]=ring.scale(ring.mul(p,p),Q(1,4))
  return V.matmul(M,C)
 HA,HB=frame('CA'),frame('CB');G=V.eye(5)
 for (i,j),x,y in zip(I['actual_group15_positions'],ev,od):G[i][j]=ring.add(G[i][j],ring.base(x),ring.mul(ring.base(y),ring.tau(0)))
 PA=[ring.add(ring.base(x),ring.mul(ring.base(y),ring.tau(0))) for x,y in pairs['CA']]
 FA=V.normalized_flag_frame(V.matmul(G,HA));FB=V.normalized_flag_frame(V.matmul(G,HB))
 def own(T):
  p=V.tensor_push(PA,T)
  return [v[0] for v in V.matmul(V.constmat(P['CA_left_inverse19']),[[p[rows[j]]] for j in P['CA_constant_left_inverse_rows19']])]
 ow=own(G);owc=own(V.matmul(FA,G))
 constraints=[FA[1][0],FA[1][3],FA[1][4],FA[0][3],FA[0][4]]+[ring.add(owc[j],ring.c(-center[j])) for j in P['CA_fiber10_fixed_coefficient_indices']]
 actual=[ring.add(ow[j],ring.c(-center[j])) for j in P['CA_free_fiber9_coefficient_indices']]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]
 kappa=[actual[j] for j in I['retained_source9_indices']];core=[kappa[j] for j in [0,1,2,6]]
 du=ring.add(ring.ONE,ring.neg(ring.mul(kappa[1],kappa[1])));pu=ring.add(ring.c(5),ring.mul(kappa[1],kappa[1]));rr=[]
 axes=[K['source9_variable_order'].index(n) for n in K['selected5_pivot_order']]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];rr.append(ring.tm(V.polynomial(p['numerator'],core),ring.power(du,-p['D_power']),ring.power(pu,-p['P_power'])))
 gs=[ring.add(kappa[j],ring.neg(x)) for j,x in zip(axes,rr)]+[V.polynomial(p,core) for p in prep['generators3']]
 return constraints,gs
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()};dump(out,'phase.json',{'phase':'INDEPENDENT_135_POSTVERIFY_SOURCE_PINNED','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['primary','oldhigh_peer','method','MS_gate'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 assert D['primary']['helper_code_sha256']==sha(ROOT/'code/ambient_B18_T_joint135_functional_tools_20261104_v63.py')
 S=D['source'];nodes=S['DAG_scalar_nodes'];unknown=D['joint']['all105_unknown_original_node_order'];eq=D['joint']['all105_equation_original_node_order']
 C20={i:sp(p) for i,p in D['cache']['all105_actual_source2_values']};phi0={i:p.get(0,Q(0)) for i,p in C20.items()};assert set(phi0)==set(unknown)
 pairnodes={f:p['complete_normalized_physical_w_coefficients'][0] for f,p in S['families'].items()};nunodes=[p['actual_tracefree_nu'] for p in S['families'].values()]
 roots=list(dict.fromkeys([i for pairs in pairnodes.values() for pair in pairs for i in pair]+nunodes))
 ad,sourceclosure,sourcepeak,_=forward(nodes,roots,{i:(p,{('physical',i):Q(1)}) for i,p in phi0.items()},{})
 ev0=[Q(p[0]) for p in D['gauge']['actual_even15_base20']];od0=[Q(p[0]) for p in D['gauge']['actual_odd15_base20']]
 Hg=[{} for _ in range(30)];Jg=[{} for _ in range(30)];Dt=[{} for _ in range(34)]
 Gbase=[list(map(Q,p[:4])) for p in D['proper']['actualG8_C80']]
 for direction in range(135):
  leaf=unknown[direction] if direction<105 else None
  gradient=lambda n:ad[n][1].get(('physical',leaf),Q(0)) if leaf is not None else Q(0)
  nu=[{1:gradient(n)} if gradient(n) else {} for n in nunodes];ring=Cover([(),(0,)],nu)
  phi={i:{j:c for j,c in [(0,phi0[i]),(1,Q(i==leaf))] if c} for i in unknown}
  pairs={f:[[{j:c for j,c in [(0,ad[x][0]),(1,gradient(x))] if c},{j:c for j,c in [(0,ad[y][0]),(1,gradient(y))] if c}] for x,y in p] for f,p in pairnodes.items()}
  ev=[{j:c for j,c in [(0,x),(1,Q(direction==105+i))] if c} for i,x in enumerate(ev0)]
  od=[{j:c for j,c in [(0,x),(1,Q(direction==120+i))] if c} for i,x in enumerate(od0)]
  con,gs=proper(ring,D,phi,pairs,ev,od)
  assert all(k[0]==1 and k[1]<2 for p in con for k in p)
  assert [[p.get((0,t),Q(0)) for t in range(4)] for p in gs]==Gbase
  for i,p in enumerate(con):
   for mask in [0,1]:
    c=p.get((1,mask),Q(0))
    if c:(Hg if leaf is not None else Jg)[15*mask+i][direction if leaf is not None else direction-105]=c
  for i,p in enumerate(gs):
   for t in range(4):
    if p.get((1,t),Q(0)):Dt[4*i+t][direction]=p[1,t]
  for i,n in enumerate(nunodes):
   if gradient(n):Dt[32+i][direction]=gradient(n)
  if direction%15==14:dump(out,'phase.json',{'phase':'INDEPENDENT_PHYSICAL_TENSOR_JACOBIAN','columns':direction+1,'seconds':time.monotonic()-st})
 assert list(map(enc,Hg))==D['J']['Jg_phi30_by105'] and list(map(enc,Jg))==D['J']['Jg_gamma30']
 assert list(map(enc,Dt))==D['J']['properG32_and_nu2_target_Jacobian34_by135']
 dump(out,'phase1-independent-J135-and-target-checks.json',{'all135_tensor_constraint_and_target_columns_verified':True,'source_dependent_monic_carry_included':True,
  'Jg_phi':list(map(enc,Hg)),'Jg_gamma':list(map(enc,Jg)),'target34_by135':list(map(enc,Dt)),'source_AD_closure':sourceclosure,'source_AD_peak':sourcepeak})
 # Rebuild all known source residuals over the tuple-indexed25-monomial
 # quotient, before any implicit cubic coefficient is supplied.
 mons=[tuple(m) for m in D['ring']['base25_monomials']];mi={m:i for i,m in enumerate(mons)}
 third=[tuple(m) for m in D['rhs']['monomials5']];assert third==[(7,7,7),(0,7,7),(1,7,7),(4,7,7),(5,7,7)]
 for m in mons:
  for i in range(len(m)):assert m[:i]+m[i+1:] in mi
 assert mons[:20]==[tuple(m) for m in D['cache']['base20_monomials_original8']]
 base=Symbolic([tuple(m.count(i) for i in range(8)) for m in mons])
 order=closure(nodes,list(dict.fromkeys(roots+eq)),set(C20));uses=collections.Counter()
 for i in order:
  if i in C20:continue
  n=nodes[i];uses.update(n[1:] if n[0] in ['a','m','n','i'] else [])
 uses.update(roots+eq);vals={};peak=0
 for i in order:
  n=nodes[i];o=n[0]
  if i in C20:v=(C20[i],{})
  elif o=='c':v=base.c(Q(n[1],n[2]))
  elif o=='a':v=base.add(vals[n[1]],vals[n[2]])
  elif o=='n':v=base.scale(vals[n[1]],-1)
  elif o=='m':v=base.mul(vals[n[1]],vals[n[2]])
  elif o=='i':v=base.inv(vals[n[1]])
  else:raise AssertionError(n)
  assert not v[1];vals[i]=v;peak=max(peak,sum(len(p) for p,g in vals.values()));assert peak<100000
  if i not in C20:
   for j in n[1:] if o in ['a','m','n','i'] else []:
    uses[j]-=1
    if not uses[j]:vals.pop(j,None)
 pairs={f:[[vals[x][0],vals[y][0]] for x,y in p] for f,p in pairnodes.items()};nu=[vals[n][0] for n in nunodes]
 ring=Cover(mons,nu);ev=[{i:Q(c) for i,c in enumerate(p) if Q(c)} for p in D['gauge']['actual_even15_base20']]
 od=[{i:Q(c) for i,c in enumerate(p) if Q(c)} for p in D['gauge']['actual_odd15_base20']]
 con,gs=proper(ring,D,C20,pairs,ev,od)
 assert all(i>=20 and mask<2 for p in con for i,mask in p)
 rhsE=[[vals[e][0].get(mi[m],Q(0)) for e in eq] for m in third]
 rhsG=[[p.get((mi[m],mask),Q(0)) for mask in [0,1] for p in con] for m in third]
 assert list(map(lambda p:list(map(str,p)),rhsE))==D['rhs']['source105_RHS3_by_monomial']
 assert list(map(lambda p:list(map(str,p)),rhsG))==D['rhs']['gauge30_RHS3_by_monomial']
 dense=lambda p:[str(p.get((i,t),Q(0))) for i in range(25) for t in range(4)]
 assert list(map(dense,gs))==D['rhs']['known_properG8_C100']
 assert [[str(p.get(i,Q(0))) for i in range(25)] for p in nu]==D['rhs']['known_nu2_base25']
 for j,e in enumerate(eq):
  target={}
  if 48<=j<57:
   v=j-48
   if v in [0,1,4,5]:target[mi[(v,)]]=Q(1)
   elif v>=7:target[mi[(7,)]]=Q(1)
  assert {i:c for i,c in vals[e][0].items() if i<20}==target
 jac=list(map(lambda p:{j:Q(c) for j,c in enumerate(p) if Q(c)},D['joint']['fixed105_Jacobian']))
 ji=list(map(lambda p:{j:Q(c) for j,c in enumerate(p) if Q(c)},D['joint']['fixed105_Jacobian_inverse']))
 jgi=list(map(sp,D['J']['Jg_gamma_inverse']))
 assert all(rowmul(p,jac)=={i:Q(1)} for i,p in enumerate(ji))
 assert all(rowmul(p,Jg)=={i:Q(1)} for i,p in enumerate(jgi))
 own3={}
 for i,p in D['joint']['all105_unknown_jets']:own3[i]=sum(Q(c) for m,c in zip(D['joint']['monomials'],p) if sum(m)==3)
 assert all(rhsE[0][i]+sum(c*own3[unknown[j]] for j,c in row.items())==0 for i,row in enumerate(jac))
 dump(out,'phase2-independent-RHS-and-response.json',{'source105_RHS3':list(map(lambda p:list(map(str,p)),rhsE)),
  'gauge30_RHS3':list(map(lambda p:list(map(str,p)),rhsG)),'all_knownG8_nu100_entries_verified':True,
  'all105_source_pure_r3_crosschecks':True,'canonical_gauge_r3_not_assumed_cached':True,'source_evaluation_peak':peak})
 del vals,pairs,con,gs,ring,ev,od
 # Independent native primitive derivatives, including reduction by monics.
 # Only N0 coefficients are frozen; omitted derivative-N3 terms multiply
 # the separately verified entire native S0*c0 vector, which is zero.
 OP=D['Op']['native_S12_by32'];N0=[[{(0,b):Q(c) for b,c in enumerate(p) if Q(c)} for p in row] for row in D['N0']['N']]
 hi=D['high'];m3=list(itertools.combinations_with_replacement(range(8),3));m3i={m:i for i,m in enumerate(m3)}
 gamma=rec(hi['Gamma4_56_times120']);grows=[{s:c for k,c in gamma.items() for s,j in [divmod(k,120)] if j==m3i[m] and s<48} for m in third]
 cases=[]
 for row in hi['all619']:cases.append({j-640:c for j,c in rec(row['actual_original768_constant_combination']).items() if j>=640})
 used=sorted({j for row in cases for j in row});prim=[[{} for j in range(128)] for m in third]
 Sn0=[[{} for j in range(32)] for i in range(12)]
 for i,j,k,sign in OP['koszul28_in_first8_rows']:Sn0[i][j]={(0,b):Q(sign)*c for b,c in enumerate(Gbase[k]) if c}
 for i,j,c in OP['free4_zero_source_columns']:Sn0[i][j]={(0,0):Q(c)}
 r0=Cover([()], [{},{}]);scols=[]
 for j in range(128):
  jn,b=divmod(j,4);vec={}
  for i in range(12):
   p=r0.mul(Sn0[i][jn],{(0,b):Q(1)})
   for (basei,mask),c in p.items():put(vec,4*i+mask,c)
  scols.append(vec)
 for case in cases:
  z={}
  for j,c in case.items():z=adds(z,scols[j],c)
  assert not z
 for var in range(34):
  r=Cover([(),(0,)], [{1:Q(1)} if var==32+i else {} for i in range(2)])
  GG=[{**{(0,b):c for b,c in enumerate(p) if c},**({(1,var%4):Q(1)} if var<32 and var//4==i else {})} for i,p in enumerate(Gbase)]
  Sn=[[{} for j in range(32)] for i in range(12)]
  for i,j,k,s in OP['koszul28_in_first8_rows']:Sn[i][j]=r.scale(GG[k],Q(s))
  for i,j,c in OP['free4_zero_source_columns']:Sn[i][j]=r.c(Q(c))
  for j in used:
   jn,b=divmod(j,4);vec={}
   for o in range(12):
    p={}
    for k in range(12):
     if Sn[k][jn] and N0[k][o]:p=r.add(p,r.mul(N0[k][o],Sn[k][jn]))
    p=r.mul(p,{(0,b):Q(1)})
    for t in range(4):
     if p.get((1,t),Q(0)):vec[4*o+t]=p[1,t]
   for m,grow in enumerate(grows):
    c=sum(a*vec.get(i,Q(0)) for i,a in grow.items())
    if c:prim[m][j][var]=c
 assert [[enc(p) for p in row] for row in prim]==D['factors']['native_source_derivative_by_monomial_and_column']
 dump(out,'phase3-independent-native-highest-factors.json',{'all619_entire_native_S0_c0_zero':True,
  'all34_native_Gnu_derivatives_replayed':True,'all5_by128_requested_primitive_covectors_match':True,
  'N3_omitted_only_after_separate_native_origin_zero':True})
 # Verify each literal adjoint; do not solve for an unknown coefficient or
 # construct the full cubic solution vector.
 req={(row['kernel619_index'],tuple(row['monomial_original8'])):row for row in D['values']['all_nonzero_requested_adjoints']}
 results=[];checked=0;direct_known=[]
 for number,(case,old) in enumerate(zip(cases,D['oldhigh']['all619'])):
  values=[];known=[]
  for m,mon in enumerate(third):
   alpha={}
   for j,c in case.items():alpha=adds(alpha,prim[m][j],-c)
   oldphi={}
   for (kind,leaf,mm),c in old['old_high_source3_sensitivity']:
    if kind=='physical' and mm==m3i[mon]:put(oldphi,unknown.index(leaf),Q(c))
   target=adds(rowmul(alpha,Dt),oldphi)
   phi={i:c for i,c in target.items() if i<105};gauge={i-105:c for i,c in target.items() if i>=105}
   if target:
    row=req.pop((number,mon));assert sp(row['native_gradient34'])==alpha and sp(row['old_physical_gradient105'])==oldphi
    rp=sp(row['rho_phi105']);rg=sp(row['rho_gauge30'])
    assert adds(rowmul(rp,jac),rowmul(rg,Hg))==phi and rowmul(rg,Jg)==gauge
    value=-sum(c*rhsE[m][i] for i,c in rp.items())-sum(c*rhsG[m][i] for i,c in rg.items())
    assert str(value)==row['actual_implicit_third_functional_value'];checked+=1
   else:value=Q(0)
   values.append(value)
   # The phase2 response is deliberately separate from the implicit value.
   # Its existence and exact source evaluation were verified above.
  assert list(map(str,values))==D['values']['all619_actual_implicit_third_functional_values_by5'][number]
  results.append(list(map(str,values)))
 assert not req and checked==152
 ep=out/'independent-selected135-source-native-and-adjoint-postverification.json'
 dump(out,ep.name,{'inputs':inputs,'all135_actual_tensor_Jacobian_columns_independently_replayed':True,
  'all_source_dependent_monic_carries_and_nonzero_odd_gauge_center_retained':True,
  'all105_source_and30_gauge_RHS3_and_knownG8_nu_response_verified':True,
  'all619_separate_native_origin_vectors_zero':True,'all_native_primitive_covectors_verified':True,
  'literal_adjoints_verified':checked,'all135_components_per_literal_adjoint_checked':True,
  'all619_implicit_contributions_by5':results,
  'formula_and_scope':'phase4 is ONLY -rho*f; direct lower target q is independently verified in phase2 and remains separate for later nonlinear native/high619 assembly',
  'no_unknown_solve_or_full105_gauge30_third_vector_by_reviewer':True,
  'no_complete619_quartics_HF4_or_exactTlength_assertion':True})
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'actual_J135_monic_carry_RHS_and152_adjoints_independently_accepted':True,
  'all619_native_S0_c0_and_implicit_contributions_accepted':True,'known_response_q_and_implicit_minus_rho_f_kept_distinct':True,
  'new_unknown_solve_or_complete619quartics_or_HF4':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-st,
  'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
