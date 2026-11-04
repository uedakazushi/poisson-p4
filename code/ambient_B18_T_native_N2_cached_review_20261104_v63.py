#!/usr/bin/env python3
"""Independent cached C44 adapter and1460 actualL2 entries/source contractions."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_T_native_S2_B11_review_20261104_v63 import B44,V

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW'
A='AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63'
PATHS={
 'report':f'reports/{A}/report.json','newR':f'reports/{A}/phase1-actual-offintersection-R12-B11.json',
 'newN':f'reports/{A}/phase2-actual-N-B11.json','values':f'reports/{A}/actual-L2-requested-values-and619-native-contributions.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'I':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json',
 'P':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json',
 'Q':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63/actual-proper-source14-equations17.json',
 'H':'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json',
 'K':'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',
 'C':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/membership-and-source-cofactors.json',
 'prep':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json',
 'E':'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json',
 'PF':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
 'B':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
 'p1':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
 'p2':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase2-gauge30-B11-correction.json',
 'p3':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'cache_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
 'firstR':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',
 'firstN':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'origin':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'S':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'A119':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'partial':'reports/AMBIENT_B18_T_NATIVE_FIRST_PARTIAL_RESIDUAL_20261104_v63/actual619-native-old-first-partial-residuals.json',
 'Gamma':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63_REVIEW/independent619-high-order-source-folds.json',
 'adjoint':'reports/AMBIENT_B18_T_NATIVE_L2_ADJOINT_PREFLIGHT_20261104_v63/all619-full-c1-native-L2-adjoints.json',
 'ring_code':'code/ambient_B18_T_native_S2_B11_review_20261104_v63.py',
 'tensor_code':'code/ambient_B18_comparison_X_origin_review_20261104_v63.py'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(x):return {int(i):Q(c) for i,c in x if Q(c)}
def enc(x):return [[list(i) if isinstance(i,tuple) else i,str(c)] for i,c in sorted(x.items()) if c]
def put(v,k,c):
 a=v.get(k,Q(0))+c
 if a:v[k]=a
 else:v.pop(k,None)
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def unpack(A,ring):return [[ring.flat(p) for p in row] for row in A]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'independent_cached_N2_inputs_saved_no_source_solve','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items() if not p.endswith('.py')}
 assert D['report']['all_passed'] and D['cache_peer']['all_passed']
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 nu=[{i:Q(c) for i,c in enumerate(p) if Q(c)} for p in D['p1']['nu_A_B_B11']]
 ring=B44(nu);ring.install();S,I,P,K,H,C,prep=[D[k] for k in ['source','I','P','K','H','C','prep']]
 first={i:rec(p) for i,p in D['B']['all105_actual_B11_root_values']}
 def base_dense(p):return ring.base({i:Q(c) for i,c in enumerate(p) if Q(c)})
 GL=D['PF']['GL18_source_parameter_positions'];positions=I['actual_group15_positions']
 def family_frame(family):
  f=S['families'][family];tau=ring.tau(0 if family=='CA' else 1)
  p=ring.add(tau,ring.scale(ring.base(first[f['p_quadratic_relation']['s']]),Q(1,2)));M=V.eye(5)
  for name,(a,b) in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma') and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);M[i][j]=ring.add(M[i][j],ring.base(first[a]),ring.mul(ring.base(first[b]),p))
  F=V.eye(5)
  if family=='CA':F[0][3]=ring.scale(p,-2);F[1][0]=p;F[1][3]=ring.neg(ring.mul(p,p))
  else:F[0][3]=ring.scale(p,Q(1,2));F[1][0]=p;F[1][3]=ring.scale(ring.mul(p,p),Q(1,4))
  return V.matmul(M,F)
 G=V.eye(5)
 for (i,j),a,b in zip(positions,D['p2']['even15_B11'],D['p2']['odd15_B11']):G[i][j]=ring.add(G[i][j],base_dense(a),ring.mul(base_dense(b),ring.tau(0)))
 assert V.normalized_flag_frame(V.matmul(G,family_frame('CA')))==V.eye(5)
 FB=V.normalized_flag_frame(V.matmul(G,family_frame('CB')))
 def physical(family):
  tau=ring.tau(0 if family=='CA' else 1)
  return [ring.add(base_dense(p),ring.mul(base_dense(q),tau)) for p,q in D['p1']['physical150_tracepairs_B11'][family]]
 PA,PB=physical('CA'),physical('CB');rows=P['zero_coordinate_rows36']
 moved=V.tensor_push(PA,G);center=list(map(Q,P['actual_CA_source19_center']))
 own=[r[0] for r in V.matmul(V.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in P['CA_constant_left_inverse_rows19']])]
 assert all(own[j]==ring.c(center[j]) for j in P['CA_fiber10_fixed_coefficient_indices'])
 actual=[ring.add(own[j],ring.c(-center[j])) for j in P['CA_free_fiber9_coefficient_indices']]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]
 kappa=[actual[j] for j in I['retained_source9_indices']]
 assert actual==[ring.flat(p) for p in D['p3']['actual14_C_B44']] and kappa==[ring.flat(p) for p in D['p3']['kappa9_C_B44']]
 targets=[18,19,28,29,60,67,72,63,64,68,73,74];F12=[ring.add(PA[j],ring.neg(PB[j])) for j in targets]
 assert F12==[ring.flat(p) for p in D['p3']['physicalF12_C_B44']] and any(F12)
 R17=[[{}]*12 for _ in range(17)];transform=V.matmul(FB,G);ann=V.constmat(P['CB_annihilator17_by36'])
 for j,t in enumerate(targets):
  x=[{}]*150;x[t]=ring.ONE
  for k,row in enumerate(D['E']['eliminated_zero_rows4']):x[row]=ring.c(D['E']['eliminated_four_coordinate_reconstruction'][k][rows.index(t)])
  pushed=V.tensor_push(x,transform);col=V.matmul(ann,[[pushed[i]] for i in rows])
  for i in range(17):R17[i][j]=col[i][0]
 assert R17==unpack(D['newR']['R17'],ring)
 F17=[[V.polynomial(p,actual)] for p in D['Q']['equations14variable_sparse']]
 assert V.matmul(R17,[[p] for p in F12])==F17
 du=ring.add(ring.ONE,ring.neg(ring.mul(kappa[1],kappa[1])));pu=ring.add(ring.c(5),ring.mul(kappa[1],kappa[1]))
 roots=[]
 for name in ['v0','v1','b','c','e']:
  p=H['five_explicit_rational_roots'][name];roots.append(ring.tm(V.polynomial(p['numerator'],kappa),ring.power(du,-p['D_denominator_power'])))
 co=[V.synthetic_prefix_coefficients(p,actual,[0,1,10,11,13],roots) for p in D['Q']['equations14variable_sparse']]
 BA=V.matmul(co[5:],V.matrix_inverse(co[:5]));R12=[[ring.add(R17[i+5][j],ring.neg(ring.add(*(ring.mul(BA[i][k],R17[k][j]) for k in range(5))))) for j in range(12)] for i in range(12)]
 for i,p in enumerate(H['all12_actual_remaining_and_divided_records']):
  unit=ring.scale(ring.power(du,p['D_denominator_power']),Q(p['constant_denominator'],p['primitive_content_unit']));R12[i]=[ring.mul(unit,q) for q in R12[i]]
 assert R12==unpack(D['newR']['R12'],ring)
 hpolys=[p['divided_primitive_polynomial'] for p in H['all12_actual_remaining_and_divided_records']];hv=[[V.polynomial(p,kappa)] for p in hpolys]
 assert V.matmul(R12,[[p] for p in F12])==[[ring.mul(kappa[1],p[0])] for p in hv]
 dump(out,'phase.json',{'phase':'independent_actual14_full17_and_R12_checked','seconds':time.monotonic()-start})
 core=[kappa[j] for j in [0,1,2,6]];axes=[K['source9_variable_order'].index(n) for n in K['selected5_pivot_order']];r2=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];r2.append(ring.tm(V.polynomial(p['numerator'],core),ring.power(du,-p['D_power']),ring.power(pu,-p['P_power'])))
 co=[V.synthetic_prefix_coefficients(p,kappa,axes,r2) for p in hpolys];perm=K['selected5_source_H_rows']+[j for j in range(12) if j not in K['selected5_source_H_rows']]
 Ai=V.matrix_inverse([co[j] for j in perm[:5]]);BA=V.matmul([co[j] for j in perm[5:]],Ai);T=[[{}]*12 for _ in range(12)]
 for i in range(5):
  for j in range(5):T[i][perm[j]]=Ai[i][j]
 for i in range(7):
  T[i+5][perm[i+5]]=ring.ONE
  for j in range(5):T[i+5][perm[j]]=ring.neg(BA[i][j])
 for i,p in enumerate(K['remaining7_actual_core_polynomials']):
  unit=ring.scale(ring.tm(ring.power(du,-p['D_denominator_power']),ring.power(pu,-p['P_denominator_power'])),Q(p['QQ_unit_content_numerator'],p['QQ_unit_content_denominator']))
  T[i+5]=[ring.mul(ring.inv(unit),x) for x in T[i+5]]
 E3=V.eye(12);E3[5][5]=ring.inv(du);E3[6][6]=ring.inv(du);r=ring.add(core[2],ring.neg(core[0]),ring.neg(core[1]))
 E3[7][5]=ring.neg(ring.tm(r,ring.inv(ring.mul(pu,du))));E3[7][7]=ring.scale(ring.inv(ring.mul(pu,du)),-2);T=V.matmul(E3,T);E4=V.eye(12)
 for j in range(4):
  for i in range(3):E4[8+j][5+i]=ring.neg(ring.tm(ring.power(du,[2,2,2,3][j]),V.polynomial(C['cofactor3_by4'][i][j],core)))
 T=V.matmul(E4,T);N=V.matmul(T,R12);Gs=[ring.add(kappa[j],ring.neg(x)) for j,x in zip(axes,r2)]+[V.polynomial(p,core) for p in prep['generators3']]
 assert Gs==[ring.flat(p) for p in D['p3']['actualG8_C_B44']]
 assert T==unpack(D['newN']['T_C_B44'],ring) and N==unpack(D['newN']['N_C_B44'],ring)
 assert V.matmul(T,hv)==[[p] for p in Gs]+[[{}] for _ in range(4)]
 assert V.matmul(N,[[p] for p in F12])==[[ring.mul(kappa[1],p)] for p in Gs]+[[{}] for _ in range(4)]
 lin=[(0,0),(1,8),(3,1),(4,2),(5,5),(6,6)]
 for mat,old in [(R12,D['firstR']['R12']),(T,D['firstN']['T']),(N,D['firstN']['N'])]:
  assert all(q.get((bi,t),Q(0))==Q(o[4*di+t]) for row,orow in zip(mat,old) for q,o in zip(row,orow) for bi,di in lin for t in range(4))
 # Independently rebuild all619 proper residuals and their L2 source adjoints.
 Ni=[[[Q(c) for c in p] for p in row] for row in D['origin']['N_inverse']]
 n0=[[[x.get((0,k),Q(0)) for k in range(4)] for x in row] for row in N]
 for L,R in [(n0,Ni),(Ni,n0)]:
  for i in range(12):
   for j in range(12):
    q=[Q(0)]*4
    for k in range(12):
     for a,c in enumerate(L[i][k]):
      for b,d in enumerate(R[k][j]):
       if not a&b:q[a|b]+=c*d
    assert q==[Q(i==j),Q(0),Q(0),Q(0)]
 ss=list(map(rec,D['S']['source741_zero_constant_combinations']));source121=list(map(rec,D['A119']['actual_selected37_original768_constant_combinations']))
 for combination in D['S']['second84_normalized_source741_combinations']:
  c={}
  for j,a in rec(combination).items():
   for k,b in ss[j].items():put(c,k,a*b)
  source121.append(c)
 assert len(source121)==121
 for c in source121:
  for row in D['X0']['full_X0_56_by768']:
   assert sum((Q(row[j])*a for j,a in c.items() if j>=640),Q(0))==0
 gamma=rec(D['Gamma']['own_Gamma4']);m2=list(itertools.combinations_with_replacement(range(8),2));i2={m:i for i,m in enumerate(m2)};m3=list(itertools.combinations_with_replacement(range(8),3))
 requests=set();sourceadj=[]
 for part,author in zip(D['partial']['all619'],D['adjoint']['all619']):
  rv=rec(part['native_partial_first_residual56_times8']);omega={}
  for index,c in rv.items():
   row,var=divmod(index,8);i,b=divmod(row,4)
   for j in range(12):
    for mask,a in enumerate(Ni[i][j]):
     if a and not mask&b:put(omega,8*(4*j+(mask|b))+var,c*a)
  assert enc(omega)==author['proper_S1c0_minus_S0Xp_a1_48_times8']
  adj={}
  for index,a in gamma.items():
   o,m=divmod(index,120)
   if o>=48:continue
   for v in set(m3[m]):
    pair=list(m3[m]);pair.remove(v);q=i2[tuple(pair)]
    for iv,b in omega.items():
     i,w=divmod(iv,8)
     if w==v:put(adj,(o,i,q),-a*b)
  assert enc(adj)==author['actual_L2_source_adjoint'];sourceadj.append(adj);requests.update(adj)
 scalar={}
 bp={(7,7):2,(0,7):7,(1,7):8,(4,7):9,(5,7):10}
 for o,i,m in requests:
  inp,b=divmod(i,4);outp,t=divmod(o,4);p=ring.mul(N[inp][outp],{(0,b):Q(1)})
  scalar[o,i,m]=p.get((bp[m2[m]],t),Q(0))
 assert len(scalar)==1460
 assert [[list(k),str(v)] for k,v in sorted(scalar.items())]==D['values']['requested_actual_scalar_L2_values']
 contributions=[sum((a*scalar[k] for k,a in adj.items()),Q(0)) for adj in sourceadj]
 assert list(map(str,contributions))==D['values']['all619_native_L2_contribution_values']
 assert sum(bool(v) for v in contributions)==27
 ev={'source_manifest':inputs,'independent_method':'sparse dictionaryC44 ring, independent exterior/symmetric tensorpush, monomialsynthetic orderedprefix division; no author worker rerun',
  'actual14_kappa9_F12_from_cached105_gauge_and_physicalpairs_reconstructed':True,
  'full17_R12_T_N_offintersection_identities_and_all_literal_matrices_checked':True,
  'all_constants_and_firsts_match_frozen_C36':True,'all121_star_L_native_origins_zero_and_N0_inverse_both_products_checked':True,
  'all619_Omega1_and_Gamma4_L2_source_adjoints_independently_checked':True,
  'all1460_actual_scalar_L2_entries':[[list(k),str(v)] for k,v in sorted(scalar.items())],
  'all619_native_L2_contributions':list(map(str,contributions)),'nonzero_native_contribution_indices':[i for i,v in enumerate(contributions) if v],
  'all_monic_nu2_carry_terms_retained':True,'no_new_source_or_gauge_solve_or_physicaljet':True,
  'these_are_not_complete619_quartic_values_or_HF4':True,'unit_inverse_checks':ring.inverse_count}
 ep=out/'independent-cached-N2-source-and1460-L2-values.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'actual_cached_N2_and_all1460_L2_entries_independently_accepted':True,
  'all619_source_adjoints_and27_nonzero_partial_native_contributions_accepted':True,
  'complete619_quartic_values_or_HF4_not_claimed':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
