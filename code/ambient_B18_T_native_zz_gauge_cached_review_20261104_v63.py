#!/usr/bin/env python3
"""Verify provided270 gauge coefficients and effective native80 source values."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time
from ambient_B18_T_native_S2_B11_review_20261104_v63 import B44,V
from ambient_B18_T_gdiff_second_univariate5_review_20261104_v63 import Symbolic,closure,rec,enc,put

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63_REVIEW'
A='AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63'
PATHS={
 'report':f'reports/{A}/report.json','physical':f'reports/{A}/phase1-actual-physical-source-base20.json',
 'gauge':f'reports/{A}/phase2-canonical270-zz-gauge.json','proper':f'reports/{A}/phase3-actual-proper-G8-base20.json',
 'values':f'reports/{A}/all619-lower-contribution-values.json',
 'cache':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
 'cache_peer':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW/report.json',
 'oldB':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW/independent-source2-transport-and619-oldB-sensitivity.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'I':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json',
 'P':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json',
 'K':'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',
 'prep':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json',
 'Op':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'J':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
 'PF':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
 'N1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'oldphysical':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
 'oldgauge':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase2-gauge30-B11-correction.json',
 'oldproper':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'oldnative':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
 'lower':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/phase4-all619-lower-knownparts-and-symbols.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
MON=[()]+[(i,) for i in [0,1,4,5,7]]+[m for m in itertools.combinations_with_replacement([0,1,4,5,7],2) if m!=(1,1)]
MI={m:i for i,m in enumerate(MON)}
def bm(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():
   k=MI.get(tuple(sorted(MON[i]+MON[j])))
   if k is not None:put(z,k,c*d)
 return z
class Cover(B44):
 def mul(self,a,b):
  out={}
  for (i,u),c in a.items():
   for (j,v),d in b.items():
    p=bm({i:c},{j:d});carry=u&v;mask=u^v
    if carry&1:p=bm(p,self.nus[0])
    if carry&2:p=bm(p,self.nus[1])
    out=self.add(out,{(k,mask):x for k,x in p.items()})
  return out

def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'provided270_and_cached_source_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()};assert D['report']['all_passed'] and D['cache_peer']['all_passed']
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 assert list(map(list,MON))==D['cache']['base20_monomials_original8']
 base=Symbolic([tuple(m.count(i) for i in range(8)) for m in MON]);source=D['source'];dag=source['DAG_scalar_nodes']
 first={i:rec(p) for i,p in D['cache']['all105_actual_source2_values']}
 roots=[i for family in source['families'].values() for pair in family['complete_normalized_physical_w_coefficients'][0] for i in pair]+[family['actual_tracefree_nu'] for family in source['families'].values()]
 order=closure(dag,roots,set(first));uses=collections.Counter()
 for i in order:
  if i in first:continue
  n=dag[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
 uses.update(roots);values={};peak=0
 for i in order:
  n=dag[i];o=n[0]
  if i in first:v=(first[i],{})
  elif o=='c':v=base.c(Q(n[1],n[2]))
  elif o=='a':v=base.add(values[n[1]],values[n[2]])
  elif o=='n':v=base.scale(values[n[1]],-1)
  elif o=='m':v=base.mul(values[n[1]],values[n[2]])
  elif o=='i':v=base.inv(values[n[1]])
  else:raise AssertionError((i,n))
  assert not v[1];values[i]=v;peak=max(peak,sum(len(a) for a,b in values.values()));assert peak<150000
  if i not in first:
   for j in n[1:3] if o in ['a','m'] else n[1:] if o in ['n','i'] else []:
    uses[j]-=1
    if not uses[j]:values.pop(j,None)
 pairs={f:[[values[a][0],values[b][0]] for a,b in data['complete_normalized_physical_w_coefficients'][0]] for f,data in source['families'].items()}
 nu=[values[data['actual_tracefree_nu']][0] for data in source['families'].values()]
 for family in ['CA','CB']:
  assert [[[str(p.get(i,Q(0))) for i in range(20)] for p in pair] for pair in pairs[family]]==D['physical']['physical150_tracepairs_base20'][family]
 assert [[str(p.get(i,Q(0))) for i in range(20)] for p in nu]==D['physical']['nu_A_B_base20']
 ring=Cover(nu);ring.install();I,P,K=[D[k] for k in ['I','P','K']];GL=D['PF']['GL18_source_parameter_positions'];rows=P['zero_coordinate_rows36'];center=list(map(Q,P['actual_CA_source19_center']))
 def frame(family):
  data=source['families'][family];p=ring.add(ring.tau(0 if family=='CA' else 1),ring.scale(ring.base(first[data['p_quadratic_relation']['s']]),Q(1,2)));M=V.eye(5)
  for name,(a,b) in data['coefficient_unknown_pairs'].items():
   if name.startswith('gamma') and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);M[i][j]=ring.add(M[i][j],ring.base(first[a]),ring.mul(ring.base(first[b]),p))
  F=V.eye(5)
  if family=='CA':F[0][3]=ring.scale(p,-2);F[1][0]=p;F[1][3]=ring.neg(ring.mul(p,p))
  else:F[0][3]=ring.scale(p,Q(1,2));F[1][0]=p;F[1][3]=ring.scale(ring.mul(p,p),Q(1,4))
  return V.matmul(M,F)
 def physical(family):return [ring.add(ring.base(a),ring.mul(ring.base(b),ring.tau(0 if family=='CA' else 1))) for a,b in pairs[family]]
 def gauge(ev,od):
  G=V.eye(5)
  for (i,j),a,b in zip(I['actual_group15_positions'],ev,od):G[i][j]=ring.add(G[i][j],ring.base(a),ring.mul(ring.base(b),ring.tau(0)))
  return G
 HA,HB=frame('CA'),frame('CB');PA,PB=physical('CA'),physical('CB')
 def own_coeffs(moved):return [r[0] for r in V.matmul(V.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in P['CA_constant_left_inverse_rows19']])]
 def constraints(G):
  FA=V.normalized_flag_frame(V.matmul(G,HA));own=own_coeffs(V.tensor_push(PA,V.matmul(FA,G)))
  return [FA[1][0],FA[1][3],FA[1][4],FA[0][3],FA[0][4]]+[ring.add(own[j],ring.c(-center[j])) for j in P['CA_fiber10_fixed_coefficient_indices']]
 bmon=[tuple(i for i,n in enumerate(m) for _ in range(n)) for m in D['B11']['B11_basis_monomials8']];oldidx=[MI[m] for m in bmon];newq=[i for i,m in enumerate(MON) if len(m)==2 and i not in oldidx]
 assert newq==D['gauge']['new_quadratic_indices9'] and len(newq)==9
 def old_embed(p):return {oldidx[i]:Q(c) for i,c in enumerate(p) if Q(c)}
 ev0=list(map(old_embed,D['oldgauge']['even15_B11']));od0=list(map(old_embed,D['oldgauge']['odd15_B11']))
 residual=constraints(gauge(ev0,od0));assert all(mask<2 and i in newq for p in residual for i,mask in p)
 rhs=[[p.get((i,t),Q(0)) for i in newq] for t in range(2) for p in residual]
 assert list(map(lambda r:list(map(str,r)),rhs))==D['gauge']['RHS30_by9']
 ev=[{i:Q(c) for i,c in enumerate(p) if Q(c)} for p in D['gauge']['actual_even15_base20']];od=[{i:Q(c) for i,c in enumerate(p) if Q(c)} for p in D['gauge']['actual_odd15_base20']]
 corr=[[Q(c) for c in p] for p in D['gauge']['new_correction30_by9']]
 for j in range(15):
  for i in range(20):
   want=corr[j][newq.index(i)] if i in newq else Q(0)
   assert ev[j].get(i,Q(0))-ev0[j].get(i,Q(0))==want
   want=corr[15+j][newq.index(i)] if i in newq else Q(0)
   assert od[j].get(i,Q(0))-od0[j].get(i,Q(0))==want
 AA=[[Q(c) for c in r] for r in D['J']['A0']];AB=[[Q(c) for c in r] for r in D['J']['A1_eta_delta']]
 JJ=[[Q(0)]*30 for _ in range(30)]
 for i in range(15):
  for j in range(15):JJ[i][j]=AA[i][j];JJ[15+i][j]=AB[i][j];JJ[15+i][15+j]=AA[i][j]
 Ji=[[Q(c) for c in r] for r in D['J']['J30_inverse']]
 assert all(sum(JJ[i][k]*Ji[k][j] for k in range(30))==sum(Ji[i][k]*JJ[k][j] for k in range(30))==Q(i==j) for i in range(30) for j in range(30))
 assert all(sum(JJ[i][j]*corr[j][q] for j in range(30))==-rhs[i][q] for i in range(30) for q in range(9))
 G=gauge(ev,od);assert not any(constraints(G));assert V.normalized_flag_frame(V.matmul(G,HA))==V.eye(5)
 FB=V.normalized_flag_frame(V.matmul(G,HB));own=own_coeffs(V.tensor_push(PA,G));assert all(own[j]==ring.c(center[j]) for j in P['CA_fiber10_fixed_coefficient_indices'])
 actual=[ring.add(own[j],ring.c(-center[j])) for j in P['CA_free_fiber9_coefficient_indices']]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j] for j in I['retained_source9_indices']]
 F12=[ring.add(PA[j],ring.neg(PB[j])) for j in [18,19,28,29,60,67,72,63,64,68,73,74]]
 du=ring.add(ring.ONE,ring.neg(ring.mul(kappa[1],kappa[1])));pu=ring.add(ring.c(5),ring.mul(kappa[1],kappa[1]));core=[kappa[j] for j in [0,1,2,6]];axes=[K['source9_variable_order'].index(n) for n in K['selected5_pivot_order']];roots2=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];roots2.append(ring.tm(V.polynomial(p['numerator'],core),ring.power(du,-p['D_power']),ring.power(pu,-p['P_power'])))
 Gs=[ring.add(kappa[j],ring.neg(p)) for j,p in zip(axes,roots2)]+[V.polynomial(p,core) for p in D['prep']['generators3']]
 def dense80(p):return [str(p.get((i,t),Q(0))) for i in range(20) for t in range(4)]
 for name,values,oldname in [('actual14_C80',actual,'actual14_C_B44'),('kappa9_C80',kappa,'kappa9_C_B44'),('physicalF12_C80',F12,'physicalF12_C_B44'),('actualG8_C80',Gs,'actualG8_C_B44')]:
  assert list(map(dense80,values))==D['proper'][name]
  assert [[str(p.get((i,t),Q(0))) for i in oldidx for t in range(4)] for p in values]==D['oldproper'][oldname]
 dump(out,'phase1-independent270-and-proper-source.json',{'provided270_J0_equations_verified_without_solving':True,'all30_canonical_constraints_zero':True,'all_actual14_kappa_F12_G8_and_B11_projections_checked':True})
 # Rebuild scalar L0+L1 independently. Drop all base-degree2 pieces only
 # here; the zero-origin source theorem, rather than an individual claim,
 # justifies the effective representative after complete combinations.
 N1=D['N1']['N'];linear_map={0:0,**{v+1:MI[(v,)] for v in [0,1,4,5,7]}};Nlin=[]
 for row in N1:
  rr=[]
  for p in row:
   v={}
   for old,new in linear_map.items():
    for t in range(4):
     if Q(p[4*old+t]):v[new,t]=Q(p[4*old+t])
   rr.append(v)
  Nlin.append(rr)
 L={}
 for j in range(48):
  inp,b=divmod(j,4)
  for o in range(12):
   p=ring.mul(Nlin[inp][o],{(0,b):Q(1)})
   for t in range(4):
    q={i:c for (i,mask),c in p.items() if mask==t and len(MON[i])<=1}
    if q:L[4*o+t,j]=q
 Sn=[[{} for _ in range(32)] for _ in range(12)]
 for i,j,g,s in D['Op']['native_S12_by32']['koszul28_in_first8_rows']:Sn[i][j]=ring.scale(Gs[g],Q(s))
 for i,j,c in D['Op']['native_S12_by32']['free4_zero_source_columns']:Sn[i][j]=ring.c(Q(c))
 parts=D['lower'];req={tuple(k) for row in parts['all619_uncomputed_native_second_source_functionals'] for k,c in row};bycol=collections.defaultdict(list)
 for key in req:bycol[key[1]].append(key)
 pairs36=list(itertools.combinations_with_replacement(range(8),2));lookup={}
 for j,keys in bycol.items():
  inp,b=divmod(j,4);sv={}
  for i in range(12):
   p=ring.mul(Sn[i][inp],{(0,b):Q(1)})
   for t in range(4):
    q={k:c for (k,mask),c in p.items() if mask==t}
    if q:sv[4*i+t]=q
  for row in {k[2] for k in keys}:
   z={}
   for k,p in sv.items():
    if (row,k) in L:
     for a,c in bm(L[row,k],p).items():put(z,a,c)
   for key in keys:
    if key[2]==row:lookup[key]=z.get(MI[pairs36[key[3]]],Q(0))
 assert set(lookup)==req
 oldnative=[{r:rec(p) for r,p in col} for col in D['oldnative']['source_columns128_rows48_sparse_B11_effective_second']]
 bmap={m:i for i,m in enumerate(bmon)}
 for (kind,j,row,m),v in lookup.items():
  if pairs36[m] in bmap:assert v==oldnative[j].get(row,{}).get(bmap[pairs36[m]],Q(0))
 newvalues=[[list(k),str(v)] for k,v in sorted(lookup.items()) if pairs36[k[3]] not in bmap]
 assert len(newvalues)==6270 and newvalues==D['values']['new_native_second_keys']
 known=list(map(Q,D['oldB']['all619_oldB_known_parts']));answers=[a+sum(Q(c)*lookup[tuple(k)] for k,c in row) for a,row in zip(known,parts['all619_uncomputed_native_second_source_functionals'])]
 assert list(map(str,answers))==D['values']['all619_actual_lower_contribution_values'] and sum(bool(a) for a in answers)==42
 ev={'source_manifest':inputs,'source_physical150_pairs_and_nu_replayed_from_verified105cache':True,
  'provided270_equations_J0_and_finalcanonicalconstraints_checked_without_solving':True,
  'new270_uniqueness':'new quadratic baseideal annihilates the maximalbaseideal; the accepted C4 origin J30 is a unit, so provided corrections are unique',
  'actual14_kappa_G8_F12_all80_entries_and_allB11_projections_checked':True,
  'effective_L0_L1_includes_only_baseconstant_and_linear_parts':True,
  'all6270_new_native_second_entries':newvalues,'all619_actual_lower_values':list(map(str,answers)),
  'nonzero_lower_indices':[i for i,a in enumerate(answers) if a],
  'no_individual_actual_nativeX2_or_N2zz_claimed':True,'whole_zero_origin_source_combinations_only':True,
  'no_new_physical_source_implicit_solve_gauge_unknown_solve_rawP4_X3_or_HF4':True,'physical_source_QQ_peak':peak}
 ep=out/'independent-provided270-and619-lower-values.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'provided270_gauge_coefficients_independently_verified_without_solving':True,
  'all6270_new_native_keys_and619_lower_values_accepted':True,'nonzero_lower_contributions':42,
  'complete619_quartic_values_or_HF4_not_claimed':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep),'phase1-independent270-and-proper-source.json':sha(out/'phase1-independent270-and-proper-source.json')},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
