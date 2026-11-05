#!/usr/bin/env python3
"""Right106 source/callback gate only; no new cubic coefficient value."""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,nqq,pivots,invq,mv
from ambient_B18_honest_nil_T_operator_20261104_v63 import fixture_h
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'cubic':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63/complete-cubic917-typed-source-and-three-NF-adjoints.json',
 'cubic_peer':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json',
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'carry':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/actual81-pivot-first-carries-and-six-H1-source-adjoints.json',
 'gradient':'reports/AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63/phase2-allB-and-P-true-honest-second-gradients.json',
 'roots2':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63/phase1-all-source-root-second-representatives.json',
 'lower2':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63/phase2-actual-native-and-m1-source-inputs.json',
 'higher':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63/actual661-first-corrections-m2-replacement-and-native-highest-sensitivities.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'proper':'reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63/phase3-actual-proper-G8-base25.json',
 'physics':'reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63/phase1-common105-and-physical-source-base25.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def put(p,k,c):
 if c:
  p[k]=p.get(k,Q(0))+c
  if not p[k]:del p[k]
def plus(p,q,c=Q(1)):
 p=dict(p)
 for k,a in q.items():put(p,k,c*a)
 return p
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def sp(v):return [[i,c]for i,c in sorted(v.items())if c]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists()
 D={k:read(ROOT/p)for k,p in PATHS.items()};assert D['cubic_peer']['all_passed']
 (out/'phase.json').write_text(json.dumps({'phase':'all fixed source inputs before first-jet span arithmetic','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()}},indent=2)+'\n')
 X,W=D['X1'],D['W'];xp=X['X0_pivot_columns'];co=[[Q(c)for c in row]for row in X['all780_X0_factorization_coefficients']];O=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[[Q(c)for c in row]for row in W['W_fixed_minor_inverse']];final=[i for i in range(762)if i not in wp]
 raw=[]
 for col in X['all780_source_columns44_sparse_dual9']:raw.append({row:rec(v)for row,v in col})
 R35=D['first']['rank_minor_columns'];assert len(R35)==35 and all(j<681 for j in R35)
 def const_source(j):
  k=final[j];w=rec(W['all_W84_columns_sparse'][k]);rhs={i:w[r]for i,r in enumerate(wr)if r in w};a=[sum((c*rhs.get(i,Q(0))for i,c in enumerate(row)if c),Q(0))for row in Wi]
  v={O[k]:Q(1)}
  for p,c in zip(xp,co[O[k]]):put(v,p,-c)
  for i,c in enumerate(a):
   if c:
    put(v,O[wp[i]],-c)
    for p,b in zip(xp,co[O[wp[i]]]):put(v,p,c*b)
  return v
 weights=[const_source(j)for j in R35];x1=[]
 for src in weights:
  v=[[Q(0)]*8 for _ in range(44)]
  for j,c in src.items():
   for row,vals in raw[j].items():
    for var in range(8):v[row][var]+=c*vals.get(var+1,Q(0))
  x1.append(v)
 m0=[[Q(c)for c in row]for row in D['m0']['m0_24_by44']]
 assert all(sum((row[j]*v[j][var]for j in range(44)),Q(0))==0 for v in x1 for row in m0 for var in range(8))
 basiscols=[[raw[p].get(j,{}).get(0,Q(0))for j in range(44)]for p in xp]+[[v[j][7]for j in range(44)]for v in x1]+[[Q(i==43)for i in range(44)]]
 basecp,baserr=pivots(basiscols);base_rank=len(basecp);assert base_rank==19
 E=D['m0']['selected_source_columns24'];EI=[[Q(c)for c in row]for row in D['m0']['selected_minor_inverse']]
 k42=[Q(i==42)for i in range(44)];a=mv(EI,[row[42]for row in m0])
 for i,c in zip(E,a):k42[i]-=c
 assert all(sum((row[j]*k42[j]for j in range(44)),Q(0))==0 for row in m0)
 basiscols.append(k42)
 cp,rr=pivots(basiscols);rank=len(cp);iv,det=invq([[basiscols[j][i]for j in cp]for i in rr]);assert rank==20 and cp==list(range(18))+[53,54]
 unknown_position=cp.index(54)
 factors=[]
 for v in x1:
  fs=[]
  for var in range(8):
   c=mv(iv,[v[i][var]for i in rr]);assert all(sum((basiscols[cp[k]][i]*c[k]for k in range(20)),Q(0))==v[i][var]for i in range(44));fs.append(c)
  factors.append(fs)
 assert all(not row[43]for row in m0)
 g3=sum(map(Q,D['C0']['cached_actual_CA_gA_all_degree3_coefficients']));assert g3==Q(125,24)
 # Ordinary quotient second source gradients and B0: every R35 raw source closes.
 GR=D['gradient'];bg={i:rec(g)for i,c,g in GR['all1239_B_root_origin_and_gradients']};pg={i:rec(g)for i,root,c,g in GR['all723_P_root_origin_and_gradients']};bv={i:rec(v)for i,v in D['roots2']['all1239_Broot_second_representatives']};Bcols=D['chain']['B48_by640_sparse_columns'];Pcols=D['radical']['actual_P8_by640_columns'];bad=[]
 for rel,src in zip(R35,weights):
  body={};grad={}
  for j,c in src.items():
   if j>=640:continue
   for row,node in Bcols[j]:
    put(body,row,c*bv[node].get(0,Q(0)))
    for leaf,a in bg[node].items():put(grad,(row,leaf),c*a)
   for row,node in enumerate(Pcols[j]):
    if node:
     for leaf,a in pg[node].items():put(grad,(48+row,leaf),c*a)
  if body or grad:bad.append([rel,sp(body),[[list(k),c]for k,c in grad.items()]])
 assert not bad,('R35_BODY_OR_HONEST_SECOND_GATE',bad)
 # All convolution duals Delta3 applied to quadraticrelation*one variable.
 m2=list(itertools.combinations_with_replacement(range(8),2));m3=list(itertools.combinations_with_replacement(range(8),3));im3={m:i for i,m in enumerate(m3)};delta=[rec(v)for v in D['cubic']['cubic_NF_dual3_by600_sparse']];duals=[]
 for f in range(3):
  for var in range(8):
   p={}
   for target in range(5):
    for q,mon in enumerate(m2):
     c=delta[f].get(120*target+im3[tuple(sorted(mon+(var,)))],Q(0))
     if c:p[36*target+q]=c
   if p:duals.append([f,var,p])
 C1=D['carry']['all81_actual_visible_C1_top2_by8_and_wedge3_r'];pairs=list(itertools.combinations(range(8),2));needs=set();rawpairs=set();effective=[]
 for f,var,pi in duals:
  adj=[[Q(0)]*8 for _ in range(81)]
  for index,c in pi.items():
   target,q=divmod(index,36);a,b=m2[q]
   if target>=2:
    assert(a,b)==(7,7)
    for k in range(81):adj[k][7]+=c*Q(C1[k][2][target-2])
   else:
    for k in range(81):
     adj[k][b]+=c*Q(C1[k][target][a])
     if a!=b:adj[k][a]+=c*Q(C1[k][target][b])
  eff=[]
  for j,h in enumerate(wr):
   for v in range(8):
    c=sum((adj[k][v]*Wi[k][j]for k in range(81)if adj[k][v]and Wi[k][j]),Q(0))
    if c:eff.append([h,v,c])
  effective.append([f,var,eff]);active={(h,v)for h,v,c in eff}
  for sec in range(3):
   for row in range(8):
    for mon in m2:
     if row==7 and mon==(7,7):continue
     ex=[0]*8
     for k in mon:ex[k]+=1
     inp=[{}for _ in range(8)];inp[row]={tuple(ex):Q(1)};hh=fixture_h(inp,[1,1,2][sec])
     relevant=False
     for h,v in active:
      if h//28!=sec:continue
      e=[0]*8;e[v]=1
      if hh[h%28].get(tuple(e),Q(0)):relevant=True;break
     if relevant:
      rawpairs.add(mon)
      if row<7:
       for a,b in [(mon[0],mon[1])]+([(mon[1],mon[0])]if mon[0]!=mon[1]else[]):
        if a==7 or a<=row:
         if any(v[j][b]for v in x1 for j in range(44)):needs.add(tuple(sorted((row,a))))
  # Direct quadratic top and mu1 term, all ordinary monomials.
  rawpairs.update(m2[q%36]for q in pi if q//36<2)
 MON=[tuple(m)for m in D['physics']['base25_monomials']];missing=sorted((rawpairs|needs)-set(MON))
 budget={'basis20_and_factor_receipts':400+35*8*20,'source_first_arrays':35*44*8,'convolution_duals_and_adjoint_scratch':8000,'cache_only35_lower_values':10000,'source_stencils':6000,'result106by3_and_rank':1000,'metadata':3000};assert sum(budget.values())<150000
 # The sole extra kernel direction is kept as TWO independent unknown
 # third-source contractions; cancel them only after ALL106 source paths.
 residual=[[factors[p][var][unknown_position]for var in range(8)]for p in range(35)]
 effmap={(f,v):{(h,w):Q(c)for h,w,c in eff}for f,v,eff in effective}
 csens=[];vsens=[]
 for syz in D['cubic']['all106_first_linear_syzygies280_sparse_columns']:
  us=[Q(0)]*3;vs=[Q(0)]*3
  for index,c in syz:
   p,var=divmod(index,8);c=Q(c)
   for f in range(3):
    if var==7:
     for n,zi in enumerate([0,4,5]):us[f]+=c*residual[p][zi]*delta[f].get(120*(n+2)+im3[(7,7,7)],Q(0))
    e=effmap.get((f,var),{})
    for i in range(7):vs[f]-=c*residual[p][i]*e.get((56+pairs.index((i,7)),7),Q(0))
  csens.append(us);vsens.append(vs)
 unknowncancel=all(not c for row in csens+vsens for c in row)
 data={'selected35_actual_first_relations':R35,'all35_constant_original780_source_combinations':list(map(sp,weights)),
  'kernel_m0_rank24_and_dimension20':True,'known_contraction_basis_rank_before_extra':base_rank,'source_contraction_basis_rank':rank,'basis55_column_labels':[['Xp0',p]for p in xp]+[['R35_X1_r',p]for p in R35]+[['literal_e43'],['k42=(I-n0m0)e42']],
  'selected20_contraction_basis_columns':cp,'fixed20_minor_rows':rr,'fixed20_minor_determinant':det,'fixed20_minor_inverse':iv,'all35_by8_X1_factorizations_into20_basis':factors,
  'all12320_first_coefficient_factorization_checks':True,'all35_Bold0_and_opaque_second_gradients_zero':True,
  'extra_kernel_k42_vector44':k42,'extra_kernel_coordinate_in_selected20':unknown_position,'all35_by8_k42_residual_coefficients':residual,
  'two_independent_retained_unknowns':['U=mu2_row7*k42','V=m1_row23_r*k42'],
  'all106_U_sensitivities3':csens,'all106_V_sensitivities3':vsens,'all_two_unknown_contractions_cancel_after106':unknowncancel,
  'selected_basis_uses_no_R35_pure_r_third_callback':True,
  'row15_m1_k42_known_from_second':'gamma_tauAB*r² readback from accepted properC100; no third callback',
  'actual_CA_gA3_r_equal_s':g3,'known_e43_contractions':{'mu2_row7':-g3,'m1_row15':Q(0),'m1_row23':Q(11,12)*g3},
  'basis_contraction_recipes':{'Xp0_mu2':'saved18=-mu1*Xp1_r-(5/2)*Xp2_scalar_row20_r2','R35_X1_r_mu2':'-mu1*X2_R35_r2-(5/2)*X3_R35_scalar_row20_r3 by exactpure-rsourcecycle','Xp0_m1_last':'-m0_last*Xp1_r','R35_X1_r_m1_last':'-m0_last*X2_R35_r2','e43':'literal -gA*tauA and gamma=(tauA,tauAB,tauB-11/12tauA); coefficients from actualCA cubic'},
  'all_nonzero_convolution_duals':[ [f,v,sp(p)]for f,v,p in duals],'all_actual_effective_C1_H2_callbacks':effective,
  'required_raw_X2_monomials':sorted(rawpairs),'required_canonical_m1_pairs':sorted(needs),'missing_pairs_outside_accepted_base25':missing,
  'scalar_X3_r3_recipe':'A1_r*Bold2_r2+A2_r2*Bold1_r plus native row20 lowerthird r3. A0=0; allR35 Bold0=0; nativehighestu5even3/nu3 constantresponses vanish for ALL681 by frozen highergate. No implicit thirdvalue assumedzero.',
  'higher_input_gate':'No actualmu2 orlast-sector higher-vector value needed: evaluate only the selected20 contractions; allR35 X1 directions factor through them.',
  'future_106_values_not_computed':True,'new_source_gauge_jet':False,'aggregate_QQ_upper':sum(budget.values()),'budget':budget}
 assert nqq(data)<150000
 p=out/'right106-full-first-span-and-cache-footprint.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'SOURCE_GATE_ONLY_right106_values_not_run','source_contraction_basis_rank':rank,'known_basis_rank':base_rank,'all35_B0_and_opaque_second_gradients_zero':True,'missing_pairs_outside_base25':missing,'all_two_unknown_third_contractions_cancel_after106':unknowncancel,'selected_basis_uses_no_R35_pure_r_third_callback':True,'cache_gate_closed':not missing and unknowncancel,'new_source_gauge_jet':False,'actual106_cubic_values_computed':False,'aggregate_QQ_upper':sum(budget.values()),'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
