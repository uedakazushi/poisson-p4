#!/usr/bin/env python3
"""Actual visible first pivot carries from immutable C44 coefficients only.
No source/gauge solve, old B2, X3, or m2 evaluation.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,nqq,polynomial,matrix_inverse,mm
from ambient_B18_T_native_S2_B11_20261104_v63 import Finite
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63'
PATHS={
 'N2':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
 'N2peer':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW/report.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'nu':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'GL':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'first2':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase2-row20-surrogates-and-true-m1row7.json',
 'first3':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'oldB1':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'd3':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_VALUES_20261104_v63/actual153-d3-quadratic-quotient-values-and-known-rank.json',
 'd3source':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63/actual153-d3-quadratic-adjoints-and12-callback-sensitivities.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'Koszul':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 budget={'C44_N_and_inverse':12672,'C44_A12_by9':4752,'C44_D_and_scratch':40000,'saved_adapters_and_carries':12000,'Wi_and_source_scratch':15000}
 assert sum(budget.values())<150000
 (out/'preflight.json').write_text(json.dumps({'bound':budget,'total_QQ_upper':sum(budget.values()),'source':'ONLY immutable N2_B11, kappa/G8 B11, oldB0 and accepted firsts','no_new_source_or_gauge_jet':True},indent=2)+'\n')
 D={k:read(ROOT/p)for k,p in PATHS.items()};assert D['N2peer']['all_passed']
 (out/'phase.json').write_text(json.dumps({'phase':'all_actual_cached_inputs_bound_before_arithmetic','sha256':{p:sha(ROOT/p)for p in PATHS.values()}},indent=2)+'\n')
 a=Finite(2,[[Q(c)for c in v]for v in D['nu']['nu_A_B_B11']]);P=D['B11'];R=D['radical']
 core=[tuple(map(Q,P['kappa9_C_B44'][i]))for i in [0,1,2,6]]
 vs=[tuple(map(Q,v))for v in P['actualG8_C_B44'][:5]]+[polynomial(a,D['GL']['forward_GL4'][s],core)for s in ['h','c','r','l']]
 N=[[tuple(map(Q,v))for v in row]for row in D['N2']['N_C_B44']]
 Ni=matrix_inverse(a,N);prod=mm(a,N,Ni)
 assert all(prod[i][j]==a.c(i==j)for i in range(12)for j in range(12))
 dd=[[polynomial(a,p,vs)for p in row]for row in R['G8_polynomial_factor_D8_by9']]
 A=mm(a,[row[:8]for row in Ni],dd)
 AF=[]
 for row in A:
  for mask in range(4):
   e=[Q(0)]*44;e[mask]=Q(1);v=a.mul(row[5],tuple(e));AF.append([v[4*i]for i in range(11)])
 assert len(AF)==48 and all(not v[0]for v in AF)
 # Each B11 coefficient is a named exact ordinary coefficient, not a derivative.
 bmons=[(),(7,),(7,7),(0,),(1,),(4,),(5,),(0,7),(1,7),(4,7),(5,7)]
 for v,w in zip(AF,D['first2']['all48_actual_adapter_row20_first9']):
  for i,m in enumerate(bmons):
   if len(m)<=1:assert v[i]==Q(w[0 if not m else m[0]+1])
 old0=[]
 for col in D['oldB1']['source_columns640_rows56_sparse_dual9']:
  old0.append({row:Q(next((c for k,c in sp if k==0),'0'))for row,sp in col if row<48})
 raw=[list(map(Q,v))for v in D['d3']['all140_actual_native_row20_z0r_z4r_z5r']]
 ix=[tuple(m)for m in D['first2']['monomials15']];needed=[ix.index((i,7))for i in [0,4,5]]
 rawold=[]
 for j,v in enumerate(D['first2']['all780_row20_second_surrogates15'][:640]):
  rawold.append([Q(v[k])+sum((c*AF[row][bi]for row,c in old0[j].items()if c),Q(0))for k,bi in zip(needed,[7,9,10])])
 raw=rawold+raw;assert len(raw)==780
 X,W,Min=D['X1'],D['W'],D['minimal'];xp=X['X0_pivot_columns'];co=[[Q(c)for c in row]for row in X['all780_X0_factorization_coefficients']]
 xs=[]
 for col in X['all780_source_columns44_sparse_dual9']:
  v={}
  for row,sp in col:
   q=[Q(0)]*9
   for k,c in sp:q[k]=Q(c)
   v[row]=q
  xs.append(v)
 Z=[Q(0)]*9;O=W['residual_original_X_columns'];Wp=W['W_pivot_residual_columns'];Wr=W['W_pivot_rows'];Wi=[[Q(c)for c in row]for row in W['W_fixed_minor_inverse']]
 Crows=[i for i in range(44)if i not in D['m0']['selected_source_columns24']];I=[Crows[i]for i in Min['fixed_minor_data']['compatible_X0_top18_rows_in_C20']]
 Ji=[[Q(c)for c in row]for row in Min['fixed_minor_data']['inverse']];hi=[[Q(c)for c in row]for row in D['first3']['all54_actual_pivot_h0_readback']];mu1=[Q(c)for c in D['first2']['true_actual_m1_row7_r_coefficients44']]
 C1=[]
 for k in Wp:
  j=O[k];nz=[(p,c)for p,c in enumerate(co[j])if c]
  x1=[[xs[j].get(i,Z)[q+1]-sum((c*xs[xp[p]].get(i,Z)[q+1]for p,c in nz),Q(0))for q in range(8)]for i in range(44)]
  a1=[[sum((row[p]*x1[I[p]][q]for p in range(18)if row[p]),Q(0))for q in range(8)]for row in Ji]
  top=[[x1[i][q]-sum((xs[xp[p]].get(i,Z)[0]*a1[p][q]for p in range(18)),Q(0))for q in range(8)]for i in [42,43]]
  bottom=[]
  for n,i in enumerate([0,4,5]):
   rr=raw[j][n]-sum((c*raw[xp[p]][n]for p,c in nz),Q(0))
   bottom.append(Q(5,2)*rr+sum((mu1[p]*x1[p][i]for p in range(44)),Q(0))-sum((hi[n][p]*a1[p][7]for p in range(18)),Q(0)))
  C1.append(top+[bottom])
 # Reconstruct precisely the previously retained twelve callbacks, including sign.
 Bg=[[Q(c)for c in row]for row in D['d3source']['actual_Bold0_residual_graph4_by81']]
 for k in range(81):
  for n,bi in enumerate([7,9,10]):
   v=Q(D['d3']['all81_pivot_C1_visible_surrogates5_by8'][k][2+n][7])
   assert C1[k][2][n]-v==Q(5,2)*sum((Bg[g][k]*AF[row][bi]for g,row in enumerate([32,36,40,44])),Q(0))
 quad=D['quad'];qm=[tuple(m)for m in quad['quadratic_monomials36']]
 adj=[[[Q(0)]*8 for _ in range(81)]for _ in range(6)]
 for f,sp in enumerate(quad['exact_quotient_duals_rows_sparse']):
  for index,c in sp:
   target,mon=divmod(index,36);p,q=qm[mon];c=Q(c)
   if target>=2:
    assert p==q==7
    for k in range(81):adj[f][k][7]+=c*C1[k][2][target-2]
   else:
    for k in range(81):
     adj[f][k][q]+=c*C1[k][target][p]
     if p!=q:adj[f][k][p]+=c*C1[k][target][q]
 eff=[];varsets=[]
 for f in range(6):
  ef=[];vv=set()
  for h,row in enumerate(Wr):
   for var in range(8):
    c=sum((adj[f][k][var]*Wi[k][h]for k in range(81)if adj[f][k][var]and Wi[k][h]),Q(0))
    if c:
     sector,pair=D['Koszul']['basis_sector_wedge'][2][row];ef.append([row,var,c]);vv.update(pair);vv.add(var)
  eff.append(ef);varsets.append(sorted(vv))
 # Independent replay of the accepted153 constant d3-zero sources with actual C1.
 projected=[]
 for col in D['Koszul']['all8_differentials_sparse_columns'][2]:
  u={}
  for row,mon,c in col:
   if sum(mon)==1 and row in Wr:
    v=mon.index(1);h=Wr.index(row)
    for k in range(81):
     if Wi[k][h]:u[k,v]=u.get((k,v),Q(0))+Q(c)*Wi[k][h]
  projected.append([-sum((adj[f][k][v]*c for (k,v),c in u.items()),Q(0))for f in range(6)])
 for comb in D['d3source']['all153_constant_d3_source_combinations']:
  assert all(sum((Q(c)*projected[j][f]for j,c in comb),Q(0))==0 for f in range(6))
 data={'actual_A12_by9_C44':A,'all48_actual_ResAT_row20_B11':AF,'all12_graph_adapter_zir_callbacks':[[i,row,AF[row][bi]]for i,bi in zip([0,4,5],[7,9,10])for row in [32,36,40,44]],
  'all81_actual_visible_C1_top2_by8_and_wedge3_r':C1,'all243_delta_identities_against_retained12_callbacks_verified':True,
  'all6_actual_C1_adjoint_81_by8':adj,'all6_actual_effective_H1_callbacks_row84_var8_sparse':eff,'callback_source_variable_unions':varsets,
  'all153_actual_zero_linear_d3_second_projections_rechecked_zero':True,
  'visible_coefficient_scope':'only top2 all8 and wedge3 r; omitted C1 directions are NOT assigned zero',
  'remaining_gate':'compose these actual H1 adjoints with full661 X-source correction; remove or bind z3/z6 callbacks before a lower-second value job',
  'new_source_gauge_jet_oldB2_X3_m2_evaluation':False}
 assert nqq(data)<150000
 p=out/'actual81-pivot-first-carries-and-six-H1-source-adjoints.json';s=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(s.encode())<32<<20;p.write_text(s)
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_cached_visible_pivot_first_carries_and_source_adjoints','all243_retained12_callback_delta_identities':True,'all153_d3_second_projections_zero':True,'source_variable_unions':varsets,'visible_QQ_support':nqq(data),'bytes':p.stat().st_size,'new_source_gauge_jet_oldB2_X3_m2_evaluation':False,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
