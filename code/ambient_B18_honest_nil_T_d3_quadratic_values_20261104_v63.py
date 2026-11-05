#!/usr/bin/env python3
"""Existing-cache six adjoints of153 actual zero-linear d3 sources only."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import pivots,invq,mv,enc,sparse,nqq,polynomial
from ambient_B18_T_native_S2_B11_20261104_v63 import Finite
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_VALUES_20261104_v63'
PATHS={
 'preflight':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63/actual153-d3-quadratic-adjoints-and12-callback-sensitivities.json',
 'preflight_report':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63/report.json',
 'first_peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW/report.json',
 'quadratic':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'phase2':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase2-row20-surrogates-and-true-m1row7.json',
 'phase3':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'nus':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'GL':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
 'Koszul':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 assert D['preflight_report']['all_passed']and D['preflight_report']['all12_adapter_second_callbacks_cancel']and D['first_peer']['all_passed']
 for p,h in D['preflight_report']['input_sha256'].items():assert sha(ROOT/p)==h
 budget={'all81_C1_partial_values':81*(2*8+3),'native140_zir3_values':420,'original780_zir3_surrogates':2340,
  'six81_by8_source_adjoints':6*81*8,'all168_projected_surrogates6':1008,'all153_actual_adjoints6':918,
  'C44_poly_and_source_scratch_upper':20000,'constant_rank_receipts_upper':2000}
 assert sum(budget.values())<150000
 (out/'preflight.json').write_text(json.dumps({'source_QQ_saved_plus_live_upper':sum(budget.values()),'budget':budget,'seconds':60,'bytes':32<<20,'new_source_X3_m2_or_gauge_solve':False},indent=2)+'\n')
 b=Finite(2,[[Q(c)for c in v]for v in D['nus']['nu_A_B_B11']])
 C=D['B11'];core=[tuple(map(Q,C['kappa9_C_B44'][i]))for i in [0,1,2,6]]
 vs=[tuple(map(Q,v))for v in C['actualG8_C_B44'][:5]]+[polynomial(b,D['GL']['forward_GL4'][s],core)for s in ['h','c','r','l']]
 R=D['radical'];native=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  pol=next((p for row,p in col if row==5),[]);v=polynomial(b,pol,vs)
  for j in range(4):
   e=[Q(0)]*44;e[j]=Q(1);w=b.mul(v,tuple(e));native.append([w[4*i]for i in [7,9,10]])
 assert len(native)==140
 mons=[tuple(v)for v in D['phase2']['monomials15']];inds=[mons.index((i,7))for i in [0,4,5]]
 raw=[[Q(v[i])for i in inds]for v in D['phase2']['all780_row20_second_surrogates15'][:640]]+native
 m1=[Q(c)for c in D['phase2']['true_actual_m1_row7_r_coefficients44']]
 X=D['X1'];W=D['W0'];Min=D['minimal'];xp=X['X0_pivot_columns'];co=[[Q(c)for c in v]for v in X['all780_X0_factorization_coefficients']]
 xs=[]
 for col in X['all780_source_columns44_sparse_dual9']:
  v={}
  for row,sp in col:
   a=[Q(0)]*9
   for k,c in sp:a[k]=Q(c)
   v[row]=a
  xs.append(v)
 Z=[Q(0)]*9;O=W['residual_original_X_columns'];Wp=W['W_pivot_residual_columns'];Wr=W['W_pivot_rows'];Wi=[[Q(c)for c in r]for r in W['W_fixed_minor_inverse']]
 M=D['m0'];Crows=[i for i in range(44)if i not in M['selected_source_columns24']]
 I=[Crows[i]for i in Min['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[[Q(c)for c in v]for v in Min['fixed_minor_data']['inverse']]
 hi0=[[Q(c)for c in v]for v in D['phase3']['all54_actual_pivot_h0_readback']]
 C1=[]
 for k in Wp:
  j=O[k];nz=[(a,c)for a,c in enumerate(co[j])if c]
  x1=[[xs[j].get(i,Z)[q+1]-sum((c*xs[xp[a]].get(i,Z)[q+1]for a,c in nz),Q(0))for q in range(8)]for i in range(44)]
  a1=[[sum((row[a]*x1[I[a]][q]for a in range(18)if row[a]),Q(0))for q in range(8)]for row in Ji]
  top=[[x1[i][q]-sum((xs[xp[a]].get(i,Z)[0]*a1[a][q]for a in range(18)),Q(0))for q in range(8)]for i in [42,43]]
  bottom=[]
  for n,i in enumerate([0,4,5]):
   r2=raw[j][n]-sum((c*raw[xp[a]][n]for a,c in nz),Q(0))
   rcoef=Q(5,2)*r2+sum((m1[a]*x1[a][i]for a in range(44)),Q(0))-sum((hi0[n][a]*a1[a][7]for a in range(18)),Q(0))
   v=[Q(0)]*8;v[7]=rcoef;bottom.append(v)
  C1.append(top+bottom)
 # Only the coefficients visible to Pi6 are assigned to C1. Other unknown
 # wedge directions do not occur in any of these six adjoints.
 QMON=[tuple(m)for m in D['quadratic']['quadratic_monomials36']]
 adj=[[[Q(0)]*8 for _ in range(81)]for _ in range(6)]
 for f,sp in enumerate(D['quadratic']['exact_quotient_duals_rows_sparse']):
  for index,c in sp:
   target,mon=divmod(index,36);p,q=QMON[mon];c=Q(c)
   if target>=2:assert p==q==7
   for k in range(81):
    adj[f][k][q]+=c*C1[k][target][p]
    if p!=q:adj[f][k][p]+=c*C1[k][target][q]
 surrogate=[]
 for source in D['Koszul']['all8_differentials_sparse_columns'][2]:
  u={}
  for row,mon,c in source:
   if sum(mon)==1 and row in Wr:
    var=mon.index(1);j=Wr.index(row)
    for k in range(81):
     if Wi[k][j]:u[k,var]=u.get((k,var),Q(0))+Q(c)*Wi[k][j]
  surrogate.append([-sum((adj[f][k][v]*c for (k,v),c in u.items()),Q(0))for f in range(6)])
 combs=D['preflight']['all153_constant_d3_source_combinations'];actual=[]
 for c in combs:
  actual.append([sum((Q(a)*surrogate[j][i]for j,a in c),Q(0))for i in range(6)])
 cp,cr=pivots(actual);rank=len(cp);Ai,det=invq([[actual[j][i]for j in cp]for i in cr])if cp else([],Q(1))
 for col in actual:
  a=mv(Ai,[col[i]for i in cr])if rank else[]
  assert all(sum((actual[cp[k]][i]*a[k]for k in range(rank)),Q(0))==col[i]for i in range(6))
 total=174+rank;rem=6-rank
 data={'all140_actual_native_row20_z0r_z4r_z5r':native,
  'all81_pivot_C1_visible_surrogates5_by8':C1,'unlisted_wedge_directions_are_not_claimed_zero':True,
  'all6_source_adjoint_rows81_by8':adj,'all168_projected_d3_second_surrogates6':surrogate,
  'individual168_surrogates_not_claimed_actual_without_callback_cancellation':True,
  'all153_actual_zero_linear_d3_quadratic_quotient_values6':actual,
  'actual_projection_rank':rank,'rank_pivot_columns':cp,'rank_pivot_rows':cr,'rank_minor_determinant':det,'rank_minor_inverse':Ai,
  'all153_rank_factorizations_verified':True,'known_total_quadratic_initial_rank':total,
  'remaining_quadratic_quotient_upper':rem,'actual_HF2_not_proved_equal_to_upper':True,
  'source_proof':'d3_min,2=−C1 Wunit0^-1 d3_pivot,1, literal kept degree2 zero; exact constant d3-only first kernel153 gives initial quadratic relations. All12 missing adapter2 responses vanish after these combinations; native B11 values are actual.',
  'no_new_source_gauge_jet_X3_m2_or_GB':True}
 assert nqq(data)<150000
 p=out/'actual153-d3-quadratic-quotient-values-and-known-rank.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_extra_known_actual_d3_quadratic_subspace_only',
  'actual_d3_quadratic_projection_rank':rank,'known_total_quadratic_initial_rank':total,'actual_HF2_upper':rem,'actual_HF2_exact':None,
  'no_new_source_gauge_jet_X3_m2_or_GB':True,'visible_QQ_support':nqq(data),'bytes':p.stat().st_size,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
