#!/usr/bin/env python3
"""Exact12 adapter-second sensitivities for153 d3-only quadratic sources.
Only frozen QQ matrices and literal Koszul entries, no coefficient values.
"""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import mv,enc,sparse,nqq
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63'
PATHS={
 'quadratic':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'quad_report':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/report.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'W0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'Koszul':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()};assert D['quad_report']['all_passed']
 A,M,W,K=(D[k]for k in ['quadratic','minimal','W0','Koszul'])
 assert A['known_quadratic_rank']==174
 ds=M['d3_only_linear'];piv=ds['pivot_columns'];rows=ds['pivot_rows'];Di=[[Q(c)for c in r]for r in ds['minor_inverse']]
 col=[]
 for sp in ds['all168_known_minimal_linear_columns40_sparse']:
  q=[Q(0)]*40
  for i,c in sp:q[i]=Q(c)
  col.append(q)
 # Exactly the full constant kernel of the d3-only first matrix.
 comb=[];labels=[]
 for j,v in enumerate(col):
  if j in piv:continue
  a=mv(Di,[v[i]for i in rows]);assert all(sum((col[piv[k]][i]*a[k]for k in range(15)),Q(0))==v[i]for i in range(40))
  c={j:Q(1)}
  for k,x in enumerate(a):
   if x:c[piv[k]]=c.get(piv[k],Q(0))-x
  comb.append(c);labels.append(j)
 assert len(comb)==153
 Wr=W['W_pivot_rows'];Wp=W['W_pivot_residual_columns'];Wi=[[Q(c)for c in r]for r in W['W_fixed_minor_inverse']];O=W['residual_original_X_columns']
 interim={j:{i:Q(c)for i,c in sp}for j,sp in M['full_first_layer_missing_data_profile']['X2_data']['constant_residual_Bold0_nonzero_list']}
 graphrows=[32,36,40,44]
 B=[[interim.get(O[k],{}).get(g,Q(0))for k in Wp]for g in graphrows]
 # Coefficient of r in the81 pivot components of literal d3.
 rho=[]
 for source in K['all8_differentials_sparse_columns'][2]:
  dr={}
  for row,mon,c in source:
   if mon==[0,0,0,0,0,0,0,1]:dr[row]=dr.get(row,Q(0))+Q(c)
  nz=[(i,dr[r])for i,r in enumerate(Wr)if r in dr]
  u=[sum((r[i]*c for i,c in nz),Q(0))for r in Wi]
  rho.append(mv(B,u))
 # Pi6 sees the wedge part ONLY at r²; top2 part is already cached X1.
 Pi=[]
 for sp in A['exact_quotient_duals_rows_sparse']:
  q=[Q(0)]*180
  for i,c in sp:q[i]=Q(c)
  assert all(not q[36*j+i]for j in [2,3,4]for i in range(35))
  Pi.append([q[36*j+35]for j in [2,3,4]])
 unknown=[];allzero=True
 for c in comb:
  r=[sum((a*rho[j][g]for j,a in c.items()),Q(0))for g in range(4)]
  # unknown(i,g)=the even coefficient ResAT(row20,graphg) at zi*r.
  sens=[[-Q(5,2)*row[i]*r[g]for i in range(3)for g in range(4)]for row in Pi]
  allzero &= all(not v for row in sens for v in row)
  unknown.append({'rho4':r,'sensitivity6_by12':sens})
 data={'actual_d3_only_first_rank':15,'all153_zero_linear_d3_source_labels':labels,
  'all153_constant_d3_source_combinations':[sparse([c.get(i,Q(0))for i in range(168)])for c in comb],
  'all153_zero_linear_combinations_verified':True,
  'potential_adapter_second_callbacks12':[[i,g]for i in [0,4,5]for g in [64,68,73,74]],
  'actual_Bold0_residual_graph4_by81':B,'actual_d3_r_source_transfer_rho4_by168':list(map(list,zip(*rho))),
  'actual_quotient_dual_wedge_r_squared_coefficients6_by3':Pi,
  'all153_unknown_callback_sensitivities':unknown,'all12_callbacks_cancel_on_all153':allzero,
  'd3_quadratic_readback':'Pi6*(-C1*Wunit0^-1*d3_pivot,1); literal kept quadratic part zero',
  'only_needed_pivot_firsts':{'top2':'all8 cached X1 directions after18 unit cancellation',
   'wedge3':'only r direction, hence only source row20 X2 at z0r,z4r,z5r plus TRUE cached m1row7 and known pivot correction'},
  'source_availability':{'native_row20_zir_second':'actual proper kappa/G8 on saved C44/B11; no new source/gauge equation',
   'first640_row20_second':'ResAT1*Bold1 plus the12 retained-or-cancelled adapter2 callbacks; ResAT0=0',
   'new_X3_m2_or_zir_squared_not_required_for_d3_only':True,
   'adapter2_if_needed':'four literal source A_i5 scalar functions at only z0r,z4r,z5r; exact proper source B11 input exists, N2 was not evaluated in earlier packets'},
  'no_actual_adapter2_C1_or_projected_quadratic_values_computed':True,
  'budget_for_source_sensitivity_only':{'upper_QQ_slots':6*153*12+4*168+4*81+153*16+81*81,'seconds':60,'address_bytes':2<<30,'bytes':32<<20}}
 assert nqq(data)<150000
 p=out/'actual153-d3-quadratic-adjoints-and12-callback-sensitivities.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_d3_quadratic_sensitivity_gate',
  'actual_zero_linear_d3_source_dimension':153,'all12_adapter_second_callbacks_cancel':allzero,
  'projected_actual_d3_quadratic_values_computed':False,'new_source_X3_m2_gauge_or_jet':False,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
