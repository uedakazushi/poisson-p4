#!/usr/bin/env python3
"""Known degree2 closure of actual35 first relations; no X3/m2 values."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import pivots,invq,mv,enc,sparse,nqq
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63'
PATHS={
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'first_report':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/report.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'Koszul':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()};A=D['first'];assert D['first_report']['all_passed']and A['actual_first_relation_rank']==35
 cols=[]
 for c in A['all849_actual_minimal_first_columns40_sparse']:
  v=[Q(0)]*40
  for i,q in c:v[i]=Q(q)
  cols.append(v)
 lp=A['rank_minor_columns'];lr=A['rank_minor_rows'];Li=[[Q(q)for q in v]for v in A['rank_minor_inverse']]
 QMON=list(itertools.combinations_with_replacement(range(8),2));QI={m:i for i,m in enumerate(QMON)}
 pre={'known_actual_linear_rank':35,'known_variable_multiples':280,'quadratic_target_dimension':180,
  'dense_product_and_quotient_and_kernel_slots_upper':180*280+180*180+849*35,
  'new_X2_X3_m1_m2_or_source_jet':False,'d3_quadratic_values_not_assumed_known':True}
 assert pre['dense_product_and_quotient_and_kernel_slots_upper']<150000
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 quad=[];labels=[]
 for j in lp:
  for var in range(8):
   v=[Q(0)]*180
   for row,c in enumerate(cols[j]):
    if c:t,k=divmod(row,8);v[36*t+QI[tuple(sorted((k,var)))]]+=c
   quad.append(v);labels.append([j,var])
 qp,qr=pivots(quad);rank=len(qp);Ai,det=invq([[quad[j][i]for j in qp]for i in qr])if qp else([],Q(1))
 complement=[i for i in range(180)if i not in qr]
 dual=[]
 for free in complement:
  a=[sum((quad[qp[k]][free]*Ai[k][j]for k in range(rank)),Q(0))for j in range(rank)]
  v=[Q(0)]*180;v[free]=Q(1)
  for i,c in zip(qr,a):v[i]-=c
  dual.append(v)
 assert all(sum((row[i]*c for i,c in enumerate(col)if c),Q(0))==0 for col in quad for row in dual)
 factorhash=hashlib.sha256()
 for c in quad:
  nz=[(j,c[i])for j,i in enumerate(qr)if c[i]]
  a=[sum((row[j]*v for j,v in nz),Q(0))for row in Ai]if rank else[]
  check=[Q(0)]*180
  for k,v in enumerate(a):
   if v:
    for i,x in enumerate(quad[qp[k]]):
     if x:check[i]+=v*x
  assert check==c
  factorhash.update((json.dumps(enc(sparse(a)),separators=(',',':'))+'\n').encode())
 # Full constant kernel of M1: every unknown actual quadratic relation not
 # from a linear coefficient is among these814 exact source combinations.
 zeros=[]
 for j,c in enumerate(cols):
  if j in lp:continue
  a=mv(Li,[c[i]for i in lr]);assert all(sum((cols[lp[k]][i]*a[k]for k in range(35)),Q(0))==c[i]for i in range(40))
  co={j:Q(1)}
  for k,v in enumerate(a):
   if v:co[lp[k]]=co.get(lp[k],Q(0))-v
  zeros.append([[k,v]for k,v in sorted(co.items())if v])
 assert len(zeros)==814
 # Profile literal d3 seconds without treating them as known initial relations.
 profile={
  'all849_M0_zero_and_M1_rank35':True,
  'every_degree2_initial_relation_sources':['variable times actual35 M1 relations','constant kernel814 combinations of the actual M2'],
  'completeness_proof':'For M=M1+M2+..., a=a0+a1+..., order1 cancellation means M1*a0=0; order2 is M2*a0+M1*a1. The814 fixed constant-kernel basis spans all a0 and the280 products span all M1*a1.',
  'd3_actual_firsts_included_in_the35_span':True,
  'd3_second_formula':'d3_keep,2 − C1*Wunit0^-1*d3_pivot,1; C0=0 on remaining five targets',
  'C1_is_firsts_of_the81_pivot_X_columns_after18_elimination':True,
  'C1_not_equal_to_final681_firsts_and_not_computed':True,
  'd3_literal_sector0_degree2_zero':'the three retained wedges lie in sector0 whose sequence is entirely linear; d3_keep,2=0 there. Top2 literal d3 entries are zero too.',
  'why_no_additional_d3_quadratics_adopted':'The correction C1*d3_pivot,1 can be nonzero and contains pivot-column source coefficients not covered by final681 sensitivity cancellation. Nonzero linear d3 columns must also be combined with actual35 M1 pivots before being initial quadratic relations.',
  'projected_missing_functionals':'Each saved quotient dual applied to (M2 times each of814 constant-kernel source combinations)',
  'source_order_frontier':{'wedge_h2_general':'orderedHadamard of (mX)_3 in sector0 row7',
    'constant_X0_zero_formula':'(mXc)_3=m0 X3c+m1 X2c+m2 X1c; only m0 row7 and m1/m2 row7 enter retained wedge readbacks before unit-pivot corrections',
    'm1_row7':'known actual C44 pure-r² data retained',
    'm2_row7':'pure-r³ coefficient of Vgamma_A(z0,r) before division, becoming the r² coefficient of m7; not provided by C44/B11 and not assumed zero',
    'possible_cancellation':'Zero M1 source combinations and homotopy boundaries may remove m2/pivot unknowns only after a separate source-specific quotient-dual proof',
    'no_source_X3_or_new_gauge_jet_authorized_or_run':True}}
 data={'actual_variable_order':A['variable_order'],'actual_target_order':A['target_order'],'quadratic_monomials36':QMON,
  'all280_known_actual_variable_multiple_columns180_sparse':[sparse(c)for c in quad],
  'all280_source_labels_relation_and_variable':labels,'known_quadratic_rank':rank,
  'fixed_rank_pivot_columns':qp,'fixed_rank_pivot_rows':qr,'fixed_minor_determinant':det,'fixed_minor_inverse':Ai,
  'all280_factorizations_verified':True,'factorization_stream_sha256':factorhash.hexdigest(),
  'quotient_dimension_after_known_products':180-rank,'quotient_coordinate_labels':[[i//36,list(QMON[i%36])]for i in complement],
  'exact_quotient_duals_rows_sparse':[sparse(v)for v in dual],
  'all814_zero_linear_source_combinations_in_original849':zeros,'all814_zero_linear_combinations_verified':True,
  'full_degree2_initial_relation_source_theorem':profile,
  'actual_HF2_or_length_not_determined':True,'no_finite_annihilator_or_GB_completeness_claim':True}
 assert nqq(data)<150000
 p=out/'known-quadratic-closure-and-complete814-source-projection.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'status':'PROVED_known_quadratic_closure_and_complete_source_projection_only',
  'known_quadratic_rank':rank,'remaining_quadratic_quotient_dimension':180-rank,
  'actual_HF2_Tnil_proved':None,'all814_actual_zero_linear_sources_complete':True,
  'd3_unknown_correction_retained':True,'new_source_X3_m2_or_gauge_jet':False,
  'visible_QQ_support':nqq(data),'bytes':p.stat().st_size,'input_sha256':{v:sha(ROOT/v)for v in PATHS.values()},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
