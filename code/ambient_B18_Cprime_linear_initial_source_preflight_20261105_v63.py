#!/usr/bin/env python3
"""Cached full-nine first V source and bounded initial-module test gate only."""
import hashlib,itertools,json,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_Cprime_ell2_nilimage_source_preflight_20261105_v63 import enc
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'value':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/report.json',
 'valuepeer':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json',
 'kappa':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase3-actual-kappa9-and-normalized-nilframe.json',
 'source33':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'table':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'radicalpeer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
 'core5':'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',
 'nilfinite':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'rank2':'reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 for k in ['value','valuepeer','radicalpeer','nilfinite','rank2']:assert D[k]['all_passed']
 assert D['value']['actual_Cprime2_length']==7
 B=D['table'];basis=B['basis33'];degree=B['basis_degrees'];lin=[i for i,d in enumerate(degree)if d==1];assert len(lin)==9
 # Formal ALL-nine first coefficient map, not a surrogate ell0/base25.
 selected=[0]+lin;index={k:i for i,k in enumerate(selected)};pr=[{index[i]:Q(1)}if i in index else{}for i in range(33)]
 def mul(v,w):
  z={}
  for i,c in v.items():
   for j,a in w.items():
    if not i or not j:z[i+j]=z.get(i+j,Q(0))+c*a
  return {i:c for i,c in z.items()if c}
 def project(v):
  z={}
  for i,c in v:
   for j,a in pr[i].items():z[j]=z.get(j,Q(0))+Q(c)*a
  return {i:c for i,c in z.items()if c}
 assert all(project(B['table_columns_for_each_left_basis'][i][j])==mul(pr[i],pr[j])for i in range(33)for j in range(33))
 assert len(D['kappa']['kappa9_C132'])==9 and all(len(v)==132 for v in D['kappa']['kappa9_C132'])
 roots={i for i,v in D['source33']['all_requested164_root_values']};assert 7790 in roots
 source_rows=D['source']['families']['CA']['actual_zero_odd_I0_generators'];assert source_rows[-1]==[63,7790]
 assert len(D['radical']['R9_generators_polynomials9'])==9
 # Source dimensions only; no first coefficient or rank is evaluated here.
 qmons=list(itertools.combinations_with_replacement(range(9),2));cmons=list(itertools.combinations_with_replacement(range(9),3));assert len(qmons)==45 and len(cmons)==165
 budget={'actual_first_trace_and_R9_input':2000,'all44_actual_linear_columns27':1188,'all396_known_quadratic_columns135_raw_upper':10692,'fixed135minor_and_inverse_upper':36450,'source_pivot_basis_live_upper':18225,'future_q_le3_and_cubic_stream':5000,'all_input_overlap_and_metadata':12000};whole=sum(budget.values());assert whole<150000
 data={'actual_R0_coordinate_order9':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],
  'actual_first_quotient_basis33_indices':selected,'actual_first_quotient_basis_labels':[basis[i]for i in selected],'all1089_formal_first_quotient_products_verified':True,
  'source_Vfirst_recipe':'Project actualkappa9_C132,actualnu_A/B33 to the full10 first-jet base. Evaluate the SAME five unit-root residuals u5 and h,c,r,l, d*R9 and36 tracebasis multiples using the exact first-trace algebra; append seven -z_i*tauA and -gA*tauA, apply fixedGamma. This gives ALL actual27x44 first columns of the full3x44 presentation.',
  'retained_trace_caution':'Proper functions have nilpotenttau constants, so evaluate their full polynomial/rational formula inC40; do not linearize proper polynomials before monictrace reduction. nu carry and unit inverses remain exact.',
  'sourcegA_first':'gA is the exact source7790 value from original10533, not guessedzero; its fullfirst functional factors the proved10dim quotient.',
  'no_new_physical_gauge_or_higher_function_input':True,'complete_first_relation_reason':'fullCprime=cokerVfull with all44 sourcecolumns and V0=0, so ALL degree1 initial relations are theconstant span of these44 actual first columns.',
  'known_quadratic_closure_shape':[135,396],'quadratic_sources':'all44 first columns timesall9variables; these are actualinitialrelations, no assertion they exhaust truequadratics',
  'conditional_cubic_compression':'Only if the known quadratic quotient dimension h2<=3: project each495 degree3monomial by a fixed quadratic divisor to27 orfewer coordinates, appendall3 ell^3targetcolumns andallconsistency I-qE. Streamall1980 known first*S2 plus3ell3, do not save adense495x1983 matrix.',
  'source_positive_test':'If all thirddegreeclasses of cokerV/(ell³) vanish using ONLY these known initialsources and the known cumulativeHF upper is9, then Nakayama givesm³=0 and length<=9. Accepted rank2 plus tors{1,e>=2} giveslength>=9, so e=2 exactly. A failed orlarger upperbound provesno exponentnegative.',
  'conditional_gate_is_not_evaluation':'No actualV1/rank/h2/cubic dimension/lengthCprime3 is computed bythispreflight. No generalGB ornewsource/gaugejet is authorized.',
  'budget':budget,'whole_saved_live_QQ_upper':whole,'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_total_bytes':32<<20,
  'scope_if_not_closed':'Keep the genuine long Cprime[ell] direction and connecting10source gate; do not substitute afinitekernel orsurrogatezero.'}
 p=out/'full-nine-first-source-and-linear-initial-completion-gate.json';p.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_GATE_ONLY_fullnine_Vfirst_and_known_initial_test','actual_first_trace_coefficients_all_cached':True,'all1089_first_functional_factors':True,'whole_saved_live_QQ_upper':whole,'actual_Vfirst_rank_or_exponent_computed':False,'new_source_or_gauge_values':False,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
