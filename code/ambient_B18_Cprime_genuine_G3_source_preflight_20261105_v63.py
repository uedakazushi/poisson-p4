#!/usr/bin/env python3
"""Source theorem/profile for cached actualCprime9 actions and genuineG3."""
import hashlib,json
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_GENUINE_G3_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'exact':'reports/AMBIENT_B18_CPRIME_EXACT_TORSION_12_20261105_v63/report.json',
 'positive_peer':'reports/AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63_REVIEW/report.json',
 'first':'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63/phase1-actual-fullnine-first-V-columns.json',
 'quadratic':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase4-all-twenty-actual-quadratic-values.json',
 'NF3':'reports/AMBIENT_B18_CPRIME_CUBIC_ELL3_SOURCE_PREFLIGHT_20261105_v63/actual-rank1-quadratic-source-and-streamed-cubic-gate.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def main():
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert D['exact']['all_passed']and D['positive_peer']['all_passed']
 F=D['first'];L=[rec(v)for v in F['all44_actual_linear_columns27_sparse']];cp=F['fixed_pivot_columns'];rows=F['fixed_pivot_rows'];iv=[[Q(c)for c in row]for row in F['fixed_minor_inverse']];free=[i for i in range(27)if i not in rows];assert free==[8,25,26]
 q1=[]
 for f in free:
  v={f:Q(1)}
  for k,r in enumerate(rows):
   c=-sum((L[cp[j]].get(f,Q(0))*iv[j][k]for j in range(24)),Q(0))
   if c:v[r]=c
  q1.append(v)
 assert all(sum((c*col.get(i,Q(0))for i,c in v.items()),Q(0))==0 for col in L for v in q1)
 q3=[rec(v)for v in D['NF3']['exact_NF3_rows_by135_sparse']];e2=D['NF3']['unit_embedding_rows135'];assert e2==[44,133,134]
 data={'coordinate_order9':F['coordinate_order9'],'Gamma_target_order':['tauA+11*tauB/12','tauA*tauB','tauB'],
  'source_L44_columns27_sparse':F['all44_actual_linear_columns27_sparse'],'source_V1_pivot24_columns':cp,'source_V1_pivot24_rows':rows,'source_V1_fixed24_inverse':iv,
  'first_NF3_rows27_sparse':[list(v.items())for v in q1],'first_basis_embedding_rows27':free,'quadratic_NF3_rows135_sparse':[list(v.items())for v in q3],'quadratic_basis_embedding_rows135':e2,
  'actual_9basis_order':'headGamma3, literalfirstfree[8,25,26], literalquadraticfree[44,133,134]',
  'all44_known_quadratic_inputs':'all44_actual_V2_pureplane_nine_values fromthepinnedquadraticartifact, ordinarycoefficients withtargets3 andpairsrr/rell/ellell',
  'head_action_formula':'linear rowe -> (0,q1(e),-pi3(V2_piv)*iv24*restriction_rows(e))',
  'first_action_formula':'literalfirstfree monomial timesvariable -> (0,0,pi3(quadraticmonomial))','quadratic_action_formula':'all9variables killquadraticbasis',
  'allorder_source_theorem':'InCprime/m³ everyvariable timesactualrelationV=L+Q isitslinearpartproduct. ThusknownK2=kerpi3 areactualrelations. ReplaceeachQbyE2*pi3Q; onlysavedpureplanevaluesareread. ThefixedLpivotreductiongivesallexact44relationsand9cyclicbasisvectors. Equalactualdimension9givesanactualmoduleisomorphism.',
  'genuine_target_proof':'Acceptedm³Cprime=ell³Cprime identifiesCprime/m³withCprime/ell³. Acceptedtorsionexponents1,2 implygenuineCprime[ell]=kerell/imell², notrawfinitekernel.',
  'future_checks':['all44L+Qrelations','all9actionscommute','strict3layerblocksupport=>m³zero','all9literalcyclicbasisvectorsI9','ellrank5','ell²rank2','ell³zero','kernel4andimage2withfixedrankcertificates','quotienttarget2withall9inducedactions','targetm²zero'],
  'future_numeric_cap':24000,'named_static_lifetime_bound':18000,'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_input_plus_output_bytes':32<<20,'hard_derived_bit_guard':2048,
  'source_bit_profile':'Inputfirst24inverse andplaneQ2/NF3 arefixedsmallQQ. Preparedworker mustpinacommon-denominator minorboundorcheckallreducedcoefficients<=2048 beforeoutput.',
  'no_new_physical_source_gauge_jet_or_connecting_values':True,'actual9actions_and_G3_values_not_computed':True}
 ep=out/'actual-nine-action-and-genuine-G3-source-interface.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_cached_actual9module_and_genuineG3','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'actual9actions_or_G3_computed':False}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
