#!/usr/bin/env python3
"""Source-only actual order-three conductor relation; no dual value evaluation."""
import hashlib,itertools,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'closure':'reports/AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63/phase2-complete-streamed-cubic-upper-certificate.json',
 'closure_report':'reports/AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63/report.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def selected(path,keys):
 text=path.read_text();d={};dec=json.JSONDecoder()
 for k in keys:
  token=json.dumps(k);assert text.count(token)==1;i=text.index(':',text.index(token))+1
  while text[i].isspace():i+=1
  d[k]=dec.raw_decode(text,i)[0]
 return d
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 C=selected(ROOT/PATHS['C0'],['actual_C0_basis8_scalar_nodes','all35_regular_coordinate_divisions','cached_actual_CA_gA_all_degree3_coefficients','C0_free9_to_H_free10_columns'])
 X=json.loads((ROOT/PATHS['closure']).read_text());assert X['known_cubic_rank']==26 and json.loads((ROOT/PATHS['C0_peer']).read_text())['all_passed']
 R=selected(ROOT/PATHS['radical'],['actual_X44_by780','actual_V3_by44_rows_after_unit1_cancellation','honest_C0_P_context'])
 free=[i for i in range(27)if i not in X['fixed_rows']];assert free==[25];v=[rec(c)for c in X['selected_original_compressed_columns']];inv=[[Q(c)for c in row]for row in X['fixed_inverse']];q={25:Q(1)}
 for k,row in enumerate(X['fixed_rows']):
  c=-sum((v[j].get(25,Q(0))*inv[j][k]for j in range(26)),Q(0))
  if c:q[row]=c
 assert all(sum((c*w.get(i,Q(0))for i,c in q.items()),Q(0))==0 for w in v)
 proj=[rec(z)for z in X['all495_source_projection_columns27_sparse']];dual={i:sum((c*z.get(k,Q(0))for k,c in q.items()),Q(0))for i,z in enumerate(proj)};dual={i:c for i,c in dual.items()if c}
 m3=list(itertools.combinations_with_replacement(range(9),3));assert all(m3[i%165][0]>=7 for i in dual);assert all(dual.get(165*k+m3.index((8,8,8)),Q(0))==0 for k in range(3))
 tail=C['all35_regular_coordinate_divisions'][0];assert tail['original_pure_CA_source_function']==7790 and tail['restriction_r_s_node']==C['actual_C0_basis8_scalar_nodes'][-1]
 original=list(map(Q,C['cached_actual_CA_gA_all_degree3_coefficients']));assert original[:6]==[Q(0)]*6;coeff=original[6:];assert coeff==[Q(3,8),Q(2),Q(13,6),Q(2,3)]
 changed=[Q(125,24),Q(5,3),Q(1,6),Q(2,375)]
 # Verify the two-variable substitution without applying the residual dual.
 own=[Q(0)]*4
 for s_power,c in enumerate(coeff):
  for ell_power in range(s_power+1):own[ell_power]+=c*Q(__import__('math').comb(s_power,ell_power),5**ell_power)
 assert own==changed
 relation={}
 for ell_power,c in enumerate(changed):
  mon=(7,)*(3-ell_power)+(8,)*ell_power;j=m3.index(mon);relation[j]=-c;relation[330+j]=Q(11,12)*c
 data={'residual_compressed_dual27_sparse':list(q.items()),'residual_actual_cubic_dual495_sparse':list(dual.items()),'residual_dual_annihilates_fixed26columns':True,
  'actual_curve_relation495_sparse':list(relation.items()),'gA_original_r_s_coefficients':coeff,'gA_changed_r_ell_coefficients':changed,
  'original_f63_node':7790,'actual_curve_gA_node':tail['restriction_r_s_node'],'all7_ordered_tail_division_source_records':tail['all7_regular_coordinate_divisions'],
  'actual_V_column_orientation':'last8 columns are -U1..-U8; replacingraw -f63*tauA bycurve -gA*tauA subtracts e_i timesnegative z_i*tauA columns. Thisis an invertibleunipotent columnchange. IfnativeU8alreadyusescurve, the curvecolumnisliteral.',
  'pure_curve_order_three_source':'gA=restriction ofrawf63 toz7=0 inpointed57unitinverse, depends ONLYr,s. Cached10coefficients certifyord3, allfirst6zero. s=r+ell/5 is theacceptedfullphysicalF9coordinatechange.',
  'Gamma_coefficient_orientation':['-gA','0','11*gA/12'],'source_actual_relation_is_not_an_assumption_about_raw_nonplane_V2':True,
  'future_value_shape':[1,8],'future_QQ_cap':8000,'future_bit_guard':2048,'future_seconds':60,'future_memory_bytes':2<<30,'future_total_bytes':32<<20,
  'actual_dual_pairing_not_evaluated':True,'new_physical_source_gauge_or_jet':False,
  'criterion':'Nonzero pairingaddsoneactualcubicinitial totheknownrank26; henceall27killed, length(Cprime/ell3)<=9 andknownlower9 givese2. Zeropairingisonlyapartialnegative.'}
 ep=out/'actual-curve-source-and-one-residual-cubic-functional.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_actual_curve_cubic_single_functional','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'dual_pairing_evaluated':False,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
