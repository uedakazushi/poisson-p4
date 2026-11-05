#!/usr/bin/env python3
"""Fixed source identities for a streamed Cprime/ell^3 cubic upper test."""
import hashlib,itertools,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CPRIME_CUBIC_ELL3_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'quadratic_report':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/report.json',
 'quadratic_values':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase4-all-twenty-actual-quadratic-values.json',
 'first':'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63/phase1-actual-fullnine-first-V-columns.json',
 'q4':'reports/AMBIENT_B18_CPRIME_PRODUCT_ANN_QUOTIENT_SOURCE_PREFLIGHT_20261105_v63/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def sp(v):return [[i,c]for i,c in sorted(v.items())if c]
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def add(v,w,c=Q(1)):
 z=dict(v)
 for i,a in w.items():
  z[i]=z.get(i,Q(0))+c*a
  if not z[i]:del z[i]
 return z
def build():
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert D['quadratic_report']['all_passed']and D['quadratic_report']['actual_complete_rank']==1
 assert D['q4']['all_passed']and D['q4']['known_rank131_preserved']
 L=[rec(v)for v in D['first']['all44_actual_linear_columns27_sparse']]
 q4=[rec(v)for v in D['q4']['literal_fixed_pi4_rows_by135_sparse']]
 E4=D['q4']['unit_NF_rows135'];assert len(E4)==4
 assert [[row.get(i,Q(0))for i in E4]for row in q4]==[[Q(i==j)for j in range(4)]for i in range(4)]
 V=D['quadratic_values'];k=V['fixed_rank_columns'][0];v=list(map(Q,V['all20_actual_quadratic_values4'][k]));p=V['fixed_rank_rows'][0];assert v[p]
 free=[i for i in range(4)if i!=p];q3=[add(q4[i],q4[p],-v[i]/v[p])for i in free];E3=[E4[i]for i in free]
 assert [[row.get(i,Q(0))for i in E3]for row in q3]==[[Q(i==j)for j in range(3)]for i in range(3)]
 m2=list(itertools.combinations_with_replacement(range(9),2));i2={m:i for i,m in enumerate(m2)}
 raw=V['all44_actual_V2_pureplane_nine_values'];K=rec(V['source_constant_kernel20'][k]);P={}
 for j,c in K.items():
  for t in range(3):
   for z,mon in enumerate([(7,7),(7,8),(8,8)]):
    a=c*Q(raw[j][3*t+z]);row=45*t+i2[mon]
    if a:P=add(P,{row:a})
 assert [sum((a*P.get(i,Q(0))for i,a in row.items()),Q(0))for row in q4]==v
 assert all(sum((a*P.get(i,Q(0))for i,a in row.items()),Q(0))==0 for row in q3)
 digest=hashlib.sha256();count=0
 for col in L:
  for var in range(9):
   z={}
   for i,c in col.items():t,x=divmod(i,9);z=add(z,{45*t+i2[tuple(sorted((x,var)))]:c})
   assert all(sum((a*z.get(i,Q(0))for i,a in row.items()),Q(0))==0 for row in q3)
   digest.update(json.dumps(enc(sp(z)),separators=(',',':')).encode()+b'\n');count+=1
 assert count==396
 m3=list(itertools.combinations_with_replacement(range(9),3));i3={m:i for i,m in enumerate(m3)};proj=[]
 for row in range(495):
  t,j=divmod(row,165);mon=m3[j];div=45*t+i2[mon[:2]];var=mon[2]
  proj.append({9*f+var:q.get(div,Q(0))for f,q in enumerate(q3)if q.get(div,Q(0))})
 return D,L,P,q3,E3,proj,m2,m3,i2,i3,{'actual_quadratic_kernel_column':k,'actual_quadratic_pivot_coordinate':p,'rank1_residue4':v,'known396_products_sha256':digest.hexdigest()}
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D,L,P,q3,E3,proj,m2,m3,i2,i3,receipt=build()
 data={**receipt,'coordinate_order9':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],
  'actual_quadratic_initial_representative135_sparse':sp(P),'exact_NF3_rows_by135_sparse':list(map(sp,q3)),
  'unit_embedding_rows135':E3,'canonical_cubic_projection_columns495_sparse':list(map(sp,proj)),
  'all396_first_products_and_actual_quadratic_killed':True,'all9_NF_unit_checks':True,
  'source_count':{'first44_times_S2':1980,'actual_quadratic_times_S1':9,'ell_cube_targets':3,'consistency':27,'total':2019},
  'compressed_target_dimension':27,'whole_saved_live_QQ_upper':24000,'derived_bit_guard':2048,'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_input_plus_output_bytes':32<<20,
  'stream_contract':'Regenerate2019 columns twice; keep at most27 selected original columns andnormalized pivots in pass1, releasepivots, then storeonlysmallinverse andfactorhash inpass2. Never retainall2019columns.',
  'source_initial_proof':'P is the pureplane part of actual V2*K18. NF4 reads onlypureplane, so its omitted part is in knownrank131 initialspan. Therefore P is itself an actualquadraticinitial relation. NF3 has kernel knownquadraticspan+P byNF4 exactness andrank1.',
  'cubic_presentation_proof':'Use canonical divisor(firsttwo sorted variables) toreduceF3 to NF3(F2)*S1. Includeall27 I-qE consistencycolumns andall1980+9 knowncubicproducts, so compressed cokernel is exactlythe knowncubic quotient. Add3ell^3 onlyforCprime/ell^3.',
  'positive_completion_criterion':'Ifall27rowsare spanned, m^3(Cprime/ell^3)=0 byNakayama; length<=3+3+3=9. IndependentDVRrank2/tors{1,e>=2} giveslength>=9, hencee=2. A nonzeroknowncubicquotient is inconclusive.',
  'quadratic_value_peer_required_before_worker':True,'no_cubic_rank_or_values_computed':True,'no_new_source_gauge_or_higherjet':True}
 ep=out/'actual-rank1-quadratic-source-and-streamed-cubic-gate.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_actual_quadratic_to_streamed_cubic_upper','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'whole_saved_live_QQ_upper':24000,'cubic_rank_or_exponent_computed':False,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':main()
