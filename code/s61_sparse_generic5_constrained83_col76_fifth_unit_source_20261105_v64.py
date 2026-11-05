#!/usr/bin/env python3
"""SOURCE only: exact fixedrow14+6row35 candidate, no fifthpivot applied."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL76_FIFTH_UNIT_SOURCE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def poly(p):return{tuple(p['shift'][j]+z[j+1]for j in range(5)):Q(p['content'])*int(z[0])for z in p['terms']}
def combine(A,B,c):
 R=dict(A)
 for e,q in B.items():
  R[e]=R.get(e,Q(0))+c*q
  if not R[e]:del R[e]
 return R
def serial(A):return[[str(q),*e]for e,q in sorted(A.items())]
begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_COEFFICIENTS_20261105_v64/report.json';assert sha(rp)=='1c6e8c908858d66248a04499b3671fc399760054088ef65b8a265d4b32f7eeba';r=read(rp);assert r['all_passed'];inputs[rp]=sha(rp);ap=rp.parent/'actual-cleared-col76-and-commonQ-representation.json';assert sha(ap)==r['evidence_sha256'][ap.name];inputs[ap]=sha(ap);api=read(ap);C={z['row']:z for z in api['actual20_col76_numerator_values']};F={}
for i in(14,35):
 z=C[i];p=rp.parent/z['file'];assert sha(p)==z['sha256'];inputs[p]=sha(p);F[i]=poly(decode(p.read_bytes()[9:]))
expected={(1,0,0,2,0):Q(50,6075),(1,0,0,1,0):Q(20,6075),(1,0,0,0,0):Q(2,6075)};assert combine(F[14],F[35],6)==expected
sp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_UNIT_COL76_SOURCE_20261105_v64/report.json';assert sha(sp)=='349d1606feb479950eb690241bfdb95c49586731465b9f456e6656d20ea467cc';s=read(sp);inputs[sp]=sha(sp);bp=sp.parent/'col76-fixed-source-operands.json';assert sha(bp)==s['evidence_sha256'][bp.name];inputs[bp]=sha(bp);B=read(bp)['source_operands'];W={}
for i in(14,35):
 assert not B['U'+str(i)]['original_literal']['terms']and not B['V'+str(i)]['original_literal']['terms'];assert B['W'+str(i)]['QQ_transform_scalar']=='1'and B['W'+str(i)]['transform_monomial_exponent5']==[0]*5;W[i]=poly(B['W'+str(i)]['original_literal'])
expected_current={(1,0,0,1,0):Q(10,6075),(1,0,0,0,0):Q(2,6075)};assert combine(W[14],W[35],6)==expected_current
v={'all_passed':True,'source_only':True,'fixed_rows':[14,35],'candidate_row_operation':'row14 <- row14+6row35, determinant1','actual_numerator_identity':'F14+6F35=(2H/6075)*(5q0+1)^2','exact_numerator_polynomial3':serial(expected),'independent_cached_W_identity':'Ui=Vi=0 for14and35, so actualS4_i76=Wi. CachedW14+6W35=(2H/6075)*(5q0+1).','exact_current_candidate_value2':serial(expected_current),'prospective_fifth_pivot_column':76,'prospective_fifth_pivot_value':'(2H/6075)*(5q0+1)','prospective_typed_inverse':'6075/(2H*(5q0+1))','no_added_divisor_over_current_Qchart':True,'Q_complement_staysOPEN':True,'fifth_pivot_NOT_applied':True,'other19_current_columns_NOT_updated':True,'current_full173_top107plus66_87extras_RHS_scope_UNCHANGED':True,'if_later_applied_shape_only':'Full172x19/top106+66,selected19square only; no such update isclaimed here.','P290onlyonce_and_no_Qcolumnoperation':True,'no_fieldGCD_ranksearch_newpoint_prime_or_otherrowsearch':True}
(out/'two-exact-identities-and-candidate-existing-chart-unit.json').write_bytes(enc(v));report={'id':ID,'all_passed':True,'source_only':True,'numerator_Qsquare_identity_checked':True,'independent_cached_W_unit_identity_checked':True,'candidate_requires_no_new_chart_factor':True,'no_fifth_pivot_or_other19_updates':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
