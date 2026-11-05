#!/usr/bin/env python3
"""Pure source corollary of exact row comparison and singular leading face."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_SOURCE_FRONTIER_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')}
def gate(stem,h):
 p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
p,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64','0ca442b4eb828e535bdb527357826c6f4cd7edb48f23735f08742b798ada9371');ep=p.parent/'whole-cached-rows11-44-exact-comparison.json';assert sha(ep)==r['evidence_sha256'][ep.name];inputs[ep]=sha(ep);v=read(ep);assert v['forced_RF_row_ratio_QQ']=='8'and v['forced_RF_row_ratio_monomial_exponent5']==[0]*5
p,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_TIED_EXACT_20261105_v64','fa0d233103696be93b8d8012186ba84dc623ad8faa76a83690dd18560109c0fb');assert r['rank_of_this_fixed_leading_face']==8 and r['kernel_dimension']==16
checks={z['column']:z for z in v['all24_complete_actual_row_comparisons']};anchor=v['anchor_row44'];assert len(anchor)==1 and not checks[49]['difference_terms'];records=[]
for col in(51,81):
 d=checks[col]['difference_terms'];assert len(d)==1;coef=-Q(anchor[0][0])*Q(d[0][0]);exp=[anchor[0][j+1]+d[0][j+1]for j in range(5)];assert coef==Q(1,27337500)and exp==[1,0,1,1,1];records.append({'rows':[11,44],'columns':[49,col],'exact_determinant_QQ_scalar':str(coef),'parameter_exponent5':exp,'identity':'det(row11,row44) = det(row11-8*row44,row44) = -delta_col*row44_col49','no_polynomial_product_expansion':True})
ev={'all_passed':True,'source_only':True,'rank_of_fixed_leading_face8_not_actual_rank':True,'two_actual_source_monomial_minors':records,'both_nonzero_in_K5':True,'unit_chart':'H*s0*q0*t0 nonzero, together with old baseunits/P290/detA/detB. Explicitly keep any added coordinate divisor in the complement obligations.','actual_cached_rows_RF_independent':True,'possible_next_method':'Use these fixed source pivots to define a two-unit Schur on the SAME fixed24 block. Every actual entry is already cached; first derive whole support/RAW bounds and keep remaining177-row compatibility. No coefficient update, alternate-row search or field pivot has been executed here.','not_complete_inverse_or_exhaustion':True}
(out/'fixed24-leading-obstruction-and-actual-two-unit-source.json').write_bytes(enc(ev));report={'id':ID,'all_passed':True,'source_only':True,'leading_rank8_actual_rank_unresolved':True,'actual_two_unit_source_minors':2,'no_new_coefficient_worker_or_rowset_search':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
