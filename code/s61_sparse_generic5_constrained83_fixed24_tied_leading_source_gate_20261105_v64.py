#!/usr/bin/env python3
"""Source-only tied-leading QQ block for the SAME fixed24, no determinant."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,factorial
import json,hashlib,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_TIED_LEADING_SOURCE_GATE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_ASSIGNMENT_20261105_v64/report.json';assert sha(rp)=='003b4201b9307abd8366140508f90eba17777df3c8827f8b05da555753e5fea4';r=read(rp);assert r['all_passed']and not r['unique_minimum_matching'];ep=rp.parent/'exact-fixed24-assignment-and-uniqueness.json';assert sha(ep)==r['evidence_sha256'][ep.name];v=read(ep)
L=[[Q(0)for _ in range(24)]for _ in range(24)];source=[]
for z in v['all_actual_leading_entries']:
 i,j=z['i'],z['j'];rc=v['all576_reduced_costs'][i][j];assert rc>=0
 if rc==0:L[i][j]=Q(z['QQ_coefficient']);source.append({'i':i,'j':j,'source':z['source'],'source_term_index':z['source_term_index'],'exponent5':z['exponent5'],'QQ_coefficient':z['QQ_coefficient']})
D=[];Z=[]
for row in L:
 d=1
 for q in row:d=d*q.denominator//gcd(d,q.denominator)
 D.append(str(d));Z.append([str((q*d).numerator)for q in row]);assert all((q*d).denominator==1 for q in row)
B=max(abs(int(z)).bit_length()for row in Z for z in row);detbits=24*B+factorial(24).bit_length()+2;prod_bits=2*detbits+2
values={'all_passed':True,'source_only':True,'same_fixed_rows24':v['fixed_rows24'],'same_fixed_columns24':v['fixed_columns24'],'same_matching24':v['matching24'],'same_Kronecker_weights5':v['weights5'],'row_LCM24':D,'integer_cleared_leading_matrix24':Z,'QQ_leading_matrix24':[[str(q)for q in row]for row in L],'all_zero_reduced_cost_sources':source,'zero_cost_nonzero_entries':len(source),'max_integer_entry_bits':B,'Hadamard_all_minor_integer_bit_upper':detbits,'Bareiss_predivision_product_bit_upper':prod_bits,'fixed24_det_or_rank_not_computed':True,'no_other_rowset_or_weight_search':True,'future_sufficient_condition':'If the fixed integer matrix has nonzero determinant, detL is nonzero and is the exact lowest-weight coefficient of the actual24det after the certified row/column cost normalization. Hence actual24det is nonzero even though minimum matchings are not unique.','if_fixed_det_zero':'Inconclusive for the actual determinant; higher weighted terms were not evaluated.','future_inverse_plan':'If detL nonzero, use exact constant L inverse to precondition the row/column t0-power normalized actual K5 matrix. Its formal leading matrix isI, so all fixed actualGaussian pivots are nonzero. Every actual entry remains in the finite RF program; no series truncation or numerical beta point.','budget_only_no_arithmetic_GO':True}
(out/'same-fixed24-tied-leading-source.json').write_bytes(enc(values));inputs={rp:sha(rp),ep:sha(ep),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')};base=sum(p.stat().st_size for p in inputs);report={'id':ID,'all_passed':True,'source_only':True,'same_fixed24':True,'zero_cost_entries':len(source),'max_integer_entry_bits':B,'all_minor_bit_upper':detbits,'proposed_fixed_det_worker_owned_upper':25000+24*24*8,'proposed_fixed_det_worker_RAW_upper':base+2*1024*1024,'proposed_fixed_det_worker_seconds':60,'no_det_rank_or_other_matching_run':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};assert report['proposed_fixed_det_worker_RAW_upper']<32*1024*1024;(out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
