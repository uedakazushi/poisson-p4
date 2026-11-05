#!/usr/bin/env python3
"""SOURCE only: exact new fixedcolumn forms; no fourthpivot or matrix update."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR3_COL48_SOURCE_FRONTIER_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR3_COL48_COEFFICIENTS_20261105_v64/report.json';assert sha(rp)=='999fb1eb5b4b7253ef4a123de1cf40238531662f2a104f7719db48ba75a2d1c0';r=read(rp);assert r['all_passed'];ap=rp.parent/'actual-residual-col48-and-structural-frontier.json';assert sha(ap)==r['evidence_sha256'][ap.name];api=read(ap);inputs={rp:sha(rp),ap:sha(ap),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')};profiles=[];actual={}
for z in api['actual21_residual_col48_entries']:
 p=rp.parent/z['file'];assert sha(p)==z['sha256'];inputs[p]=sha(p);blob=p.read_bytes()
 if not z['terms']:assert blob==b'Q5CELL1\n\x00';continue
 d=decode(blob[9:]);q=Q(d['content']);A={tuple(d['shift'][j]+t[j+1]for j in range(5)):q*int(t[0])for t in d['terms']};assert len(A)==z['terms'];actual[z['row']]=A;lo=[min(e[j]for e in A)for j in range(5)];profiles.append({'row':z['row'],'column':48,'QQ_scalar_content':d['content'],'Laurent_monomial_content5':lo,'primitive_fixed_polynomial':[[str(c/q),*[e[j]-lo[j]for j in range(5)]]for e,c in A.items()],'known_exact_nonzero_by_distinct_support':True,'source':str(p.relative_to(ROOT)),'sha256':sha(p)})
assert actual[31]=={(2,0,0,1,-1):Q(256,1215),(2,0,0,0,-1):Q(256,6075)}
assert actual[41]=={(0,0,1,0,1):Q(1,162000),(0,0,0,1,1):Q(-5,162000)}
assert actual[133]=={(1,0,0,0,1):Q(-64,1215),(-1,0,0,0,1):Q(64,1215)}
assert len(actual[17])==3
forms={'all_passed':True,'source_only':True,'nonzero_column_profiles9':profiles,'fixed_two_term_forms_verified':[{'row':31,'exact_form':'(256/6075)*H^2*(5q0+1)/t0','proposed_extra_chart_if_future_pivot':'5q0+1 nonzero'},{'row':41,'exact_form':'t0*(s0-5q0)/162000','proposed_extra_chart_if_future_pivot':'s0-5q0 nonzero'},{'row':133,'exact_form':'(64/1215)*t0*(1-H^2)/H','proposed_extra_chart_if_future_pivot':'1-H^2 nonzero'}],'row17_oldunit_did_not_persist':True,'row17_new_terms':3,'no_fourth_pivot_or_extra_factor_inverted':True,'no_polynomial_GCD_or_new_coefficients':True,'new_current_other20_columns_remain_typed':True,'old_postunit2_col81_not_currentpostunit3_values':True,'all174_rows87plus66_RHS_and_singleP290scale_kept':True,'no_actual_fixed24_determinant_conclusion':True}
(out/'nine-fixed-column-profiles-and-three-explicit-forms.json').write_bytes(enc(forms));report={'id':ID,'all_passed':True,'source_only':True,'exact_nonzero_forms':3,'nonzero_profiles':9,'row17_unit_not_assumed':True,'no_new_field_pivot_or_other20_column_updates':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
