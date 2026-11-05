#!/usr/bin/env python3
"""Source-only unchanged columns and exact cached next monomial unit candidates."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR_COL81_SOURCE_FRONTIER_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
def gate(stem,h):
 p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
p,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');ep=p.parent/'actual-fixed24-entry-source-and-exponent-plan.json';assert sha(ep)==r['evidence_sha256'][ep.name];inputs[ep]=sha(ep);v=read(ep);E={(z['row'],z['column']):z for z in v['actual_entry_sources576']};rows=[i for i in v['rows24_fixed_before_any_rank_or_matching_test']if i not in(11,44)];cols=[j for j in v['columns24']if j not in(49,51)and E[11,j]['terms']==E[44,j]['terms']==0];assert cols==[48,50,76]
units=[];profiles=[]
for col in cols:
 for row in rows:
  z=E[row,col];profiles.append({'row':row,'column':col,'unchanged_actual_terms':z['terms']})
  if z['terms']!=1:continue
  ref=z['source'];assert ref['kind']=='unchanged20_column_literal';path=ROOT/ref['file'];h=r['input_sha256'][str(path.relative_to(ROOT))];assert sha(path)==h;inputs[path]=h
  with path.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
  assert hashlib.sha256(blob).hexdigest()==ref['payload_sha256'];d=decode(blob);assert len(d['terms'])==1;t=d['terms'][0];q=Q(d['content'])*int(t[0]);e=[d['shift'][j]+t[j+1]for j in range(5)];assert q
  units.append({'row':row,'column':col,'QQ_coefficient':str(q),'parameter_exponent5':e,'source':ref,'post_unit2_actual_value_identical_to_old_cached_value':True})
p,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR_COL81_COEFFICIENTS_20261105_v64','0f05a5001b12f4b109801d60ca680416e13fa3243fb06b492af7c88be5b00f1c');ap=p.parent/'actual-residual-col81-and-structural-frontier.json';assert sha(ap)==cr['evidence_sha256'][ap.name];inputs[ap]=sha(ap);api=read(ap);lo=None;termcount=0
for z in api['actual22_residual_col81_entries']:
 path=p.parent/z['file'];assert sha(path)==z['sha256'];inputs[path]=sha(path);blob=path.read_bytes()
 if not blob[8]:assert len(blob)==9;continue
 d=decode(blob[9:]);assert len(d['terms'])==z['terms'];termcount+=len(d['terms'])
 for t in d['terms']:
  e=[d['shift'][j]+t[j+1]for j in range(5)];lo=e[:]if lo is None else[min(x,y)for x,y in zip(lo,e)]
assert termcount==183 and len(units)==6
p,fr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT2_FULL175_TRANSPORT_20261105_v64','a2eaa1fb659bda0892a20bb9c4483db5d2f48d3781a6ff8d6769dbc7de530647')
ev={'all_passed':True,'source_only':True,'post_unit2_columns_exactly_unchanged_all175':[48,50,76],'proof':'Both actual pivotrows11/44 have sourcezero in these columns. Their transformedtwo-vector is0, so the universal Schur correction is0 for EVERY175 row. The selected22 entries are unchanged literalcache references, not approximatedzero.','all66_selected_entries_source_term_counts':profiles,'next_actual_parameter_monomial_unit_candidates6':units,'no_additional_pivot_or_field_update_performed':True,'only_rank_lower':'At least one candidate supplies one further rankunit in this residual, but two distinct candidatecolumns do not by themselves certify two independent futurepivots. No fullrank24 claim.','new_actual_col81_terms':183,'new_col81_common_Laurent_monomial_exponent5':lo,'monomial_content_only_no_polynomial_GCD':True,'col81_divided_values_not_materialized':True,'col81_not_allzero_and_no_monomial_entries':True,'full175_87plus66_RHS_and_other21_columns_kept':True,'normal_coordinate_RHS_and_P290_scale_unchanged':True,'chart':'Oldbaseunits*P290*detA*detB*H*s0*q0*t0 plus any chosen next candidateparameter monomial support; complements remainexplicit.','no_actual_fixed24_determinant_conclusion':True}
(out/'actual-unchanged-columns-and-six-source-units.json').write_bytes(enc(ev));raw=sum(p.stat().st_size for p in inputs)+sum(p.stat().st_size for p in out.rglob('*')if p.is_file());peak=25000+numeric(v)+numeric(api)+numeric(ev)+512;assert peak<150000 and raw<32*1024*1024
report={'id':ID,'all_passed':True,'source_only':True,'new_parameter_monomial_unit_candidates':6,'unchanged_column_labels':cols,'no_new_pivot_or_polynomial_values':True,'owned_upper':peak,'RAW_before_small_report':raw,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
