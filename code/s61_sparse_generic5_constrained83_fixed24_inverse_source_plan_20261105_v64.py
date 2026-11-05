#!/usr/bin/env python3
"""Source-only fixed24 response availability and exponent/precision plan; no rank."""
from pathlib import Path
import json,hashlib,signal,time,resource
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def bounds(blob):
 # Primitive payload is already exact-roundtrip certified. Skip all integer
 # coefficient bytes without expanding QQ coefficients; read exponents only.
 at=0
 def u():
  nonlocal at
  n=s=0
  while True:
   assert at<len(blob);v=blob[at];at+=1;n|=(v&127)<<s
   if not v&128:return n
   s+=7
 def skip_i():
  nonlocal at
  tag=u();n=tag//2;assert at+n<=len(blob);at+=n
 skip_i();skip_i();shift=[]
 for _ in range(5):
  z=u();shift.append(z//2 if not z&1 else-(z//2)-1)
 n=u();assert n>0;lo=hi=None
 for _ in range(n):
  skip_i();e=[shift[j]+u()for j in range(5)]
  if lo is None:lo=e[:];hi=e[:]
  else:lo=[min(x,y)for x,y in zip(lo,e)];hi=[max(x,y)for x,y in zip(hi,e)]
 assert at==len(blob);return n,lo,hi
begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')}
def gate(stem,h):
 p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
pr,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_TOP_P290_BATCH7_COEFFICIENTS_20261105_v64','0fee4ac3526a15b5b322523f21f44ff16e1bd54def55eeb53d3f71c3a28627fd');ap=pr.parent/'evaluated16-cells-and-full177-frame.json';assert sha(ap)==r['evidence_sha256'][ap.name];inputs[ap]=sha(ap);api=read(ap);cells=api['all_evaluated_cells'];rows=list(dict.fromkeys(z['row']for z in cells));assert len(rows)==24 and len(cells)==96
fp,f=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FACTORED_CONSTRAINT_SCHUR_PREFLIGHT_20261105_v64','8d0c175945ac170d8275da8fb2ecd1b7673d40dfd6c3bdc035f54d58b75912bf');op=fp.parent/'actual-factored177x24-and-top111x24.json';assert sha(op)==f['evidence_sha256'][op.name];inputs[op]=sha(op);frame=read(op);cols=frame['remaining24_columns'];assert len(cols)==24
sp,s=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM49_ACTIVE_FRAGMENT_MIXED_BINARY_COMPACTION_20261105_v64','62a95473e52e2a85054586bd63e28b79ae58850bebf49d6e5a0a4daeb69e7970');ip=sp.parent/s['current_index'];mp=sp.parent/s['current_snapshot'];assert sha(ip)==s['current_index_sha256']and sha(mp)==s['current_snapshot_sha256'];inputs[ip]=sha(ip);inputs[mp]=sha(mp);ix=read(ip)
peak=numeric(ix)+numeric(api)+numeric(frame)+10000;assert peak<150000
old={(z['row'],z['column']):z for z in ix['all_literal_cells']if z['row']in rows and z['column']in cols and z['column']not in(7,13,27,33)};del ix,frame
new={(z['row'],z['column']):z for z in cells};entries=[];logical=0
with mp.open('rb')as src:
 for row in rows:
  for col in cols:
   if col in(7,13,27,33):
    z=new[row,col];p=ROOT/z['literal_cell_file'];assert sha(p)==z['sha256'];inputs[p]=sha(p);blob=p.read_bytes();assert blob[:8]==b'Q5CELL1\n'
    ref={'kind':'already_computed_P290_scaled_free','file':str(p.relative_to(ROOT)),'sha256':sha(p)}
    if blob[8]==0:assert len(blob)==9;n=0;lo=hi=None
    else:assert blob[8]==1;n,lo,hi=bounds(blob[9:]);assert n==z['terms']
   elif(row,col)in old:
    z=old[row,col];src.seek(z['payload_offset']);blob=src.read(z['payload_bytes']);assert hashlib.sha256(blob).hexdigest()==z['payload_sha256'];n,lo,hi=bounds(blob);assert n==z['terms'];ref={'kind':'unchanged20_column_literal','file':str(mp.relative_to(ROOT)),'payload_sha256':z['payload_sha256'],'offset':z['payload_offset'],'bytes':z['payload_bytes']}
   else:n=0;lo=hi=None;ref={'kind':'exact_source_zero_missing_literal_cell'}
   logical+=n;entries.append({'row':row,'column':col,'terms':n,'lower5':lo,'upper5':hi,'source':ref});peak=max(peak,numeric(entries)+numeric(old)+numeric(new)+numeric(api)+8192);assert peak<150000
nonzero=[z for z in entries if z['terms']];lo=[min(z['lower5'][j]for z in nonzero)for j in range(5)];hi=[max(z['upper5'][j]for z in nonzero)for j in range(5)];width=[24*(hi[j]-lo[j])+1 for j in range(5)];weights=[1]*5
for j in range(3,-1,-1):weights[j]=weights[j+1]*width[j+1]
zero_cols=[c for c in cols if not any(z['terms']for z in entries if z['column']==c)];zero_rows=[r for r in rows if not any(z['terms']for z in entries if z['row']==r)]
ev={'all_passed':True,'source_only':True,'rows24_fixed_before_any_rank_or_matching_test':rows,'columns24':cols,'actual_entry_sources576':entries,'logical_terms_not_resident':logical,'global_entry_exponent_lower5':lo,'global_entry_exponent_upper5':hi,'determinant_exponent_box_lower5':[24*x for x in lo],'determinant_exponent_box_upper5':[24*x for x in hi],'determinant_box_width5':width,'injective_Kronecker_weights5':weights,'formal_weight_integer_bits':max(x.bit_length()for x in weights),'structural_zero_columns':zero_cols,'structural_zero_rows':zero_rows,'unique_minimum_matching_not_tested':True,'no_rank_minor_determinant_Gaussian_or_series_worker':True,'all_actual_entries_remain_K5_functions':True,'future_nonzero_criterion':'If this fixed weighted24x24 entry-leading-cost graph has a unique minimum perfect matching, its determinant has a unique lowest-weight monomial with nonzero product coefficient. Hence determinant is nonzero overK5. This condition has NOT been checked.','conditional_inverse_instructions':'A unique matching supplies row/column order of a triangular leading matrix via the acyclic zero-reduced-cost equality graph. Formal Gaussian operations then retain every actual K5 entry; never replace an entry by its leading term. A new bounded source gate and arithmetic authorization are required.','rank24_is_only_existence_not_this_row_basis':True,'all111_evaluation_not_necessary_if_fixed24_criterion_later_passes':True,'full177_cokernel':'Full177x24 is retained; coker(full)=coker(top111x24) direct_sumK5^66. A future left inverse proves column injectivity, not deletion of66free row summands.','chart':'Old parameter monomial units*P290*detA*detB; any future nonzero24minor adds its own explicit localization. Complements remain open.'}
(out/'actual-fixed24-entry-source-and-exponent-plan.json').write_bytes(enc(ev));ledger=[{'path':str(p.relative_to(ROOT)),'sha256':h,'bytes':p.stat().st_size}for p,h in sorted(inputs.items())];raw=sum(z['bytes']for z in ledger)+sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert raw<32*1024*1024
report={'id':ID,'all_passed':True,'source_only':True,'fixed24_rows_available':True,'576_source_entries_checked':True,'logical_terms':logical,'owned_metadata_upper':peak,'fresh_RAW_bytes_before_small_report':raw,'no_fixed24_basis_or_inverse_claim':True,'no_new_coefficient_products':True,'input_sha256':{z['path']:z['sha256']for z in ledger},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
