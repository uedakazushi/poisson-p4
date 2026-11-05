#!/usr/bin/env python3
"""Source-only complete cached row comparison, anchored by an actual monomial."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,signal,resource
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def stream(blob):
 at=0
 def u():
  nonlocal at
  n=s=0
  while True:
   assert at<len(blob);v=blob[at];at+=1;n|=(v&127)<<s
   if not v&128:return n
   s+=7
 def z():
  nonlocal at
  tag=u();n=tag//2;assert at+n<=len(blob);v=int.from_bytes(blob[at:at+n],'big');at+=n;return-v if tag&1 else v
 q=Q(z(),z());shift=[]
 for _ in range(5):
  v=u();shift.append(v//2 if not v&1 else-(v//2)-1)
 n=u()
 for _ in range(n):
  c=q*z();e=tuple(shift[j]+u()for j in range(5));assert c;yield e,c
 assert at==len(blob)
def greater(e,f):return(sum(e),e)>(sum(f),f)
begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/report.json';assert sha(rp)=='a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a';r=read(rp);ep=rp.parent/'actual-fixed24-entry-source-and-exponent-plan.json';assert sha(ep)==r['evidence_sha256'][ep.name];v=read(ep);entries={(z['row'],z['column']):z for z in v['actual_entry_sources576']if z['row']in(11,44)};cols=v['columns24'];del v
inputs={rp:sha(rp),ep:sha(ep),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')};owned=25000;maxdiff=0
for z in entries.values():
 if z['terms']:
  p=ROOT/z['source']['file'];h=r['input_sha256'][str(p.relative_to(ROOT))];assert sha(p)==h;inputs[p]=h
def get(z):
 if not z['terms']:return b''
 ref=z['source'];p=ROOT/ref['file']
 if ref['kind']=='unchanged20_column_literal':
  with p.open('rb')as f:f.seek(ref['offset']);b=f.read(ref['bytes'])
  assert hashlib.sha256(b).hexdigest()==ref['payload_sha256'];return b
 b=p.read_bytes();assert sha(p)==ref['sha256']and b[:9]==b'Q5CELL1\n\x01';return b[9:]
anchor=49;assert entries[11,anchor]['terms']==entries[44,anchor]['terms']==1
a=list(stream(get(entries[11,anchor])));b=list(stream(get(entries[44,anchor])));scale=a[0][1]/b[0][1];shift=tuple(a[0][0][j]-b[0][0][j]for j in range(5));assert scale
checks=[]
for col in cols:
 za,zb=entries[11,col],entries[44,col];ba,bb=get(za),get(zb);ga=iter(stream(ba))if ba else iter(());gb=iter(stream(bb))if bb else iter(());x=next(ga,None);y=next(gb,None);diff=[];counta=countb=0
 while x is not None or y is not None:
  ey=tuple(y[0][j]+shift[j]for j in range(5))if y is not None else None
  if y is None or(x is not None and greater(x[0],ey)):
   e,c=x;counta+=1;x=next(ga,None)
  elif x is None or greater(ey,x[0]):
   e,c=ey,-scale*y[1];countb+=1;y=next(gb,None)
  else:
   assert x[0]==ey;e=x[0];c=x[1]-scale*y[1];counta+=1;countb+=1;x=next(ga,None);y=next(gb,None)
  if c:diff.append([str(c),*e]);assert len(diff)<=150000
 assert counta==za['terms']and countb==zb['terms'];owned=max(owned,25000+za['terms']+zb['terms']+6*len(diff)+2000);assert owned<=150000
 checks.append({'column':col,'row11_terms':counta,'row44_terms':countb,'difference_terms':diff,'difference_terms_count':len(diff),'all_source_terms_consumed':True});maxdiff=max(maxdiff,len(diff))
nonzero=[z for z in checks if z['difference_terms_count']];dependent=not nonzero
witness={'all_passed':True,'source_only':True,'anchor_column':anchor,'anchor_both_actual_monomials':True,'anchor_row11':[[str(c),*e]for e,c in a],'anchor_row44':[[str(c),*e]for e,c in b],'forced_RF_row_ratio_QQ':str(scale),'forced_RF_row_ratio_monomial_exponent5':shift,'all24_complete_actual_row_comparisons':checks,'all24_source_polynomial_identities':'row11 = QQscale*parameter_monomial*row44 + saved_difference','literal_RF_dependent':dependent,'RF_independent_if_nonzero_difference':not dependent,'reason':'Any RF proportionality ratio is forced by the nonzero actual monomial anchor49. A saved nonzero difference excludes that ratio exactly. No field gcd/product or new coefficients queried.','difference_nonzero_columns':[z['column']for z in nonzero],'single_monomial_difference_unit_on_old_parameter_torus':len(nonzero)==1 and len(nonzero[0]['difference_terms'])==1,'no_new_rowset_or_weight_search':True,'no_new_row_coefficient_worker_or_elimination':True,'all20_other_columns_and66_zero_rows_kept':True}
(out/'whole-cached-rows11-44-exact-comparison.json').write_bytes(enc(witness));raw=sum(p.stat().st_size for p in inputs)+sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert raw<32*1024*1024
report={'id':ID,'all_passed':True,'source_only':True,'complete_actual_columns_compared':24,'RF_dependent':dependent,'nonzero_difference_columns':witness['difference_nonzero_columns'],'single_monomial_difference_unit':witness['single_monomial_difference_unit_on_old_parameter_torus'],'owned_upper':owned,'RAW_before_small_report':raw,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
