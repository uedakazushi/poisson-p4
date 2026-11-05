#!/usr/bin/env python3
"""ONE newpurpose N6 specializedcallback; no genericzero inference or all6job."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_CLEANUP_COLUMN6_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
class GuardStop(Exception):pass
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={};prefix={};literal={};acc={};pairs=0;peak=0;current=None;completed=[]
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='1fc2e57cab9d41ee866d0b52f0fcbc59992fcfb2b13139d98213ee495c02569e';cfg=src(base+'/fresh-single-column-profiles.json')['selected_cheapest_single_column'];assert cfg['column']==6 and cfg['owned_upper']==46999;src('reports/S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_LIFT_GUARD_20261105_v63/report.json');peer=src('reports/S61_NONSQUARE_K92_WHOLE_ROW_LEADING_PEER_20261105_v63_REVIEW/report.json');assert peer['all_passed'];receipt=src('reports/S61_NONSQUARE_PREPOWER_K92_CLEANUP_COLUMN6_WORKER_GATE_20261105_v63/report.json');assert receipt['all_passed']
 rel=cfg['prefix5_polynomials_source'];p=ROOT/rel;pins[rel]=sha(p);text=p.read_text();decoder=json.JSONDecoder();pos=text.index(':',text.index('"uniform_integer_denominator"'))+1;L=int(decoder.raw_decode(text,pos+len(text[pos:])-len(text[pos:].lstrip()))[0]);wanted={'delta','w0','w1','vp0','vp1'};pos=text.index('{',text.index('"all_prefix_integer_polynomials"'))+1
 while True:
  while text[pos].isspace()or text[pos]==',':pos+=1
  if text[pos]=='}':break
  key,pos=decoder.raw_decode(text,pos)
  while text[pos].isspace():pos+=1
  assert text[pos]==':';pos+=1
  while text[pos].isspace():pos+=1
  value,pos=decoder.raw_decode(text,pos)
  if key in wanted:prefix[key]={int(g):int(v)for g,v in value}
  load_owned=sum(map(len,prefix.values()))+2*len(value)+15000;peak=max(peak,load_owned);assert load_owned<=46999
  del value
 del text
 assert set(prefix)==wanted and sum(map(len,prefix.values()))==9350
 inputs=[]
 for rec in cfg['literal_inputs']:
  p=ROOT/rec['path'];pins[rec['path']]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    if ln!=rec['line']:continue
    e=json.loads(line);assert e['entry_index']==rec['entry_index'];c=Q(e['QQ_content']);fac=c.numerator*(L//c.denominator);assert L%c.denominator==0 and fac==int(rec['cleared_factor']);name=f"{rec['block']}{rec['row']}";literal[name]={int(g):int(v)*fac for g,v in e['sparse_integer_polynomial']};inputs.append(dict(name=name,source=rec));break
 assert sum(map(len,literal.values()))==102
 dd=dict(literal['D2'])
 for g,x in literal['D0'].items():dd[g]=dd.get(g,0)-59375*x
 dd={g:x for g,x in dd.items()if x};literal['dd']=dd
 bound=int(cfg['integer_abs_upper']);terms=[('delta','dd',1),('w0','B0',-1),('w1','B1',-1),('vp0','B2',-1),('vp1','B3',-1)]
 def guard(extra=0):
  nonlocal peak
  n=sum(map(len,prefix.values()))+sum(map(len,literal.values()))+len(acc)+extra+15000;peak=max(peak,n)
  if n>46999:raise GuardStop('owned_guard')
  if time.monotonic()-start>45:raise GuardStop('45s_guard')
 try:
  for ti,(left,right,sgn)in enumerate(terms):
   pp=prefix[left];qq=literal[right];xs=sorted(pp);ys=sorted(qq);tmp={};current=dict(term_index=ti,left=left,right=right,sign=sgn,lhs_index=0,rhs_index=0,partial_product=tmp,next_pair_NOTcommitted=True)
   for ii,g in enumerate(xs):
    for jj,h in enumerate(ys):
     if pairs%2048==0:current['lhs_index']=ii;current['rhs_index']=jj;guard(len(tmp)+2)
     v=tmp.get(g+h,0)+pp[g]*qq[h];assert abs(v)<=bound
     if v:tmp[g+h]=v
     else:tmp.pop(g+h,None)
     pairs+=1
   for g,x in tmp.items():
    v=acc.get(g,0)+sgn*x;assert abs(v)<=bound
    if v:acc[g]=v
    else:acc.pop(g,None)
   guard(2*len(tmp)+2);termsha=hashlib.sha256(json.dumps(sorted(tmp.items()),separators=(',',':')).encode()).hexdigest();completed.append(dict(term_index=ti,left=left,right=right,sign=sgn,exact_product_terms=len(tmp),exact_product_canonical_sha256=termsha));current=dict(completed_terms=len(completed),next_term_index=ti+1);guard(len(tmp)+2);del tmp
 except GuardStop as err:
  signal.alarm(0);dump(out/'state-K92C6P001.json',dict(version='K92C6P001',reason=str(err),current=current,accumulator=sorted(acc.items()),completed=completed,computed_dd=sorted(dd.items()),source_pins=pins,
   immutable_prefix_and_literal_coefficients_reload_from_pins=True,no_prefix_recalculation=True,no_automatic_retry=True));dump(out/'report.json',dict(id=ID,all_passed=False,guard=str(err),state_sha256=sha(out/'state-K92C6P001.json'),owned_peak=peak,no_callback_result_claim=True));return
 assert len(acc)<=cfg['output_support_upper'];dump(out/'single-column6-specialization.json',dict(column=6,integer_polynomial=sorted(acc.items()),common_denominator=str(L**5),formal_specialization_ZERO=not acc,all5products=completed,literal_inputs=inputs,
  generic_root='typed-generic-rowcleanup-and-pivot.json N6',generic_root_RETAINED=True,no_genericzero_inference_from_formalzero=True,if_nonzero='positive specialization provesgenericRFnonzero',generic_deltaP_deltaQ_N0_chart_complement_RETAINED=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new<cfg['new_RAW_upper']and new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,ONE_newpurpose_singlecolumn6=True,all5_actual_prefix_and6literal_sources_bound=True,all5_signed_products_verified=True,formal_specialization_ZERO=not acc,integer_coefficient_terms=len(acc),
  generic_N6_callback_RETAINED=True,no_generic_RFzero_or_det11_nonzero_or_global57_claim=True,no_othercolumn_newweight_window_point_prime_or_prefix_recomputation=True,owned_scalar_peak=peak,integer_pairs=pairs,
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
