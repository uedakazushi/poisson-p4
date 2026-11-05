#!/usr/bin/env python3
"""Positive prior using independently frozen exact Hstar support; no covectors."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_SPARSE_FIRST_PRIOR_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 inputs={}
 def read(p):
  inputs[p]={'path':p,'sha256':sha(ROOT/p)}
  return json.loads((ROOT/p).read_text())
 kp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_KNOWN_INPUT_GATE_20261105_v63/pruned-known619-lower4-input-and-prior-resource-gate.json'
 d=read(kp)
 hp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_HSTAR_FIRST_20261105_v63'
 hr=read(hp+'/report.json');assert hr['all_passed'] and hr['prior_bits']==17
 h=read(hp+'/actual-cached-selected99-Hstar-r-ell-and-prior17.json')
 assert sha(ROOT/(hp+'/actual-cached-selected99-Hstar-r-ell-and-prior17.json'))==hr['evidence_sha256']['actual-cached-selected99-Hstar-r-ell-and-prior17.json']
 H=[list(map(sp,M)) for M in h['actual_Hstar_r_ell_columns99']]
 assert all(set(v)<=set([0]) for M in H for v in M)
 norm={p:sp(v) for p,v in d['B99_normalized_basis99']};rep={p:sp(v) for p,v in d['B99_normalized_selected99_representations']}
 maximum=0;stage={}
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def plus(a,b):
  if not a[1]:return mark(b)
  if not b[1]:return mark(a)
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def times(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def accum(v,k,b):
  if b[1]:v[k]=plus(v.get(k,(1,0)),b)
 gn={p:sp(v) for p,v in d['F25_normalized_basis']};gr={}
 for receipt in d['F25_fixed_receipts']:
  a={receipt['basis_column']:(1,1)}
  for p,c in receipt['subtract']:
   for j,v in gr[p].items():accum(a,j,times(ex(c),v))
  gr[receipt['pivot']]={j:times(v,ex(1/Q(receipt['divide_by']))) for j,v in a.items()}
 FF=[]
 for i in range(27):
  q={i:(1,1)};a={}
  for p in sorted(gn):
   x=q.pop(p,(1,0))
   if not x[1]:continue
   for j,c in gn[p].items():
    if j!=p:accum(q,j,times(x,ex(c)))
   for j,c in gr[p].items():accum(a,j,times(x,c))
  FF.append(a)
 stage['F25_fixed_extract']=maximum
 raw=[];specific=[]
 for chi in [sp(v) for v in d['lower_duals2_rows3']]:
  ss=[]
  for M in H:
   hb=[]
   for v in M:
    b=(1,0)
    for i,c in chi.items():b=plus(b,times(ex(c),ex(v.get(i,Q(0)))))
    hb.append(b)
   rhs={}
   for p,repc in rep.items():
    b=(1,0)
    for j,c in repc.items():b=plus(b,times(hb[j],ex(c)))
    rhs[p]=b
   s={}
   for p in reversed(sorted(norm)):
    b=rhs[p]
    for j,c in norm[p].items():
     if j!=p:b=plus(b,times(ex(c),s.get(j,(1,0))))
    if b[1]:s[p]=b
   for p,col in norm.items():
    b=(1,0)
    for j,c in col.items():b=plus(b,times(ex(c),s.get(j,(1,0))))
   ss.append(s)
  specific.append(ss)
  rr={9*i+7:b for i,b in ss[0].items()};re={9*i+7:b for i,b in ss[1].items()}
  for i,b in ss[0].items():accum(re,9*i+8,b)
  raw.extend([rr,re])
 stage['specific_row_recurrence_and_membership']=maximum
 assert not raw[0] and not raw[1]
 J=list(map(sp,d['J25_full918_columns']));Pi=list(map(sp,d['PiY0_source102_columns']));alpha=[]
 for L in raw:
  w=[]
  for col in J:
   b=(1,0)
   for i,c in col.items():b=plus(b,times(ex(c),L.get(i,(1,0))))
   w.append(b)
  projected={}
  for i,col in enumerate(FF):
   for j,b in col.items():accum(projected,i,times(w[j],b))
  a=dict(L)
  for row,col in enumerate(Pi):
   for target,c in col.items():
    for var in range(9):accum(a,9*row+var,times(ex(c),projected.get(9*target+var,(1,0))))
  alpha.append(a)
 stage['full918_cleared_covector_prior']=maximum
 for v in map(sp,d['J25_projected27_columns']):
  for out in range(25):
   b=(1,0)
   for i,c in v.items():b=plus(b,times(ex(c),FF[i].get(out,(1,0))))
 stage['F25_fixed_membership']=maximum
 numeric=count(d)+count([H,norm,rep,gn,gr,FF,raw,specific,alpha])+12000
 assert numeric<150000
 plane=sum(sum(i%9 in (7,8) for i,c in v) for v in d['G2linear_source102_columns'])
 future=count(d)+2*(sum(map(len,FF))+sum(len(c) for M in H for c in M)+4*918+3*plane)+20000
 rawbound=len(json.dumps(d).encode())+(4*918+3*plane+1000)*(2*len(str(2**maximum))+40)+2_000_000
 ready=maximum<=384 and future<150000 and rawbound<32*1024*1024
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 ep=out/'sparse-Hstar-specific-row-positive-prior.json'
 ev={'all_lowerA_raw_covectors_structurally_zero':True,'exact_Hstar_nonzero_bottom_rows':[0],'stage_positive_bound_bits':stage,'prior_bits':maximum,'declared_commit_cap_bits':384,'prior_guard_covered':ready,'next_known_source_numeric_bound':future,'next_known_source_RAW_bound':rawbound,'prior_owned_bound':numeric,'no_signed_full_covectors_or_619_values':True,'specific_row_support_upper':[[len(v) for v in M] for M in specific],'alpha_sparse_support_upper':[len(v) for v in alpha]}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'prior_guard_covered':ready,'positive_prior_bits':maximum,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'no_signed_covectors_or_values_evaluated':True,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'ready':ready,'bits':maximum,'stages':stage,'numeric':future,'RAW':rawbound,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
