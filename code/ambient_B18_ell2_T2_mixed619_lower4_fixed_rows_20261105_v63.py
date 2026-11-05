#!/usr/bin/env python3
"""Stage only fixed known rows under prior376; no full covectors or X2."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_FIXED_ROWS_20261105_v63'
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
 assert maximum==376<=384
 prior=maximum;actualmax=0
 def guard(x):
  nonlocal actualmax
  actualmax=max(actualmax,abs(x.numerator).bit_length(),x.denominator.bit_length());assert actualmax<=prior
  return x
 def axpy(v,w,c):
  a=dict(v)
  for k,x in w.items():
   y=guard(a.get(k,Q(0))+guard(c*x))
   if y:a[k]=y
   else:a.pop(k,None)
  return a
 # All arithmetic below has the positive376-bit guard established above.
 ag={}
 for receipt in d['F25_fixed_receipts']:
  a={receipt['basis_column']:Q(1)}
  for p,c in receipt['subtract']:a=axpy(a,ag[p],-Q(c))
  ag[receipt['pivot']]={j:guard(v/Q(receipt['divide_by'])) for j,v in a.items()}
 AF=[]
 for i in range(27):
  q={i:Q(1)};a={}
  for p in sorted(gn):
   x=q.get(p,Q(0))
   if x:q=axpy(q,gn[p],-x);a=axpy(a,ag[p],x)
  assert set(q)<=set([25,26])
  for j,x in a.items():
   D,U=FF[i][j];assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  AF.append(a)
 AS=[]
 for ci,chi in enumerate([sp(v) for v in d['lower_duals2_rows3']]):
  row=[]
  for vi,M in enumerate(H):
   hv=[sum((c*v.get(i,Q(0)) for i,c in chi.items()),Q(0)) for v in M]
   rhs={p:sum((hv[j]*c for j,c in rr.items()),Q(0)) for p,rr in rep.items()}
   ss={}
   for p in reversed(sorted(norm)):
    v=rhs[p]
    for j,c in norm[p].items():
     if j!=p:v=guard(v-guard(c*ss.get(j,Q(0))))
    if v:ss[p]=v
   assert not any(k in ss for k in [1,36,38])
   for p,rr in norm.items():assert sum((ss.get(j,Q(0))*c for j,c in rr.items()),Q(0))==rhs[p]
   for j,x in ss.items():
    D,U=specific[ci][vi][j];assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
   row.append(ss)
  AS.append(row)
 for j,v in enumerate(map(sp,d['J25_projected27_columns'])):
  z={}
  for i,c in v.items():z=axpy(z,AF[i],c)
  assert z=={j:Q(1)}
 # New independent positive clearing prior uses the completed known sparse
 # rows and fixedF25. It does not perform the signed918-slot clearing.
 maximum=0;actual_raw=[]
 for rr,re in AS:
  a={9*i+7:x for i,x in rr.items()};b={9*i+7:x for i,x in re.items()}
  for i,x in rr.items():b[9*i+8]=x
  actual_raw.extend([a,b])
 J=list(map(sp,d['J25_full918_columns']));Pi=list(map(sp,d['PiY0_source102_columns']));alpha=[]
 for L in actual_raw:
  w=[]
  for col in J:
   b=(1,0)
   for i,c in col.items():b=plus(b,times(ex(c),ex(L.get(i,Q(0)))))
   w.append(b)
  projected={}
  for i,col in enumerate(AF):
   for j,c in col.items():accum(projected,i,times(w[j],ex(c)))
  a={i:ex(c) for i,c in L.items()}
  for row,col in enumerate(Pi):
   for target,c in col.items():
    for var in range(9):accum(a,9*row+var,times(ex(c),projected.get(9*target+var,(1,0))))
  alpha.append(a)
 nextbits=maximum
 owned=count(d)+count([H,norm,rep,gn,gr,FF,specific,AF,AS,actual_raw,alpha])+12000
 assert owned<150000
 nextnumeric=count(d)+count([AF,AS])+2*(4*918)+20000
 rawbound=len(json.dumps(d).encode())+len(json.dumps([[[[j,str(x)] for j,x in v.items()] for v in AF],[[[[j,str(x)] for j,x in v.items()] for v in M] for M in AS]]).encode())+(4*918+1000)*(2*len(str(2**nextbits))+40)+2_000_000
 ready=nextbits<=384 and nextnumeric<150000 and rawbound<32*1024*1024
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-fixed-specific-dual-rows-and-F25-with-next-clearing-prior.json'
 ev={'actual_F25_columns27':[[[j,str(x)] for j,x in sorted(v.items())] for v in AF],'actual_specific_dual_rows2_by2':[[[[j,str(x)] for j,x in sorted(v.items())] for v in M] for M in AS],
 'completed_known_rows_prior_bits':prior,'completed_known_rows_actual_bits':actualmax,'next_clearing_prior_bits':nextbits,'next_clearing_prior_guard_covered':ready,'next_source_numeric_bound':nextnumeric,'next_source_RAW_bound':rawbound,'whole_current_owned_bound':owned,
 'no_signed_full_covectors_or_619_values':True,'row_proof':'s*normalizedBasis=chiH*selected99Representation ands vanishes onfreeE1,36,38; frozeninvertibility proves s=chiH B99 withoutconstructinguniversalB99.',
 'F25_proof':'All25 projectedsourcecolumnimagesareunitvectors; freequotient25,26retained.'}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'completed_known_row_prior_bits':prior,'actual_known_row_bits':actualmax,'next_clearing_prior_bits':nextbits,'next_clearing_ready':ready,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'no_signed_full_covectors_or_values':True,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'prior':prior,'actual':actualmax,'nextprior':nextbits,'ready':ready,'numeric':nextnumeric,'RAW':rawbound,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
