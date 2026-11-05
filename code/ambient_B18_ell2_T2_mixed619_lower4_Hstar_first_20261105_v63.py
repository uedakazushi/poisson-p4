#!/usr/bin/env python3
"""Separate cached-first source aggregation; no X2 or physical query."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,mv
from ambient_B18_ell2_T2_mixed619_lower4_known_input_20261105_v63 import count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_HSTAR_FIRST_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 ip=ROOT/'reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_KNOWN_INPUT_GATE_20261105_v63/pruned-known619-lower4-input-and-prior-resource-gate.json';d=json.loads(ip.read_text());P=list(map(sp,d['Pi3_source56_columns']));R=list(map(sp,d['original768_r_first_columns']));E=list(map(sp,d['original640_ell_first_columns']));L=list(map(sp,d['native_L0']));Se=list(map(sp,d['native_Se']));S0=list(map(sp,d['native_S0']));maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def times(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  z=math.lcm(a[0],b[0]);return mark((z,a[1]*(z//a[0])+b[1]*(z//b[0])))
 def apply(cols,v):
  a={}
  for j,b in v.items():
   for i,c in cols[j].items():a[i]=plus(a.get(i,(1,0)),times(ex(c),b))
  return a
 positive=[]
 for cd in d['selected99_original768_combinations']:
  c={j:ex(v) for j,v in cd};cn={j-640:b for j,b in c.items() if j>=640};apply(S0,cn)
  r=apply(P,apply(R,c));e=apply(E,{j:b for j,b in c.items() if j<640});n=apply(L,apply(Se,cn))
  for i,b in n.items():e[i]=plus(e.get(i,(1,0)),b)
  positive.append([r,apply(P,e)])
 assert maximum==17<384
 owned=count(d)+count([P,R,E,L,Se,S0,positive])+5000;assert owned<150000
 # Only now perform the separately17-bit-bounded known first aggregation.
 actual=[[],[]];bits=0
 for cd,(rp,ep) in zip(d['selected99_original768_combinations'],positive,strict=True):
  c=sp(cd);cn={j-640:v for j,v in c.items() if j>=640};assert not mv(S0,cn)
  r=mv(P,mv(R,c));e=add(mv(E,{j:v for j,v in c.items() if j<640}),mv(L,mv(Se,cn)));e=mv(P,e)
  for v,bound in [(r,rp),(e,ep)]:
   for i,x in v.items():
    den,up=bound[i];assert den%x.denominator==0 and abs(x.numerator)*(den//x.denominator)<=up;bits=max(bits,abs(x.numerator).bit_length(),x.denominator.bit_length())
  actual[0].append(r);actual[1].append(e)
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'actual-cached-selected99-Hstar-r-ell-and-prior17.json';path.write_text(json.dumps({'actual_Hstar_r_ell_columns99':[[[[i,str(v)] for i,v in sorted(c.items())] for c in M] for M in actual],'prior_coefficient_and_partialsum_bits':17,'actual_coefficient_bits':bits,'all99_nativeS0_zero':True,'source':'Pi3 timesactualselected99cachedr/ellfirsts;ellLeS0termcancelsONLYbecauseS0cNzero;noindividualpivotchange','no_unknown_X2_or_physical_gauge_values_queried':True},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'separate_known_first_aggregation_prior_passed_before_coefficients':True,'known_first_aggregates_computed':True,'unknown_X2_or_physical_implicit_values_queried':False,'prior_bits':17,'actual_bits':bits,'owned_numeric_bound':owned,'input_sha256':{str(ip.relative_to(ROOT)):sha(ip)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'prior':17,'actual':bits,'numeric':owned,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
