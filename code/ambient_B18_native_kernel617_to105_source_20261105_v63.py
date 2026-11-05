#!/usr/bin/env python3
"""Constant617 native source factorization through explicit105 kerS0 basis."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple,set)):return sum(map(count,x))
 if isinstance(x,str):
  try:F(x);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(x,(int,float,F)))
def select(path,names):
 s=path.read_text();de=json.JSONDecoder();i=s.index('{')+1;out={};peak=0
 while True:
  while s[i].isspace() or s[i]==',':i+=1
  if s[i]=='}':break
  k,i=de.raw_decode(s,i)
  while s[i].isspace() or s[i]==':':i+=1
  v,i=de.raw_decode(s,i);peak=max(peak,count(v))
  if k in names:out[k]=v
 assert set(out)==set(names);return out,peak
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;parsepeak=0;actual=0
 def read(p,names=None):
  nonlocal raw,parsepeak
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;assert raw<32*1024*1024
  if names is None:v=json.loads(q.read_text());parsepeak=max(parsepeak,count(v));return v
  v,peak=select(q,names);parsepeak=max(parsepeak,peak);return v
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 kp='reports/AMBIENT_B18_NATIVE_POINT_KERNEL105_SOURCE_20261105_v63/'
 kr=read(kp+'report.json');assert kr['all_passed'] and pins[kp+'report.json']=='15b6ff370efd02939a1fac290ca6c24c1237ab3046f6a82984d22b1071e8d534'
 ke='constant-kernel105-basis-unit23-minor-and-native-functional-determinacy.json';k=read(kp+ke);assert pins[kp+ke]==kr['evidence_sha256'][ke]
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json'
 s=read(sp,['original619_c0_source_combinations','kernel617_source_recipes']);c=s['original619_c0_source_combinations'];recipes=s['kernel617_source_recipes'];assert len(c)==619 and len(recipes)==617
 mp='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';S=read(mp,['native_S0_columns'])['native_S0_columns'];assert len(S)==128
 cols=k['explicit105_kernel_basis_columns'];free=[v['free_native_source'] for v in cols];assert len(set(free))==105
 # All signed calculations below involve CONSTANT source coefficients only.
 # Prior fromsourcec0, recipeweights, kernelweights, actualfull48S0; one
 # scalarfromeachclasssuffices. Path4096 dominates3*105*3*2 pointpaths.
 classes={'c0':[q for col in c for j,q in col],
  'recipeweights':['1']+[q for rec in recipes for q in rec['subtract_actual_basis563609']],
  'kernelweights':[q for col in cols for j,q in col['column']],
  'fullS0':['1']+[q for col in S for i,q in col]}
 D=1;U=4096;cl={}
 for name,v in classes.items():
  xs=list(map(F,v));d=1
  for q in xs:d=math.lcm(d,q.denominator);assert d.bit_length()<=384
  u=max([d]+[abs(q.numerator)*(d//q.denominator) for q in xs]);assert u.bit_length()<=384
  cl[name]=[str(d),str(u)];D*=d;U*=u;assert max(D.bit_length(),U.bit_length())<=384
 prior=max(D.bit_length(),U.bit_length());assert 3*prior+1<=1153
 before=count([kr,k,s,S,cols,free,classes,cl,pins])+parsepeak+20000;assert before<150000
 dump('phase-before-constant-source-arithmetic.json',{'constant_source_only':True,'positive_classes':cl,'common_D':str(D),'common_U_paths4096':str(U),'committed_prior':prior,'temporary_prior':3*prior+1,'source_upper_before_signed_constants':before})
 def ck(q):
  nonlocal actual
  actual=max(actual,abs(q.numerator).bit_length(),q.denominator.bit_length());assert actual<=prior;return q
 def put(v,i,q):
  z=ck(v.get(i,F(0))+q)
  if z:v[i]=z
  else:v.pop(i,None)
 def native(col):return {j-640:F(q) for j,q in col if j>=640 and F(q)}
 base=[native(col) for col in c];basis=[{j:F(q) for j,q in col['column']} for col in cols];Sf=[{i:F(q) for i,q in col} for col in S];results=[];maxsupport=0;total=0
 for rec in recipes:
  j=rec['free_case619'];weights=[F(q) for q in rec['subtract_actual_basis563609']];v=dict(base[j])
  for weight,b in zip(weights,[563,609]):
   for src,q in base[b].items():put(v,src,ck(-weight*q))
  z={}
  for src,q in v.items():
   for row,a in Sf[src].items():put(z,row,ck(q*a))
  assert not z # This is literalS0*cN, NEVERfromM0X cancellation.
  coords={i:v.get(src,F(0)) for i,src in enumerate(free) if v.get(src,F(0))}
  reconstruct={}
  for i,q in coords.items():
   for src,a in basis[i].items():put(reconstruct,src,ck(q*a))
  assert reconstruct==v
  assert all(v.get(u,F(0))==0 for u in k['fixed23_source_complement']) or any(v.get(u,F(0)) for u in k['fixed23_source_complement'][:7])
  # Complementmask0coefficientsmaybenonzeroandareRECONSTRUCTEDfrombasis;
  # onlyfree16coefficientshavezeroindividually.
  assert all(not v.get(u,F(0)) for u in range(112,128))
  maxsupport=max(maxsupport,len(coords));total+=len(coords)
  results.append({'free_case619':j,'subtract563609':list(map(str,weights)),
   'literal_native128_constant_combination':[[i,str(q)] for i,q in sorted(v.items())],
   'kernel105_coordinates':[[i,str(q)] for i,q in sorted(coords.items())]})
 assert len(results)==617 and {v['free_case619'] for v in results}==set(range(619))-{563,609}
 ep='all617-native-constant-sources-factor-through-fixed105-kernel.json'
 ev={'all617_literal_full48S0cN_zero_verified':True,'all617_full128_reconstructions_verified':True,
  'kernel105_free_native_indices':free,'all617_source_maps':results,'total_coordinate_nonzero_records':total,'maximum_per_source_coordinate_support':maxsupport,
  'functional_factorization_ONLY':'The nativeGammaXc3+DnonzzXc2 two-vector for eachcorrected617source is thesamefixedtwo-by105 FUNCTIONAL (unqueried) times thisconstantcoordinatecolumn. AllFIRST andindividualXp2 outerterms remainexternal.',
  'none_of_210_functional_values_or617_upperresponses_computed':True,'all_unknown_N2mixed_N3_and_honest_boundaries_still_scoped_by_existing_SOURCE':True,
  'positive_committed_prior':prior,'actual_constant_source_bits':actual}
 dump(ep,ev);owned=count([kr,k,s,S,cols,free,classes,cl,base,basis,Sf,results,ev,pins])+20000;assert owned<150000
 rawwhole=raw+sum(q.stat().st_size for q in out.iterdir());assert rawwhole<32*1024*1024
 rep={'id':STEM,'all_passed':True,'constant_SOURCE_only_no_responses':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep),'phase-before-constant-source-arithmetic.json':sha(out/'phase-before-constant-source-arithmetic.json')},'source_owned_upper':owned,'source_RAW_upper':rawwhole,'source_prior':prior,'temporary_prior':3*prior+1,'actual_constant_bits':actual,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'prior':prior,'actual':actual,'owned':owned,'coordinates':total}))
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);a=pa.parse_args();main(a.output)
