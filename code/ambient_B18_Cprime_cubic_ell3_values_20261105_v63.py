#!/usr/bin/env python3
"""Prepared two-pass fixed arithmetic; run only after separate cubic SOURCE GO."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_Cprime_cubic_ell3_source_preflight_20261105_v63 import ROOT,PATHS,STEM as GATE,build,sha,enc,sp,add
STEM='AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63'
def number_count(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(number_count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(number_count,v))
 return 0
def guard(v):
 if isinstance(v,Q):assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=2048
 elif isinstance(v,dict):
  for x in v.values():guard(x)
 elif isinstance(v,(list,tuple)):
  for x in v:guard(x)
 return v
def main(source_peer_path,quadratic_peer_path):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20))
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={};saved=0;peak=0
 def pin(path):
  p=ROOT/path;pins[path]=sha(p);return p
 def read(path):return json.loads(pin(path).read_text())
 source=read(f'reports/{GATE}/report.json');assert source['all_passed']and source['whole_saved_live_QQ_upper']==24000
 s_peer=read(source_peer_path);q_peer=read(quadratic_peer_path);assert s_peer['all_passed']and q_peer['all_passed']
 assert s_peer['input_sha256'][f'reports/{GATE}/report.json']==sha(ROOT/f'reports/{GATE}/report.json')
 assert q_peer['input_sha256'][PATHS['quadratic_report']]==sha(ROOT/PATHS['quadratic_report'])
 for p,h in source['input_sha256'].items():assert sha(pin(p))==h,p
 gate=read(f'reports/{GATE}/actual-rank1-quadratic-source-and-streamed-cubic-gate.json');assert gate['source_count']['total']==2019
 D,L,P,q3,E3,proj,m2,m3,i2,i3,receipt=build();input_count=number_count([source,s_peer,q_peer,gate,D])
 def custody(live):
  nonlocal peak
  guard(live);peak=max(peak,saved+input_count+number_count(live));assert peak<=24000,('WHOLE_NUMERIC_SLOT_GUARD',peak)
 def save(name,d,live=()):
  nonlocal saved
  custody([d,live]);blob=json.dumps(enc(d),separators=(',',':'))+'\n';raw=sum((ROOT/p).stat().st_size for p in pins);assert raw+len(blob.encode())+sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32<<20
  (out/name).write_text(blob);saved+=number_count(d);(out/'phase.json').write_text(json.dumps({'phase':name,'seconds':time.monotonic()-t,'numeric_slot_peak':peak})+'\n')
 save('preflight.json',{'input_sha256':pins,'all_input_numeric_leaves':input_count,'whole_numeric_slot_cap':24000,'source_count':2019,'compressed_rows':27,'no_new_source_gauge_or_higherjet':True})
 def project(raw):
  v={}
  for i,c in raw.items():v=add(v,proj[i],c)
  return v
 def stream():
  for j,col in enumerate(L):
   for mon in m2:
    raw={}
    for i,c in col.items():target,x=divmod(i,9);raw=add(raw,{165*target+i3[tuple(sorted(mon+(x,)))]:c})
    yield ['first44_times_S2',j,list(mon)],raw,project(raw)
  for var in range(9):
   raw={}
   for i,c in P.items():target,j=divmod(i,45);raw=add(raw,{165*target+i3[tuple(sorted(m2[j]+(var,)))]:c})
   yield ['actual_kernel18_quad_times_S1',var],raw,project(raw)
  for target in range(3):
   raw={165*target+i3[(8,8,8)]:Q(1)};yield ['ell_cube',target],raw,project(raw)
  for f,row in enumerate(E3):
   target,j=divmod(row,45)
   for var in range(9):yield ['consistency',f,var],{},add({9*f+var:Q(1)},proj[165*target+i3[tuple(sorted(m2[j]+(var,)))]],-1)
 sourcehash=hashlib.sha256();piv={};columns=[];rows=[];basis=[];labels=[];count=0
 for j,(label,raw,v)in enumerate(stream()):
  guard(v);sourcehash.update(json.dumps(enc([label,sp(raw),sp(v)]),separators=(',',':')).encode()+b'\n');count+=1;z=dict(v)
  for i,w in sorted(piv.items()):z=guard(add(z,w,-z.get(i,Q(0))))
  if z:
   i=min(z);c=z[i];piv[i]=guard({k:a/c for k,a in z.items()});columns.append(j);rows.append(i);basis.append(dict(v));labels.append(label)
  custody([L,P,q3,proj,m2,m3,piv,basis,columns,rows,labels,raw,v,z])
 assert count==2019 and len(basis)<=27;del piv
 save('phase1-streamed-fixed-pivots.json',{'source_count':count,'source_projection_stream_sha256':sourcehash.hexdigest(),'fixed_columns':columns,'fixed_rows':rows,'fixed_source_labels':labels,'selected_original_compressed_columns':list(map(sp,basis)),'rank':len(basis)},live=[L,P,q3,proj,m2,m3,basis])
 n=len(basis);minor=[[v.get(i,Q(0))for v in basis]for i in rows];aug=[row+[Q(i==j)for j in range(n)]for i,row in enumerate(minor)];det=Q(1)
 for j in range(n):
  k=next(k for k in range(j,n)if aug[k][j])
  if k!=j:aug[k],aug[j]=aug[j],aug[k];det=-det
  c=aug[j][j];det*=c;aug[j]=guard([v/c for v in aug[j]])
  for i in range(n):
   if i!=j and aug[i][j]:
    c=aug[i][j];aug[i]=guard([v-c*w for v,w in zip(aug[i],aug[j])])
  custody([L,P,q3,proj,m2,m3,basis,minor,aug])
 inverse=[row[n:]for row in aug];del aug
 assert all(sum((minor[i][k]*inverse[k][j]for k in range(n)),Q(0))==Q(i==j)for i in range(n)for j in range(n))
 assert all(sum((inverse[i][k]*minor[k][j]for k in range(n)),Q(0))==Q(i==j)for i in range(n)for j in range(n))
 factors=hashlib.sha256();backsource=hashlib.sha256();count2=0
 for label,raw,v in stream():
  f=guard([sum((c*v.get(rows[i],Q(0))for i,c in enumerate(row)),Q(0))for row in inverse]);back={}
  for col,c in zip(basis,f):back=guard(add(back,col,c))
  assert back==v;factors.update(json.dumps(enc(f),separators=(',',':')).encode()+b'\n');backsource.update(json.dumps(enc([label,sp(raw),sp(v)]),separators=(',',':')).encode()+b'\n');count2+=1
  custody([L,P,q3,proj,m2,m3,basis,minor,inverse,raw,v,back,f])
 assert count2==2019 and sourcehash.hexdigest()==backsource.hexdigest();zero=len(basis)==27
 evidence={'exact_NF3_rows_by135_sparse':list(map(sp,q3)),'unit_embedding_rows135':E3,'actual_kernel18_quadratic_representative135_sparse':sp(P),'all495_source_projection_columns27_sparse':list(map(sp,proj)),
  'all2019_source_and_projection_stream_sha256':sourcehash.hexdigest(),'fixed_columns':columns,'fixed_rows':rows,'fixed_labels':labels,'selected_original_compressed_columns':list(map(sp,basis)),
  'fixed_minor':minor,'fixed_inverse':inverse,'fixed_minor_determinant':det,'both_inverse_products_exact':True,'all2019_fixed_factorization_sha256':factors.hexdigest(),
  'source_count':count,'compressed_target_dimension':27,'known_cubic_rank':len(basis),'known_cubic_quotient_dimension':27-len(basis),'known_cubic_quotient_zero':zero,'cumulative_known_length_upper_if_zero':9 if zero else None,
  'Cprime_ell3_length9_and_long_exponent2_proved_if_zero':zero,'no_negative_exponent_inference_if_nonzero':True,'every_input_column_regenerated_twice':True}
 save('phase2-complete-streamed-cubic-upper-certificate.json',evidence,live=[basis,minor,inverse])
 report={'id':STEM,'all_passed':True,'status':'PROVED_cached_known_cubic_upper_test','known_cubic_rank':len(basis),'known_cubic_quotient_dimension':27-len(basis),'positive_Cprime_ell3_length9_and_long_exponent2':zero,'actual_Cprime_HF_through2':[3,3,3],
  'new_source_gauge_or_higherjet':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name not in['report.json','phase.json']},'numeric_slot_peak':peak,'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source-peer-report',required=True);p.add_argument('--quadratic-peer-report',required=True);a=p.parse_args();main(a.source_peer_report,a.quadratic_peer_report)
