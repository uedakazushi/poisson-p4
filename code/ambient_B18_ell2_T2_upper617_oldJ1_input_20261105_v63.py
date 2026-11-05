#!/usr/bin/env python3
"""Positive known-first OLD Jacobian over k[r,ell]/m²; no source2/3 values."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={};maximum=0;cursor={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 bp='reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63/';source=read(bp+'617-direct-old-highestthird-coherent-source-positive-input.json');assert source['all_pruned_roots_explicit_type_bound_no_missing_defaults']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63_REVIEW/report.json')['all_passed']
 fc=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json','first640');f79={i:H.parse(v) for i,v in fc['all79_atom_values_full10']};del fc
 raw=field('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json','all_requested164_root_values');f164={i:H.parse(v) for i,v in raw};del raw
 first={};classes=[]
 for i,kind,point in source['literal_boundary_origins_and_membership']:
  a=[f79[i].get(k,F(0)) for k in [0,8,9]] if i in f79 else None;b=[f164[i].get(k,F(0)) for k in [0,8,9]] if i in f164 else None
  if a is not None and b is not None:assert a==b
  v=a if a is not None else b;assert v is not None and v[0]==F(point);first[i]=v;classes.append([i,kind,list(map(str,v))])
 zero=(1,0);Z=[zero]*3
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length())
  if maximum>384:
   H.dump(out/'positive-guard-frontier.json',{'cursor':cursor,'positive_bound':list(map(str,q)),'bits':maximum,'no_signed_or_actual_source2_3_values':True});raise RuntimeError('positive384guard')
  return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else zero
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=math.lcm(a[0],b[0]);return mark((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else zero
 def vsum(a,b):return [plus(x,y) for x,y in zip(a,b,strict=True)]
 def vscale(a,c):return [mul(x,c) for x in a]
 def vmul(a,b):return [mul(a[0],b[0]),plus(mul(a[0],b[1]),mul(a[1],b[0])),plus(mul(a[0],b[2]),mul(a[2],b[0]))]
 def add(a,b,c=None):
  z=dict(a)
  for k,v in b.items():
   q=v if c is None else vmul(v,c);q=vsum(z.get(k,Z),q)
   if any(x[1] for x in q):z[k]=q
  return z
 def scale(a,c):return {i:vmul(v,c) for i,v in a.items() if any(x[1] for x in vmul(v,c))}
 one=[ex(1),zero,zero];fv={i:[ex(x) for x in v] for i,v in first.items()};co={};cg={}
 def C(i):return (co[i],cg[i]) if i in co else (fv[i],{i:one})
 for i,n in source['literal_C0_program']:
  cursor={'stage':'positive_C0_Jfirst','node':i};op=n[0]
  if op=='c':o=[ex(F(n[1],n[2])),zero,zero];g={}
  elif op=='n':o,g=C(n[1]);g=dict(g)
  elif op=='a':a,u=C(n[1]);b,v=C(n[2]);o=vsum(a,b);g=add(u,v)
  else:assert op=='m';a,u=C(n[1]);b,v=C(n[2]);o=vmul(a,b);g=add(scale(u,b),v,a)
  co[i]=o;cg[i]=g
 to={};tg={}
 for i,n in source['typed_source_program']:
  cursor={'stage':'positive_typed_Jfirst','node':i};op=n[0]
  if op=='actual_C0_scalar_function':o,g=C(n[1]);g=dict(g)
  elif op=='c':o=[ex(F(n[1],n[2])),zero,zero];g={}
  elif op=='sum':
   o=list(Z);g={}
   for j,a,b in n[1]:s=ex(F(a,b));o=vsum(o,vscale(to[j],s));g=add(g,tg[j],[s,zero,zero])
  else:
   assert op=='product';o=[ex(F(n[1],n[2])),zero,zero];g={}
   for j in n[3]:g=add(scale(g,to[j]),tg[j],o);o=vmul(o,to[j])
  to[i]=o;tg[i]=g
 gamma=[H.parse(v) for v in source['Gamma_two']];columns={}
 for j,rr in source['pruned_original_old_columns_roots']:
  roots=dict(rr);gg=[]
  for row in gamma:
   z={}
   for i,x in row.items():
    if i not in roots:continue
    tag=roots[i]
    if tag[0]=='constant':v={}
    elif tag[0]=='C0':v=C(tag[1])[1]
    else:assert tag[0]=='typed' and tag[1] in tg;v=tg[tag[1]]
    z=add(z,v,[ex(x),zero,zero])
   gg.append(z)
  columns[j]=gg
 def case(j,ch):
  z={}
  for k,x in source['original619_source768_combinations'][j]:
   if k<640:z=add(z,columns[k][ch],[ex(x),zero,zero])
  return z
 basis={j:[case(j,ch) for ch in range(2)] for j in [563,609]};common=[[zero]*3 for _ in range(2)];outentries=0
 def union(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=math.lcm(a[0],b[0]);return mark((d,max(a[1]*(d//a[0]),b[1]*(d//b[0]))))
 for r in source['case_source_recipes617']:
  cursor={'stage':'positive617_Jfirst_residual','case':r['free_case619']}
  for ch in range(2):
   z={}
   for f in r['constant_case619_source_factors']:
    j=f['case619'];v=basis[j][ch] if j in basis else case(j,ch);z=add(z,v,[ex(f['scale_operator'][1]),zero,zero])
   for v in z.values():
    for k,x in enumerate(v):common[ch][k]=union(common[ch][k],x);outentries+=bool(x[1])
 nodeentries=sum(sum(sum(bool(x[1]) for x in v) for v in d.values()) for d in list(cg.values())+list(tg.values()));colentries=sum(sum(sum(sum(bool(x[1]) for x in v) for v in d.values()) for d in pair) for pair in columns.values())
 data={k:source[k] for k in ['case_source_recipes617','original619_source768_combinations','pruned_original_old_columns_roots','Gamma_two','typed_source_program','literal_C0_program']};data.update({'literal_boundary_firsts_point_r_ell_and_membership':classes,'all96_full9_first_frame_overlaps_verified':True,'positive_common_two_by3_origin_Jr_Jell_envelopes':[[[str(x) for x in q] for q in row] for row in common],'source_committed_prior_bits':maximum,'temporary_prior_bits':3*maximum+1,'positive_node_scalar_gradient_entries':nodeentries,'positive_old_column_scalar_gradient_entries':colentries,'positive617_scalar_gradient_entries':outentries,'no_signed_Jacobian_or_actual_higher_source_response':True,'formal_only_future_scope':'J0independentthird+Jrindependentsecond_rell+Jellindependentsecond_rr ofDIRECTOLDcubicONLY;noactualhighercoeffs/native/b1/A2/M1/uppervalues.'})
 inp=count(data);future=2*(inp+1000)+4*nodeentries+4*colentries+4*outentries+16000;rawbound=2000000+future*(2*math.ceil(maximum*math.log10(2))+40);owned=count([source,first,fv,co,cg,to,tg,columns,basis,data])+25000;assert maximum<=384 and 3*maximum+1<=1153 and future<150000 and rawbound<32*1024*1024 and owned<150000
 data.update({'future_formal_numeric_bound':future,'future_formal_RAW_bound':rawbound,'metadata_reserve':1000});ep=out/'617-direct-old-known-first-formal-Jacobian-positive-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'known_first_source_and_positive_envelopes_only':True,'source_committed_bits':maximum,'source_temporary_bits':3*maximum+1,'future_formal_numeric_bound':future,'future_formal_RAW_bound':rawbound,'source_owned_bound':owned,'no_signed_Jacobian_actual_second_third_or_upper_values':True,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_committed_bits','source_temporary_bits','future_formal_numeric_bound','future_formal_RAW_bound','source_owned_bound']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
