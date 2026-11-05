#!/usr/bin/env python3
"""Source/positive first-Jacobian audit, without signed formal fold."""
import hashlib,json,math,resource,signal,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='ff4e80e5034da41bb93fe4afd9dfd83970cc7edf0d912afbb45dfeb433fb9110'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/617-direct-old-known-first-formal-Jacobian-positive-input.json';d=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 oldbase='reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63';s=read(oldbase+'/617-direct-old-highestthird-coherent-source-positive-input.json');sprev=read(oldbase+'_REVIEW/report.json');assert sprev['all_passed']
 for k in ['case_source_recipes617','original619_source768_combinations','pruned_original_old_columns_roots','Gamma_two','typed_source_program','literal_C0_program']:assert d[k]==s[k]
 # Pin actual first-slot meanings, rather than rely on a zero overlap.
 rawbase='reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63';rb=read(rawbase+'/actual-ell2-product33-and-raw45-actions.json',['actual33_basis'])['actual33_basis'];assert [rb[i] for i in [0,8,9]]==[[0,0,0],[0,1,0],[0,0,1]]
 a79={i:sp(v) for i,v in read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json',['first640'])['first640']['all79_atom_values_full10']}
 a164={i:sp(v) for i,v in read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',['all_requested164_root_values'])['all_requested164_root_values']}
 expected=[];first={}
 for i,kind,point in s['literal_boundary_origins_and_membership']:
  a=[a79[i].get(j,Q(0)) for j in [0,8,9]] if i in a79 else None;b=[a164[i].get(j,Q(0)) for j in [0,8,9]] if i in a164 else None
  assert a is not None or b is not None
  if a is not None and b is not None:assert a==b
  v=a if a is not None else b;assert v[0]==Q(point);first[i]=v;expected.append([i,kind,list(map(str,v))])
 assert expected==d['literal_boundary_firsts_point_r_ell_and_membership'] and len(first)==96
 assert all(v==['0','0','0'] for i,k,v in expected if k=='source164')
 maximum=0;zero=(1,0);Z=[zero]*3
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());return v
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else zero
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else zero
 def va(a,b):return [plus(x,y) for x,y in zip(a,b,strict=True)]
 def vs(a,c):return [mul(x,c) for x in a]
 def vm(a,b):return [mul(a[0],b[0]),plus(mul(a[0],b[1]),mul(a[1],b[0])),plus(mul(a[0],b[2]),mul(a[2],b[0]))]
 def add(a,b,c=None):
  z=dict(a)
  for i,x in b.items():
   q=va(z.get(i,Z),x if c is None else vm(x,c))
   if any(v[1] for v in q):z[i]=q
  return z
 def scale(a,c):
  z={}
  for i,x in a.items():
   q=vm(x,c)
   if any(v[1] for v in q):z[i]=q
  return z
 fv={i:[ex(x) for x in v] for i,v in first.items()};one=[ex(1),zero,zero];co={};cg={}
 def C(i):
  if i in co:return co[i],cg[i]
  assert i in fv;return fv[i],{i:one}
 for i,n in d['literal_C0_program']:
  if n[0]=='c':o=[ex(Q(n[1],n[2])),zero,zero];g={}
  elif n[0]=='n':o,g=C(n[1]);g=dict(g)
  else:
   a,u=C(n[1]);b,v=C(n[2])
   if n[0]=='a':o=va(a,b);g=add(u,v)
   else:assert n[0]=='m';o=vm(a,b);g=add(scale(u,b),v,a)
  co[i]=o;cg[i]=g
 to={};tg={}
 for i,n in d['typed_source_program']:
  if n[0]=='actual_C0_scalar_function':o,g=C(n[1]);g=dict(g)
  elif n[0]=='c':o=[ex(Q(n[1],n[2])),zero,zero];g={}
  elif n[0]=='sum':
   o=list(Z);g={}
   for j,a,b in n[1]:x=ex(Q(a,b));o=va(o,vs(to[j],x));g=add(g,tg[j],[x,zero,zero])
  else:
   assert n[0]=='product';o=[ex(Q(n[1],n[2])),zero,zero];g={}
   for j in n[3]:g=add(scale(g,to[j]),tg[j],o);o=vm(o,to[j])
  to[i]=o;tg[i]=g
 gamma=list(map(sp,d['Gamma_two']));col={}
 for j,rr in d['pruned_original_old_columns_roots']:
  roots=dict(rr);pair=[]
  for row in gamma:
   z={}
   for k,x in row.items():
    if k not in roots:continue
    tag=roots[k]
    if tag[0]=='constant':v={}
    elif tag[0]=='C0':v=C(tag[1])[1]
    else:assert tag[0]=='typed' and tag[1] in tg;v=tg[tag[1]]
    z=add(z,v,[ex(x),zero,zero])
   pair.append(z)
  col[j]=pair
 def case(j,ch):
  z={}
  for i,x in d['original619_source768_combinations'][j]:
   if i<640:z=add(z,col[i][ch],[ex(x),zero,zero])
  return z
 bas={j:[case(j,ch) for ch in range(2)] for j in [563,609]};common=[[zero]*3 for _ in range(2)];scalarout=0;uniqueout=0
 for rec in d['case_source_recipes617']:
  for ch in range(2):
   z={}
   for f in rec['constant_case619_source_factors']:
    j=f['case619'];z=add(z,bas[j][ch] if j in bas else case(j,ch),[ex(f['scale_operator'][1]),zero,zero])
   uniqueout+=len(z)
   for v in z.values():
    for k,(D,U) in enumerate(v):
     if U:scalarout+=1
     a,b=common[ch][k];dd=math.lcm(a,D);common[ch][k]=mark((dd,max(b*(dd//a),U*(dd//D)))) if b and U else (D,U) if U else (a,b)
 ns=sum(sum(sum(bool(x[1]) for x in v) for v in z.values()) for z in list(cg.values())+list(tg.values()));cs=sum(sum(sum(sum(bool(x[1]) for x in v) for v in z.values()) for z in pair) for pair in col.values());nk=sum(map(len,cg.values()))+sum(map(len,tg.values()));ckeys=sum(sum(map(len,pair)) for pair in col.values())
 assert [ns,cs,scalarout]==[d['positive_node_scalar_gradient_entries'],d['positive_old_column_scalar_gradient_entries'],d['positive617_scalar_gradient_entries']]==[632,372,646]
 assert [[[str(a),str(b)] for a,b in row] for row in common]==d['positive_common_two_by3_origin_Jr_Jell_envelopes'] and maximum==82 and 3*maximum+1==247
 assert [[str(a),str(b)] for a,b in [row[0] for row in common]]==s['positive_common_two_row_envelopes']
 # Unique key counts account for the stored zero slots in each triple, not
 # just the nonzero scalar-entry counts used by the producer forecast.
 owned=2*count(d)+4*1161+1161+4*nk+618+4*ckeys+2+16*96+2*(617+4*uniqueout+3*96+30)+4*96+96+1000+5000
 assert owned<d['future_formal_numeric_bound']==49214 and d['future_formal_RAW_bound']==6429260
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_PREPARED_CUSTODY_20261105_v63/report.json';cust=read(cp);assert cust['all_passed'];wp=ROOT/cust['worker'];assert sha(wp)==cust['worker_sha256']=='22f54bfc0d90e29031caaf5901417d6d6a1b8b7402245100df3139c65d589b13'
 for p,h in cust['literal_input_sha256'].items():assert sha(ROOT/p)==h
 z=subprocess.run([sys.executable,str(wp),'--help'],capture_output=True,text=True,check=True);assert not z.stderr and hashlib.sha256(z.stdout.encode()).hexdigest()==cust['safe_help_stdout_sha256']
 proof={'ordinary_coefficient_identity':'For highest boundary third/second symbols in directOLD coefficient r^2ell, the response is J0*delta_f_r2ell+Jr*delta_f_rell+Jell*delta_f_rr. J(u0+ur*r+ue*ell) is computed over k[r,ell]/(r,ell)^2, keeping the full product-rule SUM. No extra factorial or factor2 is introduced.',
 'independent_symbols':'All96 higher coefficients remain unassigned. Primitive source164 plane firsts are exactly0; only known79 firsts can enter Jr/Jell. Zero or nonzero formal weights say nothing yet about actual boundarysecond values or their necessity. Actualnative,knowncubic,b1,A2,M1 remain separate.',
 'descent_scope':'Only certifiedliteral164scalarR0 pureplane slots have raw33/plane9 descent. The48only79 boundaries do not inherit highercoefficients from that catalogue. This SOURCE peer accessed point/r/ell firsts only.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent96-firsts-three-component-positive-Jacobian-source-proof.json';ve.write_text(json.dumps({'all_program_roots_source617_recipes_same_as_V2':True,'raw33_indices0_8_9_are_point_r_ell_exact':True,'all96_known_first_vectors_and_overlaps_checked':True,'literal164_plane_firsts_all0':True,'scalar_support_counts':[ns,cs,scalarout],'unique_key_counts':[nk,ckeys,uniqueout],'all_positive_two_by3_envelopes_rebuilt':True,'J0_positive_projection_equals_previous_HIGH3':True,'positive_bits':82,'temporary_bits':247,'actual_worker_zero_slots_in_triples_owned_bound':owned,'whole_cap':49214,'worker_custody_and_safe_help_literal':True,'no_signed_Jacobian_or_actual_higher_response':True,'proof':proof},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'source_gate_sha256':pins[rp]['sha256'],'formal_only_source_code_custody_condition_closed':True,'actual_worker_owned_bound':owned,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)}};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owned,'unique':[nk,ckeys,uniqueout]}))
if __name__=='__main__':main()
