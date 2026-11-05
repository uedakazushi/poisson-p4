#!/usr/bin/env python3
"""Independent free-affine-polynomial replay of all617 formal OLD rows."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_FORMAL_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p):pins[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='b18d47957cf657dae893ced69de3183d3437d2db9accdc6d13fefc970ed89fd3'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/617-formal-direct-old-highestthird-symbol-residuals.json';e=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 gb='reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63';g=read(gb+'/report.json');d=read(gb+'/617-direct-old-highestthird-coherent-source-positive-input.json');peer=read(gb+'_REVIEW/report.json');assert g['all_passed'] and peer['all_passed'] and pins[gb+'/report.json']['sha256']==peer['source_gate_sha256']
 maximum=0
 def ck(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=76;return x
 def add(a,b,c=Q(1)):
  v=dict(a)
  for mon,x in b.items():
   z=ck(v.get(mon,Q(0))+ck(c*x))
   if z:v[mon]=z
   else:v.pop(mon,None)
  return v
 def scale(a,c):return {m:ck(x*c) for m,x in a.items() if x*c}
 def mul(a,b):
  v={}
  for ma,x in a.items():
   for mb,y in b.items():
    if len(ma)+len(mb)<2:v=add(v,{ma+mb:ck(x*y)})
  return v
 origin={i:Q(x) for i,kind,x in d['literal_boundary_origins_and_membership']};membership={i:kind for i,kind,x in d['literal_boundary_origins_and_membership']};co={};ty={}
 def C(i):
  if i in co:return co[i]
  assert i in origin;v={(i,):Q(1)}
  if origin[i]:v[()]=origin[i]
  return v
 for i,n in d['literal_C0_program']:
  if n[0]=='c':v={():Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='n':v=scale(C(n[1]),Q(-1))
  elif n[0]=='a':v=add(C(n[1]),C(n[2]))
  else:assert n[0]=='m';v=mul(C(n[1]),C(n[2]))
  co[i]=v
 for i,n in d['typed_source_program']:
  if n[0]=='actual_C0_scalar_function':v=dict(C(n[1]))
  elif n[0]=='c':v={():Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='sum':
   v={}
   for j,a,b in n[1]:v=add(v,ty[j],Q(a,b))
  else:
   assert n[0]=='product';v={():Q(n[1],n[2])}
   for j in n[3]:v=mul(v,ty[j])
  ty[i]=v
 gamma=list(map(sp,d['Gamma_two']));cols={}
 for j,rows in d['pruned_original_old_columns_roots']:
  rows=dict(rows);pair=[]
  for row in gamma:
   v={}
   for k,c in row.items():
    if k not in rows:continue
    tag=rows[k]
    if tag[0]=='constant':p={}
    elif tag[0]=='C0':p=C(tag[1])
    else:assert tag[0]=='typed' and tag[1] in ty;p=ty[tag[1]]
    # Only degree-one monomials enter a highest-third symbol response.
    # Do not compute irrelevant weighted origin values outside that gate.
    v=add(v,{m:x for m,x in p.items() if m},c)
   pair.append(v)
  cols[j]=pair
 source=d['original619_source768_combinations']
 def case(j,ch):
  v={}
  for i,c in source[j]:
   if i<640:v=add(v,cols[i][ch],Q(c))
  return v
 basis={j:[case(j,ch) for ch in range(2)] for j in [563,609]};answers=[];union=set();nonzero=0
 for rec,got in zip(d['case_source_recipes617'],e['all617_directOLD_high3_formal_rows'],strict=True):
  assert rec['free_case619']==got['free_case619'];pair=[]
  for ch in range(2):
   v={}
   for f in rec['constant_case619_source_factors']:
    j=f['case619'];op,s=f['scale_operator'];assert op in ['literal','negate_literal'];c=Q(s) if op=='literal' else -Q(s);v=add(v,basis[j][ch] if j in basis else case(j,ch),c)
   row={m[0]:x for m,x in v.items() if m};assert row==sp(got['two_formal_directOLD_high3_symbol_coefficients'][ch]);D,U=map(int,d['positive_common_two_row_envelopes'][ch]);assert all(D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U for x in row.values());union.update(row);pair.append(row)
  nonzero+=bool(pair[0] or pair[1]);answers.append({'free_case619':rec['free_case619'],'two_verified_formal_rows':[[[i,str(x)] for i,x in sorted(row.items())] for row in pair]})
 assert len(answers)==617 and nonzero==50 and len(union)==37 and [[i,membership[i]] for i in sorted(union)]==e['surviving_literal_boundary_symbols']
 only79=sorted(i for i in union if membership[i]=='C0_79_only');physical=sorted(i for i in union if membership[i]=='source164');assert only79==[299650,505938] and len(physical)==35
 assert e['all_scalar_third_values_unassigned'] and e['no_basis_upper_zero_fill'] and not any(i in union for i,k,x in d['literal_boundary_origins_and_membership'] if k=='C0_79_only' and i not in only79)
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent617-free-affine-polynomial-formal-symbol-proof.json';ve.write_text(json.dumps({'all617_formal_directOLD_rows_verified':answers,'nonzero_formal_case_count':nonzero,'remaining35_literal164_symbol_ids':physical,'remaining2_only79_symbol_ids':only79,'all_other_genuine_honestquotient_highestthird_symbols_cancel_in_directOLD_only':True,'known_coefficient_maximum_bits':maximum,'no_scalar_thirds_native_knowncross_b1_A2_M1_or_upper_values':True,'scope':'DIRECTOLDhighestthirdsymbolJacobianonly. No source-order necessity, actual jet values, completeupper cancellation or HF2 equality claim.'},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'primary_sha256':pins[rp]['sha256'],'all617_formal_symbol_rows_independently_verified':True,'nonzero_formal_case_count':50,'remaining_literal164_symbols':35,'remaining_only79_symbols':only79,'no_actual_higher_or_upper_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'seconds':rr['seconds']}))
if __name__=='__main__':main()
