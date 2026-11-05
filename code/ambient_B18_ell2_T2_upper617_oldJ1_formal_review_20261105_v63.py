#!/usr/bin/env python3
"""Independent two-graded free-polynomial replay of all617 J0/Jr/Je rows."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p):pins[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='01a655d356a02347b05d3ae853fe147f39e858d5eb169d0c4406b53e24d56516'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/617-formal-direct-old-third-and-second-symbol-sensitivities.json';e=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 gb='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63';g=read(gb+'/report.json');d=read(gb+'/617-direct-old-known-first-formal-Jacobian-positive-input.json');peer=read(gb+'_REVIEW/report.json');assert g['all_passed'] and peer['all_passed'] and peer['source_gate_sha256']==pins[gb+'/report.json']['sha256']
 high=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_FORMAL_20261105_v63/617-formal-direct-old-highestthird-symbol-residuals.json')['all617_directOLD_high3_formal_rows'];assert len(high)==617
 maximum=0
 def ck(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=82;return x
 def add(a,b,c=Q(1)):
  z=dict(a)
  for m,x in b.items():
   v=ck(z.get(m,Q(0))+ck(c*x))
   if v:z[m]=v
   else:z.pop(m,None)
  return z
 def scale(a,c):return {m:ck(c*x) for m,x in a.items() if c*x}
 # monomial=(base degree0/r/ell, perturbation symbol orNone).
 # Kill base degree>=2 and perturbation degree>=2, retaining r*theta/ell*theta.
 def mul(a,b):
  z={}
  for (ba,sa),x in a.items():
   for (bb,sb),y in b.items():
    if ba and bb or sa is not None and sb is not None:continue
    m=(ba or bb,sa if sa is not None else sb);z=add(z,{m:ck(x*y)})
  return z
 toy={(1,None):Q(2),(2,None):Q(3),(0,-1):Q(1)};assert mul(toy,toy)=={(1,-1):Q(4),(2,-1):Q(6)}
 first={i:list(map(Q,v)) for i,k,v in d['literal_boundary_firsts_point_r_ell_and_membership']};membership={i:k for i,k,v in d['literal_boundary_firsts_point_r_ell_and_membership']};co={};ty={}
 def C(i):
  if i in co:return co[i]
  assert i in first;v={(0,i):Q(1)}
  for base,x in enumerate(first[i]):
   if x:v[base,None]=x
  return v
 for i,n in d['literal_C0_program']:
  if n[0]=='c':v={(0,None):Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='n':v=scale(C(n[1]),Q(-1))
  elif n[0]=='a':v=add(C(n[1]),C(n[2]))
  else:assert n[0]=='m';v=mul(C(n[1]),C(n[2]))
  co[i]=v
 for i,n in d['typed_source_program']:
  if n[0]=='actual_C0_scalar_function':v=dict(C(n[1]))
  elif n[0]=='c':v={(0,None):Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='sum':
   v={}
   for j,a,b in n[1]:v=add(v,ty[j],Q(a,b))
  else:
   assert n[0]=='product';v={(0,None):Q(n[1],n[2])}
   for j in n[3]:v=mul(v,ty[j])
  ty[i]=v
 gamma=list(map(sp,d['Gamma_two']));cols={}
 for j,rr in d['pruned_original_old_columns_roots']:
  roots=dict(rr);pair=[]
  for row in gamma:
   v={}
   for i,c in row.items():
    if i not in roots:continue
    tag=roots[i]
    if tag[0]=='constant':p={}
    elif tag[0]=='C0':p=C(tag[1])
    else:assert tag[0]=='typed' and tag[1] in ty;p=ty[tag[1]]
    v=add(v,{m:x for m,x in p.items() if m[1] is not None},c)
   pair.append(v)
  cols[j]=pair
 def case(j,ch):
  v={}
  for i,c in d['original619_source768_combinations'][j]:
   if i<640:v=add(v,cols[i][ch],Q(c))
  return v
 basis={j:[case(j,ch) for ch in range(2)] for j in [563,609]};union=[set(),set(),set()];nonzero=[0,0,0];answers=[]
 for rec,got,oldhigh in zip(d['case_source_recipes617'],e['all617_formal_directOLD_J0_Jr_Jell'],high,strict=True):
  assert rec['free_case619']==got['free_case619']==oldhigh['free_case619'];pair=[];flags=[False]*3
  for ch in range(2):
   v={}
   for f in rec['constant_case619_source_factors']:
    j=f['case619'];op,x=f['scale_operator'];assert op in ['literal','negate_literal'];a=Q(x) if op=='literal' else -Q(x);v=add(v,basis[j][ch] if j in basis else case(j,ch),a)
   row={}
   for (base,i),x in v.items():
    assert i is not None
    if i not in row:row[i]=[Q(0)]*3
    row[i][base]=x
    if x:union[base].add(i);flags[base]=True
   actual={i:list(map(Q,x)) for i,x in got['two_formal_J0_Jr_Jell_symbol_rows'][ch]};assert row==actual
   j0={i:x[0] for i,x in row.items() if x[0]};assert j0==sp(oldhigh['two_formal_directOLD_high3_symbol_coefficients'][ch])
   for i,xs in row.items():
    for base,x in enumerate(xs):D,U=map(int,d['positive_common_two_by3_origin_Jr_Jell_envelopes'][ch][base]);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
   pair.append([[i,list(map(str,x))] for i,x in sorted(row.items())])
  for base,flag in enumerate(flags):nonzero[base]+=flag
  answers.append({'free_case619':rec['free_case619'],'all_two_J0_Jr_Jell_verified':pair})
 assert nonzero==[50,22,22]==e['nonzero_formal_cases_by_component'] and list(map(len,union))==[37,32,32]
 for base in range(3):assert [[i,membership[i]] for i in sorted(union[base])]==e['literal_symbols_by_component'][base]
 only79=[sorted(i for i in u if membership[i]=='C0_79_only') for u in union];assert all(u==[299650,505938] for u in only79)
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent617-two-graded-polynomial-Jacobian-and-J0-crosscheck.json';ve.write_text(json.dumps({'all617_J0_Jr_Jell_rows_verified':answers,'all617_J0_exactly_previous_HIGH3_rows':True,'nonzero_formal_case_counts':nonzero,'surviving_symbol_counts':[37,32,32],'only79_symbols_each_component':only79,'ordinary_factor2_product_fixture_verified':True,'known_coefficient_maximum_bits':maximum,'all_higher_source_coefficients_unassigned':True,'scope':'GenuinehonesthighestsecondandthirdsymbolscancelinDIRECTOLDonly. Twoonly79conductorsremainasformalvariableshere;actualzero7order>=3proofisaseparatesourcepremise. No scalarvalues/native/knownfirstcubic/b1/A2/M1/fullupperresponseorHF2equality.'},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'primary_sha256':pins[rp]['sha256'],'all617_J0_Jr_Jell_formal_coefficients_verified':True,'J0_matches_previous_ALL617_HIGH3':True,'ordinary_coefficient_product_rule_verified':True,'no_actual_higher_or_upper_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'seconds':rr['seconds']}))
if __name__=='__main__':main()
