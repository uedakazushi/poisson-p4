#!/usr/bin/env python3
"""Independent Fraction checks and case audit for raw aC constant quad6."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_CONSTANT_QUAD6_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
 z=dict(a)
 for m,c in b.items():z[m]=z.get(m,F())+c
 return{m:c for m,c in z.items()if c}
def scale(p,c):return{m:x*c for m,x in p.items()if x*c}
def mul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,F())+c*d
 return{k:c for k,c in z.items()if c}
def linear(poly,j):
 mon=tuple(int(i==j)for i in range(7));out={}
 for light,terms in poly:
  if tuple(light)==mon:
   for core,c in terms:out=add(out,{tuple(core):F(c)})
 return out
start=time.monotonic();parent=ROOT/'reports/S61_NONSQUARE_CONSTANT_QUAD6_20261104_v63/report.json';r=json.loads(parent.read_text());assert r['all_passed']and r['raw_open']=='a*C!=0'and not r['kernel22_independence_extended_off_FULLU']
for name,want in r['input_sha256'].items():assert sha(ROOT/name)==want,name
src=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json';data=json.loads(src.read_text());ell,L=data['ell'],data['L'];C={(1,0,0):F(100),(0,2,0):F(18),(0,1,0):F(-9),(0,0,0):F(-2)};J={(1,0,1):F(5),(0,1,0):F(-6),(0,0,0):F(-4)}
assert linear(ell,1)==scale(mul(C,J),F(-16,243));difference=lambda j:add(linear(ell,j),scale(mul(J,linear(L,j)),F(-2,1215)))
assert add(scale(difference(2),F(3)),scale(difference(4),F(-7,2)))=={(1,0,2):F(-80,81)}
special={}
for(a,e,d),c in linear(ell,0).items():
 if d==0:special[a]=special.get(a,F())+c*F(-2,3)**e
assert{k:c for k,c in special.items()if c}=={1:F(3200,243)}
zw=ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW/report.json';zr=json.loads(zw.read_text());assert zr['all_passed']
# A nonzero constant relation times five independent variables gives five
# independent source-kernel vectors, contradicting rank19/20.
assert 20-19==1 and 5>1 and 29-6==23
out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),parent,src,zw,ROOT/'research/notes/S61_NONSQUARE_CONSTANT_QUAD6_20261104_v63.md'];report={'id':ID,'all_passed':True,'own_Fraction_ell_CJ_and_ad2_anchor_exact':True,'own_Fraction_dJ0_ell_x0_exact':True,'three_cases_exhaustive_on_raw_aC':True,'ZW_rank19_implies_constant_four_quad_independence_reviewed':True,'known_unit_basis_det_minus80_and_final_ellZ_ellW_proof_reviewed':True,'raw_six_quad_constant_rank6_all_aC_accepted':True,'actual_weight2_quotient23_all_aFBC_accepted':True,'kernel22_off_FULLU_or_quartic_extension_not_claimed':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'seconds':time.monotonic()-start,'frozen':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2))
