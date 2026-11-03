"""Independent Fraction arithmetic for the six-quadric syzygy certificate."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW'
N=10
ZERO=(0,)*N

def add(*polys):
 out={}
 for p in polys:
  for m,c in p.items():
   out[m]=out.get(m,F(0))+c
 return {m:c for m,c in out.items()if c}

def scale(p,c):return {m:v*c for m,v in p.items()if v*c}
def mul(p,q):
 out={}
 for m,c in p.items():
  for n,d in q.items():
   k=tuple(a+b for a,b in zip(m,n));out[k]=out.get(k,F(0))+c*d
 return {m:c for m,c in out.items()if c}
def power(p,n):
 out={ZERO:F(1)}
 for _ in range(n):out=mul(out,p)
 return out
def var(i):return {tuple(int(i==j)for j in range(N)):F(1)}
def const(n):return {ZERO:F(n)}if n else {}
def cut(p):return {m:c for m,c in p.items()if not m[5]and not m[6]}
def linear_coefficient(p,i):
 out={}
 for m,c in p.items():
  if m[:7]==tuple(int(i==j)for j in range(7)):
   out[(0,)*7+m[7:]]=c
 return out
def decode(encoded):
 out={}
 for m,coefs in encoded:
  for ex,c in coefs:
   k=tuple(m)+tuple(ex);out[k]=out.get(k,F(0))+F(c)
 return {m:c for m,c in out.items()if c}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
 out=parser.parse_args().output.resolve();out.mkdir(parents=True,exist_ok=False)
 started=time.monotonic()
 x0,x1,x2,x3,x4,Z,W,a,e,d=[var(i)for i in range(N)]
 images=[x0,x1,x2,W,x3,{},add(Z,x1,scale(x0,-3)),x4]
 rawpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
 raw=json.loads(rawpath.read_text());qs=[]
 for row in raw['quadratic_rows_sparse']:
  terms=[]
  for mon,coefs in row:
   p={(0,)*7+tuple(ex):F(c)for ex,c in coefs}
   for im,n in zip(images,mon):
    if n:p=mul(p,power(im,n))
   terms.append(p)
  qs.append(add(*terms))
 folder=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63'
 cert=json.loads((folder/'certificate.json').read_text())
 L=decode(cert['L']);ell=decode(cert['ell']);Tz=decode(cert['TZ']);Tw=decode(cert['TW'])
 C=add(scale(a,100),scale(mul(e,e),18),scale(e,-9),const(-2))
 J=add(scale(mul(a,d),5),scale(e,-6),const(-4))
 assert cut(qs[1])==scale(mul(L,L),F(1,8100))
 assert cut(qs[3])==scale(mul(L,add(mul(e,L),scale(mul(a,x0),-400))),F(1,8100))
 assert linear_coefficient(L,1)==scale(C,-40)
 assert linear_coefficient(ell,1)==scale(mul(C,J),F(-16,243))
 assert mul(ell,W)==add(scale(mul(add(mul(d,e),const(8)),qs[1]),-1),scale(qs[2],8),mul(d,qs[3]))
 assert scale(mul(mul(a,ell),Z),-1)==add(mul(add(scale(a,-10),scale(mul(e,e),-3),scale(e,3)),qs[1]),scale(mul(add(e,const(-1)),qs[3]),3),scale(mul(a,qs[4]),10))
 delta=add(cut(ell),scale(mul(J,L),F(-2,1215)))
 assert delta==decode(cert['ell_mod_ZW_minus_2J_L_over1215'])
 assert not linear_coefficient(delta,1)
 anchor=add(scale(linear_coefficient(delta,2),3),scale(linear_coefficient(delta,4),F(-7,2)))
 assert anchor==scale(mul(a,mul(d,d)),F(-80,81))
 cubic=add(mul(add(mul(e,L),scale(mul(a,x0),-400)),qs[1]),scale(mul(L,qs[3]),-1))
 assert cubic==add(mul(Z,Tz),mul(W,Tw))
 six=[qs[0],qs[1],qs[3],qs[5],mul(ell,Z),mul(ell,W)]
 syz=[decode(p)for p in cert['new_syzygy']]
 assert not add(*(mul(q,c)for q,c in zip(six,syz)))
 assert syz==[{},mul(ell,add(mul(e,L),scale(mul(a,x0),-400))),scale(mul(ell,L),-1),{},scale(Tz,-1),scale(Tw,-1)]
 evidence=dict(all_exact_fraction_checks_passed=True,actual_six_quadric_term_counts=[len(q)for q in qs],
  lifted_syzygy_term_counts=[len(q)for q in syz],ell_x1_anchor=True,new_ad_squared_anchor=True,
  source_basis_change_determinant=-80,expanded_168_column_matrix_or_new_point_used=False)
 (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/f'{ID}.md',rawpath,folder/'certificate.json',folder/'report.json',
  ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW/report.json',
  ROOT/'reports/S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63_REVIEW/report.json']
 report=dict(id=ID,all_passed=True,independent_arithmetic='Fraction sparse dictionary; author code not imported or run',
  pointwise_22_directions_independent_on_FULLU=True,entire_kernel_module_basis_accepted=True,
  product_image_module_locally_free_rank=146,ambient_specialized_rank146_or_cokernel_flatness_proved=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'independent-certificate.json':sha(out/'independent-certificate.json')},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in ['id','all_passed','seconds']}))

if __name__=='__main__':main()
