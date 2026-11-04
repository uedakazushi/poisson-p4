#!/usr/bin/env python3
"""Independent tensor-resolution check accompanying the root's conceptual proof review.

No nil generators, products, or connecting-image values are evaluated here.
"""
import hashlib,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NIL_CUBE_OBSTRUCTION_REDUCTION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/(STEM+'_REVIEW');out.mkdir(exist_ok=False)
 base=ROOT/'reports'/STEM;report=json.loads((base/'report.json').read_text())
 assert sha(base/'report.json')=='64b74f5822a1f02dbb10fdf685ddc93f78e5c2e1779afb043bc8a82eee78a6c0'
 assert report['all_passed'] and report['no_actual_nil_lifts_products_or_Tor2_image_computed']
 pins={str((base/'report.json').relative_to(ROOT)):sha(base/'report.json')}
 for name,want in report['input_sha256'].items():
  p=ROOT/name;assert sha(p)==want,name;pins[name]=want
 for name,want in report['evidence_sha256'].items():
  p=base/name;assert sha(p)==want,name;pins[str(p.relative_to(ROOT))]=want
 d=json.loads((base/'exact-J-resolution-and-Tor2-input.json').read_text())
 B={int(k):[(tuple(I),j,h)for I,j,h in v]for k,v in d['basis_order'].items()}
 assert [len(B[i])for i in range(4)]==[1,10,44,112]
 zero=(0,)*9;actual_entries=0
 def mon(k,power=1):m=list(zero);m[k]=power;return tuple(m)
 # Koszul(z0,...,z6) tensor HB((r,s)^2), independently from its generators.
 for deg in [1,2,3]:
  lookup={x:i for i,x in enumerate(B[deg-1])};expected={}
  def put(row,col,m,c):
   key=(row,col);poly=expected.setdefault(key,{})
   poly[m]=poly.get(m,0)+c
  for col,(I,j,h)in enumerate(B[deg]):
   for t,k in enumerate(I):put(lookup[(I[:t]+I[t+1:],j,h)],col,mon(k),(-1)**t)
   sign=(-1)**len(I)
   if j==1:
    m=list(zero);m[7]=2-h;m[8]=h;put(lookup[(I,0,0)],col,tuple(m),sign)
   elif j==2:
    put(lookup[(I,1,h)],col,mon(8),-sign)
    put(lookup[(I,1,h+1)],col,mon(7),sign)
  expected={k:{m:c for m,c in v.items()if c}for k,v in expected.items()}
  actual={(i,j):{tuple(m):c for m,c in v}for i,j,v in d['differentials_sparse'][str(deg)]['entries']}
  assert actual==expected,deg;actual_entries+=sum(len(v)for v in actual.values())
 # Independent monomial-ideal divisibility test for the entire J^2 quotient.
 gens=[mon(i)for i in range(7)]+[mon(7,2),tuple([0]*7+[1,1]),mon(8,2)]
 squares=[tuple(a+b for a,b in zip(x,y))for x in gens for y in gens]
 standard=[]
 for total in range(5):
  for choices in itertools.combinations_with_replacement(range(9),total):
   m=[0]*9
   for k in choices:m[k]+=1
   if not any(all(a<=b for a,b in zip(g,m))for g in squares):standard.append(tuple(m))
 assert set(standard)==set(map(tuple,d['R_mod_J2_monomial_basis']))
 assert len(standard)==31 and all(sum(m)<=3 for m in standard)
 # Every degree-four monomial belongs to J^2, so no higher standard monomial exists.
 for choices in itertools.combinations_with_replacement(range(9),4):
  m=[0]*9
  for k in choices:m[k]+=1
  assert any(all(a<=b for a,b in zip(g,m))for g in squares)
 cube=json.loads((base/'cube-obstruction-tower.json').read_text())
 assert cube['total_D_generator_product_upper']==448
 assert 13*28+84==448 and [sum(len(I)==p for I,j,h in B[q])for q,p in [(1,1)]]==[7]
 note=ROOT/'research/notes'/(STEM+'_REVIEW.md')
 record=dict(id=STEM+'_REVIEW',all_passed=True,frozen=True,
  scope='Conceptual beta/gamma and connecting-map proof review; independent tensor/HB and monomial certificates. No actual nil product or lift.',
  canonical_beta_gamma_and_safe_W_quotients_accepted=True,
  beta_zero_implies_N5_zero=True,gamma_zero_implies_N4_zero=True,
  K_mod_J_factor_for_beta_NOT_assumed=True,
  Tor2_source_recipe_accepted_but_actual_image_NOT_computed=True,
  independent_tensor_coefficient_entries_checked=actual_entries,
  independent_J2_length=31,maximum_D_generator_products=448,
  original_nil_cube_still_OPEN=True,input_sha256=pins,
  code_sha256=sha(Path(__file__)),note_sha256=sha(note),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps({k:record[k]for k in ['all_passed','independent_tensor_coefficient_entries_checked','seconds']}))
if __name__=='__main__':main()
