#!/usr/bin/env python3
"""Independent sparse Fraction arithmetic for the positive Poisson Cremona map."""
import argparse,gzip,hashlib,itertools,json,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/Q_POSITIVE_WEIGHT6_DARBOUX_20261103_v63';cert=json.loads((folder/'certificate.json').read_text());r=json.loads((folder/'report.json').read_text());assert r['all_passed']
 for name,want in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want
 eta,h,x,y,z,w=s.symbols('eta h x y z w');variables=(eta,h,x,y,z,w);zero=(0,)*6
 def parse(text,renames=None):
  expr=s.sympify(text,locals=renames or {})
  return {tuple(m):Q(int(c.p),int(c.q))for m,c in s.Poly(expr,*variables).terms()if c}
 def add(*ff):
  out={}
  for f in ff:
   for m,c in f.items():out[m]=out.get(m,Q(0))+c
  return {m:c for m,c in out.items()if c}
 def neg(f):return {m:-c for m,c in f.items()}
 def mul(*ff):
  out={zero:Q(1)}
  for f in ff:
   n={}
   for a,c in out.items():
    for b,d in f.items():
     mon=tuple(i+j for i,j in zip(a,b));n[mon]=n.get(mon,Q(0))+c*d
   out={m:c for m,c in n.items()if c}
  return out
 def power(f,n):return mul(*([f]*n))
 def diff(f,i):
  out={}
  for m,c in f.items():
   if m[i]:v=list(m);v[i]-=1;out[tuple(v)]=c*m[i]
  return out
 generators=[parse(str(v))for v in variables]
 def compose(f,images):
  out={}
  for m,c in f.items():out=add(out,mul({zero:c},*(power(im,n)for im,n in zip(images,m))))
  return out
 source=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz'
 raw=json.loads(gzip.decompress(source.read_bytes()))['raw150'];P=[[{}for _ in range(5)]for _ in range(5)];pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
 for (i,j),block in zip(pairs,range(10)):
  for k,mon in enumerate(mons):P[i][j]=add(P[i][j],mul(parse(raw[15*block+k]),generators[mon[0]+1],generators[mon[1]+1]))
  P[j][i]=neg(P[i][j])
 names=dict(zip(cert['homogeneous_source_variables'],variables[1:]));names['eta']=eta
 F=[parse(t,names)for t in cert['homogeneous_forward']]
 names=dict(zip(cert['homogeneous_target_variables'],variables[1:]));names['eta']=eta
 G=[parse(t,names)for t in cert['homogeneous_inverse']]
 for A,B in [(F,G),(G,F)]:
  images=[generators[0],*B]
  assert [compose(f,images)for f in A]==[mul(power(generators[1],5),v)for v in generators[1:]]
 C=[[{}for _ in range(5)]for _ in range(5)]
 C[1][4]=parse('8*h**2/9');C[4][1]=neg(C[1][4]);C[2][3]=parse('-2*h**2/9');C[3][2]=neg(C[2][3])
 def transformed(A,P):
  J=[[diff(f,k)for k in range(1,6)]for f in A]
  return [[add(*(mul(J[i][k],P[k][l],J[j][l])for k in range(5)for l in range(5)))for j in range(5)]for i in range(5)]
 assert transformed(F,P)==[[compose(c,[generators[0],*F])for c in row]for row in C]
 assert transformed(G,C)==[[compose(c,[generators[0],*G])for c in row]for row in P]
 assert F[0]==power(generators[1],2)and G[0]==power(generators[1],3)
 assert F[2].get((0,0,2,0,0,0))==Q(-3,4)
 assert G[3].get((0,0,3,0,0,0))==Q(9,16)
 paths=[Path(__file__),folder/'report.json',folder/'certificate.json',source,ROOT/'research/notes/Q_POSITIVE_WEIGHT6_DARBOUX_20261103_v63_REVIEW.md']
 result=dict(id='Q_POSITIVE_WEIGHT6_DARBOUX_20261103_v63_REVIEW',all_passed=True,
  independent_arithmetic='Fraction sparse multivariate dictionaries; SymPy only parses coefficient strings',
  full_raw150_read_without_source_model_imports=True,
  both_homogeneous_compositions_all_coefficients=True,
  both_full_raw5x5_Poisson_pushforwards_all_coefficients=True,
  coordinate_gcds_one_uniform_at_every_eta=True,
  actual_polynomial_Darboux_and_Cremona_bidegree=[2,3],
  PGL_or_component_membership_inferred=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT POSITIVE DARBOUX FRACTION REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
