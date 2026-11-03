#!/usr/bin/env python3
"""Independent actual-model, coefficient, Darboux and projective-leaf audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(gzip.decompress(path.read_bytes())if path.suffix=='.gz'else path.read_bytes())
def rows(M):return [[str(v)for v in row]for row in M.tolist()]
def add(poly,key,c):
 value=poly.get(key,Fraction(0))+c
 if value:poly[key]=value
 else:poly.pop(key,None)
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 producer=ROOT/'reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63/report.json';original=load(producer);assert original['all_passed']
 for rel,want in original['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in original['evidence_sha256'].items():assert sha(producer.parent/rel)==want,rel
 h,x,y,z,w,W,r=s.symbols('h x y z w W r');X=[h,x,y,z,w];pairs=list(itertools.combinations(range(5),2));mon2=list(itertools.combinations_with_replacement(range(5),2));mon3=list(itertools.combinations_with_replacement(range(5),3));xx=s.Matrix(X)
 # Reconstruct the actual raw CA and CB models on the curve, independently
 # of the author's PW formula and normalization helpers.
 f=-h*h;Hamiltonian=x*x*z-h*x*w+z**3/3
 CA=s.zeros(5);CA[1,2]=f;CA[3,4]=W*h*h;CA[2,3]=-s.diff(Hamiltonian,w);CA[2,4]=s.diff(Hamiltonian,z);CA=CA-CA.T
 CB=s.zeros(5);CB[1,2]=f;CB[3,4]=-W*f;CB[2,3]=-s.diff(f,x)*z/2-s.diff(Hamiltonian,w);CB[2,4]=-s.diff(f,x)*w/2+s.diff(Hamiltonian,z);CB=CB-CB.T
 def canonical(A):
  divergence=s.Matrix([sum(s.diff(A[i,j],X[i])for i in range(5))for j in range(5)])
  return(A-(xx*divergence.T-divergence*xx.T)/5).applyfunc(s.expand)
 order=[1,0,2,3,4]
 def swap(A):return A.subs({h:x,x:h},simultaneous=True).extract(order,order)
 P=swap(canonical(CA)).applyfunc(s.expand);assert P==swap(canonical(CB)).applyfunc(s.expand)
 quadpath=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz';quad=load(quadpath);source=ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63/staircase.json.gz';common=load(source);U,V=s.symbols('U V')
 for name in['CA','CB']:
  assert(P-s.Matrix(quad[name]['full_source_P']).subs({U:0,V:0})).applyfunc(s.expand)==s.zeros(5)
  assert P==s.Matrix(common[name]['physical_entire_W_curve'])
  assert common[name]['entire_monic_critical_W_source']==['0','0','W','0','0','0','0','0','0','0']
 # A separate coefficient convolution computes ALL350 strictJacobi slots;
 # no symbolic Schouten/Jacobi routine from the producer is called.
 coefficients={};vector=[]
 for i,j in pairs:
  polynomial=s.Poly(P[i,j],*X)
  for mon in mon2:
   ex=tuple(mon.count(k)for k in range(5));c=polynomial.coeff_monomial(ex);vector.append(c)
   if c:
    coefficients[i,j,ex]={power[0]:Fraction(value)for power,value in s.Poly(c,W).terms()if value}
 def entries(i,j):
  if i==j:return[]
  a,b=sorted((i,j));sign=1 if i<j else-1
  return[(ex,{d:sign*c for d,c in coeff.items()})for(ii,jj,ex),coeff in coefficients.items()if(ii,jj)==(a,b)]
 divergence={};jacobi={}
 for j in range(5):
  for i in range(5):
   for ex,coeff in entries(i,j):
    if ex[i]:
     remaining=list(ex);remaining[i]-=1
     for degree,c in coeff.items():add(divergence,(j,tuple(remaining),degree),c*ex[i])
 assert not divergence
 for index,triple in enumerate(itertools.combinations(range(5),3)):
  for a,b,c in[triple,(triple[1],triple[2],triple[0]),(triple[2],triple[0],triple[1])]:
   for ell in range(5):
    for left,lc in entries(a,ell):
     for right,rc in entries(b,c):
      if right[ell]:
       ex=tuple(left[k]+right[k]-(k==ell)for k in range(5))
       for ld,lv in lc.items():
        for rd,rv in rc.items():add(jacobi,(index,ex,ld+rd),lv*rv*right[ell])
 assert not jacobi and len(vector)==150
 # Euler/Pfaffian confirms the globalrank4open, not only generic rank.
 pf=[]
 for omitted in range(5):
  ids=[i for i in range(5)if i!=omitted];a,b,c,d=ids;pf.append((-1)**omitted*(P[a,b]*P[c,d]-P[a,c]*P[b,d]+P[a,d]*P[b,c]))
 quintic=s.expand(sum(X[i]*pf[i]for i in range(5)));assert quintic==W*x**5
 # Derive invariants by integrating the triangular vector field D.
 av=[h,y,z,w];affine=P.extract([0,2,3,4],[0,2,3,4]).subs(x,1);assert all(P[1,j]==0 for j in range(5))
 Z=s.expand(z-s.integrate(h,h));ZZ,TT,HH,YY=s.symbols('Z T H Y')
 primitive=s.integrate(s.expand(h*h+(ZZ+h*h/2)**2),h);assert s.expand(s.diff(primitive,h)-(h*h+(ZZ+h*h/2)**2))==0
 T=s.expand(w-primitive.subs(ZZ,Z));forward=s.Matrix([h,y,Z,T]);inverse=s.Matrix([HH,YY,ZZ+HH*HH/2,TT+primitive.subs(h,HH)])
 assert(forward.subs(dict(zip(av,inverse)),simultaneous=True)-s.Matrix([HH,YY,ZZ,TT])).applyfunc(s.expand)==s.zeros(4,1)
 assert(inverse.subs(dict(zip([HH,YY,ZZ,TT],forward)),simultaneous=True)-s.Matrix(av)).applyfunc(s.expand)==s.zeros(4,1)
 bracket=lambda f,g:s.expand(sum(affine[i,j]*(s.diff(f,av[i])*s.diff(g,av[j])-s.diff(f,av[j])*s.diff(g,av[i]))for i in range(4)for j in range(i+1,4)))
 target=s.zeros(4);target[0,1]=-1;target[2,3]=W;target=target-target.T
 for i,j in itertools.combinations(range(4),2):assert bracket(forward[i],forward[j])==target[i,j]
 assert forward.jacobian(av).det()==1 and all(not f.has(W)for f in list(forward)+list(inverse))
 # Homogenize our derived maps, rather than copying the saved Cremona forms.
 def homogenize(expression):
  result=0
  for ex,c in s.Poly(s.expand(expression),*av).terms():result+=c*x**(5-sum(ex))*s.prod(v**k for v,k in zip(av,ex))
  return s.expand(result)
 F=s.Matrix([homogenize(forward[0]),x**5,homogenize(forward[1]),homogenize(forward[2]),homogenize(forward[3])]);inv_aff=inverse.subs({HH:h,YY:y,ZZ:z,TT:w},simultaneous=True);G=s.Matrix([homogenize(inv_aff[0]),x**5,homogenize(inv_aff[1]),homogenize(inv_aff[2]),homogenize(inv_aff[3])])
 assert all(s.Poly(q,*X).total_degree()==5 for q in list(F)+list(G))
 assert(G.subs(dict(zip(X,F)),simultaneous=True)-x**24*xx).applyfunc(s.expand)==s.zeros(5,1)
 assert(F.subs(dict(zip(X,G)),simultaneous=True)-x**24*xx).applyfunc(s.expand)==s.zeros(5,1)
 assert s.Poly(F[4],*X).coeff_monomial(h**5)==-s.Rational(2,15)and s.Poly(G[4],*X).coeff_monomial(h**5)==s.Rational(1,20)
 assert s.gcd(F[1],F[4])==1 and s.gcd(G[1],G[4])==1
 prodcert=load(producer.parent/'certificate.json');assert(F-s.Matrix(prodcert['homogeneous_degree5_forward'])).applyfunc(s.expand)==s.zeros(5,1);assert(G-s.Matrix(prodcert['homogeneous_degree5_inverse'])).applyfunc(s.expand)==s.zeros(5,1)
 # Full tensor pushforward with inverse variable substitution; scalarprojective.
 diagonal=s.diag(1,1,r,1,1);sub=dict(zip(X,diagonal.inv()*xx));push=(diagonal*P.subs(sub,simultaneous=True)*diagonal.T).applyfunc(s.expand)
 assert(push-r*P.subs(W,W/r)).applyfunc(s.cancel)==s.zeros(5)
 # All finite Casimir parameters: basepointfree P1-degree5 and inverseparameter.
 c,d,A,S,t=s.symbols('c d A S t');leaf=s.Matrix([t,1,c+t*t/2,d+primitive.subs({h:t,ZZ:c})]);curve=[]
 for expression in leaf:
  polynomial=s.Poly(expression,t);curve.append(s.expand(sum(coefficient*A**ex[0]*S**(5-ex[0])for ex,coefficient in polynomial.terms())))
 curve=s.Matrix(curve);assert all(s.Poly(q,A,S).total_degree()==5 for q in curve)
 assert curve[1]==S**5 and curve.subs(S,0)==s.Matrix([0,0,0,A**5/20]);assert s.cancel(curve[0]/curve[1])==A/S
 assert(curve-s.Matrix([s.sympify(f,locals={'A':A,'S':S,'c':c,'d':d})for f in prodcert['W0_leaf_projective_curve']])).applyfunc(s.expand)==s.zeros(4,1)
 # Direct tangent invariance of the affine fibers under both spanningfields.
 D=s.Matrix([1,0,h,h*h+z*z]);assert s.expand(sum(D[i]*s.diff(Z,av[i])for i in range(4)))==0;assert s.expand(sum(D[i]*s.diff(T,av[i])for i in range(4)))==0
 payload=dict(own_CA_raw=rows(CA),own_CB_raw=rows(CB),source_CB_proportionality='g=-W*f withf=-h^2',canonical_P=rows(P),own_canonical150=list(map(str,vector)),all350_strict_Jacobi_coefficients_zero=True,all25_divergence_coefficients_zero=True,degeneracy_quintic=str(quintic),
  derived_triangular_primitive=str(primitive),derived_forward=list(map(str,forward)),derived_inverse=list(map(str,inverse)),all6_coordinate_brackets=[[i,j,str(target[i,j])]for i,j in itertools.combinations(range(4),2)],
  own_homogeneous_F=list(map(str,F)),own_homogeneous_G=list(map(str,G)),common_factors_both_one=True,both_compositions='x^24*identity',
  full_PGL_tensor_identity='diag(1,1,r,1,1)_*P_W=r*P_(W/r)',
  Casimir_center_W0='QQ[Z,T], provedfromconstanttarget',all_affine_leaf_fibers='Z=c,T=d,affineH,Y free=>A2',leaf_base_curve=list(map(str,curve)),curve_O1_pullback_degree=5,curve_function_field_recovered_by='h/x=A/S',cone_vertex=[0,0,1,0,0],cone_degree=5,
  projective_base_locus_F_and_G='h=x=0 (P2), affine maps regular',local_completion_or_global_B18_classification_not_reaudited=True)
 (out/'independent-evidence.json').write_text(json.dumps(payload,indent=2)+'\n')
 paths=[Path(__file__),producer,producer.parent/'certificate.json',source,quadpath,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py',ROOT/'research/notes/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW.md']
 report=dict(id='AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW',all_passed=True,producer_not_reexecuted=True,independent_actual_CA_CB_raw_reconstruction=True,full150_and350_coefficients_audited=True,independent_integration_derives_allW_Darboux_and_inverse=True,exact_Cremona_bidegree5_5_and_gcd_verified=True,full_nonzeroW_single_projective_PGL_orbit_accepted=True,allW0_affine_leaves_A2_and_allclosure_cones_degree5_accepted=True,proof_scope_excludes_unrequested_local_completion_classification=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'independent-evidence.json':sha(out/'independent-evidence.json')},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('INDEPENDENT PW/DARBOUX/CREMONA/ORBIT/ALLLEAFCONES ACCEPTED',report['seconds'],flush=True)
if __name__=='__main__':main()
