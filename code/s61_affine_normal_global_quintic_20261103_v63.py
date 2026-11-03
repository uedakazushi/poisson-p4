#!/usr/bin/env python3
"""Global monogenic quintic and exact5x670 normal-base module recipe."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
import itertools
import numpy as np
from pathlib import Path
import time
import sympy as s
from s61_affine_fiber_macaulay import load_fiber
from s61_affine_polarization_splitter_20261103_v63 import determinant_including_zero
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(poly,mon,c):
 c=poly.get(mon,Fraction(0))+c
 if c:poly[mon]=c
 else:poly.pop(mon,None)
def padd(a,b,c=Fraction(1)):
 out=dict(a)
 for m,v in b.items():add(out,m,c*v)
 return out
def pmul(a,b):
 out={}
 for m,c in a.items():
  for n,d in b.items():add(out,tuple(sorted(m+n)),c*d)
 return out
def uadd(a,b,c=Fraction(1)):
 out={k:dict(v)for k,v in a.items()}
 for k,p in b.items():
  out[k]=padd(out.get(k,{}),p,c)
  if not out[k]:out.pop(k)
 return out
def umul(a,b):
 out={}
 for i,p in a.items():
  for j,q in b.items():out[i+j]=padd(out.get(i+j,{}),pmul(p,q))
 return {k:p for k,p in out.items()if p}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 dep=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/report.json';old=json.loads(dep.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(dep.parent/rel)==want,rel
 specpath=ROOT/'code/data/s61_affine_normal_global_quintic_20261103_v63.json';spec=json.loads(specpath.read_text());z=s.symbols('z0:3');plane=list(map(s.sympify,spec['plane']));relations=sorted(spec['relations'],key=lambda v:v['target_index']);targets=[]
 for data in relations:
  target=s.sympify(data['target']);cof=list(map(s.sympify,data['cofactors']));assert s.expand(sum(a*b for a,b in zip(cof,plane))-target)==0
  assert all(s.denom(c)%32003 for f in cof for c in s.Poly(f,*z).coeffs())
  assert max(s.Poly(f,*z).total_degree()for f in cof)<=data['degree']-2;targets.append(target)
 minor_checks=[]
 for degree in[4,5]:
  record=next(v for v in relations if v['degree']==degree);words=[m for d in range(degree+1)for m in itertools.combinations_with_replacement(range(3),d)];multipliers=[m for d in range(degree-1)for m in itertools.combinations_with_replacement(range(3),d)];monex=lambda m:tuple(m.count(k)for k in range(3));matrix=[]
  for f in plane:
   poly=s.Poly(f,*z)
   for m in multipliers:
    shift=monex(m);row=[]
    for word in words:
     targetex=monex(word);ex=tuple(a-b for a,b in zip(targetex,shift));row.append(poly.coeff_monomial(ex)if min(ex)>=0 else 0)
    matrix.append(row)
  square=[[matrix[i][j]for j in record['selected_columns']]for i in record['selected_rows']];modular=np.array([[int(s.numer(c)%32003)*pow(int(s.denom(c)%32003),-1,32003)%32003 for c in row]for row in square],dtype=np.int64);det=determinant_including_zero(modular);assert det;minor_checks.append(dict(degree=degree,rank=len(square),determinant_mod32003=det))
 p1=s.expand(z[1]-targets[0]);p2=s.expand(z[2]-targets[1]);g=s.Poly(targets[2],z[0]);assert not(p1.free_symbols|p2.free_symbols)-{z[0]};assert g.degree()==5 and g.LC()==1
 F,qids,center,original,pointpath=load_fiber();packet=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz';frame=json.loads(gzip.decompress(packet.read_bytes()));E=[[Fraction(c)for c in row]for row in frame['embedding']]
 for index,polynomial in enumerate(F[:5]):
  restriction=s.expand(sum(s.Rational(c)*s.prod(sum(s.Rational(E[j][k])*z[k]for k in range(3))for j in mon)for mon,c in polynomial.items()));assert s.expand(restriction-plane[index])==0
 # Coefficients are polynomials in n0..41,S1,S2,H. The three final
 # placeholders are NOT independent parameters: their exact normal polynomial
 # recipes are the cofactor reductions stated in the certificate and proof.
 def univariate(expr):return {int(mon[0]):{():Fraction(c)}for mon,c in s.Poly(expr,z[0]).terms()if c}
 Z=[{1:{():Fraction(1)}},uadd(univariate(p1),{0:{(42,):Fraction(1)}}),uadd(univariate(p2),{0:{(43,):Fraction(1)}})];G=univariate(g.as_expr());G=uadd(G,{0:{(44,):Fraction(-1)}});assert G[5]=={():Fraction(1)}
 Q=[]
 for i,row in enumerate(E):
  qi={}
  for k,c in enumerate(row):qi=uadd(qi,Z[k],c)
  if i>=3:qi=uadd(qi,{0:{(i-3,):Fraction(1)}})
  Q.append(qi)
 def divide(poly):
  remainder={k:dict(v)for k,v in poly.items()};quotient={}
  while remainder and max(remainder)>=5:
   degree=max(remainder);coeff=remainder[degree];term={degree-5:coeff};quotient=uadd(quotient,term);remainder=uadd(remainder,umul(term,G),Fraction(-1))
  assert uadd(umul(quotient,G),remainder)==poly
  return quotient,remainder
 products={};remainders=[];quotients=[]
 for polynomial in F:
  image={}
  for mon,c in polynomial.items():
   if mon not in products:
    term={0:{():Fraction(1)}}
    for i in mon:term=umul(term,Q[i])
    products[mon]=term
   image=uadd(image,products[mon],c)
  quotient,remainder=divide(image);assert max(remainder,default=-1)<=4;remainders.append(remainder);quotients.append(quotient)
 assert len(remainders)==134
 encode=lambda poly:[[degree,[[list(mon),str(c)]for mon,c in sorted(coeff.items())]]for degree,coeff in sorted(poly.items())]
 payload=dict(coefficient_symbols=[*[f'n{i}'for i in range(42)],'S1(n)','S2(n)','H(n)'],
  h='z0',p1=str(p1),p2=str(p2),g=str(g.as_expr()),monic_quintic_G=encode(G),all134_remainders=[encode(r)for r in remainders],all134_division_quotients=[encode(q)for q in quotients],
  coefficient_placeholders_not_free_parameters=True,normal_recipe='LetK_i(n)=Red(F_i)-f_i(z), Red(zk*nj)=P_kj(n). S_l=-Red(sum_j B_lj(z)*K_j(n));H=-Red(sum_j B_gj(z)*K_j(n)). Cofactors B are pinned code/data.',
  correction_degree_bounds=dict(S1=4,S2=4,H=5),
  full_module_presentation='C=R[h]/(g(h)-H(n)) isFREE5;actualA=C/(all134remainders). ColumnsT_h^j*r_i fori0..133,j0..4 giveactual5x670Rmodulepresentation.',
  actual_A_not_claimed_R_free=True)
 (out/'presentation.json.gz').write_bytes(gzip.compress(json.dumps(payload,separators=(',',':')).encode(),mtime=0))
 cert=dict(selected_degree45_macaulay_minors=minor_checks,all3_exact_QQ_plane_cofactor_identities=True,cofactor_terms=[sum(len(s.Poly(s.sympify(f),*z).terms())for f in d['cofactors'])for d in relations],p_integral_cofactors_at32003=True,
  all134_exact_formal_univariate_division_identities=True,all134_h_degree_atmost4=True,formal_remainder_coefficients=sum(len(p)for r in remainders for p in r.values()),
  actual_global_isomorphism='A=R[h]/(g(h)-H(n),r0(h,n),...,r133(h,n))',actual_R_module_shape_atmost=[5,670],no_cyclic_q_or_H_open_localization=True,
  generic_recipe='UseuniversalK-rationalthreeplaneframe,primitivehknownE5,anddegree4/5Macaulaycofactorinverses(selectedatfrozencenter);applyfinitecrossnormalreductionthenmonicdivision.',
  polynomial_residual_recovery='D84=e;1-D84=sum_i c_i(n)n_i byfirstnormalfactoring;Z_k(n)=sum_i c_i(n)P_ki(n),degree<=85;onB,z_k=Z_k(n).',
  free5_ambient_not_actual_freeness=True,beta_only_finiteness_or_unit_not_claimed=True,expanded_normal_corrections_not_serialized=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),specpath,dep,packet,original,pointpath,ROOT/'research/notes/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63.md']
 report=dict(id='S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63',all_passed=True,actual_global_monogenic_quintic_over_normal42=True,all134_h_reductions_saved_exactly=True,module_presentation_shape_atmost=[5,670],unknown_support_and_nilpotents_retained=True,no_cyclic_open_or_actual_freeness=True,unit_or_genericfive_not_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('GLOBAL NORMAL QUINTIC/ALL134EXACTDIVISIONS/5x670MODULE PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
