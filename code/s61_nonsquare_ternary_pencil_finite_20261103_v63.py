#!/usr/bin/env python3
"""Exact QQ finite necessary locus for a common factor in the actual ternary pair.

Every leading-coefficient case is retained. Verification uses saved polynomial
identities and good-specialization/content certificates, with no GB or fields.
"""
from pathlib import Path
import json,subprocess,sys,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'work/nonsquare-v62/quad78-pencil-pseudo-v84'
PROFILE=ROOT/'work/nonsquare-v62/quad78-small-gates-v83/pencil-profile.json'
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/ternary-pencil-finite-replay-v84'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 start=time.monotonic();t,H,v,z=s.symbols('t H v z');vars=(t,H,v)
 load=lambda terms,vs:s.Poly(sum(s.Rational(c)*s.prod(x**k for x,k in zip(vs,ex))for ex,c in terms),*vs,domain=s.QQ)
 with(OUT/'actual-source.log').open('w')as log:subprocess.run([sys.executable,str(ROOT/'code/s61_nonsquare_pair_empty_QQ_20261103_v63.py'),str(OUT/'actual-source')],check=True,stdout=log)
 source=json.loads((OUT/'actual-source/finite-source/certificate.json').read_text());rows=[[load(q,(t,v,H)).as_expr()for q in row]for row in source['natural_rows']]
 def symmetric(row):return s.Matrix([[row[0],row[1]/2,row[3]/2],[row[1]/2,row[2],row[4]/2],[row[3]/2,row[4]/2,row[5]]])
 A,B=list(map(symmetric,rows));DM=DomainMatrix.from_Matrix(A+z*B);pencil=s.Poly(DM.domain.to_sympy(DM.det()),z);assert pencil.degree()==3
 profile=json.loads(PROFILE.read_text());Ci=[load(q,(t,v,H)).reorder(t,H,v)for q in profile['normalized_pencil_coefficients']];P=H-8*t;J=5*H-2*v-4
 for i,record in enumerate(profile['normalization_records']):
  p,q,r=record['known_unit_powers_t_P_J'];actual=s.Poly(pencil.coeff_monomial(z**i),*vars,domain=s.QQ);assert actual==Ci[i]*s.Poly(s.Rational(record['QQcontent'])*t**p*P**q*J**r,*vars,domain=s.QQ)
 lead=json.loads((DATA/'leading-coefficients.json').read_text());LC,MID,CONST=[load(lead[key],(H,v))for key in ['LC','MID','CONST']];assert Ci[0]==s.Poly(LC.as_expr()*t*t+MID.as_expr()*t+CONST.as_expr(),*vars,domain=s.QQ)
 remainders=[];ab=[]
 for i in [1,2,3]:
  item=json.loads((DATA/f'remainder-{i}.json').read_text());Q=load(item['pseudo_quotient'],vars);raw=load(item['pseudo_remainder'],vars);assert s.Poly(LC.as_expr()**item['LC_power'],*vars,domain=s.QQ)*Ci[i]==Q*Ci[0]+raw
  rem=load(item['normalized_remainder'],vars);p,q,r,l=item['removed_units_t_P_J_LC'];assert raw==rem*s.Poly(s.Rational(item['QQcontent'])*t**p*P**q*J**r*LC.as_expr()**l,*vars,domain=s.QQ)
  pt=s.Poly(rem.as_expr(),t);assert pt.degree()<=1;ab.append((s.Poly(pt.coeff_monomial(t),H,v,domain=s.QQ),s.Poly(pt.coeff_monomial(1),H,v,domain=s.QQ)));remainders.append(rem)
 N=[]
 for i,j in [(0,1),(0,2)]:
  item=json.loads((DATA/f'eliminate-t-{i}-{j}.json').read_text());q=load(item['polynomial'],(H,v));p,l=item['removed_J_LC'];assert ab[i][0]*ab[j][1]-ab[j][0]*ab[i][1]==q*s.Poly(s.Rational(item['QQcontent'])*J**p*LC.as_expr()**l,H,v,domain=s.QQ);N.append(q)
 cert=json.loads((DATA/'finite-certificate.json').read_text())
 def check_coprime(data,F,G):
  assert load(data['F'],(H,v))==F and load(data['G'],(H,v))==G
  coefficients=[]
  for idx,degree in data['coefficient_labels']:coefficients.append(s.Poly(s.Poly([F,G][idx].as_expr(),H).coeff_monomial(H**degree),v,domain=s.QQ))
  weights=[load(q,(v,))for q in data['coefficient_content_Bezout']];assert sum((a*b for a,b in zip(weights,coefficients)),s.Poly(0,v,domain=s.QQ))==s.Poly(1,v,domain=s.QQ)
  value=data['special_v'];f=s.Poly(F.as_expr().subs(v,value),H,domain=s.QQ);g=s.Poly(G.as_expr().subs(v,value),H,domain=s.QQ);assert f.degree()==F.degree(H)==data['F_H_degree'] and g.degree()==G.degree(H)==data['G_H_degree'];assert f.LC() and g.LC()
  u,w=[load(q,(H,))for q in data['specialized_QQ_Bezout']];assert u*f+w*g==s.Poly(1,H,domain=s.QQ)
 check_coprime(cert['generic_LC_nonzero_projection'],*N)
 # Leading LC=0, MID!=0 branch. All divisions here are chart units, never global.
 lin=json.loads((DATA/'LC-zero-MID-nonzero.json').read_text());Rlin=load(lin['linear_pseudo_remainder'],(H,v));p,q=lin['removed_J_MID'];raw=s.Poly(s.Rational(lin['QQcontent'])*J**p*MID.as_expr()**q,*vars,domain=s.QQ)*s.Poly(Rlin.as_expr(),*vars,domain=s.QQ);quot=load(lin['quotient'],vars);degree=Ci[1].degree(t)
 assert s.Poly(MID.as_expr()**degree,*vars,domain=s.QQ)*Ci[1]==quot*s.Poly(MID.as_expr()*t+CONST.as_expr(),*vars,domain=s.QQ)+raw
 check_coprime(cert['LC_zero_MID_nonzero_projection'],LC,Rlin)
 # Leading LC=MID=0 branch; two small fixed Sylvester determinants constrain v.
 def resultant(f,g):
  ff=s.Poly(f,H);gg=s.Poly(g,H);m,n=ff.degree(),gg.degree();aa=ff.all_coeffs();bb=gg.all_coeffs();M=s.Matrix([[0]*i+aa+[0]*(n-i-1)for i in range(n)]+[[0]*i+bb+[0]*(m-i-1)for i in range(m)]);D=DomainMatrix.from_Matrix(M);return s.Poly(D.domain.to_sympy(D.det()),v,domain=s.QQ)
 last=cert['LC_MID_zero_case'];R1=load(last['resultant_R1'],(v,));R2=load(last['resultant_R2'],(v,));assert R1==resultant(LC.as_expr(),MID.as_expr());assert R2==resultant(LC.as_expr(),CONST.as_expr())
 u,w=[load(q,(v,))for q in last['resultant_gcd_Bezout']];g=load(last['resultant_gcd'],(v,));assert u*R1+w*R2==g==s.Poly((v+3)*(3*v+4)/3,v,domain=s.QQ)
 for branch in last['special_branches']:
  vv=s.Rational(branch['v']);hh=s.Rational(branch['forced_H']);leading=[s.Poly(p.as_expr().subs(v,vv),H,domain=s.QQ)for p in [LC,MID,CONST]];weights=[load(q,(H,))for q in branch['leading_Bezout']];gg=load(branch['leading_H_polynomial'],(H,));assert sum((c*q for c,q in zip(weights,leading)),s.Poly(0,H,domain=s.QQ))==gg;assert gg.degree()==1 and gg.eval(hh)==0
  higher=[s.Poly(c.as_expr().subs({v:vv,H:hh}),t,domain=s.QQ)for c in Ci[1:]];assert higher==[load(q,(t,))for q in branch['higher_t_polynomials']];cofs=[load(q,(t,))for q in branch['higher_t_Bezout']];assert sum((c*q for c,q in zip(cofs,higher)),s.Poly(0,t,domain=s.QQ))==s.Poly(1,t,domain=s.QQ)
 result=dict(all_actual_QQ_identities_passed=True,all_singular_ternary_pencil_locus_dimension_at_most_zero=True,ternary_pair_commonfactor_locus_dimension_at_most_zero=True,open='a*P*J!=0',generic_LC_nonzero_and_LC_zero_MID_unit_projections_finite=True,all_LC_MID_CONST_zero_branches_excluded=True,no_unknown_leading_coefficient_inverted_globally=True,QQ_finite_theorem_only=True,no_quad78_or_full_quartic_or_ambient_claim=True,proper_finite_arithmetic_model_not_constructed=True,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
