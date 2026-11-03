#!/usr/bin/env python3
"""Exact rank5 of the auxiliary quadrics on DB=0,a*P*J!=0; no GB.

The 804 actual coefficient terms are reconstructed before any frame operation.
All outputs go to a new work directory; prior frozen packets are read-only.
"""
from pathlib import Path
import sympy as s
import json,hashlib,time,sys,itertools
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/five-rank-DB0-certificate-v68'
OUT.mkdir(parents=True,exist_ok=True)
def main():
 start=time.monotonic();a,e,P=s.symbols('a e P');x=s.symbols('x0:6');ys=s.symbols('r U V T R');r,U,V,T,R=ys
 source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';assert hashlib.sha256(source.read_bytes()).hexdigest()=='80f7c46b1bd4ffef4a78f492ad4f7a2bb40b2ffa27cddcd23d5b5d660db4a7e0'
 data=json.loads(source.read_text());qs=[]
 for i,row in enumerate(data['cleared_qfour_sparse']):
  q=0
  for mon,terms in row:q+=sum(s.Rational(c)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,c in terms)*s.prod(v**k for v,k in zip(x,mon[:6]))
  qs.append(s.expand(q))
 J=5*P-6*e-4;DB=20*P**2+3*P*e+20*P+180*a;L=20*P+3*e+20;ad=-P*L/180;K=-130*P+168*e-5;ed=(130*P+5)/168
 Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a
 B=30*a*Rc*x[0]-9*P*Sc*x[1]+90*P*a*Tc*x[2]-18*P*Vc*x[3]-2*a*Kc*x[4];A=(B+10*a*P*J*x[4])/3
 Q2=s.cancel((qs[2]-e*qs[1])/(a*a*J));D=s.cancel((40*x[0]*qs[1]-A*Q2)/x[5]);s.Poly(D,*x,a,e,P,domain=s.QQ)
 assert s.expand(qs[1].subs(x[5],0)-A*B/6075)==0
 assert s.expand(Q2.subs(x[5],0)-40*x[0]*B/6075)==0
 assert s.expand(3*Sc-7*Vc-20*P)==0
 images=[r,V,R/a-U-6*(e+1)*V/(5*a)-(9*e+4)*T/(10*a),(T-V)/2,5*U+5*r-4*R/a+24*(e+1)*V/(5*a)+2*(9*e+4)*T/(5*a),s.Integer(0)]
 lin=s.Matrix([[s.diff(image,y)for y in ys]for image in images]);H=[]
 for q in [qs[0],qs[3],D]:
  h=(a*a*lin.T*(s.hessian(q,x)/2)*lin).applyfunc(s.cancel)
  for value in h:s.Poly(value,a,e,P,domain=s.QQ)
  H.append(h)
 expected=s.zeros(5);expected[1,4]=expected[4,1]=a*a/8;expected[2,3]=expected[3,2]=-5*a*a/8;assert H[0]==expected
 BF=s.Poly(s.expand(B.subs(dict(zip(x,images)),simultaneous=True)),*ys);b=[BF.coeff_monomial(y)for y in ys];assert s.expand(b[0]+20*a*DB)==0
 assert s.cancel(H[1][0,0].subs(a,ad))==0
 assert s.cancel(H[2][0,0].subs(a,ad)+P**5*L**4*J**2*K/s.Integer(35872267500))==0
 # Reverify the earlier exact certificate that (b,h3) has rank2 throughout DB0.
 prop=ROOT/'work/nonsquare-v62/gcd-proportional-certificate-v58/certificate.json';old=json.loads(prop.read_text());load2=lambda terms:sum(s.Rational(c)*P**ex[0]*e**ex[1]for ex,c in terms)
 G=[load2(q)for q in old['normalized_minors']];cofs=[load2(q)for q in old['Bezout_coefficients']]
 b0=[s.cancel(q.subs(a,ad))for q in b[1:]];h3=[s.cancel(H[1][0,j].subs(a,ad))for j in range(1,5)]
 for gi,ij,powers,scale in zip(G,old['corresponding_cross_minor_indices'],old['known_unit_powers_P_L_J'],old['normalization_QQ_scalars']):
  i,j=ij;factor=s.prod(v**k for v,k in zip((P,L,J),powers));assert s.cancel(s.Rational(scale)*(h3[i]*b0[j]-h3[j]*b0[i])-factor*gi)==0
 assert s.expand(sum(c*g for c,g in zip(cofs,G))-1)==0
 # On K=0, r^2 vanishes in both Q3,D. Their r-linear classes are independent modulo B.
 sub=lambda q:s.cancel(q.subs(a,ad).subs(e,ed));assert sub(H[1][0,0])==sub(H[2][0,0])==0
 cross=s.Matrix([[sub(q)for q in b[1:]],[sub(H[1][0,j])for j in range(1,5)],[sub(H[2][0,j])for j in range(1,5)]]);raw=[];primitive=[];factors=[]
 for cols in itertools.combinations(range(4),3):
  det=s.Poly(cross[:,cols].det(),P,domain=s.QQ);raw.append(det);q=det;powers=[]
  for unit in [P,L.subs(e,ed),J.subs(e,ed)]:
   up=s.Poly(unit,P,domain=s.QQ);val=0
   while not q.is_zero:
    quo,rem=q.div(up)
    if rem:break
    q=quo;val+=1
   powers.append(val)
  content,q=q.primitive();assert not q.is_zero;primitive.append(q);factors.append(dict(columns=list(cols),known_unit_powers_P_L_J=powers,QQcontent=str(content)))
 coeff=[s.Poly(1,P,domain=s.QQ)];g=primitive[0]
 for q in primitive[1:]:
  u,v,h=s.gcdex(g,q);coeff=[s.Poly(c*u,P,domain=s.QQ)for c in coeff]+[s.Poly(v,P,domain=s.QQ)];g=s.Poly(h,P,domain=s.QQ)
 assert g==s.Poly(1,P,domain=s.QQ);assert sum((c*q for c,q in zip(coeff,primitive)),s.Poly(0,P,domain=s.QQ))==g
 enc=lambda q,vs:[[list(ex),str(c)]for ex,c in s.Poly(q,*vs,domain=s.QQ).terms()]
 certificate=dict(global_H_matrices_a2_Q0_Q3_D=[[[enc(q,(a,e,P))for q in row]for row in h.tolist()]for h in H],global_B_coefficients=[enc(q,(a,e,P))for q in b],DB0_D_rr='-P^5*L^4*J^2*K/35872267500',K_zero_e=str(ed),K_zero_L=str(s.factor(L.subs(e,ed))),K_zero_J=str(s.factor(J.subs(e,ed))),K_zero_cross_matrix=[[enc(q,(P,))for q in row]for row in cross.tolist()],cross_minors=[enc(q.as_expr(),(P,))for q in raw],minor_factor_records=factors,primitive_minors=[enc(q.as_expr(),(P,))for q in primitive],Bezout_coefficients=[enc(q.as_expr(),(P,))for q in coeff])
 (OUT/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 result=dict(all_exact_QQ_checks_passed=True,actual804_terms_reconstructed=True,DB_zero_auxiliary_five_rank_exactly5=True,open='a*P*J!=0',no_FBC_needed=True,rank_proof='Q1,Q2 vanish on B=W=0 and are independent; Q0,Q3,D restricted independent by rr/r-linear cases',DB_nonzero_gate_not_proved=True,pair_gate_not_proved=True,no_global_raw24_or_product41_or_ambient_claim=True,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),prior_rank2_Bezout_source_sha256=hashlib.sha256(prop.read_bytes()).hexdigest(),seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
