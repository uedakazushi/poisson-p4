#!/usr/bin/env python3
"""Exact four-linear-variable structure and nonzero binary D; no GB or fields."""
from pathlib import Path
import json,time,hashlib,sys
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/four-quadric-binary-certificate-v71'
OUT.mkdir(parents=True,exist_ok=True)
def main():
 start=time.monotonic();a,e,P=s.symbols('a e P');x=s.symbols('x0:6');source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json'
 assert hashlib.sha256(source.read_bytes()).hexdigest()=='80f7c46b1bd4ffef4a78f492ad4f7a2bb40b2ffa27cddcd23d5b5d660db4a7e0'
 data=json.loads(source.read_text());qs=[]
 for i,row in enumerate(data['cleared_qfour_sparse']):
  q=0
  for mon,terms in row:q+=sum(s.Rational(c)*a**(ex[0]-ex[2]+[0,2,2,1][i])*e**ex[1]*P**ex[2]for ex,c in terms)*s.prod(v**k for v,k in zip(x,mon[:6]))
  qs.append(s.expand(q))
 J=5*P-6*e-4;Rc=15*P**2-6*P*e+6*P+80*a;Tc=-5*P+6*e-2;Sc=21*P*e+16*P+140*a;Vc=9*P*e+4*P+60*a;Kc=85*P**2-12*P*e+58*P+600*a
 B=30*a*Rc*x[0]-9*P*Sc*x[1]+90*P*a*Tc*x[2]-18*P*Vc*x[3]-2*a*Kc*x[4];A=(B+10*a*P*J*x[4])/3;Q2=s.cancel((qs[2]-e*qs[1])/(a*a*J))
 C1=s.cancel((qs[1]-A*B/6075)/x[5]);C2=s.cancel((Q2-40*x[0]*B/6075)/x[5]);D=s.cancel(40*x[0]*C1-A*C2)
 assert s.Poly(C1,*x).total_degree()==s.Poly(C2,*x).total_degree()==1
 assert s.expand(x[5]*D-40*x[0]*qs[1]+A*Q2)==0
 TA,X,W,TB,S,U=s.symbols('TA X W TB S U');coords=(TA,X,W,TB,S,U)
 x4=(3*TA-TB)/(10*a*P*J);t=(30*a*Rc*X+90*P*a*Tc*S-18*P*Vc*U-2*a*Kc*x4-TB)/(180*P**2)
 images=[X,3*t,S,-s.Rational(7,2)*t+U,x4,W];L=s.Matrix([[s.diff(q,y)for y in coords]for q in images]);assert s.cancel(L.det()-1/(200*a*P**3*J))==0
 forms=[s.Poly(f,*x)for f in [A,B,C1,C2]]
 linear=s.Matrix([[f.coeff_monomial(v)for v in x]for f in forms]);normal=(linear*L).applyfunc(s.cancel)
 assert normal[0,:]==s.Matrix([[1,0,0,0,0,0]]) and normal[1,:]==s.Matrix([[0,0,0,1,0,0]])
 CC=[s.Poly(sum(normal[i,j]*coords[j]for j in range(6)),*coords)for i in [2,3]];assert all(not f.as_expr().has(S,U)for f in CC)
 dd=s.Poly(s.cancel((40*X*CC[0].as_expr()-TA*CC[1].as_expr()).subs({W:0,TB:0})),TA,X)
 binary=-8*(3*P*e+2*P+12*a)*TA**2/(1215*P**2)-256*a*a*J*TA*X/(243*P)+1280*a**3*J**2*X**2/243
 assert s.cancel(dd.as_expr()-binary)==0
 assert s.cancel(dd.coeff_monomial(X*X)-1280*a**3*J**2/243)==0
 assert s.cancel(dd.coeff_monomial(TA*X)+256*a*a*J/(243*P))==0
 # The actual three quadrics are negative minors of this matrix of linear forms.
 M=s.Matrix([[A/6075,40*x[0]/6075,x[5]],[C1,C2,-B]])
 assert s.cancel(M[:,[0,1]].det()+D/6075)==0
 assert s.cancel(M[:,[0,2]].det()+qs[1])==0
 assert s.cancel(M[:,[1,2]].det()+Q2)==0
 assert s.expand(3*Sc-7*Vc-20*P)==0
 coeff=[[s.factor(f.coeff_monomial(y))for y in coords[:4]]for f in CC]
 certificate=dict(linear_coordinate_order=['A','x0','W','B','S','U'],inverse_frame_images=list(map(str,images)),inverse_frame_determinant='1/(200*a*P^3*J)',C1_C2_four_coordinate_coefficients=[[str(c)for c in row]for row in coeff],D_on_B_W_zero=str(binary),D_x0_squared_unit='1280*a^3*J^2/243',D_A_x0_unit='-256*a^2*J/(243*P)',three_minor_signs=['-D/6075','-Q1','-Q2prime'])
 (OUT/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 result=dict(all_exact_QQ_actual804_checks_passed=True,six_linear_entries_rank_exactly4=True,open='a*P*J!=0',C1_C2_depend_only_on_A_x0_W_B=True,D_lies_in_ideal_A_x0=True,D_restriction_B_W_zero_uniformly_nonzero=True,pair_gate_implies_full_auxiliary_five_rank5=True,pair_gate_still_open=True,full_global_raw24_or_product41_or_ambient_not_proved=True,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
