#!/usr/bin/env python3
"""A rank5 exact presentation of critical affine support on the quintic open.

General coefficients are specified by existing finite inverse circuits. This
producer checks the NEW original-detE norm and its complete rank5 Bezout
certificate at the frozen center; it does not assert complement exhaustion.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import sympy as s
from s61_affine_fiber_macaulay import load_fiber
from s61_nonsquare_cubic_row_basis_review_20261103_v62 import determinant_mod
ROOT=Path(__file__).resolve().parents[1];P=32003

def trim(v):
 v=[int(c)%P for c in v]
 while v and not v[-1]:v.pop()
 return v

def mul(f,g):
 h=[0]*(len(f)+len(g)-1)
 for i,a in enumerate(f):
  for j,b in enumerate(g):h[i+j]=(h[i+j]+a*b)%P
 return trim(h)

def divmodp(f,g):
 f=trim(f);g=trim(g);q=[0]*max(0,len(f)-len(g)+1);inv=pow(g[-1],-1,P)
 while len(f)>=len(g):
  d=len(f)-len(g);c=f[-1]*inv%P;q[d]=c
  for j,b in enumerate(g):f[d+j]=(f[d+j]-c*b)%P
  f=trim(f)
 return trim(q),f

def evalp(f,x):
 v=0
 for c in reversed(f):v=(v*x+c)%P
 return v

def add(f,g):
 h=[0]*max(len(f),len(g))
 for i,v in enumerate(f):h[i]+=v
 for i,v in enumerate(g):h[i]+=v
 return trim(h)

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();critical=ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62/critical-determinant.json';cyclic=ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62/report.json';quintic=ROOT/'reports/S61_AFFINE_QUINTIC_OPEN_20261103_v62/report.json';qd=json.loads(critical.read_text());cy=json.loads(cyclic.read_text());qr=json.loads(quintic.read_text());assert cy['all_passed']and qr['all_passed'];C=np.array(qd['constant_row_combination'],dtype=np.int64)
 lines=(quintic.parent/'native.txt').read_text().splitlines();g=list(map(int,next(v for v in lines if v.startswith('GCD ')).split()[1:]));roots=cy['known_five_h_values'];assert len(set(roots))==5 and all(evalp(g,h)==0 for h in roots)
 product=[1]
 for h in roots:product=mul(product,[-h,1])
 assert product==g
 original,qids,center,*_=load_fiber();knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json';known=json.loads(knownpath.read_text())['known_five_mod_32003'];nativepath=cyclic.parent/'native-certificate.txt';nl=nativepath.read_text().splitlines();coordinates={int(line.split()[1]):list(map(int,line.split()[2:]))for line in nl if line.startswith('COORDINATE ')};chi=list(map(int,next(line for line in nl if line.startswith('CHI ')).split()[1:]));assert len(coordinates)==29 and not divmodp(chi,g)[1]
 for h,point in zip(roots,known):assert [evalp(coordinates[i],h)for i in range(29)]==point+[0]*25
 plane=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';pd=json.loads(plane.read_text());vv=s.symbols('v0:3')
 def residue(v):
  v=s.Rational(v);assert int(v.q)%P;return int(v.p)*pow(int(v.q),-1,P)%P
 slopes=np.array([[residue(s.sympify(v).coeff(j))for j in vv]for v in pd['affine_plane_substitution']],dtype=np.int64);values=[];qpoints=[]
 for point in known:
  q=slopes@np.array(point[1:],dtype=np.int64)%P;qpoints.append(q.tolist());J=np.zeros((134,45),dtype=np.int64)
  for i,f in enumerate(original):
   total=0
   for mon,c in f.items():
    co=residue(c);term=co
    for j in mon:term=term*int(q[j])%P
    total=(total+term)%P
    for k,j in enumerate(mon):
     derivative=co
     for ell,v in enumerate(mon):
      if k!=ell:derivative=derivative*int(q[v])%P
     J[i,j]=(J[i,j]+derivative)%P
   assert total==0
  values.append(determinant_mod((C@J%P).tolist(),P))
 assert values==qd['known_five_determinants_mod32003']
 e=[]
 for i,h in enumerate(roots):
  basis=[1];den=1
  for j,k in enumerate(roots):
   if j!=i:basis=mul(basis,[-k,1]);den=den*(h-k)%P
  e=add(e,[(v*values[i]*pow(den,-1,P))%P for v in basis])
 assert len(e)<=5 and [evalp(e,h)for h in roots]==values
 z=s.Symbol('s');ep=s.Poly(sum(v*z**i for i,v in enumerate(e)),z,modulus=P);gp=s.Poly(sum(v*z**i for i,v in enumerate(g)),z,modulus=P);u,v,d=s.gcdex(ep,gp);assert d==s.Poly(1,z,modulus=P)
 coeff=lambda f:trim([int(f.nth(i))%P for i in range(max(0,f.degree()+1))]);eu,gv=coeff(u),coeff(v);assert add(mul(e,eu),mul(g,gv))==[1]
 M=[]
 for j in range(5):
  rem=divmodp([0]*j+e,g)[1];M.append(rem+[0]*(5-len(rem)))
 matrix=np.array(M,dtype=np.int64).T;det=determinant_mod(matrix.tolist(),P);assert det==9485
 invcols=[]
 for j in range(5):
  rem=divmodp([0]*j+eu,g)[1];invcols.append(rem+[0]*(5-len(rem)))
 inv=np.array(invcols,dtype=np.int64).T;assert np.array_equal(inv@matrix%P,np.eye(5,dtype=np.int64))
 recipe=dict(base='Existing p-integral coefficient ring restricted to admissible ss1..ss5; direction chart Y1=v1, base16 y=Yrest-vrest',localizations=['existing setup minors','Y1','cyclicDelta','quintic high-block determinant','Norm(t)'],etale_discriminant_not_needed_for_module_presentation=True,free_algebra='F=B16[s]/g5, monic g5 from existing2051 inverse',original_affine_q='reconstruct homogeneous45 q from P_Z(s),Y; divide by invertible t(s), then undo the exact affine frame',critical_element='e=rem_g(det(C*Jac_q(original134))(beta,q_affine(s)))',residuals=['rem_g(chi1031)',*['rem_g(original_homogeneous_F'+str(i)+')'for i in range(134)]],matrix_columns='All135 residuals times s^0..s^4, then e*s^0..s^4, reduced modulo g',exact_critical_module_presentation_shape=[5,680],annihilator='N_E=det(multiplication_by_e on free5 F); adjugate gives exact Bezout in(g,e)',support_criterion='Fitt0 = ideal of all5x5 minors of the5x680 matrix; this is EXACT support in this open',known_center_norm=det,cyclic_or_highblock_complement_excluded=False,full_ideal_unit_proved=False)
 (out/'finite-recipe.json').write_text(json.dumps(recipe,indent=2)+'\n');certificate=dict(prime=P,g5=g,primitive_roots=roots,original_q_points=qpoints,original_E_values=values,critical_remainder_e5=e,Bezout_e=eu,Bezout_g=gv,identity='Bezout_e*e5+Bezout_g*g5=1 modp as full polynomials',multiplication_by_e=matrix.tolist(),multiplication_inverse=inv.tolist(),norm=det,all675_original_residual_columns_zero_at_center=True)
 (out/'center-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 inputs=[Path(__file__),critical,cyclic,nativepath,quintic,quintic.parent/'native.txt',knownpath,plane,ROOT/'reports/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62/admissible-slice.json',ROOT/'code/s61_affine_fiber_macaulay.py'];result=dict(id='S61_AFFINE_CRITICAL_QUINTIC_PRESENTATION_20261103_v63',all_passed=True,new_exact_module_presentation=[5,680],direction_base_variables=16,admissible_beta_parameters=5,all_original134_retained=True,original_q_critical_E_used=True,center_prime=P,center_critical_norm=det,general_presentation_exact_finite_recipe=True,general_coefficients_expanded=False,generic_critical_unit_proved=False,global_H0_or_cyclic_complement_resolved=False,scope='Exact critical-support presentation on existing cyclic+quintic+t-unit open, with a new degree5 norm elimination. Retains all roots there; does not prove generic exhaustion.',input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in inputs},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('CRITICAL QUINTIC5x680 + ACTUAL ORIGINAL E NORM/BEZOUT PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
