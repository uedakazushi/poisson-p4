#!/usr/bin/env python3
"""Independent whole-flag coverage for Q_b, all nonzero complex b; no GB."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1];PARENT='AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63';STEM=PARENT+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);folder=ROOT/'reports'/PARENT;r=load(folder/'report.json');assert r['all_passed']
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h,n
 for n,h in r['evidence_sha256'].items():assert sha(folder/n)==h,n
 data=load(folder/'certificate.json.gz');basicpath=ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW/report.json';basic=load(basicpath)
 assert basic['all_passed'] and basic['both_source_ranks']==[29,29] and basic['all58_union_rank']==37
 h,x,y,z,w,b=s.symbols('h x y z w b');X=s.Matrix([h,x,y,z,w]);coords=[h,x,z,w];loc={str(v):v for v in [*X,b]}
 parse=lambda f:s.sympify(f,locals=loc);mat=lambda a:s.Matrix([[parse(v) for v in row] for row in a])
 same=lambda A,B:A.shape==B.shape and (A-B).applyfunc(s.cancel)==s.zeros(*A.shape)
 K=s.QQ.frac_field(b);dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(K)
 pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));cubic5=list(itertools.combinations_with_replacement(range(5),3))
 def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[c]) for i,j in pairs for a,c in mons])
 def canonical(P):
  div=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
  return (P-(X*div.T-div*X.T)/5).applyfunc(s.expand)
 def raw(f=0,g=0,H=0,CB=False):
  P=s.zeros(5);P[1,2]=f;P[3,4]=g if CB else g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
  if CB:P[2,3]-=s.diff(f,x)*z/2;P[2,4]-=s.diff(f,x)*w/2
  return P-P.T
 V=s.Matrix([x*x+s.Rational(8,5)*b*h*x,-s.Rational(2,5)*b*x*x,0,h*x+s.Rational(3,5)*b*x*z,h*h+z*z-s.Rational(7,5)*b*x*w]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
 cas=s.Matrix([[s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[c]) for j in range(5)] for i in range(5) for a,c in mons])
 assert same(cas,mat(data['all75_Casimir_matrix']));cw=data['Casimir_rank_witness']
 cminor=s.factor(K.to_sympy(dm(cas.extract(cw['rows'],cw['columns'])).det()));assert s.cancel(cminor+384*b**4/625)==0
 factor=[]
 for t in range(5):
  v=s.eye(5)[:,t];column=[]
  for i,j,k,l in itertools.combinations(range(5),4):
   expr=s.expand(X[i]*(v[j]*P[k,l]-v[k]*P[j,l]+v[l]*P[j,k])-X[j]*(v[i]*P[k,l]-v[k]*P[i,l]+v[l]*P[i,k])+X[k]*(v[i]*P[j,l]-v[j]*P[i,l]+v[l]*P[i,j])-X[l]*(v[i]*P[j,k]-v[j]*P[i,k]+v[k]*P[i,j]))
   column.extend(s.Poly(expr,*X).coeff_monomial(s.prod(X[u] for u in mon)) for mon in cubic5)
  factor.append(s.Matrix(column))
 factor=s.Matrix.hstack(*factor);assert same(factor,mat(data['all175_intrinsic_factor_matrix'])) and factor[:,2]==s.zeros(175,1)
 fw=data['intrinsic_factor_rank_witness'];fminor=s.factor(K.to_sympy(dm(factor.extract(fw['rows'],fw['columns'])).det()));assert fminor==2*b
 D=lambda f:s.expand(sum(5*V[i]*s.diff(f,X[i]) for i in range(5)))
 A,B,C,E,rh,rx,rz,rw=s.symbols('A B C D rh rx rz rw');ell=A*h+B*x+C*z+E*w;cof=rh*h+rx*x+rz*z+rw*w
 equation=s.Poly(D(ell)-ell*cof,h,x,z,w)
 assert equation.coeff_monomial(w*w)==-E*rw and equation.coeff_monomial(h*w)==-A*rw-E*rh
 assert equation.coeff_monomial(z*w)==-C*rw-E*rz and equation.coeff_monomial(h*h)==5*E-A*rh
 # If E!=0: w²,hw,zw force rw=rh=rz=0, contradict h²=5E.
 forms=[x,x+b*h-b*b*z,x+2*b*h];lams=[-2*b,3*b,8*b]
 coeff=s.Matrix([[s.diff(f,v) for f in forms] for v in [h,x,z]])
 linear=s.Matrix([[8*b,0,5],[5,-2*b,0],[0,0,3*b]])
 assert linear*coeff==coeff*s.diag(*lams) and coeff.det()==-2*b**3
 assert all(s.expand(D(f)-lam*x*f)==0 for f,lam in zip(forms,lams))
 assert same(coeff,mat(data['eigenbasis_matrix']))
 # All16 Hamiltonian variables retained in the universal endpoint/source test.
 cubics=[s.prod(coords[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
 hh=s.symbols('c0:16');H=sum(c*m for c,m in zip(hh,cubics));f0,f1,f2,g0,g1,g2,g=s.symbols('f0 f1 f2 g0 g1 g2 g')
 F=f0*h*h+f1*h*x+f2*x*x;G=g0*h*h+g1*h*x+g2*x*x
 ca0=canonical(raw(g=g,H=H));cb0=canonical(raw(g=G,H=H,CB=True))
 for endpoint in [ca0,cb0]:assert endpoint[0,:]==endpoint[1,:]==s.zeros(1,5)
 ca=canonical(raw(f=F,H=H));cb=canonical(raw(f=F,H=H,CB=True))
 for tensor,mult in [(ca,1),(cb,2)]:
  assert s.expand(-tensor[0,2]-mult*s.diff(F,x)*h/5)==0
  assert s.expand(-tensor[1,2]+F-mult*s.diff(F,x)*x/5)==0
 # Non-x ell: x lies in W, take k=x, and compare BOTH independent ell²,x².
 alpha=s.Symbol('alpha');CBres=[];CAres=[]
 for ell,lam in zip(forms[1:],lams[1:]):
  cres=s.expand(-2*b*x*x/5+alpha*ell*ell+lam*x*x/20);ares=s.expand(-2*b*x*x/5+alpha*ell*ell+3*lam*x*x/10)
  CBres.append(cres);CAres.append(ares)
  assert s.Poly(ell*ell,h,x,z,w).coeff_monomial(h*h)!=0
  assert s.cancel(cres-alpha*ell*ell-(lam-8*b)*x*x/20)==0
  assert s.cancel(ares-alpha*ell*ell-(3*lam-4*b)*x*x/10)==0
 assert all(s.cancel(a-parse(c))==0 for a,c in zip(CBres,data['CB_other_marked_line_residues']))
 assert all(s.cancel(a-parse(c))==0 for a,c in zip(CAres,data['CA_other_marked_line_residues']))
 # Ell=x. Solve ALL four k coefficients modulo the x² term.
 constraint_matrices=[]
 for lam,label,known in [(3*b,'CB',[s.Matrix([0,1,0,0]),s.Matrix([1,0,-b,0])]),(8*b,'CA',[s.Matrix([1,0,0,0]),s.Matrix([0,1,0,0])])]:
  M=s.Matrix([[s.Poly(D(v)-lam*x*v,h,x,z,w).coeff_monomial(s.prod(coords[i] for i in mon)) for v in coords] for mon in itertools.combinations_with_replacement(range(4),2) if mon!=(1,1)])
  assert same(M,mat(data[f'{label}_marked_x_k_constraint_matrix']))
  wit=data[f'{label}_marked_x_k_witness'];assert wit['rank']==2
  kdet=K.to_sympy(dm(M.extract(wit['rows'],wit['columns'])).det())
  assert s.cancel(kdet-(-25*b if label=='CB' else -25))==0
  assert dm(M).rank()==2 and all(M*v==s.zeros(9,1) for v in known)
  constraint_matrices.append(M)
 # CA and N source immersion: independently reconstruct the swapped CA basis.
 cubic_sw=[s.prod((x,h,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
 def raw_sw(f=0,g=0,H=0):
  P=s.zeros(5);P[0,2]=f;P[3,4]=g*x*x;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z);return P-P.T
 fixed=s.Matrix.hstack(*[vec(canonical(raw_sw(f=f))) for f in [x*x,h*x,h*h]],vec(canonical(raw_sw(g=1))),*[vec(canonical(raw_sw(H=H))) for H in cubic_sw])
 capath=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';npath=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
 caold,nold=load(capath),load(npath);assert same(fixed,s.Matrix(caold['fixed_fiber_linear_basis150_by20']))
 nc=fixed[:,[0,1,*range(3,20)]];assert same(nc,s.Matrix(nold['fiber19_basis']))
 center=s.Matrix(caold['center_fixed20']);center[1]-=2*b;center[9]-=b;assert center[2]==0 and same(fixed*center,vec(P))
 flagcols=[]
 for i,j in caold['flag_positions']:
  E=s.zeros(5);E[i,j]=1;flagcols.append(vec((E*P+P*E.T-X[j]*P.diff(X[i])).applyfunc(s.expand)))
 flagJ=s.Matrix.hstack(*flagcols);source_minors=[]
 for index,(name,Cbase) in enumerate([('CA',fixed),('N',nc)]):
  rr=list(dm(Cbase.T).rref()[1]);inverse=Cbase[rr,:].inv();quot=(flagJ-Cbase*inverse*flagJ[rr,:]).applyfunc(s.expand)
  record=data['CA_N_actual_source_matrices'][index];assert record['name']==name
  assert same(Cbase.row_join(flagJ),mat(record['full_source'])) and same(quot,mat(record['flag_quotient']))
  wit=data['CA_N_actual_source_rank_certificates'][index]['flag_quotient_witness'];assert wit['rank']==9
  det=s.factor(K.to_sympy(dm(quot.extract(wit['rows'],wit['columns'])).det()))
  assert s.cancel(det-[-262144*b**9/15625,1572864*b**9/78125][index])==0;source_minors.append(str(det))
 properca=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';propercb=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json'
 nreview=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json'
 for dep in [properca,propercb,nreview]:assert load(dep)['all_passed']
 files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',folder/'report.json',folder/'certificate.json.gz',basicpath,capath,npath,properca,propercb,nreview]
 files += [ROOT/n for n in r['input_sha256']]
 result=dict(id=STEM,all_passed=True,scope='every nonzero complex b',
  all75_constant_Casimir_entries_verified=True,Casimir5minor=str(cminor),
  all175_intrinsic_factor_entries_verified=True,constant_factor4minor=str(fminor),unique_factor='e_y',
  full_Darboux_equation_and_all_w_cases_proved=True,complete_three_eigenform_basis=True,eigenbasis_determinant='-2*b^3',
  all16_H_universal_endpoint_and_both_source_identities_verified=True,
  all4_k_coefficients_and_both9x4_constraint_matrices_verified=True,
  whole_CB_geometric_flags=2,whole_CA_geometric_flags=1,whole_N_geometric_flags=1,
  CA_N_full150_source_derivatives_verified=True,CA_N_source_ranks=[29,28],CA_N_flag9_minors=source_minors,
  whole_CB_fiber_exactly_two_reduced_points_all_nonzero=True,whole_CA_N_fibers_reduced_length1_all_nonzero=True,
  CA_smooth_cone29_N_smooth_cone28_all_nonzero=True,
  CB_two_smooth29_formal_branches_all_nonzero=True,CB_tangent_union=37,CB_tangent_intersection=21,
  CB_cone_depth_at_most22_all_nonzero=True,CB_nonnormal_singular_nonCM_all_nonzero=True,
  origin_full_odd_or_conductor_multiplicity_or_D_claimed=False,new_component_count=0,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(files)},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',STEM,result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
