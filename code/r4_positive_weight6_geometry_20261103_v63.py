#!/usr/bin/env python3
"""Exact geometry of the UNCONSTRAINED three-parameter99-graph Poisson family.

No modular N² condition is imposed. Component identity/count belongs to the
independent full-family/ambient-dimension packets, not this geometry producer.
"""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector,action
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X;Y=s.Matrix([x,y,z,w]);PAIRS=list(itertools.combinations(range(4),2))
def mj(M):return [[str(v)for v in row]for row in M.tolist()]
def rankcert(M):
 cc=list(M.rref()[1]);V=M[:,cc];rr=list(V.T.rref()[1]);B=V[rr,:];co=B.inv()*M[rr,:];assert V*co==M
 return dict(rank=len(cc),rows=rr,columns=cc,determinant=str(s.factor(B.det())),all_column_coefficients=mj(co),matrix=mj(M))
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();source=ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/certificate.json.gz';old=json.loads(gzip.decompress(source.read_bytes()));graph=list(map(s.sympify,old['exact99_polynomial_graph']));mons=list(map(s.sympify,old['coefficient_monomials']));P=s.Matrix([[s.sympify(v)for v in row]for row in old['positive_P']]);n0,n1,n2=s.symbols('n0 n1 n2')
 for k,(i,j)in enumerate(PAIRS):
  f=sum(c*m for c,m in zip(graph[15*k:15*k+15],mons));P[i,j]+=f;P[j,i]-=f
 P=P.applyfunc(s.expand);assert all(graph[i]==0 for i in range(90,99))
 jac=[]
 for i,j,k in itertools.combinations(range(4),3):
  f=s.expand(sum(P[i,b]*s.diff(P[j,k],Y[b])+P[j,b]*s.diff(P[k,i],Y[b])+P[k,b]*s.diff(P[i,j],Y[b])for b in range(4)));assert f==0;jac.append(str(f))
 pf=s.factor(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]);poly=s.Poly(pf,*Y);assert poly.total_degree()==2
 homog=lambda f:s.expand(sum(c*h**(2-sum(e))*s.prod(v**k for v,k in zip(Y,e))for e,c in s.Poly(f,*Y).terms()))
 Q=homog(pf);QM=s.hessian(Q,X)/2;qd=s.factor(QM.det());assert qd==s.Rational(1024,531441)*n0**3*n1*(32*n0+n1*n2)
 hminor=s.factor(QM.extract([1,2,4],[1,2,4]).det());assert hminor==-s.Rational(4,729)*n0*n1*(32*n0+n1*n2)
 vertex={h:0,x:0,y:0,z:1,w:0};gradient=[s.expand(s.diff(Q,v).subs(vertex))for v in X];assert gradient==[-s.Rational(32,27)*n0,0,0,0,0]
 raw=s.zeros(5)
 for i,j in PAIRS:raw[i+1,j+1]=homog(P[i,j]);raw[j+1,i+1]=-raw[i+1,j+1]
 actual=s.expand(h*(raw[1,2]*raw[3,4]-raw[1,3]*raw[2,4]+raw[1,4]*raw[2,3]));assert actual==s.expand(h**3*Q);Pc=canonical(raw)
 div=s.Matrix([s.expand(sum(s.diff(P[i,j],Y[i])for i in range(4)))for j in range(4)]);N=s.zeros(5);N[1:,0]=div.subs(dict.fromkeys(Y,0));N[1:,1:]=div.jacobian(Y);assert N**3==s.zeros(5)and N**2!=s.zeros(5);assert s.expand((N*N)[3,0])==4*n0*n0*n2/3
 assert s.factor(N.extract([2,3],[0,1]).det())==4*n0*n0/3;assert N.charpoly().as_expr()==s.Symbol('lambda')**5
 assert (N.T*QM+QM*N).applyfunc(s.expand)==s.zeros(5)
 assert (P*s.Matrix([s.diff(pf,v)for v in Y])+pf*div).applyfunc(s.expand)==s.zeros(4,1)
 # Exact finite modular Ga action.
 t=s.Symbol('t');G=s.eye(5)+t*N+t*t*N*N/2;Gi=s.eye(5)-t*N+t*t*N*N/2;flow=(G*raw.subs(dict(zip(X,Gi*s.Matrix(X))),simultaneous=True)*G.T-raw).applyfunc(s.expand);assert flow==s.zeros(5)
 # Compatible constant pencil and its translation realization.
 C=s.zeros(4);C[1,2]=1;C[2,1]=-1
 for i,j,k in itertools.combinations(range(4),3):assert s.expand(sum(C[i,b]*s.diff(P[j,k],Y[b])+C[j,b]*s.diff(P[k,i],Y[b])+C[k,b]*s.diff(P[i,j],Y[b])for b in range(4)))==0
 lam=s.Symbol('lambda');shift=3*lam/(4*n0);translated=P.subs(z,z-shift).applyfunc(s.cancel);assert (translated-P-lam*C).applyfunc(s.cancel)==s.zeros(4)
 # Explicit connected parameter torus and a discrete sign equivalence.
 a=s.Symbol('a',nonzero=True);parameter_transforms=[]
 for eps in [1,-1]:
  offset=(a**3-eps)/(6*n0);D=s.diag(a,a*a,a**3,eps*a);back={x:x/a,y:y/(a*a),z:(z-offset)/(a**3),w:w/(eps*a)};changed=(s.Rational(eps)*a**-2*D*P.subs(back,simultaneous=True)*D.T).applyfunc(s.cancel);expected=P.subs({n0:eps*n0,n1:eps*a*n1,n2:n2/a},simultaneous=True);assert (changed-expected).applyfunc(s.cancel)==s.zeros(4);parameter_transforms.append(dict(sign=eps,coordinate_linear=mj(D),z_translation=str(offset),tensor_scalar=str(eps*a**-2),parameter_image=[str(eps*n0),str(eps*a*n1),str(n2/a)]))
 sample={n0:1,n1:2,n2:3}
 # Generic uniqueness of the projective constant compatible pencil and compact center.
 kk=s.symbols('c0:6');CC=s.zeros(4)
 for(i,j),v in zip(PAIRS,kk):CC[i,j]=v;CC[j,i]=-v
 ce=[]
 for i,j,k in itertools.combinations(range(4),3):
  poly=s.Poly(s.expand(sum(CC[i,b]*s.diff(P[j,k],Y[b])+CC[j,b]*s.diff(P[k,i],Y[b])+CC[k,b]*s.diff(P[i,j],Y[b])for b in range(4))),*Y);ce.extend(poly.coeff_monomial(v)for v in [1,*Y])
 pencil=s.linear_eq_to_matrix(ce,kk)[0];pencilcert=rankcert(pencil.subs(sample));assert pencilcert['rank']==5
 uu,*cvv=s.symbols('u vx vy vz vw');center_equations=[]
 for i,j in PAIRS:
  poly=s.Poly(P[i,j],*Y);linear=sum(poly.coeff_monomial(v)*v for v in Y);quad=sum(co*s.prod(v**k for v,k in zip(Y,exp))for exp,co in poly.terms()if sum(exp)==2);expression=uu*linear+sum(c*s.diff(quad,v)for c,v in zip(cvv,Y));pp=s.Poly(s.expand(expression),*Y);center_equations.extend(pp.coeff_monomial(v)for v in Y)
 centermatrix=s.linear_eq_to_matrix(center_equations,[uu,*cvv])[0];centercert=rankcert(centermatrix.subs(sample));assert centercert['rank']==4
 assert centermatrix*s.Matrix([0,0,0,1,0])==s.zeros(24,1)
 Psample=Pc.subs(sample);O=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],Psample))for i,j in itertools.product(range(5),repeat=2)]);oc=rankcert(O);scalar=rankcert(O.row_join(vector(Psample)));whole=rankcert(O.row_join(vector(Psample)).row_join(s.Matrix.hstack(*[vector(Pc.diff(v).subs(sample))for v in [n0,n1,n2]])));assert oc['rank']==23 and scalar['rank']==24 and whole['rank']==26
 files=[Path(__file__),source,ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/report.json',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py'];report=dict(id='R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63',provisional_component_id='R4_POSITIVE_WEIGHT6_25',all_passed=True,unconstrained_n0_n1_n2=True,no_modular_squarezero_equation_imposed=True,affine_P=mj(P),all4_affine_Jacobi_polynomials=jac,affine_pfaffian=str(pf),actual_projective_quintic='h³*Q2',homogeneous_quadric=str(Q),quadric_matrix=mj(QM),quadric_determinant=str(qd),H_restriction_rank3_minor=str(hminor),H_tangency_vertex='[h:x:y:z:w]=[0:0:0:1:0]',quadric_gradient_at_vertex=list(map(str,gradient)),generic_geometry_open='n0*n1*n2*(32*n0+n1*n2) !=0',modular_vector_field=list(map(str,div)),homogeneous_modular_N=mj(N),N_square=mj(N*N),modular_Jordan_generic=[3,1,1],modular_ranks_generic=[2,1,0],modular_Ga_finite_action=mj(G),modular_Ga_preserves_raw_P_identically=True,modular_N_is_Q_skew=True,modular_log_Hamiltonian_identity='P*gradient(Pfaffian)=-Pfaffian*div(P)',compact_center_generic='unique [u:vx:vy:vz:vw]=[0:0:0:1:0] at infinity',compact_center_matrix=mj(centermatrix),compact_center_rank_sample=centercert,constant_pencil_generic_unique=True,constant_pencil_matrix=mj(pencil),constant_pencil_rank_sample=pencilcert,constant_pencil='C=partial_y wedge partial_z',pencil_generated_by_z_translations=True,parameter_equivalences=parameter_transforms,two_parameter_normal_form='n1=1, parameters u=n0 and v=n1*n2; simultaneous (u,v)->(-u,-v) is an explicit equivalence; full finite identifications not classified',generic_projective_orbit_dimension=23,connected_projective_stabilizer_generic='Ga generated by the modular N',generic_orbit_proof='N gives stabilizer>=1 on entire open; exact sample23rank gives lower23, hence genericexact23',saturation_projective_dimension=25,saturation_dimension_proof='two-parameter normalform gives upper23+2=25; independent exact full family rank26 cone gives lower25',orbit_sample=oc,orbit_plus_scalar_sample=scalar,family_source_sample=whole,new_component_or_lower58_claimed_here=False,scope='Geometry and dimension of the explicit3family saturation; full component recognition/count requires independent fullambient upper/exclusion packets.',input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('NEW POSITIVE3 FAMILY GEOMETRY h3Q + J3 MODULAR + TWO MODULI PASS',report['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
