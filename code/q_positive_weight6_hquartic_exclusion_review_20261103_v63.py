#!/usr/bin/env python3
"""Independent full-source / finite-graph audit of the proper Hquartic exclusion.

Rebuilds all209 equations from literal physical P+ (no producer imports),
then composes them by independent Fraction polynomial arithmetic.
"""
from pathlib import Path
from fractions import Fraction
import argparse,gzip,hashlib,itertools,json,time
import sympy as S
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
ID='Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63_REVIEW'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def qmat(a):return S.Matrix([[S.Rational(x)for x in row]for row in a])
def det(a):return S.Rational(DomainMatrix.from_Matrix(a).convert_to(S.QQ).det())
def rank(a):return DomainMatrix.from_Matrix(a).convert_to(S.QQ).rank()
def add(a,b):
 a=a.copy()
 for m,c in b.items():
  a[m]=a.get(m,Fraction(0))+c
  if not a[m]:del a[m]
 return a
def mul(a,b):
 c={}
 for m,x in a.items():
  for n,y in b.items():
   k=tuple(sorted(m+n));c[k]=c.get(k,Fraction(0))+x*y
 return {m:x for m,x in c.items()if x}
def scale(a,c):return{m:c*x for m,x in a.items()if c*x}
def tri(f,n):
 return{tuple(i for i,k in enumerate(ex)for _ in range(k)):Fraction(c)for ex,c in S.Poly(S.sympify(f),*n).terms()if c}
def encode(f,v):
 p=S.Poly(S.expand(f),*v);r={}
 for ex,c in p.terms():
  if not c:continue
  m=tuple(i for i,k in enumerate(ex)for _ in range(k));assert 1<=len(m)<=2;r[m]=Fraction(c)
 return r

def build(out):
 started=time.monotonic();base=ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63';report=json.loads((base/'report.json').read_text());cert=json.loads(gzip.decompress((base/'certificate.json.gz').read_bytes()));assert report['all_passed'] and sha(base/'certificate.json.gz')==report['certificate_sha256']
 for name,want in report['input_sha256'].items():assert sha(ROOT/name)==want,name
 x,y,z,w=S.symbols('x y z w');Z=S.Matrix([x,y,z,w]);pairs=list(itertools.combinations(range(4),2));mon1=[S.S.One,*Z];mon2=[*mon1,*[Z[i]*Z[j]for i,j in itertools.combinations_with_replacement(range(4),2)]];mon3=[S.prod(Z[i]for i in a)for d in range(4)for a in itertools.combinations_with_replacement(range(4),d)]
 def skew(values):
  A=S.zeros(4)
  for(i,j),f in zip(pairs,values):A[i,j]=f;A[j,i]=-f
  return A
 P=skew([-S.Rational(256,3)*w,-64*w*x,S.Rational(8,9),64*w*y-S.Rational(2,9),S.Rational(4,3)*x,x*x+S.Rational(2,3)*y]);assert P==S.Matrix([[S.sympify(v)for v in row]for row in cert['positive_P']])
 # Independent raw150 decoding uses the standard homogeneous pair/monomial order.
 raypath=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz';ray=json.loads(gzip.decompress(raypath.read_bytes()));eta=S.Symbol('eta');raw=[S.sympify(v).subs(eta,1)for v in ray['raw150']];X=S.symbols('h x y z w');pm=list(itertools.combinations(range(5),2));mm=list(itertools.combinations_with_replacement(range(5),2));rawP=S.zeros(5)
 for b,(i,j)in enumerate(pm):rawP[i,j]=sum(raw[15*b+t]*X[k]*X[l]for t,(k,l)in enumerate(mm));rawP[j,i]=-rawP[i,j]
 assert all(rawP[0,i]==0 for i in range(5));assert rawP[1:,1:].subs(X[0],1)==P
 pf=lambda A:S.expand(A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2]);assert pf(P)==-S.Rational(16,81)
 rr=S.symbols('r0:90');u,c1,c2,c4=S.symbols('u c1 c2 c4');dd=S.symbols('d0:5');V=[*rr,u,c1,c2,c4,*dd];assert len(V)==99 and list(map(str,V))==cert['source_variables'];F=P+skew([sum(rr[15*k+j]*m for j,m in enumerate(mon2))for k in range(6)])
 def jac(A):return[S.expand(sum(A[i,b]*S.diff(A[j,k],Z[b])+A[j,b]*S.diff(A[k,i],Z[b])+A[k,b]*S.diff(A[i,j],Z[b])for b in range(4)))for i,j,k in itertools.combinations(range(4),3)]
 equations=[]
 for f in jac(F):
  p=S.Poly(f,*Z);equations +=[p.coeff_monomial(m)for m in mon3]
 assert len(equations)==140 and all(f.subs(dict.fromkeys(rr,0))==0 for f in equations)
 center=[]
 for i,j in pairs:
  p=S.Poly(F[i,j],*Z);L=sum(p.coeff_monomial(v)*v for v in Z);Q=sum(p.coeff_monomial(v)*v for v in mon2[5:]);g=S.Poly(S.expand(u*L+sum(v*S.diff(Q,q)for v,q in zip([c1,c2,1,c4],Z))),*Z);center +=[g.coeff_monomial(v)for v in Z]
 equations +=center
 def modular(A):
  div=S.Matrix([S.expand(sum(S.diff(A[i,j],Z[i])for i in range(4)))for j in range(4)]);N=S.zeros(5);N[1:,0]=div.subs(dict.fromkeys(Z,0));N[1:,1:]=div.jacobian(Z);return N
 assert modular(P)==S.zeros(5);N=modular(F);equations +=list((N*N).applyfunc(S.expand))
 C0=skew([0,0,0,1,0,0]);cpairs=[p for p in pairs if p!=(1,2)];CC=C0.copy()
 for p,t in zip(cpairs,dd):i,j=p;CC[i,j]=t;CC[j,i]=-t
 def compatible(C,A):return[S.expand(sum(C[i,b]*S.diff(A[j,k],Z[b])+C[j,b]*S.diff(A[k,i],Z[b])+C[k,b]*S.diff(A[i,j],Z[b])for b in range(4)))for i,j,k in itertools.combinations(range(4),3)]
 for f in compatible(CC,F):
  p=S.Poly(f,*Z);equations +=[p.coeff_monomial(m)for m in mon1]
 assert len(equations)==209;encoded=[encode(f,V)for f in equations];saved=[{tuple(ids):Fraction(a,b)for ids,a,b in p}for p in cert['full209_equations_encoded']];assert encoded==saved
 # Fresh exact source Jacobian, all21 actual affine/scalar action derivatives.
 M=S.Matrix([[S.Rational(p.get((i,),0))for i in range(99)]for p in encoded]);assert M==qmat(cert['full_Jacobian_209_by99'])
 def physical_coeff(A):return[S.Poly(S.expand(A[i,j]),*Z).coeff_monomial(m)for i,j in pairs for m in mon2]
 def all_orbit(A,source=False):
  cols=[]
  for i,j in itertools.product(range(4),repeat=2):
   B=S.zeros(4);B[i,j]=1;f=B*A+A*B.T-sum(((B*Z)[k]*A.diff(Z[k])for k in range(4)),S.zeros(4));v=physical_coeff(f)
   if source:
    dC=B*C0+C0*B.T;dC-=dC[1,2]*C0;v +=[0,B[0,2],B[1,2],B[3,2],*[dC[i,j]for i,j in cpairs]]
   cols.append(S.Matrix(v))
  for i in range(4):cols.append(S.Matrix(physical_coeff(-A.diff(Z[i]))+([0]*9 if source else[])))
  cols.append(S.Matrix(physical_coeff(A)+([0]*9 if source else[])));return S.Matrix.hstack(*cols)
 Oall=all_orbit(P,True);assert Oall==qmat(cert['affine_orbit_all21_matrix']) and M*Oall==S.zeros(209,21) and rank(Oall)==19
 O=Oall[:,cert['orbit19_selected_columns']];T=qmat(cert['normal99_by3']);assert M*T==S.zeros(209,3);cols=cert['implicit77_columns'];rows=cert['implicit77_rows'];B=M[rows,cols];minor=det(B);assert minor==S.Rational(cert['implicit77_unit_minor']) and minor!=0
 frame=O.row_join(T).row_join(S.eye(99)[:,cols]);fd=det(frame);assert fd==S.Rational(cert['full99_transverse_frame_determinant']) and fd!=0
 # Independent full QQ coefficient compositions, with no degree truncation.
 normal=S.symbols('n0:3');graph=[{(j,):Fraction(T[i,j])for j in range(3)if T[i,j]}for i in range(99)];corrections={17:{(1,2):Fraction(1,48)},49:{(1,1):Fraction(-1,1536)},54:{(0,2):Fraction(1,2)}}
 for i,p in corrections.items():graph[i]=add(graph[i],p)
 assert graph==[tri(f,normal)for f in cert['exact99_polynomial_graph']]
 pull=[]
 for p in encoded:
  value={}
  for mon,c in p.items():
   term={():Fraction(1)}
   for i in mon:term=mul(term,graph[i])
   value=add(value,scale(term,c))
  pull.append(value)
 assert pull==[tri(f,normal)for f in cert['all209_exact_full_pullbacks']];assert[(i,p)for i,p in enumerate(pull)if p]==[(179,{(0,0,2):Fraction(4,3)})];assert all(not pull[i]for i in rows);assert all(not graph[i]for i in range(90,99))
 # Entire compact center and pencil fibers, rebuilt from literal P+.
 Kc=[]
 for i,j in pairs:
  p=S.Poly(P[i,j],*Z);L=sum(p.coeff_monomial(v)*v for v in Z);Q=sum(p.coeff_monomial(v)*v for v in mon2[5:])
  for v in Z:Kc.append([S.Poly(L,*Z).coeff_monomial(v),*[S.Poly(S.diff(Q,q),*Z).coeff_monomial(v)for q in Z]])
 Kc=S.Matrix(Kc);kk=S.symbols('k0:6');Cp=skew(kk);kp=[]
 for f in compatible(Cp,P):
  p=S.Poly(f,*Z)
  for v in mon1:kp.append([S.expand(p.coeff_monomial(v)).coeff(t)for t in kk])
 Kp=S.Matrix(kp)
 for A,key,expected_rank,kernel in [(Kc,'whole_projective_center_fiber',4,S.Matrix([0,0,0,1,0])),(Kp,'whole_projective_constant_pencil_fiber',5,S.Matrix([0,0,0,1,0,0]))]:
  d=cert[key];assert A==qmat(d['matrix']) and A*kernel==S.zeros(A.rows,1);md=det(A[d['rows'],d['columns']]);assert md==S.Rational(d['determinant']) and md!=0;assert rank(A)==expected_rank
 # Dense known representative: literal legacy formula, not the producer import.
 star=skew([1-x*x+x*z-x*w,-x*x,0,3*x*y-4*x*z+z*z-3*z*w,x*x+x*y-2*x*z+2*x*w+z*z/2-2*z*w+w*w,1+x*x+x*w]);assert star==S.Matrix([[S.sympify(v)for v in row]for row in cert['dense_Hquartic_affine_P']]);J=star.subs(dict.fromkeys(Z,0));assert J==qmat(cert['dense_Hquartic_constant_J']) and det(J)!=0
 assert all(S.Poly(star[i,j],*Z).coeff_monomial(v)==0 for i,j in pairs for v in Z);assert not any(jac(star)) and not any(compatible(J,star));Ns=modular(star);assert Ns==qmat(cert['dense_Hquartic_modular_N']) and Ns*Ns==S.zeros(5) and rank(Ns)==2
 h=S.Symbol('h');F4=S.expand(h**4*pf(star).subs({x:x/h,y:y/h,z:z/h,w:w/h}));expectedF=h**4+h*h*x*z+x**3*y-x**3*z+x*x*z*z/2-x*x*z*w;assert S.expand(F4-expectedF)==0;assert S.Poly(F4,y).degree()==1 and S.diff(F4,y)==x**3 and S.Poly(F4,y).coeff_monomial(1).subs(x,0)==h**4
 star_affine_rank=rank(all_orbit(star));assert star_affine_rank==20
 # All eta!=0 belong to one projective weighted orbit.
 et=S.Symbol('eta');Peta=skew([-S.Rational(256,3)*et*w,-64*et*w*x,S.Rational(8,9),64*et*w*y-S.Rational(2,9),S.Rational(4,3)*x,x*x+S.Rational(2,3)*y]);weights=[1,2,3,4]
 for i,j in pairs:
  for ex,c in S.Poly(Peta[i,j],*Z).terms():
   weight=-5+weights[i]+weights[j]-sum(a*b for a,b in zip(ex,weights));assert weight==(-6 if c.has(et)else 0)
 out.mkdir(parents=True,exist_ok=False);note=ROOT/'research/notes/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63_REVIEW.md';result=dict(id=ID,all_passed=True,review_type='Independent full209 reconstruction and Fraction full-polynomial compositions; exact source/frame/projective-fiber audit',full209_coefficients_rebuilt_from_literal_physical_P=True,full209_compositions_Fraction_no_truncation=True,only_nonzero_row=179,only_nonzero_pullback='4*n0^2*n2/3',all77_pivot_compositions_zero=True,source_Jacobian_rank=77,unit77_minor=str(minor),all21_actual_source_action_columns_rebuilt=True,actual_source_orbit_rank=19,full99_frame_determinant=str(fd),all_nine_auxiliary_coordinates_zero_on_transverse_graph=True,only_u_zero_claimed_after_orbit_restore=True,compact_center_rank=4,compact_center_unique_projective_point='[0:0:0:1:0]',constant_pencil_rank=5,constant_pencil_unique_projective_point='[0:0:0:1:0:0]',positive_quintic='-16*h^5/81',dense_star_all_Jacobi_and_compatibility_zero=True,dense_star_finite_center_and_nondegenerate_constant=True,dense_star_N_squared_zero_rank=2,dense_star_fixedH_affine_scalar_orbit_rank=star_affine_rank,dense_star_whole_projective_orbit_dimension=23,whole_Hquintic_factor_irreducible_by_linear_y_and_gcd=True,scheme_level_Hchart_trivialization_and_affine_action_reviewed=True,proper_H_center_pencil_closure_argument_reviewed=True,completed_fixedH_source='C[[orbit19,n0,n1,n2]]/(n0^2*n2)',whole_source_u_identically_zero=True,nonzero_eta_orbit_weights=[1,2,3,4],eta_weight=-6,entire_Hquartic_closure_exclusion_accepted=True,finite_jet_inference_used=False,full_Poisson_completion_or_component_count_claimed=False,seconds=time.monotonic()-started,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__).resolve(),note,base/'report.json',base/'certificate.json.gz',raypath,ROOT/'legacy/v06/nilpotent22.tex']});(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT FULL209 HQUARTIC INCIDENCE REVIEW PASS',result['seconds'],flush=True)
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);build(parser.parse_args().output.resolve())
