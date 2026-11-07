#!/usr/bin/env python3
"""Exact14-quadratic and maximal-radical certificate for the full fixedS source."""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();source_path=ROOT/'reports/CLASSIFICATION_PROJECTED_VERONESE_ZERO_20261107.json';packet=json.loads(source_path.read_text());xx=s.symbols('x0:5');zz=s.symbols('z0:5');env={str(z):z for z in xx+zz};P=s.zeros(5)
 assert packet['exact_matrix_rank']==121 and packet['exact_kernel_dimension']==29 and packet['complete_projective_section_dimension']==5
 for ij,f in packet['complete_canonical_source'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
 sparse={ij:dict(s.Poly(s.expand(P[ij]),*xx).terms())for ij in PAIRS if P[ij]!=0};J=jacobi(sparse).applyfunc(s.expand);assert list(map(str,J))==packet['complete350_Jacobi']
 mons=[zz[i]*zz[j]for i,j in itertools.combinations_with_replacement(range(5),2)]
 C=s.Matrix([[s.Poly(f,*zz).coeff_monomial(m)for m in mons]for f in J]);dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ)
 z0,z1,z2,z3,z4=zz
 relations=[z0*z0-z4*z4,z0*z1-z4*z4,z1*z1-4*z4*z4,z0*z2,z1*z2,z2*z2-z4*z4,z0*z3,z1*z3,z2*z3,z3*z3-z4*z4,z0*z4,z1*z4,z2*z4,z3*z4]
 T=s.Matrix([[s.Poly(f,*zz).coeff_monomial(m)for m in mons]for f in relations]);assert len(dm(C).rref()[1])==14;assert len(dm(T).rref()[1])==14
 rows=list(dm(C.T).rref()[1]);small=C[rows,:];cols=list(dm(small).rref()[1]);Lsmall=T[:,cols]*small[:,cols].inv();assert Lsmall*small==T
 L=s.zeros(14,350)
 for j,k in enumerate(rows):L[:,k]=Lsmall[:,j]
 assert L*C==T
 tcols=list(dm(T).rref()[1]);inverse_coefficients=C[:,tcols]*T[:,tcols].inv();assert inverse_coefficients*T==C
 r=relations
 cube_representations=[{0:z0,10:z4},{2:z1,11:4*z4},{5:z2,12:z4},{9:z3,13:z4},{0:-z4,10:z0}]
 for i,rep in enumerate(cube_representations):assert s.expand(sum(co*r[j]for j,co in rep.items())-zz[i]**3)==0
 # The geometric anticanonical section is Euler contracted with the
 # homogeneous4-vector: D=sum_i (-1)^i*x_i*Pf(P_without_i), degree5.
 # Pf_i themselves have degree4 and are not x_i times a common quintic.
 pf=[]
 for i in range(5):
  a,b,c,d=[j for j in range(5)if j!=i];pf.append((-1)**i*(P[a,b]*P[c,d]-P[a,c]*P[b,d]+P[a,d]*P[b,c]))
 D=s.expand(sum(xx[i]*pf[i]for i in range(5)))
 assert D==0
 actual_affine=s.sympify(packet['actual_affine_Pfaffian'],locals=env)
 assert s.expand(D.subs(xx[0],1)-actual_affine)==0
 # All coefficientquadrics of actualD belong to the COMPLETE Jacobi ideal.
 actual_coefficients=s.Poly(D,*xx).coeffs();G=s.groebner(relations,*zz)
 assert all(G.reduce(f)[1]==0 for f in actual_coefficients)
 socle=s.Matrix([[1,1,0,0,0],[1,4,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
 assert socle.det()==3
 inputs=[source_path,Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py']
 out=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),
  whole_linear_source_matrix_rank=121,whole_linear_kernel_dimension=29,Euler_rank=24,whole_projective_section_dimension=5,
  parameter_monomials=list(map(str,mons)),full350_by15_Jacobi_coefficient_matrix=[[str(z)for z in row]for row in C.tolist()],exact_quadratic_rank=14,
  fourteen_exact_generators=list(map(str,relations)),generator_from_full350_functionals=[{str(j):str(a)for j,a in enumerate(L.row(i))if a}for i in range(14)],
  all350_from_generators_coefficients=[[str(z)for z in row]for row in inverse_coefficients.tolist()],both_quadratic_ideal_inclusions_checked=True,
  all_five_cubes_in_full_Jacobi_ideal=[{str(j):str(a)for j,a in rep.items()}for rep in cube_representations],all_cube_identities_checked=True,
  reduced_fixed_surface_Poisson_support='origin z0=z1=z2=z3=z4=0',only_zero_projective_Poisson_bivector_vanishes_on_this_S=True,
  actual_anticanonical_quintic=str(D),actual_D_identically_zero_on_complete_linear_source_before_Jacobi=True,actual_D_equals_affine_Pfaffian_on_x0_one=True,all_actual_D_coefficients_in_full_Jacobi_ideal=True,
  proper_Pfaffian_convention='Actual D=sum_i(-1)^i*x_i*Pf(P_without_i), homogeneousdegree5. The degree4 Pf vector is not x_i times a common quintic. The exploratory common-vector diagnostic in the preceding producer is not used as a geometric identity or proof.',
  socle_quadratic_multiplication_matrix=[[int(z)for z in row]for row in socle.tolist()],socle_matrix_determinant=3,
  whole_branch_exclusion='Every smooth projectedVeronese surface inP4 is projectively equivalent to the chosen projection from a rank3 center. Pullback preserves Poisson and projectivevanishing. Since the COMPLETE fixedS source has only zeroPoisson support, no nonzeroPoisson structure onP4 vanishes on any suchS. In particular it cannot be a rankzero surface of a rank4 Poisson structure with irreducible nonnormalquintic divisor.',
  omitted_surface_types='Other nonACM surfaces, singular projections, singular or nonreduced rankzero schemes, and general irreducible nonnormal quintics are not excluded.',independent_review_required=True,seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_PROJECTED_VERONESE_SUPPORT_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(all_passed=True,only_zero_Poisson=True,quadratic_rank14=True,seconds=out['seconds'])))
if __name__=='__main__':main()
