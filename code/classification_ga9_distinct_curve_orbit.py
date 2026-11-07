#!/usr/bin/env python3
"""One fullambient/orbit and singlocus check of the actual Ga9 curve member."""
from pathlib import Path
import hashlib,itertools,json,os,resource,subprocess,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,minor_certificate
from classification_three_hyperplanes_quadric_independent import ambient_matrices
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();source_path=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_CURVE_20261107.json';source=json.loads(source_path.read_text());xx=s.symbols('h x y z w');env=dict(zip(map(str,xx),xx));P=s.zeros(5)
 for ij,f in source['nonzero_point_brackets'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
 denominator=s.ilcm(1,*[s.denom(a)for f in P for a in s.Poly(f,*xx).coeffs()]);point=denominator*P;D,V,G,F=ambient_matrices(point,xx,[s.zeros(5)]);dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);rank=lambda A:len(dm(A).rref()[1]);rD,rV,rG=map(rank,[D,V,G]);assert D*G==s.zeros(350,G.cols)
 # The canonicalEuler representative allows exact projectiveStab equations.
 E=s.Matrix(xx);div=s.Matrix([sum(s.diff(P[i,j],xx[i])for i in range(5))for j in range(5)]);canonical=(P-(E*div.T-div*E.T)/5).applyfunc(s.expand)
 def vector(Q):return s.Matrix([s.Poly(s.expand(Q[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 columns=[]
 for i,j in itertools.product(range(5),repeat=2):
  A=s.zeros(5);A[i,j]=1;columns.append(vector((A*canonical+canonical*A.T-xx[j]*canonical.diff(xx[i])).applyfunc(s.expand)))
 stabilizer_matrix=s.Matrix.hstack(*columns,-vector(canonical));kernel=dm(stabilizer_matrix).nullspace(divide_last=True).to_Matrix();generators=[]
 for row in kernel.tolist():
  A=s.Matrix(5,5,row[:25]);A-=s.trace(A)*s.eye(5)/5;kappa=row[25]
  if A==s.zeros(5):continue
  scale=s.ilcm(1,*[s.denom(a)for a in list(A)+[kappa]]);A*=scale;kappa*=scale;gcd=s.igcd(*[abs(int(a))for a in list(A)+[kappa]if a]);A/=gcd;kappa/=gcd
  if any(A==s.Matrix(r['matrix'])for r in generators):continue
  char=s.factor(A.charpoly().as_expr());powers=[rank(A**i)for i in range(1,6)];generators.append(dict(matrix=[[int(a)for a in rr]for rr in A.tolist()],projective_eigenvalue=str(kappa),characteristic_polynomial=str(char),powers_ranks=powers))
 print('full ranks',rD,rV,rG,'stab',generators,flush=True)
 quintic=s.sympify(source['nonzero_point_actual_quintic'],locals=env);gcd=s.igcd(*[abs(int(a))for a in s.Poly(quintic,*xx).coeffs()]);quintic=s.expand(quintic/gcd)
 def encode(f):
  terms=[]
  for ex,a in s.Poly(f,*xx).terms():
   co=str(a.p)if a.q==1 else '('+str(a.p)+'/'+str(a.q)+')';terms.append('*'.join([co]+[str(v)+('^'+str(k)if k!=1 else'')for v,k in zip(xx,ex)if k]))
  return '('+'+'.join(terms)+')'if terms else'0'
 work=ROOT/'work/ga9-distinct-curve-orbit';work.mkdir(parents=True,exist_ok=True);input_path=work/'singular-locus.sing';output_path=work/'singular-locus.txt'
 lines=['ring r=0,(h,x,y,z,w),dp;','poly F='+encode(quintic)+';','ideal I=diff(F,h),diff(F,x),diff(F,y),diff(F,z),diff(F,w);','ideal S=std(I);','write("'+str(output_path)+'","DIM_CONE "+string(dim(S)));','write("'+str(output_path)+'","DEGREE_TOP "+string(mult(S)));','write("'+str(output_path)+'","COMPLETE");','print("COMPLETE_SINGULAR_LOCUS");','exit;'];input_path.write_text('\n'.join(lines)+'\n')
 native_started=time.monotonic();run=subprocess.run([str(ROOT/'work/nonsquare-v62/singular-runtime/prefix/usr/bin/Singular'),'-q',str(input_path)],capture_output=True,text=True,timeout=90)
 assert run.returncode==0 and not any(line.lstrip().startswith('?')for line in run.stdout.splitlines()),run.stdout+'\n'+run.stderr
 assert output_path.exists()and output_path.read_text().splitlines()[-1]=='COMPLETE';data={line.split()[0]:int(line.split()[1])for line in output_path.read_text().splitlines()if line!='COMPLETE'}
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(source_path.relative_to(ROOT)):hashlib.sha256(source_path.read_bytes()).hexdigest(),str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'code/classification_rank3_mixed_component_review.py':hashlib.sha256((ROOT/'code/classification_rank3_mixed_component_review.py').read_bytes()).hexdigest(),'code/classification_three_hyperplanes_quadric_independent.py':hashlib.sha256((ROOT/'code/classification_three_hyperplanes_quadric_independent.py').read_bytes()).hexdigest()},producer_bytes=Path(__file__).read_text(),
  integer_scale=str(denominator),full_D_rank=rD,Euler_rank=rV,actual_group_Euler_scalar_rank=rG,actual_projective_orbit_dimension=rG-25,full_projective_tangent_dimension=125-rD,full_D_times_actual_orbit_zero=True,
  full_D_minor=minor_certificate(D,rD),full_G_minor=minor_certificate(G,rG),full_projective_stabilizer_generators=generators,canonical_divergence_free_point={str(i)+','+str(j):str(canonical[i,j])for i,j in PAIRS if canonical[i,j]!=0},
  actual_quintic_primitive=str(quintic),quintic_singular_cone_dimension=data['DIM_CONE'],quintic_projective_singular_dimension=data['DIM_CONE']-1,singular_top_degree=data['DEGREE_TOP'],quintic_is_nonnormal_if_absolutely_irreducible=data['DIM_CONE']>=3,
  component_criterion='If actualorbitdimension23, its closure contains the Ga9sourceanchor through the exactcubiccurve. The completeGa9anchor truecutdimension1 and47actualcuts give fullprojectivelocaldimension<=23. Therefore the orbitclosure itself is a full23component. WholeGa9distinctsource orbit22 lies in it by actualcurve closure and itsaccepted denseorbit theorem. No componentdistinctness or count is inferred until allprior98 are compared.',
  native_input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),native_seconds=time.monotonic()-native_started,seconds=time.monotonic()-started,independent_review_required=True)
 out=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_CURVE_ORBIT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(full_D=rD,orbit=rG-25,singdim=data['DIM_CONE']-1,seconds=result['seconds'])))
if __name__=='__main__':main()
