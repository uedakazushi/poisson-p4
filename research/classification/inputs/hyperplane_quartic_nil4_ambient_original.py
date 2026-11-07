import sympy as S
from pathlib import Path
import json,hashlib,time,math
arithmetic=Path('code/classification_two_hyperplanes_cubic_component_review.py')
raw=arithmetic.read_text();prefix=raw[:raw.index('report_path=')];exec(prefix)
input_path=Path('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json')
source=json.loads(input_path.read_text())['full_nilpotent_sources'][0]
n=S.symbols('n0:36');xx=S.symbols('x y z w');coords=xx
cases=[('delta_minus_third',{n[8]:-S.Rational(1,4),n[12]:-S.Rational(1,3),n[1]:1}),('delta_minus_half',{n[8]:-S.Rational(1,4),n[12]:-S.Rational(1,2),n[10]:1,n[6]:3})]
records=[]
for label,values in cases:
 substitution={v:values.get(v,0)for v in n}
 entries=[S.Poly(S.sympify(t).subs(substitution),*xx)for t in source['brackets']]
 integer_scale=S.ilcm(*[c.q for p in entries for c in p.coeffs()]);P={}
 for ij,p in zip(combinations(range(1,5),2),entries):
  data={}
  for exponent,c in p.terms():
   full=(2-sum(exponent),)+exponent+(0,)*(N-5);data[full]=F(int(c*integer_scale))
  P[ij]=data
 P=tensor(P);assert all(not p for p in strict_jacobi(P))
 D=[[0]*150 for _ in row_basis]
 for col,(ij,m) in enumerate(column_basis):
  variation=tensor({ij:{m+(0,)*(N-5):F(1)}})
  mixed=tensor({kl:add(P.get(kl,{}),mul(term(1,19),variation.get(kl,{})))for kl in pairs})
  for inds,p in projective_jacobi(mixed).items():
   for exponent,a in p.items():
    if exponent[19]==1:
     assert a.denominator==1 and not any(exponent[5:19]);D[row_index[inds,exponent[:5]]][col]=int(a)
 vectors=[]
 for i,j in product(range(5),repeat=2):
  A=tensor({(k,l):add(term(int(l==i),k,j),term(-int(k==i),l,j))for k,l in pairs});vectors.append(bivector_vector(A))
 for i,j in product(range(5),repeat=2):
  A=tensor({(k,l):add(scale(mul(term(1,j),der(P.get((k,l),{}),i)),-1),scale(P.get((j,l),{}),int(k==i)),scale(P.get((k,j),{}),int(l==i)))for k,l in pairs});vectors.append(bivector_vector(A))
 vectors.append(bivector_vector(P));Fmatrix=list(map(list,zip(*vectors)))
 for row in D:
  nz=[(j,a)for j,a in enumerate(row)if a];assert all(sum(a*v[j]for j,a in nz)==0 for v in vectors)
 dc=modular_minor(D);oc=modular_minor(Fmatrix)
 record=dict(source=label,point_substitutions={str(k):str(v)for k,v in values.items()},integer_scale=int(integer_scale),rank_lower_D=dc,rank_lower_orbit_scale=oc,tangent_upper=125-dc['rank_mod_prime'],orbit_lower=oc['rank_mod_prime']-25,D_times_orbit_exact_zero=True,Dshape=[350,150],D_sha256=matrix_hash(D),orbit_matrix_sha256=matrix_hash(Fmatrix))
 records.append(record);print(json.dumps({k:record[k]for k in['source','tangent_upper','orbit_lower']}),flush=True)
out=dict(task='E_reduced_nonnormal',source_bytes=Path(__file__).read_text(),input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in[arithmetic,input_path]},records=records,seconds=time.monotonic()-started,scope='Exact fullambient D/orbit diagnostics at two whole-source witnesses; no fullsource membership inference from a tangent sample.')
Path('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_AMBIENT_DIAGNOSTIC_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
