#!/usr/bin/env python3
"""Whole nil31 image-inside n7=0 sources, and the half-branch smooth25 component.

No nonlinear coordinate change is used as projective equivalence. The actual
projective permutation, Euler gauge and full150/350 dimension squeeze are
checked directly. The second branch is retained without componenthood claim.
"""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import minor_certificate,jacobi,COLUMNS,PAIRS
ROOT=Path(__file__).resolve().parents[1]

def main():
 begun=time.monotonic();xx=s.symbols('h x1:5');h,a,b,c,d=xx;v=(a,b,c,d);E=s.Matrix(xx)
 AA=s.symbols('A0:6');u,vv,t,p,q=s.symbols('u v t p q');Ah=AA[0]*h*h+AA[1]*a*h+AA[2]*d*h+AA[3]*a*a+AA[4]*a*d+AA[5]*d*d
 original_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
 original=json.loads(original_path.read_text())['cases'][0]
 H,X,Y,Z,W=oldx=s.symbols('H X Y Z W');ns=s.symbols('n0:15');oldenv={str(n):n for n in ns};oldenv.update(dict(zip(['x1','x2','x3','x4'],[X/H,Y/H,Z/H,W/H])))
 raw=s.zeros(5)
 for (i,j),f in zip(itertools.combinations(range(1,5),2),original['full_six_brackets']):
  raw[i,j]=s.cancel(H*H*s.sympify(f,locals=oldenv));raw[j,i]=-raw[i,j]
 T=s.zeros(5)
 for i,j in enumerate([3,0,1,2,4]):T[i,j]=1
 records=[]
 def sparse(P):return {ij:dict(s.Poly(s.expand(P[ij]),*xx).terms())for ij in PAIRS if P[ij]!=0}
 def encode(P):return {str(i)+','+str(j):str(s.factor(P[i,j]))for i,j in PAIRS if P[i,j]!=0}
 def digest_matrix(M):return hashlib.sha256(json.dumps(M.tolist(),default=str,separators=(',',':')).encode()).hexdigest()
 for name,r,g,scalevar in [('half',-s.Rational(1,2),p*a*a,p),('third',-s.Rational(1,3),q*h*a,q)]:
  K=3*r+t+1;S=u*a+vv*h+(t-r)*d
  entries=[r*a*c,r*h*a,0,Ah+K*(h*b-c*c/2),g+c*S,h*S]
  P=s.zeros(5)
  for (i,j),f in zip(itertools.combinations(range(1,5),2),entries):P[i,j]=f;P[j,i]=-f
  # Universally verify all four strict Jacobi identities on h=1, and the
  # actual full homogeneous Pfaffian, not only its leading quadratic part.
  J=[s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(5)))for i,j,k in itertools.combinations(range(1,5),3)]
  assert J==[0]*4
  pf=s.factor(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
  assert s.expand(pf+r*h*a*g)==0
  bind={ns[0]:AA[3],ns[1]:AA[0],ns[2]:AA[2],ns[3]:AA[5],ns[4]:AA[1],ns[8]:AA[4],ns[12]:u,ns[13]:vv,ns[14]:t,
        ns[5]:p if name=='half'else 0,ns[6]:q if name=='third'else 0,ns[7]:0,ns[9]:0,ns[10]:r,ns[11]:0}
  old=raw.subs(bind,simultaneous=True);transport=(T*old*T.T).subs(dict(zip(oldx,(a,b,c,h,d))),simultaneous=True).applyfunc(s.expand)
  gauge=s.Matrix([0,0,-r*c,-r*h,0]);assert (transport-E*gauge.T+gauge*E.T-P).applyfunc(s.expand)==s.zeros(5)
  params=list(AA)+[u,vv,t,scalevar];point=dict(zip(AA,(1,2,3,5,7,11)));point.update({u:2,vv:3,t:5,scalevar:13})
  witness=12*P.subs(point);directions=[12*P.diff(z).subs(point)for z in params]
  assert jacobi(sparse(witness))==s.zeros(350,1)
  D,Vert,G,F=ambient_matrices(witness,xx,directions);assert D*F==s.zeros(350,F.cols)
  ranks={key:len(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])for key,M in [('D',D),('vertical',Vert),('orbit_vertical_scale',G),('whole_source',F)]}
  expected={'D':100 if name=='half'else 98,'vertical':24,'orbit_vertical_scale':47 if name=='half'else 48,'whole_source':50};assert ranks==expected,(name,ranks)
  certificates={key:minor_certificate(M,ranks[key])for key,M in [('D',D),('vertical',Vert),('orbit_vertical_scale',G),('whole_source',F)]}
  records.append(dict(name='nil31_image_n7zero_r_minus_'+name,free_parameters=list(map(str,params)),source_dimension=10,
   fixed_r=str(r),old_full_source_parameter_substitution={str(z):str(a)for z,a in bind.items()},
   whole_homogeneous_h_Casimir_brackets=encode(P),all_universal_strict_Jacobi_zero=True,
   actual_affine_homogeneous_pfaffian=str(pf),actual_projective_degeneracy=str(s.factor(h*pf)),
   coordinate_permutation_new_from_old=['old_z','old_h','old_x','old_y','old_w'],coordinate_matrix_determinant=str(T.det()),
   full_universal_PGL_and_Euler_identity=True,Euler_gauge_vector=list(map(str,gauge)),
   witness_parameter_values={str(z):str(a)for z,a in point.items()},witness_tensor_scale='12',witness_brackets=encode(witness),
   all_ten_actual_parameter_derivatives=[encode(R)for R in directions],all_D_times_actual_F_zero=True,ranks_QQ=ranks,
   full_matrix_sha256={key:digest_matrix(M)for key,M in [('D',D),('vertical',Vert),('orbit_vertical_scale',G),('whole_source',F)]},certificates=certificates,
   true_full_projective_source_dimension_lower=ranks['whole_source']-25,full_projective_tangent_upper=125-ranks['D'],
   smooth_component_squeeze=name=='half',whole_source_is_full_component=name=='half',projective_component_dimension=25 if name=='half'else None,
   boundary='The zero p/q rank<=2 boundary is included in the same irreducible polynomial source closure.',
   scope='Whole source is a smooth25 component by full dimension squeeze.'if name=='half'else'Full25-dimensional source lower bound only. Two genuine normal directions remain; no componenthood or known-family membership follows.'))
 here=Path(__file__);helper=ROOT/'code/classification_three_hyperplanes_quadric_independent.py';engine=ROOT/'code/classification_rank3_mixed_component_review.py'
 report=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in [original_path,here,helper,engine]},producer_bytes=here.read_text(),branches=records,
  new_half_component_distinctness='After componenthood and dimension25: all same-dimensional known components are T01,T39,T40,positive-weight NEW25,and CN31(2,2,1). The entire new rank4 source closure has divisor H^3 L^2, a closed multiplication-image condition. It has no SNC point; positiveNEW25 has smooth rank5 residual quadric; CN31 has no triple factor at its generic three-independent-hyperplane point. Hence the half component is distinct. Other known dimensions and generic-rank2 components are excluded by dimension/closed rank.',
  new_component_count=1,second_branch_unresolved=True,seconds=time.monotonic()-begun)
 dest=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,new_half_component25=True,branches=[dict(name=r['name'],ranks=r['ranks_QQ'])for r in records],seconds=report['seconds'])))
if __name__=='__main__':main()
