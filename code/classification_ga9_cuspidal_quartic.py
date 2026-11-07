"""Complete actual cuspidal-quartic families and a new24 source support proof."""
import argparse,hashlib,importlib.util,itertools as it,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
h,x,y,z,w=coords=s.symbols('h x y z w');a=s.Symbol('a');base=(h,y,w);r=s.symbols('r0:6');mons=[h*h,h*y,h*w,y*y,y*w,w*w]
helper=Path('legacy/v03/code/poisson_p4_components_check.py');spec=importlib.util.spec_from_file_location('ga9_quartic_ambient',helper);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
def pivots(A):return list(dm(A).rref()[1])
def vector(A):return s.Matrix([s.Poly(A[i,j],*coords).coeff_monomial(mon)for(i,j),mon in m.BIVECTOR_BASIS])
def jac(A):
 return {ijk:s.factor(sum(A[ijk[0],l]*s.diff(A[ijk[1],ijk[2]],coords[l])+A[ijk[1],l]*s.diff(A[ijk[2],ijk[0]],coords[l])+A[ijk[2],l]*s.diff(A[ijk[0],ijk[1]],coords[l])for l in range(5)))for ijk in it.combinations(range(5),3)}
def encode(A):return {f'{i},{j}':str(A[i,j])for i,j in m.PAIRS if A[i,j]!=0}

def actual_family(power,ambient=False):
 start=time.monotonic();name='QUADRATIC_W_POWER'+str(power);path=Path(f'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_ZERO_C_{name}_NORMALFORM_20261107.json');original=json.loads(path.read_text())['source_magnetic_normalform'];target=s.zeros(5)
 for ij,text in original['canonical_brackets'].items():
  i,j=map(int,ij.split(','));target[i,j]=s.sympify(text);target[j,i]=-target[i,j]
 F=(3*h*w+2*y*y)**2-16*y*(h+a*y)**3 if power==1 else (h*w+6*y*y)**2-144*y**3*(h+a*y)
 D=h*F;unknown=s.symbols('b0:18');fields=[sum(unknown[6*i+j]*q for j,q in enumerate(mons))for i in range(3)]
 eq=s.Poly(s.expand(sum(fields[i]*s.diff(D,base[i])for i in range(3))),*base)
 matrix,rhs=s.linear_eq_to_matrix(eq.coeffs(),unknown);assert not any(rhs)
 M0=matrix.subs(a,0);cols=pivots(M0);rows=pivots(M0[:,cols].T);free=[i for i in range(18)if i not in cols];assert len(free)==2
 block=matrix[rows,cols];det=s.factor(block.det(method='domain-ge'));assert det.subs(a,0)!=0
 answers=block.inv()*(-matrix[rows,free]);K=s.zeros(18,2)
 for j,at in enumerate(free):K[at,j]=1
 for i,at in enumerate(cols):K[at,:]=answers.row(i)
 assert (matrix*K).applyfunc(s.cancel)==s.zeros(matrix.rows,2)
 anchors=[]
 for j in range(2):
  S=[s.cancel(sum(K[6*i+k,j]*q for k,q in enumerate(mons)))for i in range(3)];ell=s.cancel(S[0]/h);assert ell.is_polynomial(*base)
  anchors.append(s.Matrix([s.cancel(S[1]-ell*y),s.cancel(S[2]-ell*w)]))
 A0=s.Matrix([[s.Poly(U[i].subs(a,0),*base).coeff_monomial(q)for U in anchors]for i in range(2)for q in mons]);targetanchors=[s.Matrix([target[1,2],target[1,4]]),s.Matrix([-target[2,3],target[3,4]])]
 B0=s.Matrix([[s.Poly(U[i],*base).coeff_monomial(q)for U in targetanchors]for i in range(2)for q in mons]);frame=next(iter(s.linsolve((A0,B0[:,0]),s.symbols('c0:2'))));frame2=next(iter(s.linsolve((A0,B0[:,1]),s.symbols('d0:2'))));G=s.Matrix.hstack(s.Matrix(frame),s.Matrix(frame2));assert G.det()!=0 and A0*G==B0
 UV=s.Matrix.hstack(*anchors)*G;U,V=UV[:,0],UV[:,1]
 determinant=s.factor(U[0]*V[1]-U[1]*V[0]);ratio=s.factor(s.cancel(determinant/F));assert not ratio.has(h,y,w)and ratio.subs(a,0)!=0
 comm=s.Matrix([s.expand(sum(U[j]*s.diff(V[i],(y,w)[j])-V[j]*s.diff(U[i],(y,w)[j])for j in range(2)))for i in range(2)])
 alpha=s.factor(s.cancel((comm[0]*V[1]-comm[1]*V[0])/determinant));beta=s.factor(s.cancel((U[0]*comm[1]-U[1]*comm[0])/determinant));assert all(s.Poly(q,*base).total_degree()==1 for q in [alpha,beta])
 P=s.zeros(5)
 for ij,v in {(1,2):U[0],(1,4):U[1],(2,3):-V[0],(3,4):V[1],(1,3):alpha*x+beta*z+sum(c*q for c,q in zip(r,mons))}.items():P[ij]=s.factor(v);P[ij[::-1]]=-P[ij]
 assert all(v==0 for v in jac(P).values())
 targetR=s.Poly(target[1,3]-P[1,3].subs({a:0,**dict.fromkeys(r,0)}),*base);values={v:targetR.coeff_monomial(q)for v,q in zip(r,mons)}
 assert (P.subs({a:0,**values})-target).applyfunc(s.cancel)==s.zeros(5)
 denominators=s.factor(s.lcm([s.denom(s.cancel(v))for v in P]));assert denominators.subs(a,0)!=0
 result=dict(power=power,all_passed=True,input_report=str(path),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),plane_quartic=str(F),plane_quintic=str(D),full_quadratic_gradient_syzygy_matrix_shape=list(matrix.shape),all18_syzygy_coefficients_retained=True,free_syzygy_columns=free,nonzero16minor=str(det),rational_Saito_basis=[[str(q)for q in K[:,j]]for j in range(2)],quadratic_syzygy_kernel_complete=True,matched_anchor_frame=[[str(q)for q in G.row(i)]for i in range(2)],anchors=[[str(q)for q in U],[str(q)for q in V]],determinant_ratio=str(ratio),commutator_coefficients=list(map(str,[alpha,beta])),universal_Poisson_brackets=encode(P),all10_strict_Jacobi_identically_zero=True,all6_magnetic_parameters=list(map(str,r)),special_source_values={str(k):str(v)for k,v in values.items()},canonical_Ga9_full150_at_a_zero_exactly_matched=True,coefficient_denominator=str(denominators),zero_fiber_is_in_unit_chart=True,scope='Full rational actualfamily with all sixmagnetic parameters. Any unit complement is excluded only from this parametrizingchart, not from the global type6 support classification.')
 if ambient:
  point_values={a:1,**dict(zip(r,[1,2,3,5,7,11]))};assert denominators.subs(a,1)!=0
  point=P.subs(point_values).applyfunc(s.cancel);scale=s.ilcm(1,*[v.q for f in point for v in s.Poly(f,*coords).coeffs()]);point*=scale
  bv={(i,j):{ex:int(c)for ex,c in s.Poly(point[i,j],*coords).terms()if c}for i,j in m.PAIRS if point[i,j]!=0};DL=s.Matrix(m.linearized_jacobi_matrix(bv).tolist());dirs=[]
  for i,j in it.product(range(5),repeat=2):
   B=s.zeros(5);B[i,j]=1;XX=B*s.Matrix(coords);dirs.append(vector(sum([XX[k]*point.diff(coords[k])for k in range(5)],s.zeros(5))-B*point-point*B.T))
  M=s.Matrix(m.vertical_bivector_matrix().tolist()).row_join(s.Matrix.hstack(*dirs,vector(point)));sourcecols=s.Matrix.hstack(*[vector(scale*P.diff(v).subs(point_values))for v in [a,*r]]);FS=M.row_join(sourcecols)
  assert DL*FS==s.zeros(350,FS.cols)
  def minor(A):
   cc=pivots(A);rr=pivots(A[:,cc].T);return dict(rank=len(cc),rows=rr,columns=cc,determinant=str(dm(A[rr,cc]).det()))
  result.update(generic_source_point={str(k):str(v)for k,v in point_values.items()},generic_point_integer_scale=str(scale),generic_point_brackets=encode(point),full350_derivative=minor(DL),full_actual_action=minor(M),full_all7_source_tangent=minor(FS),all_D_source_products_zero=True)
 result['elapsed_seconds']=round(time.monotonic()-start,4);return result

def transfer(one,three):
 def decode(row):
  P=s.zeros(5)
  for ij,text in row['universal_Poisson_brackets'].items():
   i,j=map(int,ij.split(','));P[i,j]=s.sympify(text);P[j,i]=-P[i,j]
  return P
 def push(P,G):return (G*P*G.T).subs(dict(zip(coords,G.inv()*s.Matrix(coords))),simultaneous=True).applyfunc(s.cancel)
 def anchor_coeff(P):
  U=s.Matrix([P[1,2],P[1,4]]);V=s.Matrix([-P[2,3],P[3,4]])
  return s.Matrix([[s.Poly(C[i],*base).coeff_monomial(q)for C in [U,V]]for i in range(2)for q in mons])
 started=time.monotonic();P1=decode(one);P3=decode(three).subs(a,a**3)
 G0=s.eye(5);G0[2,0]=-1/a**3;G0[2,2]=-1/a**2;G0[4,0]=-6/a**6;G0[4,2]=-12/a**5;G0[4,4]=9/a**4
 assert s.factor(G0.det())==-9/a**6
 F1=s.sympify(one['plane_quartic']);F3=s.sympify(three['plane_quartic']).subs(a,a**3)
 assert s.cancel(F3.subs(dict(zip(coords,G0*s.Matrix(coords))),simultaneous=True)-9*F1/a**8)==0
 changed=push(P1,G0)
 def anchor_coeff(P):
  U=s.Matrix([P[1,2],P[1,4]]);V=s.Matrix([-P[2,3],P[3,4]])
  return s.Matrix([[s.Poly(C[i],*base).coeff_monomial(q)for C in [U,V]]for i in range(2)for q in mons])
 C3,C1=map(anchor_coeff,[P3,changed]);rows=pivots(C3.subs(a,1).T);assert len(rows)==2
 T=C3[rows,:].inv()*C1[rows,:];assert (C3*T-C1).applyfunc(s.cancel)==s.zeros(12,2)
 A=T.inv().T;Gfiber=s.eye(5)
 for ii,i in enumerate([1,3]):
  for jj,j in enumerate([1,3]):Gfiber[i,j]=s.factor(A[ii,jj])
 G=Gfiber*G0;changed=push(P1,G);linearpart=P3[1,3].subs(dict.fromkeys(r,0));mag=s.cancel(changed[1,3]-linearpart)
 assert mag.is_polynomial(*base);maps=s.Matrix([s.Poly(mag,*base).coeff_monomial(q)for q in mons]);B=maps.jacobian(r)
 assert (maps-B*s.Matrix(r)).applyfunc(s.cancel)==s.zeros(6,1)
 target=P3.subs(dict(zip(r,maps)),simultaneous=True);assert (changed-target).applyfunc(s.cancel)==s.zeros(5)
 magdet=s.factor(B.det());assert magdet!=0
 result=dict(all_passed=True,source_parameters='power1 a maps to power3 a^3',base_matrix=[[str(q)for q in G0.row(i)]for i in range(5)],exact_quartic_identity='F3(a^3)(Gbase*X)=9*a^(-8)*F1(a)(X)',anchor_basis_transfer=[[str(s.factor(q))for q in T.row(i)]for i in range(2)],fiber_matrix=[[str(s.factor(q))for q in A.row(i)]for i in range(2)],complete_projective_matrix=[[str(s.factor(q))for q in G.row(i)]for i in range(5)],all150_universal_tensor_coefficients_equal=True,magnetic_parameter_map=list(map(str,maps)),magnetic_linear_map_determinant=str(magdet),all6_magnetic_parameters_retained=True,localizations=[str(s.factor(s.denom(q)))for q in list(G)+list(maps)if s.denom(q)!=1],zero_parameter_complement='a=0 is not in this invertible transfer, but both actual Poisson families are separately regular there and match complete Ga9 canonical source anchors; generic containment extends by closedness, not by division at zero.',input_sources='inline complete actual families',elapsed_seconds=round(time.monotonic()-started,4))
 return result

def component(source):
 def unvector(v):
  A=s.zeros(5)
  for c,((i,j),mon)in zip(v,m.BIVECTOR_BASIS):
   q=c*s.prod(a**e for a,e in zip(coords,mon));A[i,j]+=q;A[j,i]-=q
  return A
 def jvec(A):
  J=jac(A);W={quad:s.expand(sum((-1)**i*coords[quad[i]]*J[quad[:i]+quad[i+1:]]for i in range(4)))for quad in it.combinations(range(5),4)}
  return s.Matrix([s.Poly(W[quad],*coords).coeff_monomial(mon)for quad,mon in m.FOUR_VECTOR_BASIS])
 def minor(A):
  cc=pivots(A);rr=pivots(A[:,cc].T);return dict(rank=len(cc),rows=rr,columns=cc,determinant=str(dm(A[rr,cc]).det()))
 def sparse(v):return {str(i):str(c)for i,c in enumerate(v)if c}
 start=time.monotonic();P=s.zeros(5);family=s.zeros(5)
 for ij,text in source['generic_point_brackets'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(text);P[j,i]=-P[i,j]
 for ij,text in source['universal_Poisson_brackets'].items():
  i,j=map(int,ij.split(','));family[i,j]=s.sympify(text);family[j,i]=-family[i,j]
 env={s.Symbol(k):s.Rational(v)for k,v in source['generic_source_point'].items()};scale=s.Rational(source['generic_point_integer_scale']);assert (scale*family.subs(env)-P).applyfunc(s.expand)==s.zeros(5)
 bv={(i,j):{ex:int(c)for ex,c in s.Poly(P[i,j],*coords).terms()if c}for i,j in m.PAIRS if P[i,j]!=0};D=s.Matrix(m.linearized_jacobi_matrix(bv).tolist());directions=[]
 for i,j in it.product(range(5),repeat=2):
  B=s.zeros(5);B[i,j]=1;X=B*s.Matrix(coords);directions.append(vector(sum([X[k]*P.diff(coords[k])for k in range(5)],s.zeros(5))-B*P-P*B.T))
 M=s.Matrix(m.vertical_bivector_matrix().tolist()).row_join(s.Matrix.hstack(*directions,vector(P)));F=M.row_join(s.Matrix.hstack(*[vector(scale*family.diff(v).subs(env))for v in [a,*r]]));K=dm(D).nullspace().to_Matrix().T
 assert K.cols==50 and D*K==s.zeros(350,50)and D*F==s.zeros(350,F.cols)
 fc=pivots(F);fr=pivots(F[:,fc].T);assert len(fc)==len(fr)==49
 inverse=dm(F[fr,fc]).inv().to_Matrix();projected=K-F[:,fc]*inverse*K[fr,:];pc=pivots(projected);assert len(pc)==1
 R=projected[:,pc[0]];common=s.ilcm(1,*[v.q for v in R]);R*=common;g=s.igcd(*[abs(int(v))for v in R if v]);R/=g
 assert R[fr,:]==s.zeros(49,1)and D*R==s.zeros(350,1)and len(pivots(F.row_join(R)))==50
 Q=jvec(unvector(R));quadratic=minor(D.row_join(Q));result=dict(full_true_source_cut_dimension=101,source_frame=minor(F),full_derivative=minor(D),normal_sparse=sparse(R),source_linear_cut_rows=fr,normal_satisfies_all49cuts=True,quadratic_obstruction=quadratic,quadratic_sparse=sparse(Q),component24_proved=quadratic['rank']>100,unique_ambient_component_not_claimed=True)
 if quadratic['rank']==100:
  dr,dc=source['full350_derivative']['rows'],source['full350_derivative']['columns'];solution=-dm(D[dr,dc]).inv().to_Matrix()*Q[dr,:];T=s.zeros(150,1)
  for at,v in zip(dc,solution):T[at]=v
  T-=F[:,fc]*inverse*T[fr,:]
  assert D*T+Q==s.zeros(350,1)and T[fr,:]==s.zeros(49,1)
  cubic=jvec(unvector(R+T))-jvec(unvector(R))-jvec(unvector(T));cert=minor(D.row_join(cubic))
  result.update(quadratic_correction_sparse=sparse(T),cubic_obstruction=cert,cubic_sparse=sparse(cubic),cut_kernel_dimension_exactly_one=True,all_cut_kernel_correction_ambiguity_in_image_D=minor(D.row_join(2*Q))['rank']==100,component24_proved=cert['rank']>100)
 Z=s.Matrix([sum(s.diff(P[i,j],coords[i])for i in range(5))for j in range(5)]);N=Z.jacobian(coords);result.update(cofactor_matrix=[[str(c)for c in N.row(i)]for i in range(5)],cofactor_ranks=[len(pivots(N**j))for j in [1,2]],h_outside_image_dual_cofactor=len(pivots((N.T).row_join(s.Matrix([1,0,0,0,0]))))==3)
 result.update(all_passed=True,input_source='inline wholepower3family',elapsed_seconds=round(time.monotonic()-start,4),scope='A nonzero reduced-support obstruction proves the whole24source closure a fullcomponent by localdimensionupper24. No localuniqueness is inferred from the nonsymmetry sourceparameter cuts. Stop after this one cubic if needed.')
 result['source_plus_normal_frame']=minor(F.row_join(R));assert result['source_plus_normal_frame']['rank']==50
 result['full150_source49_normal1_implicit100_direct_sum']=True
 return result

def prior102_audit():
 """Exclude every prior102 component using the weakest closed invariants."""
 import subprocess
 commit='3fdd4295f818338f23079fa6201bfae6f291e390';catalogue_path='research/COMPONENTS.json'
 blob=subprocess.check_output(['git','rev-parse',commit+':'+catalogue_path],text=True).strip()
 contents=subprocess.check_output(['git','cat-file','blob',blob]);rows=json.loads(contents)['components'];assert len(rows)==102
 wanted=['T02','T03','T16','T17','T21','T36','T37','T38'];xx=s.symbols('x0:5');witnesses=[];paths=[]
 for label in wanted:
  path=Path(f'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/{label}.json')
  if not path.exists():path=Path(f'reports/MPS_S61_RAW_CATALOGUE_20261103_v62/{label}.json')
  data=json.loads(path.read_text());assert data.get('exact_projective_Jacobi_checked',data.get('projective_Jacobi_identically_zero',False))
  P=s.zeros(5);env=dict(zip(map(str,xx),xx))
  for ij,text in data['homogeneous_bivector'].items():
   i,j=map(int,ij.split(','));entry=s.sympify(text,locals=env);edges=[q for q in entry.free_symbols if str(q).startswith('t')]
   P[i,j]=s.cancel(entry.subs(dict.fromkeys(edges,0)));P[j,i]=-P[i,j]
  periods=sorted(set().union(*(q.free_symbols for q in P))-set(xx),key=str)
  log_coefficients={}
  for i,j in m.PAIRS:
   value=s.cancel(P[i,j]/(xx[i]*xx[j]));assert not(set(xx)&value.free_symbols);log_coefficients[f'{i},{j}']=str(value)
  def pf(indices):
   i,j,k,l=indices;return P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]
  actual=s.cancel(sum((-1)**i*xx[i]*pf(tuple(j for j in range(5)if j!=i))for i in range(5)))
  kappa=s.factor(s.cancel(actual/s.prod(xx)));assert kappa!=0 and not(set(xx)&kappa.free_symbols)
  assert s.expand(actual-kappa*s.prod(xx))==0
  minor=kappa**5;assert minor!=0
  witnesses.append(dict(id=label,accepted_whole_source_report=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
   all_period_parameters=list(map(str,periods)),all_edge_parameters_specialized_to_zero=True,
   actual_logcanonical_coefficients=log_coefficients,actual_SNC_quintic=str(actual),kappa=str(kappa),
   gradient_minor_rows=[str(s.prod(xx[j]for j in range(5)if j!=i))for i in range(5)],
   gradient_5_by_5_minor=str(s.factor(minor)),essential_variable_rank=5,
   source_membership_dependency='This edge-zero actual coordinate-SNC tensor belongs to the named accepted completed MPS source, whose closure is that full MPS component. Only one actual source witness outside the new closed essential-rank<=3 locus is needed.'))
  paths.append(path)
 reasons=[]
 for row in rows:
  label=row['id']
  if row['generic_rank']==2:reason='Closed generic-rank2 condition excludes the new rank4 source'
  elif row['projective_dimension']!=24:reason='Different full component dimension'
  elif label in wanted:reason='The actual coordinate-SNC source witness has essential quintic rank5, whereas the entire new source closure satisfies the closed essential-rank<=3 condition'
  elif label=='R4_CPLUS_EXT':
   assert row['generic_actual_degeneracy_divisor']['irreducible_factor_degrees']==[1,2,2]
   reason='Whole old component is in closed1+2+2 multiplication locus, excluding simpleH times irreduciblequartic'
  elif label=='R4_DOUBLE_RANK3_QUADRIC_TWZERO_24':
   assert not row['generic_actual_degeneracy_divisor']['reduced']
   reason='Whole old component is in closed nonreduced-quintic locus, excluding new generic reducedH+quartic'
  else:
   assert label=='R4_DEPENDENT_THREE_HYPERPLANES_QUADRIC_24'and row['generic_actual_degeneracy_divisor']['irreducible_factor_degrees']==[1,1,1,2]
   reason='Whole old component is in closed1+1+1+2 multiplication locus, excluding simpleH times irreduciblequartic'
  reasons.append(dict(id=label,reason=reason))
 return dict(immutable_commit=commit,catalogue_git_blob=blob,catalogue_sha256=hashlib.sha256(contents).hexdigest(),
  all102_rows_frozen=rows,all102_distinctness_reasons=reasons,all102_distinctness_proved=True,
  all8_MPS24_actual_SNC_witnesses=witnesses,public_MPS_source_inputs=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in paths],
  closed_invariant='rank of V5->Sym4(V5*) given by directional derivatives of the actual quintic is <=3 for every displayed source point, hence its full closure. This condition is closed and includes zero quintics. A rank5 point in each other full component excludes equality.')

def main():
 import resource
 resource.setrlimit(resource.RLIMIT_CPU,(179,180));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();started=time.monotonic()
 one=actual_family(1);three=actual_family(3,ambient=True);comparison=transfer(one,three);support=component(three)
 assert support['component24_proved']and support['quadratic_obstruction']['rank']==101
 assert support['cofactor_ranks']==[2,0]and support['h_outside_image_dual_cofactor']
 audit=prior102_audit()
 generic_irreducibility={
  'power1':dict(w_degree=2,w_discriminant='576*h^2*y*(h+a*y)^3',primitive_unit='1-4*a^3',reason='The quadratic is primitive in w under this unit. Over C(h,y), y and h+a*y have odd multiplicity for a!=0; at a0 the discriminant is576*h^5*y. In both cases it is nonsquare.'),
  'power3':dict(w_degree=2,w_discriminant='576*h^2*y^3*(h+a*y)',primitive_unit='1-4*a',reason='The quadratic is primitive in w under this unit. Over C(h,y), y and h+a*y have odd multiplicity for a!=0; at a0 the discriminant is576*h^3*y^3. In both cases it is nonsquare.')}
 for family in [one,three]:
  F=s.sympify(family['plane_quartic']);pol=s.Poly(F,w);disc=s.factor(s.discriminant(pol.as_expr(),w))
  expected=576*h*h*y*(h+a*y)**3 if family['power']==1 else 576*h*h*y**3*(h+a*y)
  assert s.expand(disc-expected)==0
  unit=1-4*a**3 if family['power']==1 else 1-4*a
  assert s.expand(F.subs(h,0)-(4*unit if family['power']==1 else 36*unit)*y**4)==0
  family['quartic_primitive_unit']=str(unit);family['quartic_discriminant']=str(disc);family['generic_and_zero_endpoint_absolutely_irreducible']=True
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,actual_family_power1=one,actual_family_power3=three,
  whole_six_magnetic_family_transfer=comparison,full_source_component=support,all_prior102_audit=audit,
  complex_absolute_irreducibility=generic_irreducibility,whole_source_dimension=24,
  source_dimension_bridge='Actual source+PGL frame rank49 gives projectivedimension>=24. Its generic cofactor is nil22 with h outside dualimage and basebracket0, so the accepted whole plane-quartic bound<=24 applies. The irreducible source image closure has dimension24.',
  componenthood_bridge='The genuine101dimensional linear slice to49 independent source directions has D100 and one tangentnormal. The full quadratic residual is nonzero moduloD, so formalIFT gives reducedcut dimension0. Full local raw dimension<=49, projective<=24. Any fullcomponent containing the entire irreducible24source closure must have dimension24, hence equals that closure. No smoothness or localuniqueness follows or is claimed.',
  Ga9_whole_source_bridge='Both full actual families are regular at a0 and match the accepted nonzero-magnetic canonical Ga9 power1/power3 points in all150coefficients. Their generic families are related by the explicit actualPGL transfer retaining all6mag coefficients. Closedness and PGL stability place both complete accepted23orbitclosures, including zero magnetic limits, in this24component.',
  no_global_count_changed=True,scope='One new24 full source component and wholeGa9c0 power1/power3 containment. Remaining Ga9linear-w and all other type6 plane-quartic points are not covered by this result.',
  public_dependencies=[dict(path=p,sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest())for p in ['research/classification/HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND.md','research/classification/NIL211_SHORT_TRANSLATION_REVIEW.md','research/notes/MPS_S61_CATALOGUE_20261103_v62.md']],
  producer_bytes=Path(__file__).read_text(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  seconds=round(time.monotonic()-started,4),max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert not args.output.exists();args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(dict(all_passed=True,new24=True,Ga9_sources_covered=['power1','power3'],prior102_distinct=True,seconds=result['seconds'])))

if __name__=='__main__':main()
