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

def d5_source_component(output):
 """Isolated whole D5 marked-source anchor; old functions are unchanged."""
 import ast,resource
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(119,120));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
 original=Path('reports/CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json');old=json.loads(original.read_text())['producer_bytes']
 tree=lambda text:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(text).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
 before,after=tree(old),tree(Path(__file__).read_text());assert all(after[k]==v for k,v in before.items())
 F=h*h*y*w+h**4+y**4;R=sum(c*q for c,q in zip(r,mons));family=s.zeros(5)
 entries={(1,2):y*y,(1,4):4*h*h+3*y*w,(2,3):-4*h*h-y*w,(3,4):-16*y*y+3*w*w,(1,3):3*w*x-y*z+R}
 for ij,v in entries.items():family[ij]=v;family[ij[::-1]]=-v
 assert all(v==0 for v in jac(family).values())
 env=dict(zip(r,[1,0,0,0,0,0]));P=family.subs(env)
 bv={(i,j):{ex:int(c)for ex,c in s.Poly(P[i,j],*coords).terms()if c}for i,j in m.PAIRS if P[i,j]!=0}
 D=s.Matrix(m.linearized_jacobi_matrix(bv).tolist());action=[]
 for i,j in it.product(range(5),repeat=2):
  B=s.zeros(5);B[i,j]=1;X=B*s.Matrix(coords);action.append(vector(sum([X[k]*P.diff(coords[k])for k in range(5)],s.zeros(5))-B*P-P*B.T))
 V=s.Matrix(m.vertical_bivector_matrix().tolist());M=V.row_join(s.Matrix.hstack(*action,vector(P)));FS=M.row_join(s.Matrix.hstack(*[vector(family.diff(v))for v in r]))
 assert D*FS==s.zeros(350,FS.cols)
 def cert(A):
  cc=pivots(A);rr=pivots(A[:,cc].T);return dict(shape=list(A.shape),rank=len(cc),rows=rr,columns=cc,determinant=str(dm(A[rr,cc]).det()))
 def sp(v):return {str(i):str(c)for i,c in enumerate(v)if c}
 def unvec(v):
  A=s.zeros(5)
  for c,((i,j),ex)in zip(v,m.BIVECTOR_BASIS):
   q=c*s.prod(a**b for a,b in zip(coords,ex));A[i,j]+=q;A[j,i]-=q
  return A
 def jvec(A):
  J=jac(A);W={quad:s.expand(sum((-1)**i*coords[quad[i]]*J[quad[:i]+quad[i+1:]]for i in range(4)))for quad in it.combinations(range(5),4)}
  return s.Matrix([s.Poly(W[quad],*coords).coeff_monomial(mon)for quad,mon in m.FOUR_VECTOR_BASIS])
 assert jvec(P)==s.zeros(350,1)
 dc,vc,mc,fc=map(cert,[D,V,M,FS]);assert vc['rank']==24 and mc['rank']==fc['rank']==48
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,whole_source_proof='research/classification/PLANE_HYPERPLANE_QUARTIC_D5_TANGENT_SOURCE.md',
  whole_marked_quartic='H=h,F=h^2*y*w+h^4+y^4',all_six_magnetic_family=encode(family),canonical_source_values={str(k):str(v)for k,v in env.items()},canonical_brackets=encode(P),
  actual_projective_degeneracy_quintic=str(16*h*F),all10_strict_source_Jacobi_zero=True,all350_anchor_Jacobi_zero=True,
  full_D=dc,vertical=vc,full51_action=mc,all6_actual_source=fc,all_D_times_actual_directions_zero=True,
  original_all_functions_AST_preserved=True,original_function_AST_sha256={k:hashlib.sha256(v.encode()).hexdigest()for k,v in before.items()},
  actual_projective_orbit_dimension=23,canonical_normal_count=150-dc['rank']-48,independent_review_status='pending',new_component_counted=False)
 if dc['rank']==102:
  result.update(status='SMOOTH_UNIQUE_23_ORBIT_COMPONENT_CANDIDATE',component23_proved_by_dimension=True,local_scheme_smooth=True,
   dimension_argument='Fullraw tangent48 equals the actual48 Euler/PGL/scale orbit. Thus projectivepoint smooth23/uniquecomponent, and its23orbit closure is that fullcomponent. Whole6source closure equals this orbitclosure by explicit pureproof normalizer.')
 elif dc['rank']==101:
  cols=mc['columns'];rows=mc['rows'];F48=M[:,cols];inverse=dm(F48[rows,:]).inv().to_Matrix();K=dm(D).nullspace().to_Matrix().T
  reduced=(K-F48*inverse*K[rows,:]).applyfunc(s.cancel);pc=pivots(reduced);assert len(pc)==1
  N=reduced[:,pc[0]];N=N/next(c for c in N if c);assert D*N==s.zeros(350,1)and N[rows,:]==s.zeros(48,1)
  I=s.eye(150)[:,dc['columns']];U=I-F48*inverse*I[rows,:];frame=F48.row_join(N).row_join(U);direct=cert(frame);assert direct['rank']==150
  Q=jvec(unvec(N));quad=cert(D.row_join(Q));result.update(true_cut_linear_coordinate_rows=rows,true_cut_normal=sp(N),true48plus1plus101_direct_sum=direct,quadratic_Jacobi_normal=sp(Q),quadratic_augmented_D=quad)
  if quad['rank']==102:
   dr,dd=dc['rows'],dc['columns'];di=dm(D[dr,dd]).inv().to_Matrix();res=Q-D[:,dd]*di*Q[dr,:];bad=next(i for i,v in enumerate(res)if v)
   lam=s.zeros(1,350);lam[bad]=1
   for j,k in enumerate(dr):lam[k]-=(D[bad,dd]*di)[j]
   assert lam*D==s.zeros(1,150)and(lam*Q)[0]!=0
   result.update(status='ONE_NORMAL_NONZERO_QUADRATIC_23_COMPONENT_CANDIDATE',component23_proved_by_dimension=True,local_scheme_smooth=False,
    all150_D_annihilating_functional=sp(lam),nonzero_normal_quadratic_value=str((lam*Q)[0]),
    dimension_argument='Actual48symmetryspace plus true102cut is full150directsum; D101 and1normal, nonzero quadratic cokernel forces reducedcut support0. Krull givesrawdimension<=48 matching actualorbit48. The irreducible23orbitclosure is a fullcomponent; no higherjet or exactthickness needed.')
  else:result.update(status='ONE_NORMAL_QUADRATIC_ZERO_STOP_NO_HIGHER_JET',component23_proved_by_dimension=False)
 else:result.update(status='TANGENT_NORMAL_AT_LEAST_TWO_STOP_NO_HIGHER_JET',component23_proved_by_dimension=False)
 inputs=[Path(__file__),original,helper,Path(result['whole_source_proof']),Path('reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_D5_SOURCE_20261107.json')]
 result.update(inputs=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in inputs],seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
  scope='One justified canonicalpoint of the complete D5+tangent/full6mag irreduciblesource. Distinctness against allknowncomponents and independentverification still required before counting.')
 assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['all_passed','status','canonical_normal_count','component23_proved_by_dimension','seconds','max_rss_kib']}))

def mu6_marked_sources(output):
 """One finite batch: complete eleven marked candidates plus two controls."""
 import ast,resource
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(119,120));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
 jobpath=Path('research/classification/jobs/plane_hyperplane_quartic_mu6_sources.json');job=json.loads(jobpath.read_text())
 functions={n.name:hashlib.sha256(ast.dump(n,include_attributes=False).encode()).hexdigest()for n in ast.parse(Path(__file__).read_text()).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
 assert all(functions[k]==v for k,v in job['old_function_AST_sha256'].items())
 inputpath=Path('reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_MARKED_LIST_20261107.json');packet=json.loads(inputpath.read_text());cases=packet['all11_candidates']+packet['controls'];assert len(cases)==13
 env={'h':h,'y':y,'w':w};base=(h,y,w);degree6=[h**(6-i-j)*y**i*w**j for i in range(7)for j in range(7-i)];assert len(degree6)==28
 def cert(A):
  cc=pivots(A);rr=pivots(A[:,cc].T);return dict(shape=list(A.shape),rank=len(cc),rows=rr,columns=cc,determinant=str(dm(A[rr,cc]).det()))
 def mat(A):return [[str(v)for v in A.row(i)]for i in range(A.rows)]
 records=[]
 for case in cases:
  F_old=s.sympify(case['F'],locals=env);H=s.sympify(case['H'],locals=env);first=s.Matrix([[s.Poly(H,*base).coeff_monomial(v)for v in base]])
  G=first
  for i in range(3):
   trial=G.col_join(s.eye(3)[i,:])
   if trial.rank()>G.rank():G=trial
   if G.rows==3:break
  assert G.det()!=0;back=G.inv()*s.Matrix(base);F=s.expand(F_old.subs(dict(zip(base,back)),simultaneous=True));D=h*F
  assert s.expand(H.subs(dict(zip(base,back)),simultaneous=True)-h)==0
  unknown=s.symbols('b0:18');fields=[sum(unknown[6*i+j]*q for j,q in enumerate(mons))for i in range(3)]
  polynomial=s.Poly(s.expand(sum(fields[i]*s.diff(D,base[i])for i in range(3))),*base)
  A,rhs=s.linear_eq_to_matrix([polynomial.coeff_monomial(q)for q in degree6],unknown);assert rhs==s.zeros(28,1)
  K=dm(A).nullspace().to_Matrix().T;assert A*K==s.zeros(28,K.cols)and 18-cert(A)['rank']==K.cols
  anchors=[]
  for j in range(K.cols):
   field=s.Matrix([sum(K[6*i+k,j]*q for k,q in enumerate(mons))for i in range(3)]);ell=s.cancel(field[0]/h);assert ell.is_polynomial(*base)
   anchors.append(s.Matrix([s.expand(field[1]-ell*y),s.expand(field[2]-ell*w)]))
  pairs=[]
  for i,j in it.combinations(range(len(anchors)),2):
   determinant=s.expand(anchors[i][0]*anchors[j][1]-anchors[i][1]*anchors[j][0]);ratio=s.cancel(determinant/F)
   if ratio!=0 and not ratio.free_symbols:pairs.append((i,j,ratio))
   else:assert determinant==0,'Nonzero determinant without Saito divisibility is a defect.'
  row=dict(name=case['name'],input_quartic=case['F'],input_marked_H=case['H'],actual_H_to_h_base_matrix=mat(G),actual_inverse_base_matrix=mat(G.inv()),base_determinant=str(G.det()),transformed_irreducible_quartic=str(F),full28x18_gradient_matrix=mat(A),gradient_rank=cert(A),full_gradient_kernel=mat(K),full_kernel_dimension=K.cols,all_kernel_to_quadratic_h_Casimir_anchors=[list(map(str,U))for U in anchors])
  if not pairs:
   row.update(status='NO_FULLRANK_QUADRATIC_SAITO_FRAME',all_anchor_pair_determinants_zero=True);records.append(row);continue
  assert K.cols==2 and len(pairs)==1
  ii,jj,ratio=pairs[0];U,V=anchors[ii],anchors[jj];det=s.expand(U[0]*V[1]-U[1]*V[0]);assert s.expand(det-ratio*F)==0
  bracket=s.Matrix([s.expand(sum(U[j]*s.diff(V[i],(y,w)[j])-V[j]*s.diff(U[i],(y,w)[j])for j in range(2)))for i in range(2)])
  alpha=s.cancel((bracket[0]*V[1]-bracket[1]*V[0])/det);beta=s.cancel((U[0]*bracket[1]-U[1]*bracket[0])/det)
  assert all(v.is_polynomial(*base)and s.Poly(v,*base).total_degree()<=1 for v in [alpha,beta])
  assert (bracket-alpha*U-beta*V).applyfunc(s.expand)==s.zeros(2,1)
  P=s.zeros(5);R=sum(c*q for c,q in zip(r,mons))
  for ij,v in {(1,2):U[0],(1,4):U[1],(2,3):-V[0],(3,4):V[1],(1,3):alpha*x+beta*z+R}.items():P[ij]=v;P[ij[::-1]]=-v
  strict=jac(P);assert all(v==0 for v in strict.values())
  projective=s.Matrix([s.Poly(s.expand(sum((-1)**i*coords[quad[i]]*strict[quad[:i]+quad[i+1:]]for i in range(4))),*coords).coeff_monomial(mon)for quad,mon in m.FOUR_VECTOR_BASIS]);assert projective==s.zeros(350,1)
  pf=s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3]);assert s.expand(pf-ratio*F)==0
  logarithms=[s.cancel(sum(T[j]*s.diff(F,(y,w)[j])for j in range(2))/F)for T in [U,V]];assert all(v.is_polynomial(*base)for v in logarithms)
  Z=s.Matrix([s.expand(sum(s.diff(P[i,j],coords[i])for i in range(5)))for j in range(5)]);N=Z.jacobian(coords);assert Z==s.Matrix([0,-logarithms[0],0,-logarithms[1],0]);assert N*N==s.zeros(5)
  spatial=(x,y,z,w);AZ=Z[1:,:].subs(h,1);AA=AZ.jacobian(spatial);cc=AZ.subs(dict.fromkeys(spatial,0));nr=N.rank();ar=AA.rank();cin=AA.row_join(cc).rank()==ar;h_inside=N.T.row_join(s.Matrix([1,0,0,0,0])).rank()==nr
  effects=[s.expand(-alpha*l-sum(V[j]*s.diff(l,(y,w)[j])for j in range(2)))for l in [h,y,w]]+[s.expand(-beta*l+sum(U[j]*s.diff(l,(y,w)[j])for j in range(2)))for l in [h,y,w]]
  T=s.Matrix([[s.Poly(v,*base).coeff_monomial(q)for v in effects]for q in mons]);Tcert=cert(T)
  category='nil22_image_type6'if nr==ar==2 and cin else 'nil211_translation_zero_covered_type8'if nr==ar==1 and cin else 'nil211_short_translation_type9'if nr==2 and ar==1 and not cin else 'OTHER_KEEP_EXPLICIT'
  row.update(status='ACTUAL_WHOLE_SIX_MAGNETIC_SOURCE_VERIFIED',anchor_U=list(map(str,U)),anchor_V=list(map(str,V)),Saito_determinant_ratio=str(ratio),commutator=[str(alpha),str(beta)],all6_magnetic_parameters=list(map(str,r)),all150_source_brackets=encode(P),all10_strict_Jacobi_zero=True,all350_projective_Jacobi_zero=True,actual_quintic=str(s.expand(ratio*h*F)),logarithmic_cofactors=list(map(str,logarithms)),canonical_modular_field=list(map(str,Z)),homogeneous_cofactor_matrix=mat(N),homogeneous_cofactor_power_ranks=[(N**j).rank()for j in range(1,6)],distinguished_h_in_dual_image=h_inside,affine_linear_matrix=mat(AA),affine_constant=list(map(str,cc)),affine_power_ranks=[(AA**j).rank()for j in range(1,5)],affine_constant_in_linear_image=cin,canonical_type=category,full_magnetic_translation_matrix=mat(T),full_magnetic_translation_rank=Tcert,magnetic_cokernel_basis=[list(map(str,v))for v in T.T.nullspace()],no_source_membership_or_component_claim=True)
  records.append(row)
 report=dict(task_id='E_reduced_nonnormal',all_passed=True,whole_geometry_input=str(inputpath),input_sha256=hashlib.sha256(inputpath.read_bytes()).hexdigest(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),old_all_FunctionDef_AST_preserved=True,all11_whole_candidates_retained=True,all2_controls_retained=True,cases=records,seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='Complete finite markedgeometry batch with actual H/inverse, full28x18kernel/Saito/6mag/Jac/Pf/cofactor. No additionalambient anchors, no Eparent/count claim; magnetic normalization/component recognition separate.')
 assert not output.exists();output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,cases=[dict(name=v['name'],status=v['status'],type=v.get('canonical_type'),magnetic_rank=v.get('full_magnetic_translation_rank',{}).get('rank'))for v in records],seconds=report['seconds'],max_rss_kib=report['max_rss_kib'])))

if __name__=='__main__':
 import sys
 if '--mu6-marked-sources' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--mu6-marked-sources',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();mu6_marked_sources(args.output)
 elif '--d5-source' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--d5-source',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();d5_source_component(args.output)
 else:main()
