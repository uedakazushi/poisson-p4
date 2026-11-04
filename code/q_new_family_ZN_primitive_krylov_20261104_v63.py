#!/usr/bin/env python3
"""ONE oldpoint/prime finite-free10 Krylov nonzero certificate; NO new beta point."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63';P=32003
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mod(z):
 q=Q(str(z));assert q.denominator%P;return q.numerator*pow(q.denominator,-1,P)%P
def mm(A,B):return [[sum(A[i][k]*B[k][j]for k in range(4))%P for j in range(4)]for i in range(4)]
def determinant(A):
 A=[r[:]for r in A];value=1;rank=0
 for k in range(len(A)):
  j=next((i for i in range(k,len(A))if A[i][k]%P),None)
  if j is None:return 0,rank
  if j!=k:A[j],A[k]=A[k],A[j];value=-value
  pivot=A[k][k]%P;value=value*pivot%P;inverse=pow(pivot,-1,P);rank+=1
  for i in range(k+1,len(A)):
   q=A[i][k]*inverse%P
   if q:A[i]=[(a-q*b)%P for a,b in zip(A[i],A[k])]
 return value%P,rank
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 pointpath=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/certificate.json.gz'
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 Gpath=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63/primitive.txt'
 Cpath=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_PRS_CONTENT_PREFLIGHT_20261104_v63/C1-primitive.txt.gz'
 shape=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json'
 gates=[ROOT/'reports'/stem/'report.json'for stem in['Q_NEW_FAMILY_TWO_SHAPE_CONTENT_PRIMITIVITY_20261104_v63_REVIEW','Q_NEW_FAMILY_TWO_SHAPE_COPRIME_DEGREES_20261104_v63_REVIEW','Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in gates)
 point,spin,data=map(read,[pointpath,spinpath,shape]);assert point['same_existing_point_index']==0
 certificate=dict(prime=P,oldpoint_index=0,point_provenance=str(pointpath.relative_to(ROOT)),same_saved_source_point=point['same_source_point'],beta16_QQ=point['beta16'])
 stage='oldpoint_loaded';positive=False;delta=None;Glead=None;error=None;checks_passed=True
 try:
  beta=[mod(z)for z in point['beta16']];mon=spin['quadratic_monomials'];products=[beta[i]*beta[j]%P for i,j in mon];mats=[]
  for key in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
   v=[sum(mod(c)*z for c,z in zip(row,products))%P for row in spin[key]];mats.append([v[i:i+4]for i in range(0,16,4)])
  C,D=mats;M=[[(-C[i][j]+2*D[i][j])%P for j in range(4)]for i in range(4)];M2=mm(M,M);F4=sum(M2[i][i]for i in range(4))%P
  pencil=[0]*5
  for perm in itertools.permutations(range(4)):
   sign=(-1)**sum(perm[i]>perm[j]for i in range(4)for j in range(i+1,4));cur=[sign]
   for i,j in enumerate(perm):
    nxt=[0]*(len(cur)+1)
    for k,c in enumerate(cur):nxt[k]=(nxt[k]+c*C[i][j])%P;nxt[k+1]=(nxt[k+1]+c*D[i][j])%P
    cur=nxt
   pencil=[(a+b)%P for a,b in zip(pencil,cur)]
  cc,cd,dd=mm(C,C),mm(C,D),mm(D,D)
  assert sum(cc[i][i]for i in range(4))%P==F4*pow(5,-1,P)%P
  assert sum(cd[i][i]for i in range(4))%P==-F4*pow(10,-1,P)%P
  assert sum(dd[i][i]for i in range(4))%P==F4*pow(10,-1,P)%P
  assert(F4*F4-400*pencil[0]-1600*pencil[4])%P==0
  certificate.update(beta16_modp=beta,C4_modp=C,D4_modp=D,F4_modp=F4,pencil_p0_to_p4_modp=pencil,quartic_and_pencil_normalization_checks=True)
  if not F4 or not pencil[0]:stage='old_F4_or_p0_unit_gate_failed';raise ArithmeticError(stage)
  invF2=pow(F4*F4%P,-1,P);base=dict(R=(pencil[3]+4*pencil[4])*pow(pencil[0],-1,P)%P,Z0=pencil[0]*invF2%P,Z1=pencil[1]*invF2%P,Z2=pencil[2]*invF2%P);certificate['coefficient_base_modp']=base
  # Evaluate only the FOUR coefficient parameters, leaving r as the algebra generator.
  G=[0]*11
  with Gpath.open()as f:
   assert f.readline().split()==['6','891']
   for line in f:
    c,*e=map(int,line.split());assert e[0]==0 and e[1]<=10
    z=c%P
    for name,k in zip(['R','Z0','Z1','Z2'],e[2:]):z=z*pow(base[name],k,P)%P
    G[e[1]]=(G[e[1]]+z)%P
  Glead=G[10];certificate['specialized_G_coefficients_ascending']=G;certificate['specialized_G_leading_coefficient']=Glead
  if not Glead:stage='G10_degree_leading_unit_gate_failed';raise ArithmeticError(stage)
  monic=[c*pow(Glead,-1,P)%P for c in G];certificate['monic_G_coefficients_ascending']=monic
  def reduce(v):
   v=list(v)+[0]*max(0,10-len(v));v=[c%P for c in v]
   for k in range(len(v)-1,9,-1):
    z=v[k]
    if z:
     for j in range(10):v[k-10+j]=(v[k-10+j]-z*monic[j])%P
    v[k]=0
   return v[:10]
  one=[1]+[0]*9;rr=[0,1]+[0]*8
  def add(*terms):return [sum(c*v[i]for c,v in terms)%P for i in range(10)]
  def mul(a,b):
   v=[0]*19
   for i,x in enumerate(a):
    for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%P
   return reduce(v)
  def pw(a,n):
   b=one
   for _ in range(n):b=mul(b,a)
   return b
  C1=[]
  with gzip.open(Cpath,'rt')as f:
   assert f.readline().split()==['5','2']
   for index in range(2):
    n=int(f.readline());coef=[0]*32
    for _ in range(n):
     c,z2,z1,z0,re,er=map(int,f.readline().split());assert er<=31
     z=c%P
     for name,k in[('Z2',z2),('Z1',z1),('Z0',z0),('R',re)]:z=z*pow(base[name],k,P)%P
     coef[er]=(coef[er]+z)%P
    C1.append(coef)
   assert f.read().strip()==''
  B=reduce(C1[0]);A=reduce(C1[1]);certificate.update(C1_evaluated_coefficients_constant_then_linear=C1,C1_reduced_constant_B=B,C1_reduced_linear_A=A)
  rvar,svar,Rvar,Zvar,avar,bvar=S.symbols('r s R Z0 A B');v3=(rvar,svar,Rvar)
  def original(rec):return S.Poly.from_dict({tuple(e):S.Rational(n,d)for n,d,*e in rec['coefficients']},v3,domain=S.QQ)
  N1rec=data['rational_coefficient_pairs']['E1']['numerator'];N2rec=data['rational_coefficient_pairs']['E2']['numerator'];Fnrec=data['factors']['Fn'];Jrec=data['factors']['J'];N1,N2,Fn,J=map(original,[N1rec,N2rec,Fnrec,Jrec])
  assert data['rational_coefficient_pairs']['E1']['constant']=='1/20'and data['rational_coefficient_pairs']['E2']['constant']=='5/4'
  assert N1.degree(svar)==3 and N2.degree(svar)==2 and Fn.degree(svar)==2 and J.degree(svar)==1
  # Exact symbolic balanced-degree identity, without inserting the LARGE C1 polynomials.
  def hom(poly,degree):return sum(co*rvar**e[0]*Rvar**e[2]*avar**(degree-e[1])*(-bvar)**e[1]for e,co in poly.terms())
  dc=rvar*(rvar-1)**2;JH=hom(J,1);FH=hom(Fn,2);NH1=hom(N1,3);NH2=hom(N2,2)
  Usymbol=Zvar*dc*JH*NH1+25*Zvar**2*FH*NH2;Vsymbol=20*dc*FH*JH**2
  origU=S.Poly(Zvar*dc*J.as_expr()*N1.as_expr()+25*Zvar**2*Fn.as_expr()*N2.as_expr(),rvar,svar,Rvar,Zvar,domain=S.QQ)
  origV=S.Poly(20*dc*Fn.as_expr()*J.as_expr()**2,rvar,svar,Rvar,Zvar,domain=S.QQ)
  def hom4(poly):return sum(co*rvar**e[0]*Rvar**e[2]*Zvar**e[3]*avar**(4-e[1])*(-bvar)**e[1]for e,co in poly.terms())
  assert S.Poly(hom4(origU)-Usymbol,rvar,Rvar,Zvar,avar,bvar,domain=S.QQ).is_zero
  assert S.Poly(hom4(origV)-Vsymbol,rvar,Rvar,Zvar,avar,bvar,domain=S.QQ).is_zero
  # Evaluate SAME recipe in finite-free A_p, with NO division by A or V.
  def hom_mod(rec,degree):
   terms=[]
   for n,d,er,es,re in rec['coefficients']:
    z=mod(Q(int(n),int(d)))*pow(base['R'],re,P)%P
    value=mul(pw(rr,er),mul(pw(A,degree-es),pw(add((-1,B)),es)));terms.append((z,value))
   return add(*terms)
  jh,fh,nh1,nh2=[hom_mod(rec,k)for rec,k in[(Jrec,1),(Fnrec,2),(N1rec,3),(N2rec,2)]]
  dcr=mul(rr,pw(add((1,rr),(-1,one)),2));U=add((base['Z0'],mul(dcr,mul(jh,nh1))),(25*base['Z0']**2,mul(fh,nh2)));V=add((20,mul(dcr,mul(fh,pw(jh,2)))))
  pu=[one];pv=[one]
  for i in range(9):pu.append(mul(pu[-1],U));pv.append(mul(pv[-1],V))
  columns=[mul(pu[j],pv[9-j])for j in range(10)];matrix=[[columns[j][i]for j in range(10)]for i in range(10)];delta,rank=determinant(matrix);positive=bool(delta)
  certificate.update(symbolic_degree4_balanced_homogeneity_identity_verified=True,homogeneous_pair_recipe=dict(U='Z0*r*(r-1)^2*JH*NH1+25*Z0^2*FH*NH2',V='20*r*(r-1)^2*FH*JH^2'),
   intermediate_residues=dict(JH=jh,FH=fh,NH1=nh1,NH2=nh2,r_times_rminus1_squared=dcr),U=U,V=V,
   Krylov_column_order='U^j V^(9-j),j0..9',basis_order='1,r,...,r^9',Krylov_matrix_rows=matrix,Krylov_det_mod32003=delta,Krylov_rank_mod32003=rank,
   no_A_or_V_inverse_in_special_algebra=True,special_V_unit_not_assumed_or_required=True)
  stage='ZN_primitive_degree10_by_nonzero_finitefree_coefficient_determinant'if positive else 'one_oldpoint_Krylov_det_zero_inconclusive_STOP'
 except ArithmeticError as exc:error=str(exc)
 except Exception as exc:error=repr(exc);checks_passed=False;stage='implementation_or_input_check_failed_STOP'
 signal.alarm(0)
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n');assert(out/'certificate.json').stat().st_size<32*1024*1024
 report=dict(id=STEM,all_passed=checks_passed,status=stage,gate_stop_reason=error,prime=P,old_saved_beta_point_index=0,number_of_oldpoints_used=1,number_of_newpoints_used=0,
  G10_leading_unit=bool(Glead),Krylov_determinant_mod32003=delta,ZN_primitive_degree10_proved=positive,
  finitefree_model='localization of Z_(p)[R,Z0,Z1,Z2] at oldpoint coefficient maximalideal and Gleading unit; monicG gives free rank10',
  generic_V_nonzero_reason='actual source chart Fn,J,r,r-1 units and birational linearC1 recovery; fieldG10 irreducible',
  full_NEW_target_index_one_NOT_claimed_in_this_packet=True,
  no_full_symbolic_N_norm_elimination_newInv14_source26_or_automatic_point_prime_retry=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),pointpath,spinpath,Gpath,Cpath,shape,*gates]},
  evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','status','G10_leading_unit','Krylov_determinant_mod32003','ZN_primitive_degree10_proved','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
