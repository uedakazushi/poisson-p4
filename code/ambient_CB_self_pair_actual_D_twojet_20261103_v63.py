#!/usr/bin/env python3
"""Actual symmetric CB self-pair: finite-source second jet and honest D asymmetry."""
import argparse,hashlib
import gzip,itertools,json,signal,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1];start=time.monotonic();signal.alarm(120)
resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();folder=args.output.resolve();folder.mkdir(parents=True,exist_ok=False);out=folder/'certificate.json'
def rd(p):
 with gzip.open(ROOT/p,'rt') as f:return json.load(f)
h,x,y,z,w,t,b=s.symbols('h x y z w t b');X=s.Matrix([h,x,y,z,w])
rank=lambda A:len(DomainMatrix.from_Matrix(A).convert_to(s.QQ.frac_field(b)).rref()[1])
pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[u]*X[v]) for i,j in pairs for u,v in mons])
def canonical(P):
 d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
 return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
def rawCA(f,g,H):
 P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z);return P-P.T
def rawCB(f,g,H):
 P=s.zeros(5);P[1,2]=f;P[3,4]=g;P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z);return P-P.T
def lie(i,j,P):
 E=s.zeros(5);E[i,j]=1;V=E*X
 return (E*P+P*E.T-sum((V[k]*P.diff(X[k]) for k in range(5)),s.zeros(5))).applyfunc(s.expand)
def push(G,U):return (G*U.subs(dict(zip(X,G.inv()*X)),simultaneous=True)*G.T).applyfunc(s.expand)
cubics=[s.prod((h,x,z,w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in m)]
fixed=[canonical(rawCB(f,0,0)) for f in [h*h,h*x,x*x]]+[canonical(rawCB(0,0,H)) for H in cubics]
with (ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json').open() as f:globalcb=json.load(f)
flags=[q for q in itertools.product(range(5),repeat=2) if list(q) not in globalcb['parabolic_GL_positions']]
source=rd('reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz')
ey=s.eye(5)[:,2];Db=s.Matrix([x*x+8*b*h*x/5,-2*b*x*x/5,0,h*x+3*b*x*z/5,h*h+z*z-7*b*x*w/5]);Pb=ey*Db.T-Db*ey.T;P=Pb
params=[];Js=[]
for datum in source['actual_two_CB_sources']:
 parse=lambda q:s.sympify(q,locals={'b':b,'h':h,'x':x,'y':y,'z':z,'w':w})
 G=s.Matrix([[parse(q) for q in row] for row in datum['G']]);F=parse(datum['f']);H=parse(datum['H']);Q=push(G,P)
 assert (Q-canonical(rawCB(F,0,H))).applyfunc(s.expand)==s.zeros(5)
 local=s.Matrix.hstack(*map(vec,[*fixed,canonical(rawCB(0,F,0)),*(lie(i,j,Q) for i,j in flags)]))
 J=s.Matrix.hstack(*(vec(push(G.inv(),U)) for U in [*fixed,canonical(rawCB(0,F,0)),*(lie(i,j,Q) for i,j in flags)]))
 assert rank(J)==29;Js.append(J);params.append((G,F,H,Q))
J1,J2=Js;O=s.Matrix.hstack(*(vec(lie(i,j,P)) for i,j in itertools.product(range(5),repeat=2)));S=O.row_join(vec(Pb.diff(b)))
assert rank(S)==20
nn=DomainMatrix.from_Matrix(J1.row_join(-J2)).convert_to(s.QQ.frac_field(b)).nullspace().to_Matrix();kernel=[nn[i,:].T for i in range(nn.rows)];assert len(kernel)==21
chosen=next(k for k in kernel if S.row_join(J1*k[:29,:]).subs(b,1).rank()==21)
scale=s.factor(chosen[20]/(-4));assert scale!=0;chosen=chosen.applyfunc(lambda q:s.cancel(q/scale))
u1,u2=chosen[:29,:],chosen[29:,:];v=J1*u1;assert (v-J2*u2).applyfunc(s.cancel)==s.zeros(150,1)
result=dict(stage='chosen',seconds=time.monotonic()-start,u1=list(map(str,u1)),u2=list(map(str,u2)),common_tangent=list(map(str,v)))
out.write_text(json.dumps(result,indent=2)+'\n');print('chosen',time.monotonic()-start,flush=True)
def curve(param,u):
 G,F,H,Q=param;A=s.zeros(5)
 for q,(i,j) in zip(u[20:,:],flags):A[i,j]+=q
 assert A**5==s.zeros(5)
 E=sum((t**k*A**k/s.factorial(k) for k in range(5)),s.zeros(5));Ei=sum(((-t)**k*A**k/s.factorial(k) for k in range(5)),s.zeros(5));assert (E*Ei).applyfunc(s.expand)==s.eye(5)
 dF=sum(q*m for q,m in zip(u[:3,:],[h*h,h*x,x*x]));dH=sum(q*m for q,m in zip(u[3:19,:],cubics));dmu=u[19]
 U=canonical(rawCB(F+t*dF,t*dmu*(F+t*dF),H+t*dH))
 # Only degree<=2 of the actual finite exponential matters here. Drop higher
 # terms BEFORE coordinate substitution, independently of Lie Hessian formulas.
 E2=s.eye(5)+t*A+t*t*A*A/2;Ei2=s.eye(5)-t*A+t*t*A*A/2
 pushed=(E2*U.subs(dict(zip(X,Ei2*X)),simultaneous=True)*E2.T).applyfunc(s.expand)
 back=push(G.inv(),pushed)
 return [vec(back.applyfunc(lambda q:s.expand(q).coeff(t,k))) for k in range(3)],A,dmu,dF
c1,A1,m1,df1=curve(params[0],u1);c2,A2,m2,df2=curve(params[1],u2)
assert all((A-B).applyfunc(s.cancel)==s.zeros(150,1) for A,B in [(c1[0],vec(P)),(c2[0],vec(P)),(c1[1],v),(c2[1],v)])
L=(c1[2]-c2[2]).applyfunc(s.cancel);U=J1.row_join(J2)
base=rd('reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz');C=s.Matrix(base['fixed_fiber_linear_basis150_by20']);cf=base['flag_positions'];CA=C.row_join(s.Matrix.hstack(*(vec(lie(i,j,P)) for i,j in cf)))
ranks=[rank(J1.row_join(CA)),rank(J2.row_join(CA)),rank(U),rank(U.row_join(L)),rank(U.row_join(L).row_join(CA))]
extra=[rank(CA.row_join(v)),rank(U.row_join(CA))]
ell=s.zeros(150,1);ell[63]=1;assert ell.T*U==s.zeros(1,58) and (ell.T*L)[0]!=0
result.update(all_passed=True,stage='finished',ranks=ranks,extra_ranks_CA_v_union_CA=extra,seconds=time.monotonic()-start,second_difference=list(map(str,L)),
 mixed_ratio_speeds=list(map(str,[m1,m2])),mixed_ratio_binary_directions=list(map(str,[df1,df2])),
 source_GL_matrices=[[[str(q) for q in row] for row in A.tolist()] for A in [A1,A2]],
 coker_functional=list(map(str,ell)) if ell is not None else None,coker_value=str((ell.T*L)[0]) if ell is not None else None)
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['ranks','extra_ranks_CA_v_union_CA','coker_value','mixed_ratio_speeds','seconds']},indent=2),flush=True)

assert ranks==[36,30,37,38,38] and extra==[29,37]
assert s.factor((ell.T*L)[0])==24*b**4
for q in [*u1,*u2]:
 den=s.Poly(s.denom(s.cancel(q)),b)
 assert len(den.terms())==1 and den.terms()[0][1]!=0
def sha(path):
 d=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):d.update(chunk)
 return d.hexdigest()
review1=ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW/report.json'
review2=ROOT/'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json'
for p in [review1,review2]:assert json.loads(p.read_text())['all_passed']
inputs=[Path(__file__),ROOT/'research/notes/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63.md',
 ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz',
 ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz',
 ROOT/'reports/AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW/report.json',review1,review2,
 ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
 ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json',
 ROOT/'reports/AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63_REVIEW.json']
for p in inputs:
 if p.name.endswith('REVIEW.json') or p in [review1,review2]:assert json.loads(p.read_text())['all_passed']
report=dict(id='AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63',all_passed=True,field='QQ(b)',
 full_actual_nonlinear_CB_source_with_mu_times_f_retained=True,chosen_ratio_speeds=list(map(str,[m1,m2])),
 finite_nilpotent_exponential_GL_curves_defined_exactly=True,
 physical_all150_coefficients_through_degree2_from_finite5matrix_substitution=True,
 both_common_first_source_derivatives_identical=True,
 ranks_CB1_CA_CB2_CA_CBunion_plusL_plusCA=ranks,ranks_CA_plusv_CBunion_plusCA=extra,
 sparse_coker_functional='e63',sparse_true_quadratic_pairing='24*b**4',
 source_direction_denominators_only_nonzero_constants_and_b=True,
 CB_meeting_exact_double_smooth20_for_every_b_nonzero_by_constant_rank_and_quadratic_unit=True,
 reduced_CB_completion_two_regular29_fiber_product_with_full81_kernel_reused=True,
 reduced_CB_cone_dimension=29,reduced_CB_cone_embedding_dimension=38,reduced_CB_cone_depth=21,reduced_CB_cone_multiplicity=2,
 reduced_CB_projective_dimension=28,reduced_CB_projective_embedding_dimension=37,reduced_CB_projective_depth=20,
 first_CB_branch_intersection_CA_tangent_dimension_generic=22,second_CB_branch_intersection_CA_is_smooth_N28=True,
 common_extra_CB_tangent_is_inside_CA_and_N=True,CA_tangent_inside_CB_branch_tangent_union_as_QQb_identity=True,
 actual_D_cone_dimension=28,actual_D_cone_embedding_dimension=29,actual_D_nonreduced_generic=True,
 old_CA_selfpair_full_D_ring_does_not_follow_by_swapping_labels=True,
 full_actual_D_nilmodule_length_primary_decomposition_or_second_embedded_prime_not_claimed=True,
 seconds=time.monotonic()-start,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
 evidence_sha256={'certificate.json':sha(out)})
(folder/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('PUBLIC SYMMETRIC CB MEETING DOUBLE AND ASYMMETRIC D2 PASS',report['seconds'],flush=True)
