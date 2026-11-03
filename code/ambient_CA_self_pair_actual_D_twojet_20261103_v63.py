import gzip,itertools,json,sys,time,hashlib,argparse
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
used_inputs=set()
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
STEM='AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63'
out=args.output or ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action
from sympy.polys.matrices import DomainMatrix
from signal import alarm
alarm(120)
import resource
resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
def rank(A):return len(DomainMatrix.from_Matrix(A).convert_to(s.QQ.frac_field(a)).rref()[1])
st=time.monotonic();pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));a=s.Symbol('a')
def rd(stem):
 p=ROOT/'reports'/stem/'certificate.json.gz';used_inputs.add(p);return json.loads(gzip.decompress(p.read_bytes()))
base=rd('AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63');C=s.Matrix(base['fixed_fiber_linear_basis150_by20']);flags=base['flag_positions']
def tensor(v):
 P=s.zeros(5)
 for k,(i,j)in enumerate(pairs):P[i,j]=sum(v[15*k+l]*X[u]*X[w]for l,(u,w)in enumerate(mons));P[j,i]=-P[i,j]
 return P
def trans(g):
 xx=g.inv()*s.Matrix(X);E=s.Matrix([[g[i,k]*g[j,l]-g[i,l]*g[j,k]for k,l in pairs]for i,j in pairs]);M=s.Matrix([[s.Poly(s.expand(xx[k]*xx[l]),*X).coeff_monomial(X[i]*X[j])for k,l in mons]for i,j in mons]);return s.kronecker_product(E,M)
source=rd('AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63');pv=s.Matrix([sum(s.Rational(c)*a**power for power,c in v)for v in source['first_source_physical150']]);qv=s.Matrix([sum(s.Rational(c)*a**power for power,c in v)for v in source['second_source_physical150']]);P=tensor(pv);Q=tensor(qv)
def flagJ(P):return s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flags])
M=s.zeros(5);M[1,0]=1;M[0,3]=-2;back=trans(s.eye(5)-a*M);assert back*qv==pv
J1=C.row_join(flagJ(P));J2local=C.row_join(flagJ(Q));J2=back*J2local
orbit=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in itertools.product(range(5),repeat=2)])
dp=s.Matrix([sum(s.Rational(c)*power*a**(power-1)for power,c in v if power)for v in source['first_source_physical150']]);S=orbit.row_join(dp);assert rank(S)==20
nn=DomainMatrix.from_Matrix(J1.row_join(-J2)).convert_to(s.QQ.frac_field(a)).nullspace().to_Matrix();null=[nn[k,:].T for k in range(nn.rows)];assert len(null)==21
chosen=next(k for k in null if S.row_join(J1*k[:29,:]).subs(a,1).rank()==21)
u1=chosen[:29,:];u2=chosen[29:,:];v=J1*u1;assert (v-J2*u2).applyfunc(s.cancel)==s.zeros(150,1)
def quadratic(P,u):
 A=sum((c*s.eye(5)[:,i]*s.eye(5)[j,:]for c,(i,j)in zip(u[20:,:],flags)),s.zeros(5));delta=tensor(C*u[:20,:])
 return vector(action(A,delta)+action(A,action(A,P))/2)
L=(quadratic(P,u1)-back*quadratic(Q,u2)).applyfunc(s.cancel)
G=s.eye(5);G[0,0]=-3*a;G[0,1]=-1
cb=rd('AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63');B=trans(G.inv())*s.Matrix(cb['full150_by29_actual_CB_source_derivative'])
union=J1.row_join(J2);extended=union.row_join(L)
result=dict(rank_CA1_CB=rank(J1.row_join(B)),rank_CA2_CB=rank(J2.row_join(B)),rank_CA1_CA2=rank(union),rank_plus_quadratic_difference=rank(extended),rank_plus_CB=rank(extended.row_join(B)),seconds=time.monotonic()-st)
print(json.dumps(result,indent=2),flush=True)
(out/'certificate.json').write_text(json.dumps(dict(result=result,u1=list(map(str,u1)),u2=list(map(str,u2)),common_tangent=list(map(str,v)),second_difference=list(map(str,L))),indent=2)+'\n')

assert list(result[k]for k in ['rank_CA1_CB','rank_CA2_CB','rank_CA1_CA2','rank_plus_quadratic_difference','rank_plus_CB'])==[30,38,37,38,38]
review_names=['AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json','AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json','AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63_REVIEW.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name in review_names:
 p=ROOT/'reports'/name;d=json.loads(p.read_text());assert d['all_passed'];used_inputs.add(p)
 for fn,h in d['input_sha256'].items():assert sha(ROOT/fn)==h,fn
used_inputs.update([Path(__file__).resolve(),ROOT/'research/notes'/f'{STEM}.md',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py'])
report=dict(id=STEM,all_passed=True,field='QQ(a); generic nonzero parameter',whole150_source_derivatives_and_true_second_difference=True,CA_branch_CB_clean_N_and_smooth20_meetings=True,whole_CA_tangent_contains_all_CB_tangent_by_true_twojet_identity=True,actual_D_cone_dimension=28,actual_D_cone_embedding_dimension=29,actual_D_reduced_support_smooth_N=True,actual_D_nonreduced=True,full_actual_D_ring_or_nil_length_or_Ass_proved=False,scope='Nonempty open of the actual P_a curve and its G saturation; not B18 origin or all special parameters.',input_sha256={str(p.relative_to(ROOT)):sha(p)for p in sorted(used_inputs)},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-st)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('ACTUAL D TRUE TWOJET NONREDUCED GENERIC Pa PASS',report['seconds'])
