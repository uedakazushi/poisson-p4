#!/usr/bin/env python3
"""Actual D quadratic restriction to the two-normal CB1/CA plane."""
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
K=s.QQ.frac_field(b)
piv=lambda A:list(DomainMatrix.from_Matrix(A).convert_to(K).rref()[1])
rank=lambda A:len(piv(A))
invmat=lambda A:DomainMatrix.from_Matrix(A).convert_to(K).inv().to_Matrix()
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

JCB=Js[0]
base=rd('reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz');C=s.Matrix(base['fixed_fiber_linear_basis150_by20']);caflags=base['flag_positions'];JCA=C.row_join(s.Matrix.hstack(*(vec(lie(i,j,P)) for i,j in caflags)))
A=JCB.row_join(-JCA);pc=piv(A);assert len(pc)==36
pr=piv(A[:,pc].T);inv=invmat(A[pr,pc])
O=s.Matrix.hstack(*(vec(lie(i,j,P)) for i,j in itertools.product(range(5),repeat=2)));S=O.row_join(vec(Pb.diff(b)));sc=piv(S);assert len(sc)==20;S=S[:,sc]
def left(J):
 rows=piv(J.T);return invmat(J[rows,:]),rows
lcb,rcb=left(JCB);lca,rca=left(JCA)
KS=(lcb*S[rcb,:]).col_join(lca*S[rca,:]);assert (A*KS).applyfunc(s.cancel)==s.zeros(150,20)
nn=DomainMatrix.from_Matrix(A).convert_to(K).nullspace().to_Matrix();ker=[nn[i,:].T for i in range(nn.rows)];assert len(ker)==22
for i,k in enumerate(ker):
 scale=next(q for q in reversed(k) if q!=0);ker[i]=k.applyfunc(lambda q:s.cancel(q/scale))
normals=[];span=S
for k in ker:
 v=JCB*k[:29,:]
 if span.row_join(v).subs(b,1).rank()>span.cols:normals.append(k);span=span.row_join(v)
 if len(normals)==2:break
assert len(normals)==2
F=KS.row_join(s.eye(58)[:,pc]).row_join(s.Matrix.hstack(*normals));detF=F.subs(b,1).det();assert detF!=0
result=dict(stage='frame',seconds=time.monotonic()-start,frame58_by58=[[str(q)for q in row]for row in F.tolist()],frame_determinant=str(detF),implicit36_target_rows=pr,implicit36_source_columns=pc,implicit36_inverse=[[str(q)for q in row]for row in inv.tolist()],normal58_by2=[[str(q)for q in row]for row in s.Matrix.hstack(*normals).tolist()])
out.write_text(json.dumps(result,indent=2)+'\n');print('frame',time.monotonic()-start,flush=True)


D=json.loads((ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63/certificate.json').read_text());parse=lambda q:s.sympify(q,locals={'b':b})
L=s.Matrix(list(map(parse,D['second_difference'])));E=Js[0].row_join(Js[1]).row_join(L);ec=piv(E);assert len(ec)==38;er=piv(E[:,ec].T);Einv=invmat(E[er,ec])
Q=rd('reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63/generic-quadratic.json.gz');raw=s.Matrix([[parse(q)for q in row]for row in Q['raw150_by3']]);R=(raw-E[:,ec]*Einv*raw[er,:]).applyfunc(s.cancel);rr=piv(R.T)
M=json.loads((ROOT/'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63/quotient-map.json').read_text());c0,c1=map(parse,M['normal2']);eta2=s.Matrix([[c1*c1,-2*c0*c1,c0*c0]])
result=dict(all_passed=True,seconds=time.monotonic()-start,whole_CB_tangent_rank=38,actual_D_quadratic_rank=len(rr),basis_rows=rr,basis_coefficients=[[str(q)for q in R[i,:]]for i in rr],eta_square_coefficients=list(map(str,eta2)),rank_with_eta_squared=rank(R.col_join(eta2)),all150_actual_D_quad_coefficients=[[str(q)for q in row]for row in R.tolist()])
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in['actual_D_quadratic_rank','basis_rows','basis_coefficients','rank_with_eta_squared','seconds']},indent=2),flush=True)

assert len(rr)==0 and R==s.zeros(150,3) and result['rank_with_eta_squared']==1
def sha(path):
 d=hashlib.sha256()
 with path.open('rb')as f:
  for q in iter(lambda:f.read(1<<20),b''):d.update(q)
 return d.hexdigest()
paths=[Path(__file__),ROOT/'research/notes/AMBIENT_CB_ACTUAL_D_QUADRATIC_GATE_20261103_v63.md',
 ROOT/'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63/generic-quadratic.json.gz',
 ROOT/'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63/quotient-map.json',
 ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63/certificate.json',
 ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz',
 ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz',
 ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json']
for p in[paths[2],paths[5]]:assert json.loads(p.read_text())['all_passed']
report=dict(id='AMBIENT_CB_ACTUAL_D_QUADRATIC_GATE_20261103_v63',all_passed=True,field='QQ(b)',
 full_physical_CB1_CB2_CA_sources_rebuilt=True,actual_whole_CB_tangent_rank=38,
 full150_by3_true_CB1_minus_CA_quadratics_reused_from_completed_generic_run=True,
 exact38block_quotient_of_all150_normal_quadratics_rank=0,
 actual_D_quadratic_restriction_to_two_normal_plane_zero=True,
 eta_squared_is_nonzero_quadratic_and_cannot_be_actual_D_consequence_here=True,
 all_base20_directions_preserved_by_intrinsic_G_times_b_argument=True,
 other7_CA_transverse_variables_and_their_quadratics_not_settled_by_zero_block=True,
 whole_D_quadratic_ideal_zero_not_claimed=True,
 full_actual_D_ring_nilmodule_rank3versus4_or_embedded_primary_type_not_claimed=True,
 next_gate='actual cubic restriction: possible eta^2*rho socle kill remains open',
 seconds=time.monotonic()-start,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out)})
(folder/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('ACTUAL TWO-NORMAL D QUADRATIC ZERO GATE PASS',report['seconds'],flush=True)
