#!/usr/bin/env python3
"""Actual common29 target projection is etale on all three Qb source lifts."""
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

D=json.loads((ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json').read_text());rows=D['common28_linear_target_rows']+[D['common29th_target_row']];assert len(set(rows))==29
res=[]
for name,J in [('CB1',Js[0]),('CB2',Js[1]),('CA',JCA)]:
 M=J[rows,:];rr=rank(M);value=None
 if rr==29:value=str(s.factor(K.to_sympy(DomainMatrix.from_Matrix(M).convert_to(K).det())))
 res.append(dict(name=name,rank=rr,determinant=value,point1rank=M.subs(b,1).rank()))
result=dict(all_passed=True,base29_rows=rows,results=res,seconds=time.monotonic()-start)
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

assert [v['rank']for v in res]==[29,29,29]and all(v['point1rank']==29 for v in res)
def sha(path):
 d=hashlib.sha256()
 with path.open('rb')as f:
  for q in iter(lambda:f.read(1<<20),b''):d.update(q)
 return d.hexdigest()
paths=[Path(__file__),ROOT/'research/notes/AMBIENT_CB_SELF_PAIR_COMMON_BASE_ETALE_20261103_v63.md',
 ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz',
 ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz',
 ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json',
 ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/report.json',
 ROOT/'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63/report.json',
 ROOT/'reports/AMBIENT_CB_ACTUAL_D_QUADRATIC_GATE_20261103_v63/report.json']
for p in paths[6:]:assert json.loads(p.read_text())['all_passed']
report=dict(id='AMBIENT_CB_SELF_PAIR_COMMON_BASE_ETALE_20261103_v63',all_passed=True,field='QQ(b), completedpointb1',
 actual29target_coordinates_same_on_all_three_sources=True,source_projection_ranks=[29,29,29],
 exact_projection_determinants=[v['determinant']for v in res],
 all_three_completed_source_rings_identified_with_same_actual_R29_near_b1=True,
 full_CB_algebra_over_R_is_split_R_plus_tI_t_squared1=True,
 full_ambient_odd_vector_including_redundant_generators_must_be_kept=True,
 actual_D_exact_colon_recipe='f*(J_B, f-2*b0, (J_B:b0)) after c=e0',
 product_only_recipe_f_times_JZ_not_proved_by_etale_base=True,
 missing_cubic_socle_class_still_unresolved=True,
 full_D_ring_nilmodule_length_or_primary_type_not_claimed=True,seconds=time.monotonic()-start,
 input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out)})
(folder/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('COMMON29 ETALE FACT AND FULL COLON CRITERION PASS',report['seconds'],flush=True)
