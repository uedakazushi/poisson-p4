from pathlib import Path
from fractions import Fraction
import gzip,hashlib,json,resource,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
resource.setrlimit(resource.RLIMIT_CPU,(60,61))
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
begun=time.monotonic()
path=Path('reports/CLASSIFICATION_NORMAL_LINEAR_OBSTRUCTION_20261107.json')
r=json.loads(path.read_text()); prime=r['prime']; field=s.GF(prime)
source=Path(r['source_cross_path'])
assert hashlib.sha256(source.read_bytes()).hexdigest()==r['source_cross_sha256']
lines=gzip.decompress(source.read_bytes()).decode().splitlines()
assert lines[0]=='129 126 945' and len(lines)==259
def conv(v):
 q=Fraction(v); assert q.denominator%prime
 return q.numerator*pow(q.denominator,-1,prime)%prime
C=[[conv(x) for x in line.split()] for line in lines[1:130]]
D=[[conv(x) for x in line.split()[:42]] for line in lines[130:]]
selected=r['selected_cross_rows']; assert len(set(selected))==126
CS=DomainMatrix.from_list([C[i] for i in selected],field)
DS=DomainMatrix.from_list([D[i] for i in selected],field)
Cdet=int(CS.det())%prime; assert Cdet==19216
inverse=CS.inv()
I=DomainMatrix.eye(126,field).to_dense()
assert CS*inverse==I and inverse*CS==I
T=-(inverse*DS)
Ts=[T.extract(range(42*h,42*h+42),range(42)) for h in range(3)]
def ints(A):return [[int(c)%prime for c in row] for row in A.to_list()]
assert [ints(A)for A in Ts]==r['operators_T3_mod_p']
actualL=[[0]*129 for _ in range(126)]
for i,row in enumerate(ints(inverse)):
 for j,k in enumerate(selected):actualL[i][k]=row[j]
assert actualL==r['left_inverse_mod_p']
comm=[Ts[i]*Ts[j]-Ts[j]*Ts[i]for i,j in [(0,1),(0,2),(1,2)]]
rows=[row for A in comm for row in ints(A)]
small=DomainMatrix.from_list([rows[i]for i in r['selected_commutator_rows']],field)
value=int(small.det())%prime
assert value==21866!=0
assert r['selected_commutator_rows']==list(range(42))
universal_path=Path('reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json')
point_path=Path('reports/S61_KNOWN_FIBER_AUDIT_20261031.json')
frame_path=Path('reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz')
section_path=Path('reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json')
universal=json.loads(universal_path.read_text())
point=next(x for x in json.loads(point_path.read_text())['points'] if x['component_id']=='T01')['S61_coordinates_u1_through_u61']
subs={s.Symbol('u'+str(i+1)):s.Rational(v)for i,v in enumerate(point)}
qbar=s.Matrix([s.sympify(x).subs(subs)for x in universal['q_center']])
N=s.Matrix([[s.sympify(x).subs(subs)for x in row]for row in universal['q_linear_matrix_N']])
kern=s.Matrix([[s.Rational(v)for v in row]for row in universal['generic_plane_recipe']['sample_kernel_T']])
section=json.loads(section_path.read_text());ids=universal['weight2_u_indices_zero_based']
q52=s.Matrix([s.Rational(section['q52_S61_coordinates'][i])for i in ids])
q51=s.Matrix([s.Rational(section['q51_S61_coordinates'][i])for i in ids])
assert 3*qbar==2*q51+q52
W=(N*kern).row_join(q52-qbar)
frame_minor=W[:3,:].det(); assert frame_minor==s.Rational(162221,576240)
E=W*W[:3,:].inv()
saved=s.Matrix([[s.Rational(x)for x in row]for row in json.loads(gzip.decompress(frame_path.read_bytes()))['embedding']])
assert E==saved
shift=s.Matrix([s.Rational(421,48),-s.Rational(767,84),s.Rational(499,84)])
assert qbar+E*shift==s.Matrix([s.Rational(point[i])for i in ids])
result=dict(task_id='G',all_passed=True,independent_backend='SymPy DomainMatrix finite-field exact inverse/determinant; producer used NumPy elimination',source_cross_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),producer_report_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),cross_determinant_mod32003=Cdet,first_commutator_T0_T1_determinant_mod32003=value,all_three_operators_and_selected_left_inverse_match=True,rank3_branch_excluded_over_QQ=True,generic_transport='Requires the rational original-equation family and nonzero specialization; reviewed separately in NORMAL_DIRECTION_INCIDENCE.md',generic_exhaustion_proved=False,seconds=time.monotonic()-begun)
result['rational_global_center45_match']=True
result['canonical_frame45x3_match']=True
result['sample_global_frame_minor']=str(frame_minor)
result['additional_source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [universal_path,point_path,frame_path,section_path]}
result['independent_checker_source']=Path(__file__).read_text()
result['independent_checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
out=Path('reports/CLASSIFICATION_NORMAL_LINEAR_OBSTRUCTION_REVIEW_20261107.json')
with out.open('x')as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
print({k:v for k,v in result.items()if k not in ['independent_checker_source']})
