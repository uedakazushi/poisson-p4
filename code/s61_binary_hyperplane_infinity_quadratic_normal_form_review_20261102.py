"""Independent polynomial coefficient audit of the supplied infinity quadrics."""
import argparse, hashlib, json, time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
import s61_binary_hyperplane_infinity_quadratic_normal_form_20261102 as primary
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_HYPERPLANE_INFINITY_QUADRATIC_NORMAL_FORM_REVIEW_20261102'

def verify():
    start=time.monotonic()
    p='reports/S61_BINARY_HYPERPLANE_INFINITY_QUADRATIC_NORMAL_FORM_20261102.json'
    old=json.loads((ROOT/p).read_text()); fresh=primary.verify()
    assert {k:v for k,v in old.items() if k!='seconds'}=={k:v for k,v in fresh.items() if k!='seconds'}
    a,e=s.symbols('a e'); field=QQ.frac_field(a,e)
    K,*x=ring('x0,x1,x2,x3,x4,x5,x6,x7',field)
    coefficient=lambda v:field.from_sympy(s.sympify(v))
    data_path='code/certificates/s61_binary_hyperplane_infinity_quadrics_20261102.json'
    data=json.loads((ROOT/data_path).read_text())
    f=[K.from_dict({tuple(m):coefficient(c) for m,c in row['terms']}) for row in data['rows']]
    T=[[coefficient(c) for c in row]for row in old['coordinate_matrix8']]
    Ti=[[coefficient(c) for c in row]for row in old['coordinate_inverse8']]
    R=[[coefficient(c) for c in row]for row in old['generator_matrix']]
    Ri=[[coefficient(c) for c in row]for row in old['generator_inverse']]
    def product(A,B):
        return [[sum((A[i][k]*B[k][j] for k in range(len(B))),field.zero) for j in range(len(B[0]))]for i in range(len(A))]
    for A,B in [(T,Ti),(Ti,T),(R,Ri),(Ri,R)]:
        assert product(A,B)==[[field.one if i==j else field.zero for j in range(len(A))]for i in range(len(A))]
    y=[sum((c*z for c,z in zip(row,x)),K.zero)for row in T]
    L,W,C,P,E,H,r,ss=y
    g=[L*L,L*W,L*(E-C),W*W-W*C-L*P,P*E,E*C]
    for z,row in zip(g,R): assert z==sum((c*ff for c,ff in zip(row,f)),K.zero)
    for z,row in zip(f,Ri): assert z==sum((c*gg for c,gg in zip(row,g)),K.zero)
    # Distinct determinant implementation over exact rational-function DomainMatrix.
    def determinant(A):
        from sympy.polys.matrices import DomainMatrix
        return DomainMatrix(A,(len(A),len(A)),field).det()
    assert determinant(T)==coefficient('2278125*(50*a-6*e-1)/(4*a)')
    assert determinant(R)==coefficient('91125/(32768*a**3)')
    paths=[str(Path(__file__).relative_to(ROOT)),p,data_path,
           'research/notes/'+STEM+'.md',*old['input_sha256']]
    hashes={q:hashlib.sha256((ROOT/q).read_bytes()).hexdigest()for q in dict.fromkeys(paths)}
    for q,h in old['input_sha256'].items(): assert hashes[q]==h,q
    return dict(id=STEM,all_passed=True,primary_replay_seconds=fresh['seconds'],
                independent_sparse_QQae_two_way_generator_equalities=True,
                independent_exact_matrix_products_and_determinants=True,
                proof_accepted_on_open='a*(50*a-6*e-1)!=0',
                original135_independently_rederived=False,
                higher_degree_or_ambient_claim=False,input_sha256=hashes,
                seconds=time.monotonic()-start)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    r=verify()
    with args.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
