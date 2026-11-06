import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
resource.setrlimit(resource.RLIMIT_CPU,(45,46));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();path=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_CURVE_20261107.json');src=json.loads(path.read_text());xx=s.symbols('x0:5');t=s.Symbol('t');P=s.zeros(5)
for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(src['family'].get(str(ij),'0'));P[ij[::-1]]=-P[ij]
div=s.Matrix([sum(s.diff(P[i,j],xx[i]) for i in range(5)) for j in range(5)]);E=s.Matrix(xx);U=(P-(E*div.T-div*E.T)/5).applyfunc(s.expand)
assert all(s.expand(sum(s.diff(U[i,j],xx[i]) for i in range(5)))==0 for j in range(5))
q=-t*xx[2]**2+10*xx[0]*xx[2]+5*xx[1]**2-5*xx[3]**2
X=[]
for j in range(5):
    value=s.expand(sum(s.diff(q,xx[i])*U[i,j] for i in range(5)));a,r=s.div(value,q,*xx);assert r==0;X.append(s.expand(a))
A=s.Matrix(X).jacobian(xx);lam=s.Symbol('L');char=s.factor(A.charpoly(lam).as_expr());A1=A.subs(t,1);eig=A1.eigenvals()
D=s.sympify(src['Euler_contracted_Pfaffian']);c=s.cancel(D/q);assert s.denom(c)==1 or s.denom(c).is_number
out={'task_id':'E_reduced_nonnormal','source_report':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'unimodular_brackets':{str(ij):str(U[ij]) for ij in itertools.combinations(range(5),2) if U[ij]},'canonical_cofactor':list(map(str,X)),'canonical_cofactor_matrix':[[str(v) for v in row] for row in A.tolist()],'characteristic_polynomial':str(char),'at_t_one_eigenvalues':{str(k):int(v) for k,v in eig.items()},'q_cofactor_exact':True,'unimodular_exact':True,'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start}
Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_CURVE_RECOGNITION_20261107.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['canonical_cofactor','characteristic_polynomial','at_t_one_eigenvalues','seconds']},indent=2))
