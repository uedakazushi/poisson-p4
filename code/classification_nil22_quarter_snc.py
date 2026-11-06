#!/usr/bin/env python3
"""Universal actual Poisson curve and an exact SNC orbit degeneration."""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
import classification_rank3_mixed_component_review as engine
ROOT=Path(__file__).resolve().parents[1]

def main():
    started=time.monotonic();h,x,y,z,w=s.symbols('h x1 x2 x3 x4');xx=(h,x,y,z,w);t=s.Symbol('t')
    path=ROOT/'reports/CLASSIFICATION_NIL22_QUARTER_WEIGHTED_SLICE_20261107.json'
    packet=json.loads(path.read_text());tau=s.Symbol('tau0')
    env={str(v):v for v in xx};env['tau0']=tau
    P=s.zeros(5)
    for key,value in packet['full_curve_brackets'].items():
        i,j=map(int,key.split(','));P[i,j]=s.expand(s.sympify(value,locals=env).subs(tau,135*t)/45)
    P-=P.T
    def sparse(M):return {ij:dict(s.Poly(M[ij],*xx).terms()) for ij in engine.PAIRS if M[ij]!=0}
    assert engine.jacobi(sparse(P))==s.zeros(350,1)
    M=s.Matrix([[1,0,0,0,0],[0,1,1,0,0],[0,-1,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    inv=M.inv()*s.Matrix(xx)
    T=(M*P.subs(t,1)*M.T).subs(dict(zip(xx,inv)),simultaneous=True).applyfunc(s.expand)
    expected={(0,3):-2*x*w,(1,2):4*h*(x-y)-4*x*y,(1,3):2*x*w-x*z-h*z,
              (2,3):(-2*x+4*y)*w+y*z-h*z,(1,4):(x+h)*w,(2,4):(h-y)*w,
              (3,4):28*h*h+3*z*w+w*w}
    assert all(s.expand(T[i,j]-expected.get((i,j),0))==0 for i,j in engine.PAIRS)
    weights=[1,0,0,-1,1];limit=s.zeros(5);terms=[]
    for i,j in engine.PAIRS:
        for e,c in s.Poly(T[i,j],*xx).terms():
            if not c:continue
            weight=sum(a*b for a,b in zip(weights,e))-weights[i]-weights[j]
            assert weight>=0
            terms.append({'pair':[i,j],'monomial':e,'coefficient':str(c),'weight':weight})
            if weight==0:limit[i,j]+=c*s.prod(a**b for a,b in zip(xx,e))
    limit-=limit.T
    assert all(limit[0,j]==0 for j in range(5))
    pf=s.expand(h*(limit[1,2]*limit[3,4]-limit[1,3]*limit[2,4]+limit[1,4]*limit[2,3]))
    assert pf==-12*h*x*y*z*w
    assert engine.jacobi(sparse(limit))==s.zeros(350,1)
    original=ROOT/'reports/CLASSIFICATION_NIL22_QUARTER_CANONICAL_20261107.json'
    input_paths=[path,original,Path(__file__),Path(engine.__file__)]
    result={'task_id':'E_multiple_hyperplane','all_passed':True,'method':'Entire actual curve, invertible coordinate change, every coefficient weight, rank-four SNC limit; no local-dimension assumption needed for coverage.',
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths},
            'universal_curve_brackets':{str(i)+','+str(j):str(P[i,j]) for i,j in engine.PAIRS if P[i,j]!=0},
            'whole350_projective_Jacobi_identically_zero':True,
            'coordinate_matrix':[[str(v) for v in row] for row in M.tolist()],
            'coordinate_matrix_determinant':str(M.det()),
            'transformed_parameter_one_brackets':{str(i)+','+str(j):str(T[i,j]) for i,j in engine.PAIRS if T[i,j]!=0},
            'degeneration_weights':weights,'every_coefficient_term':terms,
            'limit_brackets':{str(i)+','+str(j):str(limit[i,j]) for i,j in engine.PAIRS if limit[i,j]!=0},
            'limit_pfaffian_section':str(pf),'entire_nil22_A_one_quarter_source_MPS_covered':True,
            'component_label_not_identified':True,'seconds':time.monotonic()-started}
    out=ROOT/'reports/CLASSIFICATION_NIL22_QUARTER_SNC_20261107.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'limit_pfaffian':str(pf)}))

HISTORICAL_QUARTER_PRODUCER_BYTES = '#!/usr/bin/env python3\n"""Universal actual Poisson curve and an exact SNC orbit degeneration."""\nimport hashlib,itertools,json,time\nfrom pathlib import Path\nimport sympy as s\nimport classification_rank3_mixed_component_review as engine\nROOT=Path(__file__).resolve().parents[1]\n\ndef main():\n    started=time.monotonic();h,x,y,z,w=s.symbols(\'h x1 x2 x3 x4\');xx=(h,x,y,z,w);t=s.Symbol(\'t\')\n    path=ROOT/\'reports/CLASSIFICATION_NIL22_QUARTER_WEIGHTED_SLICE_20261107.json\'\n    packet=json.loads(path.read_text());tau=s.Symbol(\'tau0\')\n    env={str(v):v for v in xx};env[\'tau0\']=tau\n    P=s.zeros(5)\n    for key,value in packet[\'full_curve_brackets\'].items():\n        i,j=map(int,key.split(\',\'));P[i,j]=s.expand(s.sympify(value,locals=env).subs(tau,135*t)/45)\n    P-=P.T\n    def sparse(M):return {ij:dict(s.Poly(M[ij],*xx).terms()) for ij in engine.PAIRS if M[ij]!=0}\n    assert engine.jacobi(sparse(P))==s.zeros(350,1)\n    M=s.Matrix([[1,0,0,0,0],[0,1,1,0,0],[0,-1,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])\n    inv=M.inv()*s.Matrix(xx)\n    T=(M*P.subs(t,1)*M.T).subs(dict(zip(xx,inv)),simultaneous=True).applyfunc(s.expand)\n    expected={(0,3):-2*x*w,(1,2):4*h*(x-y)-4*x*y,(1,3):2*x*w-x*z-h*z,\n              (2,3):(-2*x+4*y)*w+y*z-h*z,(1,4):(x+h)*w,(2,4):(h-y)*w,\n              (3,4):28*h*h+3*z*w+w*w}\n    assert all(s.expand(T[i,j]-expected.get((i,j),0))==0 for i,j in engine.PAIRS)\n    weights=[1,0,0,-1,1];limit=s.zeros(5);terms=[]\n    for i,j in engine.PAIRS:\n        for e,c in s.Poly(T[i,j],*xx).terms():\n            if not c:continue\n            weight=sum(a*b for a,b in zip(weights,e))-weights[i]-weights[j]\n            assert weight>=0\n            terms.append({\'pair\':[i,j],\'monomial\':e,\'coefficient\':str(c),\'weight\':weight})\n            if weight==0:limit[i,j]+=c*s.prod(a**b for a,b in zip(xx,e))\n    limit-=limit.T\n    assert all(limit[0,j]==0 for j in range(5))\n    pf=s.expand(h*(limit[1,2]*limit[3,4]-limit[1,3]*limit[2,4]+limit[1,4]*limit[2,3]))\n    assert pf==-12*h*x*y*z*w\n    assert engine.jacobi(sparse(limit))==s.zeros(350,1)\n    original=ROOT/\'reports/CLASSIFICATION_NIL22_QUARTER_CANONICAL_20261107.json\'\n    input_paths=[path,original,Path(__file__),Path(engine.__file__)]\n    result={\'task_id\':\'E_multiple_hyperplane\',\'all_passed\':True,\'method\':\'Entire actual curve, invertible coordinate change, every coefficient weight, rank-four SNC limit; no local-dimension assumption needed for coverage.\',\n            \'input_sha256\':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths},\n            \'universal_curve_brackets\':{str(i)+\',\'+str(j):str(P[i,j]) for i,j in engine.PAIRS if P[i,j]!=0},\n            \'whole350_projective_Jacobi_identically_zero\':True,\n            \'coordinate_matrix\':[[str(v) for v in row] for row in M.tolist()],\n            \'coordinate_matrix_determinant\':str(M.det()),\n            \'transformed_parameter_one_brackets\':{str(i)+\',\'+str(j):str(T[i,j]) for i,j in engine.PAIRS if T[i,j]!=0},\n            \'degeneration_weights\':weights,\'every_coefficient_term\':terms,\n            \'limit_brackets\':{str(i)+\',\'+str(j):str(limit[i,j]) for i,j in engine.PAIRS if limit[i,j]!=0},\n            \'limit_pfaffian_section\':str(pf),\'entire_nil22_A_one_quarter_source_MPS_covered\':True,\n            \'component_label_not_identified\':True,\'seconds\':time.monotonic()-started}\n    out=ROOT/\'reports/CLASSIFICATION_NIL22_QUARTER_SNC_20261107.json\';out.write_text(json.dumps(result,indent=2)+\'\\n\')\n    print(json.dumps({\'all_passed\':True,\'seconds\':result[\'seconds\'],\'limit_pfaffian\':str(pf)}))\nif __name__==\'__main__\':main()\n'


def minus_two():
    started=time.monotonic();xx=s.symbols('h x1 x2 x3 x4');h,x,y,z,w=xx;t=s.Symbol('s');a,b=s.symbols('tau0 tau1')
    path=ROOT/'reports/CLASSIFICATION_NIL22_MINUS_TWO_WEIGHTED_SLICE_20261107.json'
    source=ROOT/'reports/CLASSIFICATION_NIL22_MINUS_TWO_CANONICAL_20261107.json'
    packet=json.loads(path.read_text());env={str(v):v for v in xx+(a,b)}
    P=s.zeros(5)
    for key,value in packet['full_curve_brackets'].items():
        i,j=map(int,key.split(','));P[i,j]=s.expand(s.sympify(value,locals=env).subs({a:-4*t,b:3*t},simultaneous=True))
    P-=P.T
    ps={ij:dict(s.Poly(P[ij],*xx).terms()) for ij in engine.PAIRS if P[ij]!=0}
    assert engine.jacobi(ps)==s.zeros(350,1)
    weights=packet['coordinate_weights'];pw=packet['P_weight']
    for i,j in engine.PAIRS:
        for mon,c in s.Poly(P[i,j],*(xx+(t,))).terms():
            if c:assert weights[i]+weights[j]-sum(u*v for u,v in zip(weights,mon[:5]))-pw==-mon[5]
    for key,value in packet['canonical_brackets'].items():
        i,j=map(int,key.split(','));assert s.expand(P[i,j].subs(t,0)-s.sympify(value,locals=env))==0
    section=0
    for k in range(5):
        i,j,l,m=[v for v in range(5) if v!=k]
        section+=(-1)**k*xx[k]*(P[i,j]*P[l,m]-P[i,l]*P[j,m]+P[i,m]*P[j,l])
    actual=s.factor(section.subs(t,1));content,factors=s.factor_list(actual)
    assert len(factors)==2 and all(e==1 for _,e in factors)
    by_degree={s.Poly(f,*xx).total_degree():f for f,_ in factors};assert set(by_degree)=={2,3}
    Q=by_degree[2];C=by_degree[3];Qmatrix=s.hessian(Q,xx)/2
    assert Qmatrix.det()!=0
    assert s.Poly(C,w).degree()==1 and s.diff(C,z)==0
    A=s.diff(C,w);g=C.subs(w,0);Amatrix=s.hessian(A,(h,x,y))/2
    assert Amatrix.det()!=0
    at={h:0,x:3,y:1};av=s.expand(A.subs(at));gv=s.expand(g.subs(at))
    assert av==0 and gv!=0
    old=json.loads((ROOT/'reports/CLASSIFICATION_NIL22_QUARTER_SNC_20261107.json').read_text())
    old_hash=hashlib.sha256(HISTORICAL_QUARTER_PRODUCER_BYTES.encode()).hexdigest()
    assert old_hash==old['input_sha256']['code/classification_nil22_quarter_snc.py']
    paths=[path,source,Path(__file__),Path(engine.__file__)]
    result={'task_id':'E_multiple_hyperplane','all_passed':True,
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'producer_bytes':Path(__file__).read_text(),
            'historical_quarter_producer_sha256':old_hash,'historical_quarter_producer_bytes':HISTORICAL_QUARTER_PRODUCER_BYTES,
            'curve_parameter_substitution':{'tau0':'-4*s','tau1':'3*s'},
            'full_curve_brackets':{str(i)+','+str(j):str(s.factor(P[i,j])) for i,j in engine.PAIRS if P[i,j]!=0},
            'all350_universal_projective_Jacobi_zero':True,'endpoint_equals_canonical_frozen_witness':True,
            'all_nonzero_curve_parameters_same_projective_orbit':True,'equivariant_coordinate_weights':weights,
            'parameter_one_pfaffian_section':str(actual),'factor_content':str(content),
            'smooth_quadric':str(Q),'smooth_quadric_matrix_determinant':str(Qmatrix.det()),
            'irreducible_cubic':str(C),'cubic_linear_variable':'x4','linear_variable_coefficient':str(A),
            'coefficient_quadric_rank':3,'coefficient_quadric_matrix_determinant':str(Amatrix.det()),
            'coprimality_point':{str(k):str(v) for k,v in at.items()},'coefficient_value':str(av),'cubic_remainder_value':str(gv),
            'irreducibility_reason':'The coefficient A has rank3 and is absolutely irreducible. A(point)=0,g(point)!=0 proves gcd(A,g)=1. Primitive degree-one polynomial A*w+g is irreducible over C by Gauss.',
            'coverage_dependency':'SMOOTH_QUADRIC_CUBIC_SUPPORT.md: every actual smooth quadric plus irreducible cubic divisor is covered by MPS.',
            'whole_nil22_A_minus_two_source_MPS_covered':True,'seconds':time.monotonic()-started}
    out=ROOT/'reports/CLASSIFICATION_NIL22_MINUS_TWO_MPS_20261107.json'
    if out.exists():raise FileExistsError(str(out))
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'Qdet':str(Qmatrix.det()),'coprimality_nonzero':str(gv)}))

def half_main():
    started=time.monotonic();h,x,y,z,w=xx=s.symbols('h x1 x2 x3 x4');B,tau=s.symbols('B tau')
    source=ROOT/'reports/CLASSIFICATION_NIL22_HALF_MAIN_CANONICAL_20261107.json'
    packet=json.loads(source.read_text());env={str(v):v for v in xx};env['n14']=B
    P=s.zeros(5)
    for key,value in packet['canonical_brackets'].items():
        i,j=map(int,key.split(','));P[i,j]=s.sympify(value,locals=env)
    P-=P.T;R=s.zeros(5)
    for ij,value in {(0,3):2*y*z/B,(0,4):2*y*w/B,(1,4):-(2*B+1)*x*z/B,(3,4):z*z}.items():R[ij]=value
    R-=R.T;curve=P+tau*R
    def sparse(M):return {ij:dict(s.Poly(s.expand(M[ij]),*xx).terms()) for ij in engine.PAIRS if M[ij]!=0}
    assert engine.jacobi(sparse(curve)).applyfunc(s.cancel)==s.zeros(350,1)
    # Canonical stabilizer makes all nonzero pencil members one orbit.
    canonical_weights=packet['weight']
    for M,expected in ((P,1),(R,-1)):
        for i,j in engine.PAIRS:
            for mon,c in s.Poly(M[i,j],*xx).terms():
                if c:assert canonical_weights[i]+canonical_weights[j]-sum(a*b for a,b in zip(canonical_weights,mon))==expected
    point=curve.subs(tau,-B/2).applyfunc(s.cancel)
    weights=[1,0,1,0,0];limit=s.zeros(5);terms=[]
    for i,j in engine.PAIRS:
        for mon,c in s.Poly(point[i,j],*xx).terms():
            if not c:continue
            weight=sum(a*b for a,b in zip(weights,mon))-weights[i]-weights[j]
            assert weight>=0
            terms.append(dict(pair=[i,j],monomial=mon,coefficient=str(c),weight=weight))
            if weight==0:limit[i,j]+=c*s.prod(a**b for a,b in zip(xx,mon))
    limit-=limit.T
    # New coordinate order X=x,U=h-y,V=h+y,Z=2w-z,W=2w+z.
    M=s.Matrix([[0,1,0,0,0],[1,0,-1,0,0],[1,0,1,0,0],[0,0,0,-1,2],[0,0,0,1,2]])
    inverse=M.inv()*s.Matrix(xx)
    toric=(M*limit*M.T).subs(dict(zip(xx,inverse)),simultaneous=True).applyfunc(s.expand)
    gamma={}
    for i,j in engine.PAIRS:
        value=s.cancel(toric[i,j]/(xx[i]*xx[j]))
        assert not any(value.has(v) for v in xx)
        if value:gamma[f'{i},{j}']=str(value)
    section=0
    for k in range(5):
        i,j,l,m=[v for v in range(5) if v!=k]
        section+=(-1)**k*xx[k]*(toric[i,j]*toric[l,m]-toric[i,l]*toric[j,m]+toric[i,m]*toric[j,l])
    section=s.factor(section);factor=s.cancel(section/s.prod(xx));assert factor!=0 and not any(factor.has(v) for v in xx)
    assert s.factor(factor/(2*B+1))!=0
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,
                source_identity_report=str(source.relative_to(ROOT)),source_identity_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),
                universal_endpoint_brackets={f'{i},{j}':str(P[i,j]) for i,j in engine.PAIRS if P[i,j]!=0},
                universal_pencil_direction={f'{i},{j}':str(R[i,j]) for i,j in engine.PAIRS if R[i,j]!=0},
                full350_universal_projective_Jacobi_zero=True,pencil_endpoint_weights=[1,-1],
                selected_nonzero_parameter='-B/2',degeneration_weights=weights,every_coefficient_term=terms,
                coordinate_matrix=[[str(v) for v in row] for row in M.tolist()],coordinate_determinant=str(M.det()),
                toric_log_coefficients=gamma,toric_pfaffian_section=str(section),
                canonical_open='B*(2*B+1)!=0',remaining_B_zero_covered_by_explicit_source_density=True,
                entire_half_main_source_MPS_covered=True,seconds=time.monotonic()-started)
    out=ROOT/'reports/CLASSIFICATION_NIL22_HALF_MAIN_MPS_20261107.json'
    if out.exists():raise FileExistsError(str(out))
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'],SNC_factor=str(factor))))


if __name__=='__main__':
    import sys
    half_main() if '--half-main' in sys.argv else minus_two() if '--minus-two' in sys.argv else main()
