#!/usr/bin/env python3
"""Exact orbit normal form of the entire open four-parameter nil22 source."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();sourcepath=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_COMPONENTS_20261107.json');source=json.loads(sourcepath.read_text());r=next(c for c in source['cases'] if c['source']=='nil22_nonzero');xx=s.symbols('x0:5');n0,n2,n3,n6=s.symbols('n0 n2 n3 n6');lam=s.Symbol('lam',nonzero=True);P=s.zeros(5)
for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(r['universal_rational_family'].get(str(ij),'0'))/125;P[ij[::-1]]=-P[ij]
P0=P.subs({n0:1,n2:0,n3:0,n6:0});k=3*n3/5
G0=s.diag(1,1,lam,lam,1);Gz=s.eye(5);Gz[4,0]=3*(n2-k*n6)/10;Gz[4,2]=3*n6/5;G3=s.eye(5);G3[2,0]=-k;G3[3,1]=-k;G=G3*Gz*G0
push=(G*P0.subs(dict(zip(xx,G.inv()*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.expand)
assert (push-P.subs(n0,lam**-2)).applyfunc(s.expand)==s.zeros(5)
weights=[2,0,1,-1,-2];pweight=-2
for ij in itertools.combinations(range(5),2):
    for mon,coef in s.Poly(P0[ij],*xx).terms():
        if coef:assert weights[ij[0]]+weights[ij[1]]-sum(a*b for a,b in zip(weights,mon))==pweight
record={'task_id':'E_reduced_nonnormal','entire_four_parameter_open_source_is_one_orbit':True,'open_condition':'n0!=0','surjective_parameterization_n0':'lam**(-2), lam!=0','composition_order':'G3*Gz*diag(1,1,lam,lam,1)','coordinate_change':[[str(v) for v in row] for row in G.tolist()],'full_source_equals_pushforward_exact':True,'canonical_brackets':{','.join(map(str,ij)):str(P0[ij]).replace('x0','h') for ij in itertools.combinations(range(5),2) if P0[ij]},'branch':'q4_nil22_nonzero','weight':weights,'P_weight':pweight,'expected_full_D_rank':100,'expected_actual_group_rank':47,'source_report':str(sourcepath),'source_sha256':hashlib.sha256(sourcepath.read_bytes()).hexdigest(),'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CANONICAL_20261107.json')
with out.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps({k:record[k] for k in ['entire_four_parameter_open_source_is_one_orbit','weight','P_weight','seconds']},indent=2))
