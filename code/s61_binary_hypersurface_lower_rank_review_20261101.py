#!/usr/bin/env python3
"""Review all remaining outer-block defects on the binary K=0 hypersurface."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.time()
    primary_path=ROOT/'reports/S61_BINARY_HYPERSURFACE_LOWER_RANK_REDUCTION_20261101.json'
    cover_path=ROOT/'reports/S61_BINARY_RANK_DROP_COVER_20261101.json'
    rank6_path=ROOT/'reports/S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_REVIEW_20261101.json'
    c3_path=ROOT/'reports/S61_BINARY_MIXED_STRATUM_REMOVAL_REVIEW_20261101.json'
    primary=json.loads(primary_path.read_text());cover=json.loads(cover_path.read_text())
    for path in [primary_path,cover_path,rank6_path,c3_path]:
        assert json.loads(path.read_text())['all_passed']
    for path,digest in primary['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    u=s.symbols('u1:62');a,b,c,d,e,f=[u[i] for i in [3,13,34,35,48,49]]
    core=[a,b,c,d,e,f];L=5*b-e
    A,B,C=[s.Matrix(x).applyfunc(s.sympify) for x in cover['block_matrices']]
    for label,M in [('A',A),('C_linear',C)]:
        for record in primary['row_column_weight_certificates'][label]:
            weights=record['core'];rw=list(map(s.Rational,record['rows']));cw=list(map(s.Rational,record['columns']))
            for i in range(M.rows):
                for j in range(M.cols):
                    for mon,coef in s.Poly(M[i,j],*core).terms():
                        if coef:assert sum(k*w for k,w in zip(mon,weights))==rw[i]+cw[j]
    A=A.subs(f,-5*c);C=C.subs(f,-5*c)
    assert s.expand(C[[0,1,2,5,8],:].det(method='bareiss')-a*L**4/625)==0
    D=C[[0,1,2,6],:].col_join(C[5,:]-C[8,:])
    assert s.expand(D.det(method='bareiss')-(25*b-14*e)*L**4/3125)==0
    assert C[:,1].subs({a:0,b:14*e/25})==s.zeros(12,1)
    compatibility=list(map(s.sympify,cover['first70_t_compatibility']))
    assert s.expand(10*compatibility[8]-20*compatibility[28]-L**2)==0
    # Build a fresh ideal script from the reviewed block, rather than execute
    # the new producer's saved script unchanged.
    variables=s.symbols('a b c d e');AA=A.subs(dict(zip(core[:5],variables)))
    expr=lambda x:str(x).replace('**','^')
    matrix='matrix{'+','.join('{'+','.join(map(expr,row))+'}' for row in AA.tolist())+'}'
    program='R=QQ[a,b,c,d,e];\nM='+matrix+''';
I=minors(4,M);
P=ideal(c,d,5*b-2*e);
Q=ideal(5*b-e,5*a*d-6*c*e,20*c^2+3*d*e,50*a*c+9*e^2);
J=intersect(P,Q);
assert (numgens J==6);
assert (dim(R/I)==2 and dim(R/J)==2);
assert isSubset(I,J);
assert all(flatten entries gens J,x->(x^2)%I==0);
print("REVIEW_PASS: dimension2; six square containments; ideal containment");
'''
    env=os.environ.copy();env['DYLD_LIBRARY_PATH']=':'.join([str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib','/opt/homebrew/opt/givaro/lib'])
    with tempfile.TemporaryDirectory(dir=ROOT/'work',prefix='outer-review-') as td:
        path=Path(td)/'review.m2';path.write_text(program)
        run=subprocess.run([primary['m2_binary'],'--script',str(path)],env=env,capture_output=True,text=True,timeout=30)
    assert run.returncode==0 and 'REVIEW_PASS' in run.stdout,(run.stdout,run.stderr)
    cases=primary['remaining_general_beta_K_zero_cases']
    assert [(r['A_rank'],r['B_rank'],r['C_linear_rank'],r['leading_rank'],r['full_rank']) for r in cases]==[(4,8,5,33,34),(4,7,5,32,33)]
    assert all(16+r['A_rank']+r['B_rank']+r['C_linear_rank']==r['leading_rank'] and r['full_rank']==r['leading_rank']+1 for r in cases)
    paths=[Path(__file__),primary_path,cover_path,rank6_path,c3_path,
           ROOT/'research/notes/S61_BINARY_HYPERSURFACE_LOWER_RANK_REVIEW_20261101.md']
    return dict(id='S61_BINARY_HYPERSURFACE_LOWER_RANK_REVIEW_20261101',all_passed=True,accepted=True,field='QQ',
        independent_A_ideal_script=program,independent_A_ideal_output=run.stdout,
        C_two_minor_identities_rechecked=True,t_compatibility_L_squared_rechecked=True,
        stabilizer_weights_rechecked=True,remaining_general_beta_cases=cases,
        limitations=['Reuses the previously reviewed triangular block decomposition and primitive stabilizer representation.',
                    'Classifies remaining linear ranks only; higher obstructions are not exhausted.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','seconds']}))
