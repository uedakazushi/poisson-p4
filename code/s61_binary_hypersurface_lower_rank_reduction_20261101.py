#!/usr/bin/env python3
"""Complete outer-block defects on K=0 and remaining general-beta ranks33/34."""
import argparse, hashlib, json, os, subprocess, tempfile, time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]


def weight_certificate(matrix,core,weights):
    rw=s.symbols('r0:'+str(matrix.rows));cw=s.symbols('c0:'+str(matrix.cols))
    equations=[]
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            if matrix[i,j]:
                terms=s.Poly(matrix[i,j],*core).terms()
                observed={sum(mi*wi for mi,wi in zip(mon,weights)) for mon,_ in terms}
                assert len(observed)==1
                equations.append(rw[i]+cw[j]-observed.pop())
    sol=list(s.linsolve(equations,rw+cw))[0]
    free=set().union(*(v.free_symbols for v in sol))
    sol=[v.subs(dict.fromkeys(free,s.S.Zero)) for v in sol]
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            deriv=sum(w*v*s.diff(matrix[i,j],v) for w,v in zip(weights,core))
            assert s.expand(deriv-(sol[i]+sol[matrix.rows+j])*matrix[i,j])==0
    return dict(core=list(weights),rows=list(map(str,sol[:matrix.rows])),columns=list(map(str,sol[matrix.rows:])))


def verify(m2_binary):
    start=time.time()
    names=['S61_BINARY_RANK_DROP_COVER','S61_BINARY_DEEPER_REDUCTION',
        'S61_BINARY_DEEPER_RANK_SIX_EXCLUSION','S61_BINARY_MIXED_STRATUM_REMOVAL']
    paths=[ROOT/'reports'/(name+'_20261101.json') for name in names]
    inputs=[json.loads(p.read_text()) for p in paths]
    assert all(r['all_passed'] for r in inputs)
    cover=inputs[0]
    assert inputs[2]['all_K_zero_B_rank_at_most_six_avoided_by_general_beta']
    assert inputs[3]['theorem_conclusions']['entire_general_beta_C3_stratum_removed']
    u=s.symbols('u1:62');core=[u[i] for i in [3,13,34,35,48,49]]
    a,b,c,d,e,f=core;K,L=5*c+f,5*b-e
    A0,B0,C0=[s.Matrix(v).applyfunc(s.sympify) for v in cover['block_matrices']]
    A=A0.subs(f,-5*c);C=C0.subs(f,-5*c)
    minor_a=s.factor(C[[0,1,2,5,8],:].det(method='domain-ge'))
    comb=C[[0,1,2,6],:].col_join(C[5,:]-C[8,:])
    minor_b=s.factor(comb.det(method='domain-ge'))
    assert s.expand(minor_a-a*L**4/625)==0
    assert s.expand(minor_b-L**4*(25*b-14*e)/3125)==0
    assert C[:,1].subs({a:0,b:14*e/25})==s.zeros(12,1)
    # Thus C-rank<5 on L!=0 is exactly the previously removed C3 locus.
    weights=[(-3,-1,1,3,-1,1),(2,1,0,-1,1,0)]
    weighted={label:[weight_certificate(M,core,w) for w in weights]
        for label,M in [('A',A0),('C_linear',C0)]}
    # Exact small ideal certificate for the A-rank defect. Radical need not be guessed.
    aa,bb,cc,dd,ee=s.symbols('a b c d e')
    sym=A.subs(dict(zip([a,b,c,d,e],[aa,bb,cc,dd,ee])))
    ex=lambda v:str(v).replace('**','^')
    matrix='matrix{'+','.join('{'+','.join(ex(v) for v in row)+'}' for row in sym.tolist())+'}'
    program='R=QQ[a,b,c,d,e];\nA='+matrix+''';
I=minors(4,A);
P1=ideal(d,c,5*b-2*e);
P2=ideal(5*b-e,5*a*d-6*c*e,20*c^2+3*d*e,50*a*c+9*e^2);
J=intersect(P1,P2);
assert isSubset(I,J);
assert all(flatten entries gens J,g->g^2%I==0);
assert (dim(R/J)==2);
assert (dim(R/I)==2);
print("A_RANK_DEFECT_CORE_DIMENSION=2;ALL_SIX_GENERATOR_POWERS=2");
print("PASS");
'''
    env=os.environ.copy()
    env['DYLD_LIBRARY_PATH']=':'.join([str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib','/opt/homebrew/opt/givaro/lib'])
    with tempfile.TemporaryDirectory(prefix='binary-outer-blocks-',dir=ROOT/'work') as temp:
        source=Path(temp)/'certificate.m2';source.write_text(program)
        run=subprocess.run([str(m2_binary),'--script',str(source)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=30)
    assert run.returncode==0 and 'PASS' in run.stdout,(run.stdout,run.stderr)
    paths += [Path(__file__),ROOT/'research/notes/S61_BINARY_HYPERSURFACE_LOWER_RANK_REDUCTION_20261101.md']
    return dict(id='S61_BINARY_HYPERSURFACE_LOWER_RANK_REDUCTION_20261101',all_passed=True,field='QQ',
        A_rank_defect_core_dimension=2,A_rank_defect_beta_image_dimension_upper_bound=15,
        A_membership_certificate=program,A_membership_output=run.stdout,
        C_minor_a=str(minor_a),C_combination_minor_b=str(minor_b),
        C_rank_defect_on_L_nonzero=['u4','25*u14-14*u49'],
        C_rank_defect_entire_general_beta_stratum_removed=True,
        row_column_weight_certificates=weighted,
        remaining_general_beta_K_zero_cases=[
            dict(condition='Qstar != 0',A_rank=4,B_rank=8,C_linear_rank=5,leading_rank=33,full_rank=34),
            dict(condition='Qstar = 0',A_rank=4,B_rank=7,C_linear_rank=5,leading_rank=32,full_rank=33)],
        interpretation='After excluding proper beta images and the already saturated-away C3 stratum.',
        proof_note='research/notes/S61_BINARY_HYPERSURFACE_LOWER_RANK_REDUCTION_20261101.md',
        limitations=['This classifies the remaining linear ranks on K=0 over sufficiently general beta.',
            'It does not identify higher obstructions or prove free65 at every remaining point.',
            'Other rank-drop presentations and square/flag endpoints are separate.'],
        m2_binary=str(m2_binary),m2_binary_sha256=hashlib.sha256(m2_binary.read_bytes()).hexdigest(),
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--m2-binary',type=Path,default=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary'))
    args=parser.parse_args();report=verify(args.m2_binary)
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','seconds']}))
