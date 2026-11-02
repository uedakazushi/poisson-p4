#!/usr/bin/env python3
"""Identify gr of the actual flat quotient and bound the actual kernel over QQ[[v0,v1]]."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_zero12_initial_module_20261101 import module_program

ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/S61_FLAG_ZERO12_ACTUAL_KERNEL_FINITE_20261102.md'


def verify():
    started=time.monotonic()
    primary=ROOT/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json'
    squarezero=ROOT/'reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json'
    a=json.loads(primary.read_text());d=json.loads(squarezero.read_text())
    assert a['all_passed'] and d['all_passed']
    assert a['module_certificate']['expected_primary_ideal']=='(v4,v7)+(v2,v3+2*v6,v5,v8)^2'
    assert a['module_certificate']['embedded_free_base_variables']==[0,1]
    point=s.Matrix(a['point61']).applyfunc(s.Rational)
    model=S61();local=calculate(model,True,point)
    T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    N=s.Matrix(d['squarezero_normal_directions']).applyfunc(s.sympify)
    N=N.subs({s.Symbol(f'u{i+1}'):point[i] for i in range(61)})
    NL=N[local['chart_columns'],:].col_join(s.zeros(1,4))
    rows=list(T.T.to_DM().rref()[1])
    X=T[rows,:].inv()*NL[rows,:]
    assert T*X==NL
    C=s.Matrix(a['geometry']['actual_flag_tangent']).applyfunc(s.Rational)
    L=s.zeros(2,9);L[0,4]=1;L[1,7]=1
    E=s.zeros(4,9);E[0,2]=1;E[1,3]=1;E[1,6]=2;E[2,5]=1;E[3,8]=1
    assert L*C==s.zeros(2,3) and L*X==s.zeros(2,4)
    assert E*C==s.zeros(4,3) and (E*X).det()==s.Rational(4,63)
    assert C.row_join(X).rank()==7
    # The primary module proof is replayed; no equality with full gr(S) is used.
    program,certificate=module_program(a['initial_relations'])
    env=os.environ.copy()
    env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
    with tempfile.TemporaryDirectory(prefix='zero12-actual-kernel-') as temp:
        source=Path(temp)/'module.m2';source.write_text(program)
        run=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(source)],
                           cwd=ROOT,env=env,capture_output=True,text=True,timeout=30)
        assert run.returncode==0 and 'EMBEDDED FREE111 PASS' in run.stdout,run.stdout+run.stderr
    enc=lambda M:[[str(c) for c in row] for row in M.tolist()]
    paths=[Path(__file__),primary,squarezero,ROOT/NOTE,ROOT/'code/s61_model.py',
           ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_flag_zero12_initial_module_20261101.py']
    return dict(id='S61_FLAG_ZERO12_ACTUAL_KERNEL_FINITE_20261102',all_passed=True,field='QQ',
                point61=list(map(str,point)),actual_C_tangent=enc(C),squarezero_normal_tangent=enc(X),
                vanishing_linear_forms=enc(L),nilpotent_linear_forms=enc(E),normal_determinant='4/63',
                graded_actual_quotient_ideal='(v4,v7)+(v2,v3+2*v6,v5,v8)^2',
                actual_kernel_finite_base_variables=[0,1],actual_kernel_generator_bound=111,
                associated_graded_kernel_is_quotient_of_certified_free111_module=True,
                associated_graded_kernel_equality_claimed=False,actual_kernel_freeness_claimed=False,
                ambient_t_nilpotence_claimed=False,
                scope='Specified fixed-beta point; actual kernel M, not just initial E, is finite over QQ[[v0,v1]].',
                module_certificate=certificate,module_log=run.stdout+run.stderr,
                proof_note=NOTE,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as out:
        json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
