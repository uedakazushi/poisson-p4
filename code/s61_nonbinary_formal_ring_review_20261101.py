#!/usr/bin/env python3
"""Independent coefficient/rank audit and formal-module proof review.

Rebuilds canonical leading derivatives from immutable S61 equations, checks the
full canonical beta locus, saved rational kernels and minors, and normal-square
obstruction spaces. The completed-ring argument is audited in the note.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from datetime import datetime,timezone
from pathlib import Path
import sympy as s
from s61_model import S61
ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.monotonic();model=S61();u=s.symbols('u1:62')
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    position={old:new for new,old in enumerate(two)};fixed=40
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra=json.loads(extra_path.read_text())['universal_new_equation']
    encoded=model.reduced['equations']+[extra]
    values=[u[i] if i in one else (-s.S.One if i==fixed else s.S.Zero) for i in range(61)]
    equations=[];derivative=s.zeros(135,45)
    for row,equation in enumerate(encoded):
        total=0
        for mon,numerator,denominator in equation:
            if sum(i in position for i in mon)!=2:continue
            coefficient=s.Rational(numerator,denominator)
            total+=coefficient*s.prod(values[i] for i in mon)
            for k,index in enumerate(mon):
                if index in position:
                    derivative[row,position[index]]+=coefficient*s.prod(values[j] for h,j in enumerate(mon) if h!=k)
        equations.append(s.expand(total))
    linear=s.Matrix(equations[:134]).jacobian([u[i] for i in one])
    expected=s.Matrix([u[3],u[4],u[11],u[12],u[49]+5*u[34],u[50]-5*u[25]]).jacobian([u[i] for i in one])
    assert linear.rref()[0][:6,:]==expected.rref()[0]
    assert not any(linear.rref()[0][6:,:])
    sub={u[3]:0,u[4]:0,u[11]:0,u[12]:0,u[49]:-5*u[34],u[50]:5*u[25]}
    assert all(s.expand(f.subs(sub))==0 for f in equations)
    derivative=derivative.subs(sub).applyfunc(s.expand)
    formal_path=ROOT/'reports/S61_NONBINARY_FORMAL_RING_20261101.json'
    formal=json.loads(formal_path.read_text());certificates={}
    for name,count,expected_rank in [('primary134',134,41),('augmented135',135,42)]:
        record=formal[name];kernel=s.Matrix(record['constant_kernel']);constant=derivative[:70,:]
        assert len(constant.rref()[1])==24
        assert constant*kernel==s.zeros(70,21)
        assert len(kernel.rref()[1])==21
        residual=derivative[70:count,:]*kernel
        assert residual==s.Matrix(record['residual_matrix'])
        relations=s.Matrix(record['independent_rational_kernel'])
        assert (residual*relations).applyfunc(s.cancel)==s.zeros(count-70,45-expected_rank)
        assert len(relations.rref()[1])==45-expected_rank
        sample={s.Symbol(k):s.Rational(v) for k,v in record['residual_sample'].items()}
        witness=record['residual_minor'];sample_matrix=residual.subs(sample)
        minor=sample_matrix.extract(witness['rows_zero_based'],witness['columns_zero_based']).det()
        assert minor==s.Rational(witness['determinant']) and minor!=0
        assert 24+witness['rank']==expected_rank
        certificates[name]={'generic_rank':expected_rank,'kernel_dimension':45-expected_rank,
                            'nonzero_residual_minor':str(minor)}
    local_path=ROOT/'reports/S61_NONBINARY_LOCAL_DIMENSION_20261101.json'
    local=json.loads(local_path.read_text())
    for name,n in [('primary134',3),('augmented135',2)]:
        record=local[name];tangent=s.Matrix(record['tangent_basis']);jac=s.Matrix(record['jacobian'])
        obstruction=s.Matrix(record['quadratic_obstruction_basis']);radical=s.Matrix(record['quadratic_common_radical'])
        assert jac*tangent==s.zeros(jac.rows,n) and tangent[-1,:]==s.zeros(1,n)
        assert len(tangent.rref()[1])==n and len(obstruction.rref()[1])==n*(n-1)//2
        assert radical.shape==(n,1) and radical!=s.zeros(n,1)
        # Explicitly check every quadratic is independent of the common-radical direction.
        monomials=[tuple(mon) for mon in record['quadratic_monomials']]
        for row in obstruction.tolist():
            coefficients=dict(zip(monomials,row,strict=True))
            hessian=s.Matrix(n,n,lambda i,j:coefficients[tuple(sorted((i,j)))]*(2 if i==j else 1))
            assert hessian*radical==s.zeros(n,1)
    paths=[Path(__file__),ROOT/'code/s61_model.py',extra_path,formal_path,local_path,
        ROOT/'reports/S61_R0_INFINITY_20261101.json',
        ROOT/'research/notes/S61_NONBINARY_FORMAL_RING_20261101.md',
        ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md',
        model.source/'symplectic_slice_reduced_equations.json']
    return {'id':'S61_NONBINARY_FORMAL_RING_REVIEW_20261101','all_passed':True,
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'canonical_beta_locus_exactly_six_linear_conditions':True,
        'leading_derivatives_rebuilt_from_immutable_equations':True,'rank_certificates':certificates,
        'normal_quadrics_equal_full_Sym2_of_normal_space':True,
        'formal_module_proof_accepted':True,
        'conclusions':{'134':'C[[u,e1,e2]]/(e1,e2)^2','135':'C[[u,e]]/(e^2)','homogenizer':0},
        'scope':'The specified nonbinary point and the open rank locus in its family; not all R0 or full-fiber exhaustion.',
        'proof_note':'research/notes/S61_NONBINARY_FORMAL_RING_REVIEW_20261101.md',
        'elapsed_seconds':round(time.monotonic()-started,3),
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as handle:json.dump(report,handle,indent=2);handle.write('\n')
    print(json.dumps({'all_passed':True,'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
