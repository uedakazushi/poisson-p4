#!/usr/bin/env python3
"""Independent coefficient and proof audit of the flag completed local rings."""
from __future__ import annotations
import argparse,hashlib,json,time
from datetime import datetime,timezone
from pathlib import Path
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.monotonic();model=S61();u=s.symbols('u1:62')
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    position={old:new for new,old in enumerate(two)};fixed=37
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra=json.loads(extra_path.read_text())['universal_new_equation']
    values=[u[i] if i in one else (s.S.One if i==fixed else s.S.Zero) for i in range(61)]
    equations=[];derivative=s.zeros(135,45)
    for row,equation in enumerate(model.reduced['equations']+[extra]):
        total=0
        for mon,numerator,denominator in equation:
            if sum(i in position for i in mon)!=2:continue
            coefficient=s.Rational(numerator,denominator)
            total+=coefficient*s.prod(values[i] for i in mon)
            for k,index in enumerate(mon):
                if index in position:
                    derivative[row,position[index]]+=coefficient*s.prod(values[j] for h,j in enumerate(mon) if h!=k)
        equations.append(s.expand(total))
    linear=s.Matrix(equations[:134]).jacobian([u[i] for i in one]);reduced,pivots=linear.rref()
    assert len(pivots)==1 and reduced[0,one.index(3)]==1 and sum(abs(v) for v in reduced[0,:])==1
    assert all(s.expand(f.subs(u[3],0))==0 for f in equations)
    derivative=derivative.subs(u[3],0).applyfunc(s.expand)
    main_path=ROOT/'reports/S61_FLAG_FORMAL_RING_20261101.json';main=json.loads(main_path.read_text())
    results=[]
    for record,expected in zip(main['generic_ranks'],[36,37],strict=True):
        count=record['equation_count'];kernel=s.Matrix(record['constant_kernel']);constant=derivative[:70,:]
        assert len(constant.rref()[1])==15 and constant*kernel==s.zeros(70,30)
        assert len(kernel.rref()[1])==30
        residual=derivative[70:count,:]*kernel
        assert residual==s.Matrix(record['residual'])
        relation=s.Matrix(record['rational_kernel'])
        assert len(relation.rref()[1])==45-expected
        assert (residual*relation).applyfunc(s.cancel)==s.zeros(count-70,45-expected)
        witness=record['residual_minor'];sample={s.Symbol(k):s.Rational(v) for k,v in record['sample'].items()}
        minor=residual.subs(sample).extract(witness['rows_zero_based'],witness['columns_zero_based']).det()
        assert minor==s.Rational(witness['determinant']) and minor!=0
        assert 15+witness['rank']==expected
        results.append({'equations':count,'generic_rank':expected,'residual_minor':str(minor)})
    local_path=ROOT/'reports/S61_FLAG_LOCAL_DIMENSION_20261101.json';old=json.loads(local_path.read_text())
    records=[(main['local134']['jacobian'],main['local134']['tangent'],main['local134']['quadratic_obstructions'],main['local134']['common_radical'],main['local134']['monomials']),
             (old['homogenized_jacobian'],old['tangent_basis'],old['quadratic_obstruction_basis'],old['quadratic_common_radical'],old['quadratic_monomials'])]
    for (jj,tt,oo,rr,mm),n in zip(records,[8,7],strict=True):
        jac=s.Matrix(jj);tangent=s.Matrix(tt);obstruction=s.Matrix(oo);radical=s.Matrix(rr)
        assert jac*tangent==s.zeros(jac.rows,n) and tangent[-1,:]==s.zeros(1,n)
        assert len(tangent.rref()[1])==n and len(radical.rref()[1])==3
        assert len(obstruction.rref()[1])==(n-3)*(n-2)//2
        monomials=[tuple(mon) for mon in mm]
        for row in obstruction.tolist():
            coefficient=dict(zip(monomials,row,strict=True))
            hessian=s.Matrix(n,n,lambda i,j:coefficient[tuple(sorted((i,j)))]*(2 if i==j else 1))
            assert hessian*radical==s.zeros(n,3)
    z,uu,beta=primitive_beta(model);a=s.Symbol('a');v=s.Matrix([0,-1,0,0]);w=s.Matrix([0,0,a,-1])
    incidence=s.expand((v.T*beta.subs(dict(zip(z,v)))*w)[0])
    sample=dict(zip(uu,range(2,18)));sample[u[3]]=0
    transverse_derivative=s.diff(incidence,a).subs(sample)
    assert transverse_derivative in (4,-4)
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',extra_path,main_path,local_path,
        ROOT/'research/notes/S61_FLAG_FORMAL_RING_20261101.md',ROOT/'research/notes/S61_NONBINARY_FORMAL_RING_20261101.md',
        ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md',model.source/'symplectic_slice_reduced_equations.json']
    return {'id':'S61_FLAG_FORMAL_RING_REVIEW_20261101','all_passed':True,
        'generated_utc':datetime.now(timezone.utc).isoformat(),'canonical_beta_condition_exactly_u4_zero':True,
        'independent_generic_ranks':results,'normal_quadrics_span_full_Sym2':True,
        'flag_incidence_transverse_derivative':str(transverse_derivative),
        'formal_module_proof_over_three_variable_base_accepted':True,
        'complete_rings':{'134':'C[[u1,u2,u3]] plus a free rank5 square-zero ideal','135':'C[[u1,u2,u3]] plus a free rank4 square-zero ideal'},
        'homogenizer_zero_in_both':True,
        'scope':'The saved flag and open incidence rank locus; no all-infinity or full-affine-fiber claim.',
        'proof_note':'research/notes/S61_FLAG_FORMAL_RING_REVIEW_20261101.md','elapsed_seconds':round(time.monotonic()-started,3),
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as handle:json.dump(report,handle,indent=2);handle.write('\n')
    print(json.dumps({'all_passed':True,'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
