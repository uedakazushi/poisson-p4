#!/usr/bin/env python3
"""Exact implicit-source certificates: ALL Qa in seventeen more actual MPS components.

The Jacobian used here belongs to the actual toric-period parameter source.
Its invertible minor integrates by the algebraic implicit function theorem;
it is not a Poisson tangent-vector integration assumption.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as sp

from q_mps_toric_arc_memberships_20261103_v63 import PAIRS4,PERMS,annihilator,period_permutations

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    source=ROOT/'code/data/q_mps_implicit_toric_arc_witnesses_20261103_v63.json'
    witnesses=json.loads(source.read_text())
    catalogue=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json'
    records=json.loads(catalogue.read_text())['diagrams']
    maps=period_permutations()
    alpha,beta,z,x,y,rho=sp.symbols('alpha beta z x y rho')
    variables=[alpha,beta,z,x,y];di=[0,1,z,x,y]
    gamma=sp.Matrix([alpha*(di[i]*di[j]**2-di[i]**2*di[j])
        +beta*(di[i]*di[j]**3-di[i]**3*di[j])
        +rho*(di[i]*di[j]**4-di[i]**4*di[j])
        +di[i]**2*di[j]**3-di[i]**3*di[j]**2 for i,j in PAIRS4])
    vandermonde=sp.prod(di[i]for i in range(1,5))*sp.prod(di[j]-di[i]for i in range(1,5)for j in range(i+1,5))
    assert sp.expand(gamma[0]*gamma[5]-gamma[1]*gamma[4]+gamma[2]*gamma[3]-rho*vandermonde)==0
    expected=[f'T{i:02d}'for i in range(4,16)]+[f'T{i:02d}'for i in range(18,22)]+['T28']
    assert [w['component']for w in witnesses]==expected
    evidence=[]
    for w in witnesses:
        rec=next(r for r in records if r['id']==w['component'])
        M=sp.Matrix(maps[PERMS.index(tuple(w['permutation']))].tolist())
        L=sp.Matrix(annihilator(rec).tolist());period=M*gamma;F=L*period
        value={alpha:sp.Rational(w['alpha']),beta:sp.Rational(w['beta']),
            z:sp.Integer(2),x:sp.Integer(3),y:sp.Integer(4),rho:sp.Rational(w['rho'])}
        assert F.subs(value)==sp.zeros(F.rows,1)
        assert value[rho] not in[0,1,-1]
        assert sp.Rational(w['a'])==8*value[rho]/(9*(value[rho]**2-1))
        basis=sp.Matrix(rec['period_basis']).applyfunc(sp.Rational).T
        prows=list(basis.T.rref()[1]);coordinates=basis[prows,:].inv()*period.subs(value)[prows,:]
        assert basis*coordinates==period.subs(value)
        selected=w['Jacobian_selected_columns'];assert len(selected)==F.rows
        derivative=F.jacobian(variables);J=derivative.subs(value)
        H=J[:,selected];det=H.det()
        assert det==sp.Rational(w['Jacobian_minor'])!=0
        # Exact first and second rho derivatives of the local algebraic curve,
        # with the unselected nonrho variables fixed.  Existence to all orders
        # follows from this invertible derivative, not from these finite jets.
        tangent=-H.inv()*F.diff(rho).subs(value)
        velocity=[sp.S.Zero]*5+[sp.S.One]
        for k,c in enumerate(selected):velocity[c]=tangent[k]
        allvariables=variables+[rho]
        second_rhs=sp.Matrix([sum(velocity[i]*velocity[j]*sp.diff(f,allvariables[i],allvariables[j]).subs(value)
            for i in range(6)for j in range(6))for f in F])
        acceleration=-H.inv()*second_rhs
        assert H*tangent+F.diff(rho).subs(value)==sp.zeros(F.rows,1)
        assert H*acceleration+second_rhs==sp.zeros(F.rows,1)
        # The period membership itself also integrates: F=0 is exactly the
        # condition that the six period coordinates belong to this linear span.
        assert basis.rank()==6-F.rows and L*basis==sp.zeros(F.rows,basis.cols)
        evidence.append(dict(component=w['component'],MPS_internal_component_id=True,
            projective_component_dimension=rec['projective_component_dimension'],
            original_period_basis=rec['period_basis'],permutation=w['permutation'],
            source_point={str(k):str(v)for k,v in value.items()},
            exact_period_equations=list(map(lambda f:str(sp.expand(f)),F)),
            period_coordinates_at_source_point=list(map(str,coordinates)),
            exact_nonrho_Jacobian=[[str(v)for v in row]for row in J.tolist()],
            selected_columns=selected,selected_determinant=str(det),
            rho_tangent=list(map(str,tangent)),rho_second_derivative=list(map(str,acceleration)),
            eigen_distinctness_product=str(vandermonde.subs(value)),
            bivector_Pfaffian=str((rho*vandermonde).subs(value)),
            all_period_equations_exact_zero=True,full_nonrho_Jacobian_rank=True,
            formal_rho_curve_exists_by_algebraic_implicit_function_theorem=True,
            all_finite_Qa_membership_by_dominant_rho_and_closedness=True))
    (args.output/'implicit-certificates.json').write_text(json.dumps(evidence,indent=2)+'\n')
    inputs=[Path(__file__),source,catalogue,
        ROOT/'code/q_mps_toric_arc_memberships_20261103_v63.py',
        ROOT/'code/q_mps_extended_toric_arcs_20261103_v63.py',
        ROOT/'research/notes/Q_MPS_IMPLICIT_TORIC_ARCS_20261103_v63.md',
        ROOT/'reports/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63/report.json']
    report=dict(id='Q_MPS_IMPLICIT_TORIC_ARCS_20261103_v63',all_passed=True,
        all_finite_complex_Qa_in_actual_MPS_components=expected,
        exact_QQ_source_points_and_nonrho_derivative_minors=True,
        first_and_second_rho_jets_checked_but_not_used_as_all_order_proof=True,
        actual_period_source_implicit_function_theorem=True,
        proof_uses_no_modular_specialization=True,includes_a_zero_by_closedness=True,
        component_or_local_branch_exhaustion_claimed=False,known_global_component_lower_bound_unchanged=57,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact QQ implicit source; ALL Qa in',','.join(expected),';',report['seconds'])


if __name__=='__main__':main()
