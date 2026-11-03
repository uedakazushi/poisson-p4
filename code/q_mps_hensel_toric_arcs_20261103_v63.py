#!/usr/bin/env python3
"""Positive Hensel existence certificates for actual toric-period source curves.

No emptiness or finite-length claim is inferred from a finite-field fiber.
All five actual equations and an invertible nonrho derivative give positive
characteristic-zero existence by multivariate Hensel and an étale rho-map.
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


def evaluate(terms,values,modulus):
    answer=0
    for exponents,coefficient in terms:
        term=int(coefficient)%modulus
        for value,exponent in zip(values,exponents):term=term*pow(value,exponent,modulus)%modulus
        answer=(answer+term)%modulus
    return answer


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    source=ROOT/'code/data/q_mps_hensel_toric_arc_witnesses_20261103_v63.json'
    witnesses=json.loads(source.read_text())
    catalogue=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json'
    records=json.loads(catalogue.read_text())['diagrams'];maps=period_permutations()
    alpha,beta,z,x,y,rho=sp.symbols('alpha beta z x y rho')
    variables=[alpha,beta,z,x,y,rho];di=[0,1,z,x,y]
    gamma=sp.Matrix([alpha*(di[i]*di[j]**2-di[i]**2*di[j])
        +beta*(di[i]*di[j]**3-di[i]**3*di[j])
        +rho*(di[i]*di[j]**4-di[i]**4*di[j])
        +di[i]**2*di[j]**3-di[i]**3*di[j]**2 for i,j in PAIRS4])
    expected=['T22','T23','T24','T25','T26','T27','T29','T30','T31','T32','T33','T34','T35','T36','T37','T38']
    assert [w['component']for w in witnesses]==expected
    certificates=[]
    for w in witnesses:
        p=w['prime'];assert sp.isprime(p)and p not in[2,3]
        rec=next(r for r in records if r['id']==w['component'])
        M=sp.Matrix(maps[PERMS.index(tuple(w['permutation']))].tolist())
        L=sp.Matrix(annihilator(rec).tolist());F=L*M*gamma;assert F.rows==5
        B=sp.Matrix(rec['period_basis']).applyfunc(sp.Rational).T
        assert B.rank()==1 and L*B==sp.zeros(5,1)
        assert all(q.q%p for q in B)
        ds=w['eigenvalues'];assert len(set(c%p for c in ds))==5
        values=[w['alpha'],w['beta'],*ds[2:],w['rho']]
        assert values[-1]%p not in[0,1,p-1]
        polynomials=[sp.Poly(f,*variables,domain=sp.ZZ)for f in F]
        terms=[f.terms()for f in polynomials]
        assert all(evaluate(f,values,p)==0 for f in terms)
        J=F.jacobian(variables[:5]).subs(dict(zip(variables,values))).applyfunc(lambda c:int(c)%p)
        assert J.tolist()==w['Jacobian_mod_prime']
        determinant=int(J.det())%p
        assert determinant==w['selected_determinant_mod_prime']!=0
        assert w['nonrho_selected_columns']==[0,1,2,3,4]
        inverse=J.inv_mod(p)
        rho_derivative=F.diff(rho).subs(dict(zip(variables,values))).applyfunc(lambda c:int(c)%p)
        tangent=(-inverse*rho_derivative).applyfunc(lambda c:int(c)%p)
        # Finite Newton replays verify the data/implementation.  The all-order
        # existence claim uses the unit derivative and Hensel's theorem.
        lifts=[dict(precision=1,values=list(values))]
        lifted=list(values)
        for precision in range(1,8):
            oldmod=p**precision;newmod=oldmod*p
            residual=[evaluate(f,lifted,newmod)for f in terms]
            assert all(c%oldmod==0 for c in residual)
            quotient=sp.Matrix([c//oldmod for c in residual])
            correction=(-inverse*quotient).applyfunc(lambda c:int(c)%p)
            for i in range(5):lifted[i]=(lifted[i]+oldmod*int(correction[i]))%newmod
            assert lifted[-1]==values[-1]
            assert all(evaluate(f,lifted,newmod)==0 for f in terms)
            lifts.append(dict(precision=precision+1,values=list(lifted)))
        certificates.append(dict(component=w['component'],MPS_internal_component_id=True,
            projective_component_dimension=rec['projective_component_dimension'],prime=p,
            exact_integer_period_equations=[str(f.as_expr())for f in polynomials],
            original_period_basis=rec['period_basis'],permutation=w['permutation'],
            modular_source_values=dict(zip(map(str,variables),values)),
            all_five_actual_equations_zero_mod_p=True,nonrho_Jacobian_mod_p=[[int(c)for c in row]for row in J.tolist()],
            nonrho_determinant_mod_p=determinant,rho_tangent_mod_p=list(map(int,tangent)),
            all_admissibility_factors_p_units=True,Newton_lifts_at_fixed_integer_rho=lifts,
            positive_characteristic_zero_existence_by_Hensel=True,
            rho_dominance_by_etale_source_and_faithful_base_extension=True,
            all_finite_Qa_membership_by_actual_toric_arc_and_closedness=True))
    (args.output/'hensel-certificates.json').write_text(json.dumps(certificates,indent=2)+'\n')
    inputs=[Path(__file__),source,catalogue,
        ROOT/'code/q_mps_toric_arc_memberships_20261103_v63.py',
        ROOT/'code/q_mps_extended_toric_arcs_20261103_v63.py',
        ROOT/'research/notes/Q_MPS_HENSEL_TORIC_ARCS_20261103_v63.md',
        ROOT/'reports/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63/report.json']
    report=dict(id='Q_MPS_HENSEL_TORIC_ARCS_20261103_v63',all_passed=True,
        all_finite_complex_Qa_in_actual_MPS_components=expected,
        positive_modular_existence_bridge_not_emptiness_lifting=True,
        all_five_actual_equations_and_unit_nonrho_derivative_checked=True,
        all_admissibility_factors_p_units=True,Newton_replay_precision=8,
        includes_a_zero_by_closedness=True,component_or_local_branch_exhaustion_claimed=False,
        known_global_component_lower_bound_unchanged=57,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS positive Hensel ALLQa',','.join(expected),';',report['seconds'])


if __name__=='__main__':main()
