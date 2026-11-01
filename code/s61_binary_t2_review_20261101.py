#!/usr/bin/env python3
"""Independent source/proof audit with explicitly labelled producer replays."""
import argparse, hashlib, json, time
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify(formal_replay=None, double_replay=None):
    started = time.time()
    names = ['CONNECTION', 'FORMAL', 'RELATIVE_DOUBLE_LINE', 'STRATUM_REMOVAL']
    paths = [ROOT / 'reports' / ('S61_BINARY_T2_' + n + '_20261101.json') for n in names]
    connection, formal, double, removal = [json.loads(p.read_text()) for p in paths]
    for r in [connection, formal, double, removal]:
        assert r['all_passed']
        for path, digest in r['input_sha256'].items():
            assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, path

    # Independently check the irreducible linear source, stabilizer invariance,
    # and the orbit dimension/dominance calculation.
    family_path = ROOT / 'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    rank_path = ROOT / 'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    family, rank = [json.loads(p.read_text()) for p in [family_path, rank_path]]
    u = list(map(s.Symbol, family['sample']))
    cons = s.Matrix([s.sympify(v) for v in connection['source_equations']])
    assert cons.jacobian(u).rank() == 5
    stab = s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    gb = s.groebner(list(cons), *u, domain=s.QQ)
    assert all(gb.reduce(s.expand(v))[1] == 0 for v in cons.jacobian(u) * stab)
    action = s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    sample = {s.Symbol(k): s.Rational(v) for k, v in connection['point']['sample'].items()}
    normal = (cons.jacobian(u) * action).subs(sample)
    assert normal == s.Matrix(connection['beta_normal_matrix']).applyfunc(s.Rational)
    assert normal[:, connection['beta_minor_columns']].det() == s.Rational(2248, 5)
    normal2 = s.Matrix(double['beta_normal_matrix']).applyfunc(s.Rational)
    assert normal2[:, double['beta_minor_columns']].det() == -s.Rational(2248, 25)
    assert 11 + 10 - 5 == connection['joint_image_dimension'] == 16

    if formal_replay is None:
        from s61_binary_t2_formal_20261101 import verify as run_formal
        current_formal = run_formal()
    else:
        current_formal = json.loads(Path(formal_replay).read_text())
    if double_replay is None:
        from s61_binary_t2_relative_double_line_20261101 import verify as run_double
        current_double = run_double()
    else:
        current_double = json.loads(Path(double_replay).read_text())
    for old, new in [(formal, current_formal), (double, current_double)]:
        assert new['all_passed']
        assert {k: v for k, v in old.items() if k != 'seconds'} == {
            k: v for k, v in new.items() if k != 'seconds'}
    assert formal['actual_binary_kernel_torsion_length_upper_bound'] == 1
    assert formal['difference_module_concentrated_degree'] == 3
    assert removal['relative_t_zero_boundary_dimension_upper_bound'] == 8 < 9
    assert removal['entire_general_beta_T2_connection_stratum_removed']
    assert removal['ambient_specified_point_t_nilpotence_asserted'] is False
    paths += [Path(__file__), family_path, rank_path,
              ROOT / 'research/notes/FINITE_BASE_TORSION_SPECIALIZATION_20261101.md',
              ROOT / 'research/notes/S61_BINARY_T2_FORMAL_20261101.md',
              ROOT / 'research/notes/S61_BINARY_T2_STRATUM_REMOVAL_20261101.md',
              ROOT / 'research/notes/S61_BINARY_T2_REVIEW_20261101.md']
    return dict(
        id='S61_BINARY_T2_REVIEW_20261101', all_passed=True, accepted=True,
        verification_kind='Independent source geometry and mathematical proof review; exact producer replays',
        source_stabilizer_invariance_verified=True, source_dominant_image_dimension=16,
        fixed_beta_formal_producer_replay_seconds=current_formal['seconds'],
        relative_double_line_producer_replay_seconds=current_double['seconds'],
        proof_checks=[
            'All lower-degree initial relations are present before the length-one cubic socle argument.',
            'Actual double-line quotient is strict; its kernel has at most one finite-length torsion element.',
            'Equal-rank specialization makes relative torsion finite over the nine-dimensional beta base.',
            'S/ann(T), rather than T[1/t], is the finite base module controlling the saturated support.',
            'The t=0 boundary has dimension at most eight and cannot contain the nine-dimensional section.',
            'Smooth equivariant base change and stable colons justify descent from the relative slice.',
            'The sixteen-dimensional irreducible stratum then has exceptional beta image dimension at most fifteen.'],
        limitations=[
            'The exact jet/initial-ideal calculations are replays, not an independent implementation.',
            'The relative argument does not assert ambient t^2=0 at the specified point.',
            'Other binary and flag strata and affine exhaustion remain unresolved.'],
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time() - started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--formal-replay', type=Path)
    parser.add_argument('--double-replay', type=Path)
    args = parser.parse_args()
    result = verify(args.formal_replay, args.double_replay)
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({k: result[k] for k in ['all_passed', 'seconds']}))
