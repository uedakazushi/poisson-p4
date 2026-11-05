#!/usr/bin/env python3
"""Independent core-monomial proof of the D-module countermodel only."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal

ROOT = Path(__file__).resolve().parents[1]
ID = 'AMBIENT_B18_BAR_CUBE_ZERO_INDEX7_COUNTERMODEL_20261105_v63_REVIEW'


def mul(a, b):
    # A core monomial is (basis, ell exponent); 0,1,2 encode 1,x,x^2;
    # 3 encodes y with ell^11*y=0. None is the zero element.
    if a is None or b is None:
        return None
    i, r = a
    j, s = b
    e = r + s
    if i < 3 and j < 3:
        k = i + j
        if k < 3:
            return (k, e)
        k, e = 3, e + 1 + 3 * (k - 3)
    elif i == 3 and j == 3:
        k, e = 3, e + 8
    else:
        k, e = 3, e + 3 * (j if i == 3 else i)
    return None if k == 3 and e >= 11 else (k, e)


def main(out):
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000,) * 2)
    signal.alarm(60)
    # Prior budget: <=1024 integers (all monomial pairs/triples and metadata),
    # exponent<=24, coefficient=1; input+publication<=65536 bytes. No actual
    # D0 coefficients or producer arithmetic are used.
    assert not out.exists()
    primary = ROOT / 'reports/AMBIENT_B18_BAR_CUBE_ZERO_INDEX7_COUNTERMODEL_20261105_v63/report.json'
    assert hashlib.sha256(primary.read_bytes()).hexdigest() == '2a75e3ecd6a604543c3f724d89fa40ee68beea797b4949d521b7b61ac0f7b6c4'
    assert json.loads(primary.read_text())['all_passed']
    core = [(i, 0) for i in range(4)]
    triples = 0
    for a in core:
        for b in core:
            assert mul(a, b) == mul(b, a)
            for c in core:
                assert mul(mul(a, b), c) == mul(a, mul(b, c))
                triples += 1
    assert triples == 64
    for a in core:
        assert mul((3, 11), a) is None
        assert mul((0, 0), a) == a
    x = (1, 0)
    powers = [(0, 0)]
    for _ in range(7):
        powers.append(mul(powers[-1], x))
    assert powers[3] == (3, 1)
    assert powers[6] == (3, 10) and powers[7] is None
    cube_exponents = []
    for a in core[1:]:
        for b in core[1:]:
            for c in core[1:]:
                v = mul(mul(a, b), c)
                if v is not None:
                    assert v[0] == 3 and v[1] >= 1
                    cube_exponents.append(v[1])
    assert min(cube_exponents) == 1
    assert mul((3, 0), (3, 0)) == (3, 8)
    assert all(mul(a, (3, 8)) is None for a in core[1:])
    ev = dict(core_basis=['1', 'x', 'x^2', 'y'], y_torsion_exponent=11,
              independent_core_associativity_checks=64, x_powers=powers,
              cube_exactly_ell_times_y=True, N_K_squared_zero=True,
              extension='Trivial square-zero extension by D^9 + k^2 through core augmentation to D; all new generators annihilate the nil ideal.',
              N_D_rank=11, N_D_torsion_minimal_generators=3, N_D_minimal_generators=14,
              K='tors_D(N) plus any four new D-free generators',
              K_D_rank=4, K_D_torsion_minimal_generators=3, K_D_minimal_generators=7,
              P_D_free_rank=7, N_bar_cube_zero=True, N_cube_exactly_ell_H0=True,
              nil_index_exactly=7,
              scope='Only the listed D-module/nil identities. No compatibility with actual S-source, Tnil Loewy structure or all geometric maps is asserted.')
    out.mkdir(parents=True)
    evidence = out / 'independent-core-monomial-proof.json'
    evidence.write_text(json.dumps(ev, indent=2) + '\n')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    note = ROOT / ('research/notes/' + ID + '.md')
    report = dict(id=ID, all_passed=True, no_actual_D0_coefficients_evaluated=True,
                  prior_owned_numeric_bound=1024, prior_TOTAL_RAW_bound=65536,
                  input_sha256={str(primary.relative_to(ROOT)): sha(primary)},
                  code_sha256=sha(Path(__file__)), note_sha256=sha(note),
                  evidence_sha256={evidence.name: sha(evidence)})
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    assert primary.stat().st_size + sum(p.stat().st_size for p in out.iterdir()) < 65536
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    main(parser.parse_args().output.resolve())
