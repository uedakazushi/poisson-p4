#!/usr/bin/env python3
"""Record the proved one-step DVR reduction, with actual valuations unresolved."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_REDUCED_UNION_DVR_CONTACT_TARGET_20261103_v63'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    paths = [
        Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/report.json',
        ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/report.json',
        ROOT / 'reports/AMBIENT_B18_NILRADICAL_LINEAR_VISIBILITY_20261103_v63_REVIEW.json',
    ]
    for p in paths[2:]:
        assert json.loads(p.read_text())['all_passed']
    # Both alternatives really occur in a reduced union with a DVR branch:
    # S=C[[x,y,z,w]], I=(z,w,y²-x³), K=(y,z-x^d,w).
    # The curve parameter is x=t; nonlinear restriction is -t³.
    toy = []
    for d in [2, 4]:
        k = min(d, 3)
        h = min(d + 1, 3)
        assert h - k in [0, 1]
        toy.append(dict(d=d, k=k, h=h, union_embedding_dimension=2+h-k))
    assert [q['union_embedding_dimension'] for q in toy] == [3, 2]
    result = dict(
        id=STEM, all_passed=True,
        proof_status='Proved DVR ideal-valuation reduction, actual k and h not computed.',
        actual_common_transverse_S_dimension=20, actual_N_embedding_dimension=18,
        actual_R_transverse_branch='C[[t]]', actual_R_tangent_contained_in_N_tangent=True,
        J='image(I_N)=(t^k), k>=2', H='I_N intersect m_S^2',
        J_H='image(H)=(t^h)', h_minus_k_in=[0, 1],
        reduced_union_transverse_embedding_dimension='18+(h-k)',
        actual_nilradical_linear_image_dimension='2-(h-k)',
        actual_k_h_or_exact_reduced_embedding_dimension_claimed=False,
        exact_depth_contact_length_full_D_ring_or_uniform_jet_bound_claimed=False,
        two_toy_alternatives=toy,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    assert not out.exists()
    out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'all_passed': True, 'actual_valuations_computed': False}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())
