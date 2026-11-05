#!/usr/bin/env python3
"""Pure consequences of accepted DVR results; no product coefficient evaluation."""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'AMBIENT_B18_NIL_CUBE_280_AND_ELL_DETECTION_20261105_v63'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main(out):
    started = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    stems = [
        'AMBIENT_B18_NIL_CUBE_OBSTRUCTION_REDUCTION_20261104_v63',
        'AMBIENT_B18_NIL_CUBE_OBSTRUCTION_REDUCTION_20261104_v63_REVIEW',
        'AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63',
        'AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW',
        'AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63',
        'AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63_REVIEW',
    ]
    paths = [ROOT / 'reports' / s / 'report.json' for s in stems]
    data = [json.loads(p.read_text()) for p in paths]
    assert all(d['all_passed'] for d in data)
    old, _, exact, _, bounds, _ = data
    assert old['canonical_beta_gamma_reduction']
    assert exact['actual_T_length'] == 11
    assert exact['actual_mu_D_K'] == 7
    assert exact['actual_mu_D_N'] == 14
    assert exact['actual_mu_D_torsion_K'] == 3
    assert bounds['actual_T2_length_interval'] == [17, 20]
    sym2p, sym3p = math.comb(8, 2), math.comb(9, 3)
    beta_count, gamma_count = 7 * sym2p, sym3p
    assert (beta_count, gamma_count, beta_count + gamma_count) == (196, 84, 280)
    evidence = {
        'unconditional': {
            'D': 'k[[ell]], ell=5(sA-rA), characteristic zero',
            'accepted_inputs': 'N^3 subset H=tors_D N=tors_D K; NK^2=0; P=N/K is D-free of rank7',
            'mu_D_K': 7,
            'mu_D_N': 14,
            'canonical_obstructions': [
                'beta: K tensor_D Sym^2(P) -> N^3, image Hbeta=KN^2',
                'gamma: Sym^3(P) -> N^3/Hbeta',
            ],
            'beta_generator_products_at_most': beta_count,
            'gamma_generator_products_at_most': gamma_count,
            'total_generator_products_at_most': beta_count + gamma_count,
            'criterion': 'N^3=0 iff beta=0 and gamma=0',
            'finite_reduction': 'N/ell*N is a 14-dimensional nonunital k-algebra',
            'nonzero_detection': 'Any nonzero triple product in N/ell*N proves N^3 nonzero',
        },
        'conditional': {
            'additional_hypothesis_not_established': 'actual length(T2)=17, equivalently ell*H=0 and length(H)=3',
            'module_identity': 'H intersect ell*N = ell*H, by the DVR structure theorem for the finite D-module N',
            'consequence_under_hypothesis': 'H injects into N/ell*N; thus N^3=0 iff Sym^3_k(N/ell*N) -> N/ell*N is zero',
            'at_most_triples_on_a_basis': math.comb(16, 3),
            'hypothesis_is_not_a_consequence_of_mu6_or_HF66': True,
        },
        'not_computed': [
            'actual seven K generators or seven lifts of a P basis',
            'actual fourteen N/ell*N basis lifts and their products',
            'beta, gamma, any triple product, actual T2 total length or H length',
        ],
        'original_nil_index_interval_unchanged': [3, 7],
        'no_new_source_jet_gauge_or_matrix_job': True,
    }
    ep = out / 'formal-280-products-and-conditional-ell-detection.json'
    ep.write_text(json.dumps(evidence, ensure_ascii=False, indent=2) + '\n')
    note = ROOT / 'research/notes' / (ID + '.md')
    report = {
        'id': ID, 'all_passed': True, 'formal_corollary_only': True,
        'maximum_D_generator_products': 280,
        'ell_fiber_zero_detection_requires_unproved_H0_annihilation': True,
        'actual_nil_cube_or_H0_or_T2_total_decided': False,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'code_sha256': sha(__file__), 'note_sha256': sha(note),
        'evidence_sha256': {ep.name: sha(ep)}, 'seconds': time.monotonic() - started,
    }
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output)
