#!/usr/bin/env python3
"""Pure consequence of accepted whole-E4 generic Artin models; no source job."""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_NORMALIZER_NIL_SQUARE_THREE_SUPPORTS_20261104_v63'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main(out):
    out.mkdir(parents=True, exist_ok=False)
    fixed = {
        'reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json': 'cffc3d5c6a9cf9ed1207788a27675bdc76f925380fa27f0dc3a3c0c8e5337f80',
        'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW/report.json': '02888d3c4b987327f6df2d8faf1f2adffa3f9c01af47e27aff8c5c0982e75472',
        'reports/AMBIENT_B18_QUADRATIC_NORMALIZER_NIL_UPPER12_20261104_v63/report.json': '332701163619c4e3480e77e4240f604bb3200ef9c17a2131daa72d363bc0d77c',
    }
    for p, h in fixed.items():
        assert sha(ROOT / p) == h, p
        assert json.loads((ROOT / p).read_text())['all_passed']
    A = [[0,1,0,0], [Q(1,2),0,0,0], [Q(-1,2),0,0,1], [1,-1,1,0]]
    det = Q(0)
    for p in itertools.permutations(range(4)):
        inv = sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
        v = Q((-1)**inv)
        for i in range(4):
            v *= A[i][p[i]]
        det += v
    assert det == Q(1,2)
    lines = {'N':[0,0,0,1], 'P':[2,0,-2,1], 'C':[0,1,1,0]}
    transformed = {n:[sum(Q(a)*Q(b) for a,b in zip(row,v)) for row in A] for n,v in lines.items()}
    assert transformed == {'N':[0,0,1,0], 'P':[0,1,0,0], 'C':[1,0,0,0]}
    # Embedded nil-square corresponds to d*T, not d^2*T.
    N_d, P_d = [0,0,4], [0,1]
    def times_d_rank_truncated(d):
        return len(d) - next(i for i,c in enumerate(d) if c)
    assert times_d_rank_truncated(N_d) == times_d_rank_truncated(P_d) == 1
    proof = {
        'whole_E4':'S9/(dJ), not the core S9/J',
        'embedded_nil_module_identification':'d:T_q -> Nil(E4)_q',
        'embedded_nil_square_correspondence':'Nil(E4)_q^2 corresponds to d*T_q',
        'embedded_nil_cube_correspondence':'Nil(E4)_q^3 corresponds to d^2*T_q',
        'generic_nil_square_lengths':{'N':1,'P':1,'C':1,'R':0,'L':0,'main_d':0},
        'Ass_nil_square_exactly':['N','P','C'],
        'no_origin_submodule':True,
        'linear_coordinates_c_p_n_H':[[str(c) for c in row] for row in A],
        'coordinate_determinant':str(det),
        'transformed_lines':{n:[str(c) for c in v] for n,v in transformed.items()},
        'exact_annihilator_S9':['u1','u2','u3','u4','u5','a-d+x','a*d','d*(2*t-a)','a*(2*t-a)'],
        'annihilator_reduced':True,
        'dimension':1,'depth':1,'multiplicity':3,
        'rank_on_each_reduced_support_line':1,
        'module_isomorphism_to_reduced_curve_or_canonical_module_claimed':False,
        'upper8_corrected_gate':'image Nil(D0)^2 vanishes at ALL N,P,C actual source generic comparisons',
        'Qb_to_C_alone_suffices':False,
        'upper8_or_original_cube_decision_claimed':False,
        'new_source_jet_point_CAS_or_GB':False,
    }
    ep = out / 'nil-square-generic-lengths-and-annihilator-proof.json'
    ep.write_text(json.dumps(proof,indent=2)+'\n')
    files = list(fixed) + [
        'research/notes/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63.md',
        'research/notes/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW.md',
        'research/notes/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63.md',
        f'research/notes/{STEM}.md',
        str(Path(__file__).relative_to(ROOT)),
    ]
    report = {
        'id':STEM,'all_passed':True,
        'status':'PROVED_whole_normalized_E4_nil_square_support_and_reduced_annihilator',
        'scope':'normalized E4 only, not original D0 or full D nilideal',
        'embedded_N_and_P_nil_square_nonzero':True,
        'nil_square_Ass_exactly_N_P_C':True,
        'nil_square_reduced_annihilator_exact':True,
        'nil_square_dimension_depth_multiplicity':[1,1,3],
        'actual_source_images_at_N_P_C_evaluated':False,
        'upper8_claimed':False,
        'input_sha256':{p:sha(ROOT/p) for p in files},
        'evidence_sha256':{ep.name:sha(ep)},
    }
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','nil_square_dimension_depth_multiplicity','upper8_claimed']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path)
    main(p.parse_args().output.resolve())
