#!/usr/bin/env python3
"""Length65 local opens on both codimension-two binary rank-drop charts."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_binary_rank_drop_point_engine_20261101 import verify_point
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time();cover_path=ROOT/'reports/S61_BINARY_RANK_DROP_COVER_20261101.json'
    cover=json.loads(cover_path.read_text());assert cover['all_passed'];cases=[]
    for source in cover['orbit_strata']:
        if source['name']=='hypersurface':continue
        linear=source['name']=='linear';qr=13 if linear else 14;params=[0,1,4] if linear else [0,1,3]
        result=verify_point(source['sample'],qr,params)
        # The source constraint is smooth at the stated sample, and its group
        # directions map onto all four normal directions to the base slice.
        normal=s.Matrix(source['beta_normal_action']).applyfunc(s.Rational);cols=source['beta_minor_columns']
        assert normal.shape==(4,10) and normal[:,cols].det()==s.Rational(source['beta_minor'])!=0
        result.update(name=source['name'],base_generators=source['generators'],base_dimension=12,joint_image_dimension=17,
            beta_dominance_minor=source['beta_minor'],relative_known_support_dimension=15,
            relative_free_rank=65,ambient_nilpotence_on_nonempty_open=True,whole_generic_beta_curve_removed=False)
        cases.append(result)
    assert len(cases)==2
    paths=[Path(__file__),ROOT/'code/s61_binary_rank_drop_point_engine_20261101.py',cover_path,
        ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
        ROOT/'research/notes/S61_BINARY_UNIFORM_CUBIC_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_BINARY_CODIM_TWO_FORMAL_20261101',all_passed=True,prime=32003,cases=cases,
        proof_note='research/notes/S61_BINARY_CODIM_TWO_FORMAL_20261101.md',
        limitations=['Specified points and resulting nonempty opens only; no removal of all generic-beta curves.',
                    'The complement of each17D open may contain16D strata dominating W16.',
                    'No full binary or affine-fiber exhaustion.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,cases=[{k:c[k] for k in ['name','exact_transverse_quadratic_rank','quintic_initial_minor_mod_prime','relation_minor_mod_prime','augmented_homogenizer_minor_mod_prime','transverse_hilbert_function']} for c in r['cases']],seconds=r['seconds'])))

if __name__=='__main__':main()
