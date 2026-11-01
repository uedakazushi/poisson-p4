#!/usr/bin/env python3
"""Pin the accepted all-type reduced infinity assembly and corrected inputs."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEMS=[
    'S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101','S61_NONNILPOTENT_INFINITY_REVIEW_20261101',
    'S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101','S61_SMALL_NILPOTENT_INFINITY_REVIEW_20261101',
    'S61_NIL22_INFINITY_20261101','S61_NIL22_INFINITY_REVIEW_20261101',
    'S61_ZERO_MODULAR_INFINITY_REVIEW_20261101','S61_R0_GENERIC_EXHAUSTION_20261101',
    'S61_R0_EXHAUSTION_REVIEW_20261101','S61_R0_UNIFORM_VERTEX_20261101',
    'S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101',
    'S61_MINIMAL_VERTEX_EXHAUSTION_20261101','S61_MINIMAL_VERTEX_REVIEW_20261101',
    'S61_UNIVERSAL_WEIGHT6_REVIEW_20261031']


def run():
    paths=[Path(__file__),ROOT/'research/notes/S61_GENERIC_REDUCED_INFINITY_EXHAUSTION_20261101.md',
           ROOT/'research/notes/S61_GENERIC_REDUCED_INFINITY_REVIEW_20261101.md']
    reports={}
    for stem in STEMS:
        path=ROOT/'reports'/f'{stem}.json'
        r=json.loads(path.read_text());assert r['all_passed'] is True
        reports[stem]=r;paths.append(path)
        note=ROOT/'research/notes'/f'{stem}.md'
        assert note.exists();paths.append(note)
    boundary=reports['S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101']['G2_A_zero_boundary']
    assert boundary['generic_rank']==12 and len(boundary['rank_drop_points'])==4
    assert all(v['rank']==11 for v in boundary['rank_drop_points'])
    families=reports['S61_MINIMAL_VERTEX_EXHAUSTION_20261101']['families']
    expected={'M0':{('-1/3','g0'),('-2/3','g1')},'M1':{('-4/9','g0'),('-5/9','g1')}}
    for f in families:
        assert {(v['k'],v['parameter']) for v in f['baseline_free_lines']}==expected[f['label']]
    return dict(id='S61_GENERIC_REDUCED_INFINITY_REVIEW_20261101',all_passed=True,
        review_kind='Geometric assembly and frozen-dependency audit; individual exact certificates are not all rerun here.',
        corrected_boundary_included=True,corrected_baseline_free_lines=expected_to_json(expected),
        modular_cases=['zero','regular nilpotent (4)','nilpotent (2,2)',
                       'minimal nilpotent (2,1,1): open stratum and M0/M1','nonnilpotent'],
        accepted_conclusion='On a nonempty open of W16, every reduced point of unsaturated135 infinity is covered by R0 or the corrected M0/M1 charts.',
        limitations=['No assertion about saturation of every retained point.','No affine-finiteness or five-point-exhaustion assertion.',
                     'Archived complete fixed-trace classifications remain explicit dependencies.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


def expected_to_json(d):return {k:[list(v) for v in sorted(values)] for k,values in d.items()}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();r=run()
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'dependencies':len(STEMS),'output':str(args.output)}))
