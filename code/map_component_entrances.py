#!/usr/bin/env python3
"""Map current component records using recorded actual divisors, never ID names.

This is a navigation table, not an independent component or source verification.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    components=json.loads((ROOT/'research/COMPONENTS.json').read_text())['components']
    mps={r['component_id']:r for r in json.loads((ROOT/'code/data/mps_s61_coordinate_catalogue_20261103_v62.json').read_text())['records']}
    crosswalk={r['id']:r for r in json.loads((ROOT/'code/data/component_manuscript_catalogue.json').read_text())['records']}
    rows=[]
    for index,record in enumerate(components):
        cid=record['id']
        proofs=[p.get('path',p.get('url',''))+(' # '+p['anchor'] if 'anchor' in p else '') for p in record['proof_sources']]
        if record['generic_rank']==2:
            entrance='R2'; divisor='Pfaffian identically zero; no actual quintic divisor'
            evidence=f'research/COMPONENTS.json:components[{index}].generic_rank'
        else:
            d=record.get('generic_actual_degeneracy_divisor')
            if d is not None:
                degrees=d['irreducible_factor_degrees']; multiplicities=d['corresponding_multiplicities']
                evidence=f'research/COMPONENTS.json:components[{index}].generic_actual_degeneracy_divisor'
            else:
                # These are the generic smoothed source divisors, not the toric boundary witness.
                degrees=mps[cid]['source_diagram']['generic_degeneracy_degree_partition']
                multiplicities=[1]*len(degrees)
                evidence=f'code/data/mps_s61_coordinate_catalogue_20261103_v62.json:{cid}.source_diagram.generic_degeneracy_degree_partition; research/notes/MPS40_CROSSWALK_20261031.md'
            assert sum(d*m for d,m in zip(degrees,multiplicities))==5
            entrance='A' if 1 in degrees else 'B' if len(degrees)>1 else 'C'
            divisor=' + '.join((str(m)+'*' if m!=1 else '')+f'D{d}' for d,m in zip(degrees,multiplicities))
            if d and d.get('description'):
                divisor+='; '+d['description']
        whole={'A':'KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md §§3-4; GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md',
               'B':'KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md §3.4 (all Q ranks)',
               'C':'KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md §2 (normal and nonnormal)',
               'R2':'RANK_TWO_SUPPORT.md; KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md §6'}[entrance]
        mentions=[crosswalk[cid]['main_label'],crosswalk[cid]['construction_label']]
        rows.append({'id':cid,'generic_actual_divisor':divisor,'entrance':entrance,
                     'divisor_evidence':evidence,'whole_support_input':'research/classification/'+whole,
                     'componenthood_and_distinctness_sources':' | '.join(proofs),
                     'main_manuscript_locator':'; '.join(mentions+['prop:global-exact147'])})
    assert len(rows)==147 and len({r['id'] for r in rows})==147
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print(json.dumps({'rows':len(rows),'generic_entrances':dict(collections.Counter(r['entrance'] for r in rows)),
        'scope':'Navigation of recorded generic divisors and proof sources; pointwise entrances are not a partition of component closures.',
        'output':str(args.output)},ensure_ascii=False))


if __name__=='__main__':
    main()
