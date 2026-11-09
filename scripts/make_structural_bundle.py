#!/usr/bin/env python3
"""Build the current research bundle without nested archives or old run caches.

Includes the complete current main manuscript and its shared proofs, original
textual inputs, and standalone checks for this milestone. It does not package
all historical finite certificates or independently certify the 147 components.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    files={
        'README.md','AGENTS.md','RESEARCH_POLICY.md','requirements.txt',
        'research/STATE.md','research/NEXT_TASKS.md','research/REFERENCES.md',
        'research/COMPONENTS.json','research/CLAIMS.json',
        'research/conceptual/GLOBAL_PLAN.md','research/conceptual/NOTES.md','research/conceptual/FINAL_TARGET.md',
        'research/conceptual/component_entrances.csv','research/classification/CURRENT.md',
        'papers/poisson_p4_classification_ja.tex','papers/poisson_p4_classification_ja.pdf',
        'papers/poisson_divisor_symmetries_ja.tex','papers/poisson_divisor_symmetries_ja.pdf',
        'papers/marked_hyperplane_reconstruction.tex','papers/poisson_lifting_examples.tex',
        'code/check_global_structure.py','code/check_marked_reconstruction.py',
        'code/check_quadric_residue_reconstruction.py','code/data/quadric_cofactor_sources.json',
        'code/check_conductor_conormal_reconstruction.py','code/check_normalization_residue_reconstruction.py',
        'code/check_rank_two_reconstruction.py','reports/rank_two_reconstruction_checks.json',
        'reports/rank_two_component_recheck.json',
        'reports/normalization_residue_checks.json','reports/normalization_residue_validation.json',
        'reports/conductor_conormal_validation.json','reports/conductor_conormal_checks.json',
        'code/map_component_entrances.py','code/data/mps_s61_coordinate_catalogue_20261103_v62.json',
        'research/notes/MPS40_CROSSWALK_20261031.md',
        'scripts/build_structural_manuscripts.sh','scripts/make_structural_bundle.py',
        'reports/quadric_residue_validation.json','reports/quadric_residue_checks.json','reports/marked_reconstruction_checks.json',
        'reports/marked_reconstruction_base_identity_checks.json',
    }
    # Preserve all current textual classification proofs. The large finite
    # certificate collection remains indexed in the unchanged original ledger.
    files.update(str(p.relative_to(ROOT)) for p in (ROOT/'research/classification').glob('*.md'))
    # Keep the original textual component inputs, including materialized legacy
    # TeX, instead of putting previous-version ZIPs inside the new bundle.
    records=json.loads((ROOT/'research/COMPONENTS.json').read_text())['components']
    for record in records:
        for source in record['proof_sources']:
            name=source.get('path','')
            if name and Path(name).suffix in {'.md','.tex','.json'} and (ROOT/name).is_file():
                files.add(name)
    files.update(str(p.relative_to(ROOT)) for p in (ROOT/'papers').glob('*.tex'))
    # Expand the manuscript includes and any explicit materialized legacy TeX
    # dependencies mentioned by the complete main text.
    maintext=(ROOT/'papers/poisson_p4_classification_ja.tex').read_text()
    for name in re.findall(r'legacy/[A-Za-z0-9_./-]+\.(?:tex|json)',maintext):
        if (ROOT/name).is_file():
            files.add(name)
    missing=[name for name in files if not (ROOT/name).is_file()]
    if missing:
        raise FileNotFoundError(missing)
    forbidden={'.zip','.gz','.xz','.bz2','.tar','.log','.aux','.out','.pyc'}
    assert all(Path(name).suffix not in forbidden for name in files)
    manifest={'version':'v72-normalization-residue-integration','date':'2026-11-09',
        'repository':'https://github.com/uedakazushi/poisson-p4',
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'scope':'Complete current main and structural manuscripts, shared proofs, original textual classification/component inputs, and standalone checks of this milestone.',
        'limits':['Not the full transitive closure of all 147 historical finite computation certificates.',
                  'Not an independent verification of the 147-component classification.',
                  'Legacy original TeX is evidence; rebuild the current self-contained manuscripts, not every legacy draft.'],
        'library_save':'Not performed by this script. The separate Page upload receipt records Library persistence.',
        'reproduction':['python code/check_global_structure.py --output work/global-structure.json',
            'python code/check_marked_reconstruction.py --output work/marked-reconstruction.json',
            'python code/check_quadric_residue_reconstruction.py --output work/quadric-residue-checks.json',
            'python code/check_conductor_conormal_reconstruction.py --output work/conductor-conormal-checks.json',
            'python code/check_normalization_residue_reconstruction.py --output work/normalization-residue-checks.json',
            'python code/check_rank_two_reconstruction.py --output work/rank-two-checks.json',
            'python code/map_component_entrances.py --output work/component_entrances.csv',
            'bash scripts/build_structural_manuscripts.sh'],
        'files':[]}
    output=args.output
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    with zipfile.ZipFile(output,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for name in sorted(files):
            data=(ROOT/name).read_bytes()
            archive.writestr('poisson-p4/'+name,data)
            manifest['files'].append({'path':name,'bytes':len(data),'sha256':digest(data)})
        archive.writestr('poisson-p4/BUNDLE_MANIFEST.json',json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    with zipfile.ZipFile(output) as archive:
        assert archive.testzip() is None
        for record in manifest['files']:
            assert digest(archive.read('poisson-p4/'+record['path']))==record['sha256']
    print(json.dumps({'output':str(output),'files':len(files),'bytes':output.stat().st_size,
        'sha256':digest(output.read_bytes()),'zip_crc':'PASS','all_member_sha256':'PASS',
        'scope':manifest['scope'],'limits':manifest['limits'],'library_save':manifest['library_save']},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
