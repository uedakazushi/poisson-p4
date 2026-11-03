#!/usr/bin/env python3
"""Audit complete coverage and unchanged archived evidence for the28 RNC4 families.

No old symbolic calculation is rerun; this checks hashes, family correspondence,
and a complete explicit assignment of the proof's closed-point exit conditions.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def audit(output):
    inputs={};base=ROOT/'legacy/v16'
    manifest=dict((line.split(maxsplit=1)[1].lstrip('*'),line.split(maxsplit=1)[0])
                  for line in (base/'SHA256SUMS.txt').read_text().splitlines() if line.strip())
    def record(p):
        value=hashlib.sha256(p.read_bytes()).hexdigest();inputs[str(p.relative_to(ROOT))]=value
        if p.is_relative_to(base):assert manifest[str(p.relative_to(base))]==value,p
        return p
    for folder in ['sections','code/geometry','code/semisimple','code/nilpotent','data/reproduced/reports']:
        for p in sorted((base/folder).rglob('*')):
            if p.is_file() and p.suffix in ['.py','.tex','.json']:record(p)
    summary=json.loads(record(base/'data/reproduced/verification_summary.json').read_text())
    assert summary['all_passed'] and len(summary['checks'])==10
    for check in summary['checks']:
        assert check['passed']
        assert hashlib.sha256(record(base/check['script']).read_bytes()).hexdigest()==check['script_sha256']
    original=base/'previous_v15/code/normal_geometry/rnc4_slices.json'
    reduction=original.with_name('rnc4_equation_reduction.json')
    candidate=base/'code/geometry/candidate_parameterizations.json'
    for name,p in [('v15_rnc4_slices',original),('v15_equation_reduction',reduction),('v16_candidates',candidate)]:
        assert hashlib.sha256(record(p).read_bytes()).hexdigest()==summary['input_sha256'][name]
    data=json.loads(candidate.read_text());cases=data['cases'];assert len(cases)==28
    reports=base/'data/reproduced/reports'
    cover=json.loads((reports/'conormal_cover_28.json').read_text())
    assert cover['all_passed'] and cover['candidate_families_verified']==28 and cover['all_branch_bases_verified_complete']
    semi=json.loads((reports/'semisimple_21.json').read_text());nil=json.loads((reports/'nilpotent_7.json').read_text())
    assert semi['all_21_passed'] and nil['all_passed']
    assert {c['id'] for c in semi['checks']}=={c['id'] for c in cases if c['residue']=='semisimple'}
    assert {c['case'] for c in nil['cases']}=={c['id'] for c in cases if c['residue']=='nilpotent'}
    coverage=[]
    for case in cases:
        d=case['d'];h=case['h_exponents'];branch=case['branch'];residue=case['residue']
        if residue=='semisimple':
            if d==2 and h==[1,1] and branch=='diagonal_zero':reason='v0 or explicit reducible quintic FB'
            elif d==3 and branch=='diagonal_zero':reason='v(0,1)=0 contradicts primitive degree3'
            elif d==3:reason='primitive implies ab=2/3; entireCstar orbit has explicit singular surface'
            else:reason='finite powers of every v coefficient in Jacobi ideal imply v0'
        else:
            if d==2 and branch=='diagonal_zero':reason='quartic k roots:2 give v0;1 nonreducedF;1 entireGa orbit singular surface'
            elif d==2:reason='Q=0 contradicts required Qnonzero; boundary retained in diagonalzero family'
            elif d==3:reason='v(0,1)=0 contradicts primitive degree3'
            else:reason='squares of every v coefficient in Jacobi ideal imply v0'
        coverage.append(dict(id=case['id'],residue=residue,d=d,h_exponents=h,branch=branch,
                             v_dimension=len(case['v_coefficient_basis']),exit=reason,
                             closed_point_only=True))
    assert Counter(c['residue'] for c in coverage)=={'semisimple':21,'nilpotent':7}
    inputs[str(Path(__file__).relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report=dict(id='NORMAL_ACTUAL_RNC4_COVER_AUDIT_20261103_v62',all_passed=True,
                theorem='No normal actual degeneration divisor has reduced singular support an irreducible rational normal quartic.',
                candidate_count=28,residue_counts=dict(Counter(c['residue'] for c in coverage)),families=coverage,
                complete_cover_dependencies=['nonzeroZ normalized byPGL2+scalar toAh orAe',
                                             'full36variable affine lift, with no imposed[A,P]=0',
                                             'conormalE=O(-6)^3; nonzero rank1 symmetric sexticB',
                                             'primitiveB=hvvT with0<=d<=3,zero(h)subsetzero(Z)',
                                             'rankL33,rankT9,TL0,andB0inkernelT imply exact affineimagekernelT',
                                             'Cartan identity splits Q0 versus linear branch; all kernels complete',
                                             'all28originalJacobi substitutions retain3kernelvariables'],
                auxiliary_fixed_symmetry_results_used=False,
                archived_verification_groups_preserved=10,
                essential_unchanged_groups=['v15_original_Jacobi','conormal_cover_28','semisimple_21',
                                            'semisimple_nonnormal','nilpotent_7','nilpotent_nonnormal'],
                additional_independent_hard_nilpotent_preserved=True,
                old_certificates_rerun=False,unresolved_finite_problem=None,
                limits=['Human proof coverage plus preserved existing exact certificates; not a new symbolic execution.',
                        'Primitive-factor exit conditions concern complex closed points, not absence of nilpotents in coefficient schemes.',
                        'No classification of the full36variable schemes or of nonnormalD is claimed.'],
                input_sha256=inputs)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    r=audit(p.parse_args().output);print('PASS',r['candidate_count'],r['residue_counts'],len(r['input_sha256']))
