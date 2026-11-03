#!/usr/bin/env python3
"""Preserve the dependency/coverage audit of the normal-actual inactive-line proof.

This verifies archived input integrity and the proof dependency graph, not the
geometric theorems by executable algebra.  Unchanged old certificates are read,
not rerun.  The accompanying note gives the independent mathematical argument.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]

NODES=[
 ('canonical_BE',[],['generic_rank4','normal_actual_D'],
  'canonical quadratic lift; b(P)=dF; P*dF=0; exact height3 BE presentation'),
 ('line_triple',['canonical_BE'],['inactive_embedded_line'],
  'P vanishes on coneL, hence F in I_L^3 globally'),
 ('local_five_generators',['canonical_BE','line_triple'],[],
  'at eta_L exact presentation modulo m has zero relation matrix: mu(J_eta)=5'),
 ('multiplicity_3_or_4',['local_five_generators'],[],
  'degree5 excludes m>5; m5 and dependent quartic pencil contradict five minimal generators'),
 ('all_point_affine_nilpotence',['canonical_BE','line_triple'],[],
  'formal semisimple Casimir lemma at every nonzero coneL point, including all intersections'),
 ('m3_fixed_cube_and_flag',['all_point_affine_nilpotence'],[],
  'nilpotent rank4 center gives cube; primitive binary factorization 2<3 fixes w; sl2 pencil fixes flag'),
 ('m3_rank1',['m3_fixed_cube_and_flag'],[],
  'b nonzero: exact full Casimir and divergence imply w divides F'),
 ('m3_zero_block',['m3_fixed_cube_and_flag'],[],
  'b=0 and lower block zero: one-point Heisenberg/normal Jacobi imply w divides F'),
 ('m3_squarefree',['m3_fixed_cube_and_flag','m3_zero_block'],[],
  'b=0, c squarefree: two1:-4 resonances exchanged by s,t; coefficient certificate48gamma^2x'),
 ('m3_double',['m3_fixed_cube_and_flag','m3_zero_block'],[],
  'b=0, c double: two nonzero-delta branches plus delta0 exception; singular plane in all3'),
 ('m4_fixed_flag',['multiplicity_3_or_4','all_point_affine_nilpotence'],[],
  'B=0; cross(Mz,Nz)=0; exceptional3x3 pencil has no bidegree(1,4) invariant; M,N notboth0'),
 ('m4_trichotomy',['m4_fixed_flag'],[],
  'constant regular direction, constant common image, or constant proportional pair with moving direction'),
 ('m4_regular',['m4_trichotomy'],[],
  'invertibleV; rank1 nonzero trace; rank1 zero trace; retain h2 shear residual d until closedness kills it'),
 ('m4_fixed_image',['m4_trichotomy'],[],
  'gcd/syzygy for all e=0,1,2, including C=0 orZ=0, implies F in planeideal^2'),
 ('m4_three_moving',['m4_trichotomy'],[],
  'independent ell,q; dependent nonzero ell,q and independentp; ell0 with independentp,q; fullV transforms'),
 ('no_inactive_line',['multiplicity_3_or_4','m3_rank1','m3_zero_block','m3_squarefree','m3_double',
                     'm4_regular','m4_fixed_image','m4_three_moving'],[],
  'all closed-point alternatives force factorization or a two-dimensional singular plane'),
 ('smooth_inactive_component',['no_inactive_line'],['v16_componentwise_smooth_curve_to_line'],
  'no smooth inactive irreducible component even inside reducible C'),
 ('v20_singular_component_connection',['smooth_inactive_component'],['v20_secant_directrix_lemma'],
  'singular inactive component: spanP4 gives singular join surface, spanP3 factor, plane quadratic vanishing'),
]


def verify(output):
    start=time.monotonic();inputs={};archive_checks=[]
    def record(path):
        p=ROOT/path;digest=hashlib.sha256(p.read_bytes()).hexdigest();inputs[path]=digest
        if path.startswith('legacy/'):
            version=Path(path).parts[1];local='/'.join(Path(path).parts[2:])
            manifest=ROOT/'legacy'/version/'SHA256SUMS.txt'
            hashes=dict((line.split(maxsplit=1)[1].lstrip('*'),line.split(maxsplit=1)[0])
                        for line in manifest.read_text().splitlines() if line.strip())
            assert hashes[local]==digest,path
            archive_checks.append(path)
        return p
    sources=['legacy/v14/sections/normal_quintics.tex','legacy/v16/sections/setup.tex',
             'legacy/v16/sections/smooth_curves.tex','legacy/v17/sections/setup.tex',
             'legacy/v17/sections/nilpotent_isotropy.tex','legacy/v17/sections/jacobian.tex',
             'legacy/v17/sections/line_algebra.tex','legacy/v18/sections/setup.tex',
             'legacy/v18/sections/m3_restrictions.tex','legacy/v18/sections/m4_flag.tex',
             'legacy/v18/sections/m4_regular_direction.tex','legacy/v19/sections/setup.tex',
             'legacy/v19/sections/conclusions.tex','legacy/v19/sections/verification.tex',
             'legacy/v20/sections/setup.tex','legacy/v20/sections/inactive.tex',
             'legacy/v20/sections/secant.tex']
    sources += [str(p.relative_to(ROOT)) for p in sorted((ROOT/'legacy/v19/sections').glob('m*.tex'))]
    for p in sources:record(p)
    for version in [18,19]:
        path=f'legacy/v{version}/code/run_checks.py';tree=ast.parse(record(path).read_text())
        jobs=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
                  and any(isinstance(t,ast.Name) and t.id=='JOBS' for t in n.targets))
        summary=json.loads(record(f'legacy/v{version}/data/reproduced/verification_summary.json').read_text())
        assert summary['all_passed'] and summary['check_count']==len(jobs)
        for name,sub,script,reportname in jobs:
            record(f'legacy/v{version}/proofs/{sub}/{script}')
            archived=json.loads(record(f'legacy/v{version}/data/reproduced/reports/{name}.json').read_text())
            assert archived['all_passed'],name
    ids={node[0] for node in NODES};visited=set()
    for identifier,dependencies,hypotheses,conclusion in NODES:
        assert set(dependencies)<=visited,(identifier,dependencies)
        visited.add(identifier)
    assert visited==ids
    forbidden=['literal_local_length40','whole_Y_finite_flat_degree40_over_L',
               'global_splitting_for_support_L','radical_J_equals_I_L','tau_mu_numerical_bound']
    assert not ids.intersection(forbidden)
    record('audits/2026-10-30/geometry_strategy.txt');record(str(Path(__file__).relative_to(ROOT)))
    report=dict(id='NORMAL_ACTUAL_INACTIVE_LINE_SCOPE_AUDIT_20261103_v62',all_passed=True,
                human_proof_scope_verdict='No C=L-only result is used inside the reducible-C inactive-line exclusion.',
                theorem_hypotheses=['complex closed point','generic rank4 Poisson onP4',
                                    'normal actual divisor D, not only D_red','inactive embedded projective line L'],
                dependency_nodes=[dict(id=i,depends_on=d,extra_hypotheses=h,conclusion=c) for i,d,h,c in NODES],
                excluded_whole_support_dependencies=forbidden,
                archived_v18_certificate_groups=5,archived_v19_certificate_groups=9,
                archived_checksum_matches=len(archive_checks),
                archived_certificates_rerun=False,
                routine_ambiguities_resolved=['five-generator proof uses right exactness only',
                                             'all-point normality bound allows intersecting C components',
                                             'constant cube direction extends across zeros without division there',
                                             'm4 moving rank1 image is distinct from fixed image',
                                             'h2 shear residual retained and killed by -5alpha*d',
                                             'xi^2 ideal membership is used only for complex closed points'],
                unresolved_finite_problem_for_inactive_line=None,
                limits=['This is a human proof/coverage audit plus archived-certificate integrity, not a fresh execution of old algebra.',
                        'The entire normal-D classification still requires the separate E6/E8 and eight-support chain.',
                        'No conclusion for nonnormal or nonreduced actual divisors; no whole Jacobi component count.'],
                input_sha256=inputs,seconds=time.monotonic()-start)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    r=verify(p.parse_args().output);print('PASS',r['id'],r['archived_checksum_matches'])
