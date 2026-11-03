#!/usr/bin/env python3
"""Record the independent proof audit of v14 normal-quintic recognition.

The universal and relative assertions are proved in the accompanying note.
Here unchanged archived sources and saved sample certificates are identified;
the examples are not treated as proof of the universal statement or rerun.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def audit(output):
    base=ROOT/'legacy/v14';inputs={}
    manifest=dict((line.split(maxsplit=1)[1].lstrip('*'),line.split(maxsplit=1)[0])
                  for line in (base/'SHA256SUMS.txt').read_text().splitlines() if line.strip())
    def record(path):
        p=ROOT/path;value=hashlib.sha256(p.read_bytes()).hexdigest();inputs[path]=value
        if p.is_relative_to(base):assert manifest[str(p.relative_to(base))]==value,p
        return p
    for path in ['legacy/v14/sections/normal_quintics.tex','legacy/v14/proofs/normal_quintic_full_scheme_audit.txt',
                 'legacy/v14/data/reproduced/verification_summary.json']:
        record(path)
    for p in sorted((base/'code/normal').rglob('*')):
        if p.is_file() and p.suffix in ['.py','.json']:record(str(p.relative_to(ROOT)))
    summary=json.loads((base/'data/reproduced/verification_summary.json').read_text())
    for name in ['normal_quintic_symbolic','normal_quintic_independent_integer']:
        row=next(c for c in summary['checks'] if c['name']==name);assert row['passed'] and row['exit_code']==0
        p=record('legacy/v14/code/'+row['script']);assert hashlib.sha256(p.read_bytes()).hexdigest()==row['source_sha256']
        data=json.loads(record(f'legacy/v14/data/reproduced/outputs/{name}.json').read_text())
        if name=='normal_quintic_symbolic':
            assert data['Fisher_gradient_equals_signed_Pfaffians'] and data['Fisher_divergence_zero'] and data['Fisher_Jacobi_zero']
            assert data['Fisher_annihilator_rank']==149
        else:
            assert data['exact_Pfaffian_gradient_identity'] and data['exact_Jacobi_zero'] and data['exact_divergence_zero']
    record(str(Path(__file__).relative_to(ROOT)))
    report=dict(id='NORMAL_QUINTIC_RECOGNITION_AUDIT_20261103_v62',all_passed=True,
                necessary_direction='genericrank4 projectivePoisson and normalACTUALD imply canonicalaffinePoissonP,gradF=Pfaffians,height3,exactminimalBE(0,4,6,10)',
                closed_point_recognition_equivalences=['D is actualdegeneracydivisor ofrank4Poisson',
                                                       'heightJ3,GorensteinS/J,fiveindependentquarticpartials',
                                                       'centralitykernelK_F containsrank4quadraticalternatingmatrix',
                                                       'K_F isaline withnonzeroPfaffianvector'],
                proof_obligations_audited=['canonicaldivergencefree splitting and affineJacobi fromprojectiveJacobi viaEuler-divergence identity',
                                          'BV andCartan fixactualF gradient andfactor5 onaffinechart',
                                          'normalitygivesheight>=3,Pfaffianheight<=3,regularringheight=grade forBE',
                                          'Gorensteinselfdualgradedresolution giveslastshift10 andquadraticmatrix',
                                          'GL5congruence realizesanyquarticgeneratorbasisoverC',
                                          'normalF gradientgcd1 convertsfractionfieldproportionality toconstant',
                                          'no linearJacobiideal syzygies implydivP0 andBV impliesJacobi',
                                          'exactBEdegreecomparison provesK_F uniqueness'],
                full_scheme_recognition=True,
                relative_scope='geometricallynormalhypersurfaceFIBRES onthenormalopen; totalDovernonreducedbase isnotassertednormal',
                relative_proof_obligations_audited=['actual150minors plusunit149minor forcezeroSchurcolumn overANYring,notjustresiduefields',
                                                   'Pfaffiancomplex ispolynomialcomplex beforeexactness',
                                                   'eachhomogeneousstrand finitefree; residueacyclicpositive => splitcontractiblepairs andfreeH0',
                                                   'BE appliedonlyoverresiduefields,notassumedoverarbitrarynonCMbase',
                                                   'relativeexactnessdegree10 givesgradF=cPfaffians withcunit',
                                                   'degree5 hasnosyzygy givesdivP0 exactlyincludingnilpotents',
                                                   'projectivescalingglueswithoutasquareroot; inverseconstructionscommutewithallbasechange'],
                unresolved_proof_gap=None,
                existing_example_certificates_rerun=False,
                limits=['These are proof-audit conclusions; the two saved examples are not evidence of universal coverage by themselves.',
                        'Recognition is independent of v22closedpointclassification and doesnot implythePoissonschemeisreduced.',
                        'No claim of arbitrary-base orbitclassification or an explicitPGL5transport toq51.'],
                input_sha256=inputs)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    r=audit(p.parse_args().output);print('PASS',len(r['input_sha256']),r['full_scheme_recognition'])
