#!/usr/bin/env python3
"""Coverage audit of all eight v22 supports and all localized certificate cells.

Archived algebra certificates are not rerun.  This independently enumerates the
finite rational spectrum steps and checks the complete certificate manifest.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def spectrum(h,m):return (Q(1),-Q(h,2),-Q(m,3),Q(h,2)+Q(m,3)-1)


def audit(output):
    base=ROOT/'legacy/v22';inputs={}
    manifest=dict((line.split(maxsplit=1)[1].lstrip('*'),line.split(maxsplit=1)[0])
                  for line in (base/'SHA256SUMS.txt').read_text().splitlines() if line.strip())
    def record(p):
        value=hashlib.sha256(p.read_bytes()).hexdigest();inputs[str(p.relative_to(ROOT))]=value
        if p.is_relative_to(base):assert manifest[str(p.relative_to(base))]==value,p
        return p
    for directory in ['sections','code','data','audits','proofs/independent']:
        for p in sorted((base/directory).rglob('*')):
            if p.is_file() and p.suffix in ['.tex','.py','.json','.txt']:record(p)
    summary=json.loads((base/'audits/summary.json').read_text())
    assert summary['all_passed'] and summary['group_count']==10 and summary['checks']==3229
    for group in summary['groups']:
        saved=json.loads((base/'audits'/f"{group['name']}.json").read_text())
        assert saved.get('all_passed',saved.get('all_checks_passed',False)) or (
            isinstance(saved.get('checks'),dict) and all(saved['checks'].values())),group['name']
    jets=json.loads((base/'data/e8_jets/summary.json').read_text())
    names={m['name'] for m in jets['models']};assert len(names)==8 and jets['all_structural_checks_passed']
    for name in names:
        original=ROOT/'legacy/v21/data/tree_jacobi_models'/f'{name}.json'
        inherited=base/'previous_v21/data/tree_jacobi_models'/f'{name}.json'
        assert original.read_bytes()==record(inherited).read_bytes(),name
    coverage={
     '1111_chain':'two distinct marked nodes force simplezeros, but line-line theorem forces doublezero',
     '1111_triple_tail':'central line has a separate node andtriplepoint; node doublezero exceeds totaldegree2',
     '1111_concurrent':'four independentaxes force allmodularslopes0; allP independentofvertex; BE excludescone',
     '211_conic_center':'conic-plane4jet gives nowherezero bracket; both intersections doublezeros ononeP1',
     '211_line_center':'central line getsdoublezero atconic intersection plusanotherzero atlineintersection',
     '211_triple':'allslopes0; four nonzero normalizations; rankone and polynomial3term certificate forces z20=0',
     '22_node':'bothplaneHessian coefficientsnonzero; equalslopes; SSunit, NNgenericunit, lambda-1spectralunit',
     '31_external':'BTnode0; mixedslopesexcluded;7SS cells and NN k3/IA/IB/II withallexceptions'}
    assert set(coverage)==names
    line_matches=[]
    for h,m,hh,mm in it.product([1,2],range(3),[1,2],range(3)):
        for ratio in set(spectrum(h,m))-{Q(0)}:
            if Counter(spectrum(h,m))==Counter(ratio*x for x in spectrum(hh,mm)):
                line_matches.append(dict(h=h,m=m,h2=hh,m2=mm,ratio=str(ratio)))
    unequal=[c for c in line_matches if c['ratio']!='1']
    assert unequal==[dict(h=2,m=m,h2=2,m2=m,ratio='-1') for m in range(3)]
    strong=[]
    for m in range(3):
        for ratio in set(Q(1)/x for x in spectrum(1,m) if x):
            if Counter([Q(1),ratio,-ratio/2,-ratio/2-1])!=Counter(ratio*x for x in spectrum(1,m)):continue
            for k in range(1,4):
                if -Q(k+1,2) in [ratio,-ratio/2]:strong.append((str(ratio),k,m))
    assert sorted(strong)==[('-2',3,0),('-3/2',2,2)]
    ss=json.loads((base/'audits/31_ss_stdlib.json').read_text())
    assert len(ss['certificates'])==7 and not ss['parameter_denominators']
    assert sum(not c['unit'] for c in ss['certificates'])==1
    nn=json.loads((base/'audits/31_nn_stdlib.json').read_text())
    outcomes={c['name']:c for c in nn['certificates']}
    assert set(outcomes['31_NN_k1_IA_certificate']['denominator_factors'])=={'z17','69*z17+28'}
    assert set(outcomes['31_NN_k1_IB_certificate']['denominator_factors'])=={'z17','3*z17+1','621*z17**2+87*z17+52'}
    assert set(outcomes['31_NN_k1_II_certificate']['denominator_factors'])=={'z21'}
    assert outcomes['31_NN_k1_IA_rational_certificate']['outcome']=='unit'
    assert outcomes['31_NN_k1_IB_quadratic_field_certificate']['outcome']=='unit'
    assert outcomes['31_NN_k1_IB_rational_certificate']['outcome']=='actual nonnormal quintic family with singular plane'
    assert outcomes['31_NN_partition_certificate']['target_count']==12
    for label,p in [('model',base/'previous_v21/data/tree_jacobi_models/31_external.json'),
                    ('jets',base/'data/e8_jets/31_external.json')]:
        assert nn['source_sha256'][label]['sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
    two=json.loads((base/'audits/22_stdlib.json').read_text())
    assert len(two['certificates'])==4
    primary=ROOT/'work/cloud-d0-20261103/gualtieri-pym-1203.4293.pdf'
    record(Path(__file__))
    report=dict(id='NORMAL_ACTUAL_EIGHT_SUPPORT_AUDIT_20261103_v62',all_passed=True,
                support_coverage=coverage,
                independent_finite_line_spectrum_matches=line_matches,
                independent_strong_TC_spectrum_survivors=strong,
                certificates_22=two['certificates'],certificates_TC_SS=ss['certificates'],certificates_TC_NN=nn['certificates'],
                archived_exact_groups=10,archived_exact_checks=3229,unchanged_old_certificates_rerun=False,
                primary_GP_source=dict(url='https://arxiv.org/pdf/1203.4293',
                                       sha256=hashlib.sha256(primary.read_bytes()).hexdigest(),
                                       read='Section8 pp28-30: Theorem8.1 normality,Theorem8.3 actualdivisorequality,Proposition8.5 nonvanishingresidue'),
                geometric_audit_findings=['fourthjet conic obstruction is coordinate/unit invariant afterMorse elimination',
                                         'projective cubic identity applied onlyafterlinearcoordinatechanges',
                                         'allresiduenormalizations use verified nonzero coefficients andcomplete-linear-system supportautomorphisms',
                                         'SS finitecases require isotropy derivation andkernelweights, notweak spectralmatchingalone',
                                         'allNN generic denominators include everyintermediatecombination/substitution; allspecialvaluescovered',
                                         'actual nonnormalSS andNN families are retained, notcalledemptyideals'],
                unresolved_finite_problem=None,
                theorem_scope='Together with the independently audited v14-v21 inputchain, validates normal ACTUAL divisor classification atcomplexclosedpoints.',
                limits=['Human proof/coverage audit with archived exact algebra, not a fresh execution of3229 checks.',
                        'No full moduli-scheme reducedness or total Poisson component count; nonnormal/nonreduced actual divisors remain separate.'],
                input_sha256=inputs)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    r=audit(p.parse_args().output);print('PASS',len(r['support_coverage']),len(r['input_sha256']))
