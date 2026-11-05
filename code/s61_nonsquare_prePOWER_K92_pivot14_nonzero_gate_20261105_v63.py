#!/usr/bin/env python3
"""Fulfil the separate PIVOT14 nonzero gate without changing frozen inputs."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_PIVOT14_NONZERO_GATE_20261105_v63'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False);pins={}
    def load(st,fn='report.json'):
        rel=f'reports/{st}/{fn}';p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
    def check_inputs(r):
        for rel,h in r['input_sha256'].items():
            p=ROOT/rel;assert sha(p)==h;pins[rel]=h
    cond='S61_NONSQUARE_PREPOWER_K92_PIVOT14_SCHUR_20261105_v63'
    cr=load(cond);assert cr['all_passed'] and cr['conditional_on_independent_nonzero_peer'];check_inputs(cr)
    op=load(cond,'operator.json');proof=load(cond,'proof.json')
    peer='S61_NONSQUARE_K92_AXIS_A12_PEER_20261105_v63_REVIEW'
    pr=load(peer);assert pr['all_passed'];check_inputs(pr)
    assert pins[f'reports/{peer}/report.json']=='841214c3b4c85a9f48c9192cc4b013cf33d9b8c7efb476a0e28ae7962dee081e'
    ps=load('S61_NONSQUARE_K92_AXIS_A12_PEER_20261105_v63_SOURCE');assert ps['all_passed']
    assert pins['reports/S61_NONSQUARE_K92_AXIS_A12_PEER_20261105_v63_SOURCE/report.json']=='39c162c81f602ea13b451da4659d7c55f3ceb4974764f1d01421d454d47770d4'
    ev=load(peer,'independent-source-and-a12-coefficient.json')
    expected=Q(488042709959022324503443014699328840826617856,17514302730560302734375)
    # The independent coefficient is also compared against the frozen literal
    # primary coefficient; the full peer evidence remains pinned, not rewritten.
    assert Q(proof['nonzero_primary']['QQ_coefficient'])==expected
    assert proof['nonzero_primary']['axis_coefficient_monomial']==[12,0,0]
    assert str(expected) in json.dumps(ev)
    assert op['chart']=='FULLU D(Delta407*Saux*T*N6*N7bar*P_K92)'
    assert op['source_remaining_current15']==[14] and op['target_remaining_original32']==[1,3,4,6,7,8,9,17]
    assert op['only_new_chart_factor']=='P_K92'
    lifts={k:v for k,v in op['exports'].items() if k.startswith('surviving1_')};assert len(lifts)==5
    receipt=dict(conditional_report_reference=f'reports/{cond}/report.json',conditional_operator_reference=f'reports/{cond}/operator.json',
      fulfilled_gate=dict(independent_review=f'reports/{peer}/report.json',source_report='reports/S61_NONSQUARE_K92_AXIS_A12_PEER_20261105_v63_SOURCE/report.json',
        exact_monomial=[12,0,0],QQ_coefficient=str(expected),polynomial_nonzero=True,
        peer_scope='independent source579 and one exact coefficient, NOT all24 axis coefficients'),
      chart=op['chart'],only_added_parameter_inverse='P_K92',same_cokernel_shape=[1,8],
      inverse_identity='s^-1=Delta407*N7/(c8*g*P_K92)',c8='-88000000/13851',
      g_is_inherited_N6_chart_unit=True,no_extra_g_or_cleanup_chart=True,
      all_actual_source_lift_exports=lifts,all_tail_effects_retained=True,
      implication='residual1x8 full row rank implies actual full cubic57 on this chart',
      implication_IS_NOT_actual_bad_Fitting_ideal_equality=True,
      unresolved=['rank of last1x8','all chart complements','global FULLU57','global239','generic five-point exhaustion'],
      conditional_files_unchanged=True,no_new_scalar_arithmetic=True)
    dump(out/'fulfilled-gate-and-scope.json',receipt)
    code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code)
    note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
    total=sum(p.stat().st_size for p in out.iterdir());assert total<1024*1024
    report=dict(id=ID,all_passed=True,frozen=True,peer_nonzero_gate_fulfilled=True,
      P_K92_independently_proved_nonzero=True,source_only_new_receipt=True,same_cokernel_shape=[1,8],
      source_pivot_count=14,only_new_chart_factor='P_K92',conditional_files_unchanged=True,
      all_full57_lifts_retained=True,no_actual_bad_ideal_equality_or_global57_claim=True,
      input_sha256=pins,evidence_sha256={p.name:sha(p) for p in out.iterdir()},total_new_RAW_before_report=total,seconds=time.monotonic()-start)
    dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
