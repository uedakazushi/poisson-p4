#!/usr/bin/env python3
"""Logical gluing ONLY of independently accepted H15 complement and full finite cover."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_WHOLE_QUARTIC150_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out,peer,expected):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists();phase=read(out/'phase.json')
 old=ROOT/phase['old_exact_H15_primary'];new=ROOT/phase['new_finite_full239_primary'];op=ROOT/'reports/S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63_REVIEW/report.json'
 assert sha(old)==phase['old_primary_SHA256']and sha(new)==phase['new_primary_SHA256']
 assert sha(op)=='2b4abe5029157c3aa08dbdd9ba26f256cb9b49b03bfa77146f88cf91d1e48bcd'
 assert sha(peer)==expected
 a,b,c,d=map(read,[old,new,op,peer]);assert a['complete']and a['exact_QQ_gcd_degree']==15
 assert b['all_passed']and b['positive_full239_on_finite_H15']and b['native_summary']['completed_full98_factors']==7
 assert c['all_passed']and c['old_five_minors_exact_common_gcd_H15_verified']
 assert d['all_passed']and d.get('frozen')is True
 assert c['input_sha256'][str(old.relative_to(ROOT))]==sha(old)
 assert d['primary_report_sha256']==sha(new)
 assert d['positive_full239_on_entire_H15_only']
 assert d['independent_arithmetic']['all319_catalogue_retained']
 assert d['independent_arithmetic']['section_uses_actual_delta7_L_and_corrected_original112_tails']
 assert d['native_fixed_rank_audit']['CRT_full_Hbar_row_surjectivity_gate']
 nb=read(new.with_name('native-summary.json'));assert nb['all_passed']
 for i in range(7):
  f=read(new.with_name(f'factor{i}.json'));assert f['rank']==98 and f['positive_full98']and any(f['selected98_minor_determinant'])
 # Mathematical lemma, not a new coefficient identity or minor computation.
 cert=dict(ring='T100=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))]',
  actual_map_shape=[319,239],Gamma_shape=[319,98],global_removed_unit_target_directions=141,
  finite_chosen_cover_degree15=True,
  proof_steps=[
   'Accepted global NF117, H2 unit1 and H1 section23 identify coker(original319x239) with coker(Gamma319x98) overT100.',
   'Consequently Fitt89(coker)=I150(original)=I9(Gamma).',
   'The five actual9-minors have ideal(H), with only known T100-unit row LCDs; hence H lies in I9(Gamma).',
   'Independently accepted all-seven-factor positive certificates and finite V-Nakayama give Gamma modH surjective, hence I98(Gamma)+(H)=T100.',
   'Since I98(Gamma) is contained in I9(Gamma), the preceding two statements imply I9(Gamma)=T100.',
   'Thus the actual quartic matrix has rank at least150 at EVERY permitted complex P100 point. OnH=0 it has rank239.'
  ],whole_Fitting_I150_equals_T100=True,no_fixed150_minor_or_expanded_unit_cofactors_computed=True,
  finiteH15_full239=True,wholeP100_pointwise_lower150=True,wholeP100_full239_NOT_proved=True,
  no_fullFULLU_full239_generic_exhaustion_or_global_component_upper_claim=True,
  old_peer= str(op.relative_to(ROOT)),finite_peer=str(peer.relative_to(ROOT)),
  no_new_coefficient_determinant_rank_GB_point_prime=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 phase.update(phase='both_independent_final_gates_passed_logical_combination_complete',old_peer_SHA256=sha(op),finite_peer_SHA256=sha(peer))
 (out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
 context=ROOT/'reports/S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63/operator.json'
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',old,new,op,peer,context,
  ROOT/'reports/S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/report.json']
 result=dict(id=ID,all_passed=True,frozen=True,whole_P100_pointwise_quartic_lower150=True,
  whole_P100_Fitting_I150_is_unit_ideal=True,finite_H15_quartic_rank239=True,
  whole_P100_quartic239_NOT_proved=True,independent_gate_reports_final=True,
  no_new_arithmetic_rank_minor_source_point_or_prime=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},certificate_sha256=sha(out/'certificate.json'),phase_sha256=sha(out/'phase.json'),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--finite-peer',required=True,type=Path);p.add_argument('--finite-peer-sha',required=True);a=p.parse_args();main(a.output.resolve(),a.finite_peer.resolve(),a.finite_peer_sha)
