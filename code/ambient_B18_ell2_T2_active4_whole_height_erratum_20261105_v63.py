#!/usr/bin/env python3
"""Immutable withdrawal receipt for the old nativePivotE height bound."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_ACTIVE4_WHOLE_HEIGHT_ERRATUM_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=ROOT/'reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_INPUT_GATE_20261105_v63/report.json';r=ROOT/'reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_INPUT_GATE_20261105_v63_REVIEW/report.json'
 assert sha(p)=='7b6a485ad19ffee2002c2c099133e56e4a4b7609f0179c3bd2c3f180874bb88b' and sha(r)=='afa5180917535d254fb9ec3a726af8deec95e274461bfbe3b793f1d061b6f815'
 # A union-max controls choosing either summand, not their actual sum.
 assert max(1,1)<1+1
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 d={'id':STEM,'all_passed':True,'erratum':True,'old_whole_height_gate_and_peer_withdrawn_as_value_GO':True,
 'withdrawn_primary':{'path':str(p.relative_to(ROOT)),'sha256':sha(p)},'withdrawn_review':{'path':str(r.relative_to(ROOT)),'sha256':sha(r)},
 'exact_error':'nativePivotE is the actual sum L0*Se+Le*S0; positive envelopes must use LCM/SUM, not union-max. The prior independent review incorrectly reproduced the union choice.',
 'proof':'If both terms have bounds(D1,U1),(D2,U2), their sum has D=lcm(D1,D2), U=U1*(D/D1)+U2*(D/D2). A max envelope is insufficient, as1+1>max(1,1). All downstream pivot,psi,coordinate,rank,ownership/RAW envelopes must be revalidated in a distinct V2 gate.',
 'old_source_custody_frame_sign_checks_not_retracted':True,'old_948_height_and_value_authorization_not_accepted':True,
 'no_actual_psi_dot_or_rank_values_run_before_detection':True,'all_old_bytes_preserved':True,'fresh_corrected_SUM_gate_peer_and_root_custody_REQUIRED':True,
 'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md')}
 (out/'report.json').write_text(json.dumps(d,separators=(',',':'))+'\n');print(sha(out/'report.json'))
if __name__=='__main__':main()
