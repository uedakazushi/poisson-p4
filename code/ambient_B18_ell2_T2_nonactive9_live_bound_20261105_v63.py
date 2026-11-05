#!/usr/bin/env python3
"""Metadata/count proof only: no lower product, minor, or rank evaluation."""
import argparse,json,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_LIVE_BOUND_20261105_v63'
def main(out,peer=None):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/'reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63/report.json';g=json.loads(gp.read_text());ip=gp.parent/'exact171-decomposition-162-zero-9-lower-first-source-input.json';d=json.loads(ip.read_text());assert H.sha(ip)==g['evidence_sha256'][ip.name]
 rows=[
 ['temporary_raw_JSON_copy_of_source_input',count(d),'One temporary read_text/loads overlap; count all source numeric records again.'],
 ['parsed_Hstar_r_and_ell',2*99*3*2,'Two lists, at most3 sparse entries per selected99 column; index and coefficient counted.'],
 ['current_parsed_b1_r_and_ell',2*99*2,'Only the current source is parsed; two sparse99 vectors.'],
 ['literal_and_integer_matrices_plus_all_rank_arrays',2*(2*9*4+3*4*4+9*4),'Two9-column matrices; original basis, normalized basis and source representations rank<=4; all9 factors.'],
 ['all_source_kernel_vectors',2*9*5,'At most9 free vectors, each one free entry and at most4 selected coefficients.'],
 ['rank_receipts_and_partial_image_accumulators',200,'At most4 receipts/6 prior pivots; two4-row image accumulators and axpy copy.'],
 ['gate_peer_custody_and_all_small_metadata',900,'Gate+peer numerical metadata must count<=900 before arithmetic. No evidence matrix is loaded from the peer.'],
 ['scalar_registers_loops_and_determinant_scratch',320,'At most24 determinant terms streamed; no determinant table; current dot accumulators.'],
 ['complete_output_object_and_serialization_overlap',1500,'Includes all36+36 matrix entries, selected factors/receipts/kernels and duplicated record encoding, report/phase metadata.']]
 scratch=sum(r[1] for r in rows);assert scratch<9000
 decoded=count(d);whole=decoded+scratch;assert whole<=decoded+9000<=g['future_lower_owned_bound']==11408
 # Input added exactly three numeric bound fields after the source producer's
 # count. They are now included literally in decoded, not removed or ignored.
 assert all(k in d for k in ['future_lower_value_owned_bound','source_owned_bound','future_lower_value_RAW_bound'])
 metadata=count(g);peercheck=None
 if peer:
  pp=ROOT/peer;p=json.loads(pp.read_text());assert p['all_passed'] and H.sha(gp) in json.dumps(p);metadata+=count(p);assert metadata<=900;peercheck={'path':peer,'sha256':H.sha(pp),'gate_plus_peer_numeric':metadata}
 report={'id':STEM,'all_passed':True,'bound_only_no_signed_values':True,'decoded_input_including_all_three_metadata_fields':decoded,'scratch_categories':rows,'proven_scratch_bound':scratch,'worker_scratch_allowance':9000,'proven_whole_live_bound':whole,'unchanged_authorized_whole':11408,'metadata_gate_count':count(g),'mandatory_metadata_count_cap_before_values':900,'actual_peer_count_check':peercheck,'parent_value_GO_and_FINAL_peer_still_required':peer is None,
  'input_sha256':{str(gp.relative_to(ROOT)):H.sha(gp),str(ip.relative_to(ROOT)):H.sha(ip)},'code_sha256':H.sha(Path(__file__)),'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['decoded_input_including_all_three_metadata_fields','proven_scratch_bound','proven_whole_live_bound','unchanged_authorized_whole']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer');a=p.parse_args();main(a.output.resolve(),a.source_peer)
