"""Metadata-only clarification of frozen direct peer height and RAW bounds."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bits(v):
 if isinstance(v,list):return max(map(bits,v),default=0)
 if isinstance(v,dict):return max(map(bits,v.values()),default=0)
 q=F(v);return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def main():
 peer=ROOT/'reports/AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json'
 assert sha(peer)=='10844d8720465054d6b9f7a7aa408dd95dd613a2e857eb52e50cd5a71a847e92'
 r=json.loads(peer.read_text());pins={str(peer.relative_to(ROOT)):sha(peer),**r['input_sha256']}
 for p,h in pins.items():assert sha(ROOT/p)==h,p
 src=[p for p in pins if p.endswith('factored-native23-Gamma-pivot-plane2-source-and-full-L2S0.json')];assert len(src)==1
 d=json.loads((ROOT/src[0]).read_text());names=['cached_N_12_by12_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE','canonical_G8_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE','canonical_nu_two_plane5','Gamma_two']
 input_bits={k:bits(d[k])for k in names};assert max(input_bits.values())==152
 source_peak=280;assert max(input_bits.values(),default=0)<source_peak and r['actual_max_bits']==112
 producer=ROOT/'code/ambient_B18_native92_ROOT_direct_review_20261105_v63.py';assert sha(producer)==r['code_sha256'];pins[str(producer.relative_to(ROOT))]=sha(producer)
 raw_read=sum((ROOT/p).stat().st_size for p in r['input_sha256'])+producer.stat().st_size
 # All VALUE/guard output contains <1k fractions, each <=280 bits, plus fixed
 # source IDs/cursors and four input hashes: <1MiB even without cancellations.
 raw_prior=raw_read+1024*1024;assert raw_prior<4_882_000
 out=ROOT/'reports/AMBIENT_B18_NATIVE92_ROOT_METADATA_20261105_v63';out.mkdir(exist_ok=False)
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 rep=dict(id=out.name,all_passed=True,metadata_only_no_arithmetic_response_rerun=True,frozen_direct_peer_unchanged=True,
 actual_max_bits_field_means_new_arithmetic_intermediates_ONLY=112,loaded_mathematical_input_peak_bits=152,input_field_bits=input_bits,
 whole_mathematical_input_and_intermediate_peak=152,unchanged_source_committed_bits=280,unchanged_temporary_bits=841,unchanged_numeric_prior=12432,
 literal_input_and_code_RAW_bytes=raw_read,whole_RAW_source_prior=raw_prior,unchanged_old_RAW_limit=4_882_000,
 source_output_bound='92+322 table or guard prefix, <1000 fractions at <=280bits plus fixed cursor/hash metadata; all generated files <1MiB',input_sha256=pins)
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':main()
