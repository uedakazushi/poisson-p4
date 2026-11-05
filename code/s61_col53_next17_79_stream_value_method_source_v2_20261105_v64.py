"""SOURCE proof/forecasts only for streamed coefficient evaluation. No signed values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_STREAM_VALUE_METHOD_SOURCE_V2_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
SOURCE='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64'
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
if __name__=='__main__':
 start=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def read(f,expected=None):
  p=ROOT/f;h=sha(p)
  if expected:assert h==expected
  pins[f]=h;return json.loads(p.read_text())
 r=read(SOURCE+'/report.json','8cba44d77933886756700a63fb250058256db3310272b6ca626936116c785d60');assert r['all_passed']and r['SOURCE_support_only_no_signed_values']
 f=read(FINAL+'/report.json','548ebf5f4078707d6ee2bb6dd1867c7161c2d32bd553ebef4ce454600c7db877')
 D=read(FINAL+'/compiled-full-SLP-and-fixed-projection-plan.json',f['evidence_sha256']['compiled-full-SLP-and-fixed-projection-plan.json'])
 countsfile=SOURCE+'/all-projected-SOURCE-node-counts.json';C=dict(read(countsfile,r['evidence_sha256']['all-projected-SOURCE-node-counts.json']))
 nodes=D['nodes'];order=sorted(C);N=len(C);assert N==1884 and sum(C.values())==71584 and max(C.values())==129
 weights=[];maxM=0;rawbits=0;failures=[];maxqnum=0;maxqden=0
 for k in order:
  n=nodes[k];qi=Q(n['q']);Mi=int(n['M']);args=n['args'];maxM=max(maxM,Mi.bit_length());maxqnum=max(maxqnum,abs(qi.numerator).bit_length());maxqden=max(maxqden,qi.denominator.bit_length())
  w=[];bound=0
  if not qi or not C[k]:pass
  elif n['op']=='input':bound=Mi
  elif n['op']=='copy':
   a=args[0];z=Q(nodes[a]['q'])*Q(n['scalar'])/qi;assert z.denominator==1;w=[z.numerator];bound=abs(z.numerator)*int(nodes[a]['M'])
  elif n['op']=='multiply':
   a,b=args;z=Q(nodes[a]['q'])*Q(nodes[b]['q'])*Q(n['scalar'])/qi;assert z.denominator==1;w=[z.numerator];raw=int(nodes[a]['M'])*int(nodes[b]['M'])*abs(z.numerator);rawbits=max(rawbits,raw.bit_length());bound=min(C[a],C[b])*raw
  else:
   assert n['op']=='add'
   for j,a in enumerate(args):
    z=Q(nodes[a]['q'])*(Q(n['scalar'])if j else Q(1))/qi;assert z.denominator==1;w.append(z.numerator)
   bound=sum(abs(z)*int(nodes[a]['M'])for z,a in zip(w,args))
  if bound>Mi:failures.append([k,bound.bit_length(),Mi.bit_length()])
  weights.append({'node':k,'integer_weights':list(map(str,w)),'positive_partial_bound_le_M':bound<=Mi})
 assert not failures,failures[:5]
 observed_needed_maxima=[maxM,maxqnum,maxqden]
 assert maxM<=824 and maxqnum<=28 and maxqden<=164
 # Preserve the certified global2998 limits; the needed1884 subset is smaller.
 maxM,maxqnum,maxqden=824,28,164
 # Exact unique metadata slots with node-cell segments only, no full G/content
 # resident table; qStored/qNode is checked each time a child is loaded.
 meta_parts={'runtime':25000,'graph_numeric':65701,'parsed_mode_numeric_keys':2*N,'order':N,'source_count_dict':2*N,'thin_segment_index_key_offset_bytes_terms':4*N,'ROOT_gate_records_and_boundaries':4096}
 meta_exact=sum(meta_parts.values());assert meta_exact==111753
 meta=118000
 forecasts=[];peak=0;at=None
 for k in order:
  n=nodes[k];children=list(dict.fromkeys(n['args']if n['op']in['copy','multiply']else D['active_tied_HL_children'][str(k)]))
  childterms=sum(C[c]for c in children);ni=C[k]
  evalphase=meta+7*(childterms+ni)+4096
  loadphase=meta+7*childterms+8*max([C[c]for c in children]+[0])+4096
  inputphase=meta+16152+7*ni+4096 if n['op']=='input'and n['q']!='0'else 0
  codecphase=meta+7*(childterms+ni)+3*ni+4096
  whole=max(evalphase,loadphase,inputphase,codecphase)
  assert whole<150000
  if whole>peak:peak=whole;at=k
  forecasts.append({'node':k,'source_terms':ni,'unique_child_terms':childterms,'whole_numeric_upper':whole})
 assert peak==138535,(peak,at)
 # Existing Q5 cell codec computes G=gcd(abs(ci)), shift=min coordinates. It
 # changes stored primitive content to qi*G; fixed source qi never changes.
 # G<=Mi (824 bits). Decode primitive di*G -> fixed-qi ci; ratio is integral.
 codec={'fixed_q_NOT_assumed_primitive':True,'G_gcd_abs_fixed_q_coefficients_bits_le':maxM,'qi_G_numerator_bits_le':maxqnum+maxM,'qi_G_denominator_bits_le':maxqden,'naive_valuehash_qi_di_G_numerator_bits_le':maxqnum+2*maxM,'decoded_di_times_G_bits_before_final_bound':2*maxM,'qStored_div_qi_naive_cross_bits_le':maxqnum+maxM+maxqden,'all_temporary_bits_le_1921_with_hard2048':True,'fixed_q_restore_rule':'G=qStored/qi positive integer; ci=di*G; verify abs(ci)<=Mi; do NOT use qStored for recurrences','empty_cell_rule':'empty map stores flag0/content0; no division by qi=0'}
 assert max(codec['naive_valuehash_qi_di_G_numerator_bits_le'],codec['decoded_di_times_G_bits_before_final_bound'],codec['qStored_div_qi_naive_cross_bits_le'])<=1921<2048
 onegen=512*N+116*sum(C.values());assert onegen==9268352
 storage={'one_codec_generation_bytes_upper':onegen,'TWO_byte_identical_generations_bytes_upper':2*onegen,'manifest_line_RAW_max':2048,'manifest_RAW_upper':2048*N,'root_report_and_control_RAW_reserve':524288,'stop_partial_codec_and_stream_tail_RAW_reserve':1048576,'ALL_generated_RAW_upper':2*onegen+2048*N+524288+1048576,'temporary_serializations_preserved_under_work_AND_exactbytes_in_immutable_stream':True,'all_generated_versions_counted_not_gzip':True}
 assert storage['ALL_generated_RAW_upper']==23968000
 for p in[Path(__file__),ROOT/'code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'reports/S61_COL53_NEXT17_79_STREAM_VALUE_METHOD_SOURCE_20261105_v64/failure.json',ROOT/'reports/S61_COL53_NEXT17_79_STREAM_VALUE_METHOD_SOURCE_20261105_v64/failed-source.py']:pins[str(p.relative_to(ROOT))]=sha(p)
 data={'SOURCE_only_no_signed_coefficient_job':True,'observed_needed_field_maxima':observed_needed_maxima,'old_full_vs_subset_metadata_assertion_preserved':True,'integer_q_normalization_and_partial_majorant_checks':weights,'metadata_exact_parts':meta_parts,'metadata_exact_total':meta_exact,'future_fixed_metadata_base':meta,'all1884_unique_child_and_codec_phase_forecasts':forecasts,'whole_VALUE_forecast':peak,'peak_node':at,'codec_preconditions_and_bounds':codec,'storage':storage,'full5keys_plus_coeff_7_per_mapterm':True,'allresident_186939_method_NOT_authorized':True,'next_immutable_VALUE_worker_code_and_ROOT_review_REQUIRED':True}
 (out/'streamed-residency-height-codec-storage-source-proof.json').write_bytes(enc(data))
 report={'id':ID,'all_passed':True,'SOURCE_METHOD_ONLY_no_signed_values':True,'nodes':N,'positive_term_instances':sum(C.values()),'max_single_terms':max(C.values()),'whole_VALUE_forecast':peak,'generated_RAW_forecast':storage['ALL_generated_RAW_upper'],'all_integer_weights_and_partial_majorants_checked':True,'codec_fixed_q_to_primitive_G_and_restore_proof':True,'immutable_VALUE_worker_and_fresh_ROOT_GO_STILL_REQUIRED':True,'input_sha256':pins,'evidence_sha256':{'streamed-residency-height-codec-storage-source-proof.json':sha(out/'streamed-residency-height-codec-storage-source-proof.json')},'seconds':time.monotonic()-start}
 (out/'report.json').write_bytes(enc(report));print('report_sha256',sha(out/'report.json'));print('whole/generated',peak,storage['ALL_generated_RAW_upper'])
