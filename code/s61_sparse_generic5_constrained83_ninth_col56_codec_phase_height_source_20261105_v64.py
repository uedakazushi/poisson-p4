#!/usr/bin/env python3
"""SOURCE-only codec-phase integer envelopes and unchanged RAW caps."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,hashlib,time,signal,resource,re
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_CODEC_PHASE_HEIGHT_SOURCE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def scalars(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(scalars,v))
 if isinstance(v,dict):return sum(map(scalars,v.values()))
 return 0
def ub(n):
 assert n>=0;return max(1,(n.bit_length()+6)//7)
def bb(n):return(n.bit_length()+7)//8
def digits(n):return 1 if n==0 else(n.bit_length()*30103+99999)//100000+1

def main(a):
 assert a.freeze
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')}
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
 def evidence(p,r,name):
  f=p.parent/name;assert sha(f)==r['evidence_sha256'][name];inputs[f]=sha(f);return read(f)
 pp,pr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64','e332bea912ec182d8c6f7ed614a1cb7a25cc0dda4ab13f5d6033c85c9200b5ab');program=evidence(pp,pr,'staged-literal-col56-program-and-positive-profiles.json');P=program['positive_support_and_height_profiles']
 gp,gr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_IMPLEMENTATION_GATE_V2_20261105_v64','c1dab9549187c3161a46393b909d426e331285a695c56da4727bddfe26a9a7fa');proof=evidence(gp,gr,'implementation-custody-streaming-lifetime-proof.json');assert proof['scalar_content_unreduced_temporary_bit_upper']==174 and proof['all_integer_coefficient_positive_bound_bits']==272
 worker=ROOT/'code/s61_sparse_generic5_constrained83_ninth_col56_staged_coefficients_v2_20261105_v64.py';codec=ROOT/'code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py';assert sha(worker)==proof['worker_sha256']and sha(codec)==proof['stream_codec_sha256'];inputs[worker]=sha(worker);inputs[codec]=sha(codec)
 cg,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STREAM_CODEC_SOURCE_GATE_20261105_v64','c9921d8f8db77f5a4cf2739b26d88a4728cb871c4e58b4ed9aa748a6d09e2d6f');assert cr['all4160_existing_coefficients_wire_identical']and cr['all16_existing_and2_fixtures_stream_readback']
 allq=[(k,Q(z['known_QQ_content']))for k,z in P.items()];allnum=max(abs(q.numerator).bit_length()for k,q in allq);allden=max(q.denominator.bit_length()for k,q in allq);allM=max(int(z['primitive_integer_abs_upper']).bit_length()for z in P.values());assert max(allnum,allden)<=174 and allM==272
 names=[('shared35',program['literal_shared_N8_35_56_node'])]+[(str(z['row']),z['node'])for z in program['row15_output_nodes']];assert len(names)==16;rows=[];codecmax=0;headmax=0;indexmax=0;coordinate_max=0;hashbytes=0
 for row,name in names:
  p=P[name];q=Q(p['known_QQ_content']);M=int(p['primitive_integer_abs_upper']);n=p['terms_upper'];assert q and M>0 and n>0;A=abs(q.numerator);D=q.denominator;Gmax=M;Cmax=M
  # Writer q*G; reader header Fraction; writer hash q*z; reader LEFT-TO-RIGHT
  # q*c*G. M bounds integers in the resident map; primitive c and G<=M.
  qGnum=A*Gmax;hash1=A*Cmax;hash2=hash1*Gmax;cross=qGnum*D;CG=Cmax*Gmax;nodebits=max(qGnum.bit_length(),hash2.bit_length(),cross.bit_length(),CG.bit_length(),D.bit_length(),M.bit_length());codecmax=max(codecmax,nodebits)
  box=p['exponent_box5'];shift=[a for a,b in box];spans=[b-a for a,b in box];assert all(v>=0 for v in spans);zig=[2*a if a>=0 else-2*a-1 for a in shift];coordabs=[max(abs(a),abs(b))for a,b in box];sumabs=sum(coordabs);heapindex=2*n+1;indexbits=heapindex.bit_length();coordbits=max([v.bit_length()for v in coordabs+zig+spans]+[sumabs.bit_length()]);indexmax=max(indexmax,indexbits);coordinate_max=max(coordinate_max,coordbits);nodebits=max(nodebits,indexbits,coordbits,8)
  Lnum=bb(qGnum);Lden=bb(D);Lcoef=bb(M);header=9+ub(2*Lnum+1)+Lnum+ub(2*Lden+1)+Lden+sum(map(ub,zig))+ub(n);expbytes=sum(map(ub,spans));term=Lcoef+ub(2*Lcoef+1)+expbytes;cap=proof['output_caps'][row];sourcecap=2048+Lden+n*(Lnum+expbytes+16);wire=header+n*term
  assert sourcecap==cap and header<=2048+Lden and term<=Lnum+expbytes+16 and wire<=cap
  # Normalization/GCD never enlarges the numerator or denominator. Even if
  # Fraction equality used cross multiplication it stays under cross above.
  vhbytes=digits(hash2)+1+digits(D)+1+sum(digits(v)+1 for v in coordabs)+6;hashbytes=max(hashbytes,vhbytes);headmax=max(headmax,header)
  assert nodebits<2048 and vhbytes<2048
  rows.append({'row':row,'node':name,'map_terms_upper':n,'map_integer_bound_bits_unchanged':M.bit_length(),'known_q_numerator_bits':A.bit_length(),'known_q_denominator_bits':D.bit_length(),'G_and_primitive_c_bits_upper':M.bit_length(),'write_q_times_G_raw_numerator_bits_upper':qGnum.bit_length(),'write_q_times_map_c_raw_numerator_bits_upper':hash1.bit_length(),'validate_LEFT_TO_RIGHT_q_times_c_times_G_raw_numerator_bits_upper':hash2.bit_length(),'c_times_G_integer_bits_upper':CG.bit_length(),'normalization_or_comparison_cross_bits_upper':cross.bit_length(),'all_codec_integer_temporary_bits_upper':nodebits,'header_bytes_upper':header,'coefficient_record_bytes_upper':term,'full_wire_bytes_upper':wire,'unchanged_gate_RAW_cap':cap,'header_absorbed_by_existing2048_plus_den_bytes':True,'valuehash_text_buffer_bytes_upper':vhbytes,'heap_index_bits_upper':indexbits,'signed_exponent_and_degree_sum_bits_upper':coordbits})
 exactEncodedNum=max(z['known_q_numerator_bits']for z in rows);exactEncodedDen=max(z['known_q_denominator_bits']for z in rows);assert exactEncodedNum==1 and exactEncodedDen==82
 broad={'DAG_content_unreduced_bound_unchanged':174,'integer_map_bound_unchanged':272,'write_qG_unreduced_num_bits':174+272,'validate_qcG_naive_unreduced_num_bits':174+272+272,'Fraction_comparison_naive_cross_bits':174+272+174,'normalization_denominator_bits':174,'uniform_CODEC_phase_envelope_bits':2048};assert max(broad['write_qG_unreduced_num_bits'],broad['validate_qcG_naive_unreduced_num_bits'],broad['Fraction_comparison_naive_cross_bits'])<2048
 # Separate existing-data SOURCE fixture is not the future map272 assertion.
 fixtureQ=2**293+1;fixtureC=2**1651+1;fixtureG=1;fixturebits=(fixtureQ*fixtureC*fixtureG).bit_length();assert fixturebits<2048
 operation_proof=[
  'Integer gcd for primitive G uses only map integers bounded by M; every Euclidean remainder and G<=M.',
  'write_cell q*G has numerator<=abs(q.num)*M, denominator dividing q.den. Fraction normalization/GCD cannot increase either bound.',
  'sint/uint abs, sign, byte length and canonical parsing stay inside each field bound; uint masks/shifts do not enlarge integer bit length. Header numerator, denominator and primitive coefficient signed tags fit the explicit lengths.',
  'write valuehash uses q*map_coefficient. validate valuehash evaluates (q*c)*G left-to-right; both raw multiplication envelopes are recorded, even before Fraction cross reduction.',
  'validate constructor Q(num,den), reduced header checks and equality normalization use integers within header bounds; even naive numerator/denominator cross products fit the recorded comparison envelope.',
  'Primitive c=map_coefficient//G is exact, c*G reconstructs the map integer. A deliberately looser independent c<=M,G<=M envelope also covers their temporary product.',
  'Alias-key heapsort indices are<=2n+1; each comparison computes only two sums of five signed coordinates. Zigzag shifts and relative exponents are bounded by the recorded boxes and spans.',
  'For actual frozen writer files, incremental canonical readback reconstructs exactly these bounded fields. Empty output uses the nine-byte zero path and never constructs q*G.',
  'No term/value/decode arrays or new scalar-coordinate maps are introduced; the 20k metadata/130791 scalar-count bound is unchanged. Bit width is a distinct phase bound.'
 ]
 evidence={'all_passed':True,'source_only':True,'all16_codec_phase_rows':rows,'all700_stored_q_exact_max_bits':{'numerator':allnum,'denominator':allden},'encoded16_stored_q_exact_max_bits':{'numerator':exactEncodedNum,'denominator':exactEncodedDen},'all_codec_raw_integer_envelope_bits_exact_profile_max':codecmax,'uniform_broad_envelope':broad,'max_header_bytes_upper':headmax,'all16_original_header_RAW_caps_unchanged':True,'max_heap_index_bits':indexmax,'max_coordinate_degree_bits':coordinate_max,'max_valuehash_text_buffer_bytes':hashbytes,'separate_existing_SOURCE_wire_fixture_raw_hash_bits':fixturebits,'operation_by_operation_proof':operation_proof,'worker_codec_gate_and_all_map_DAG_scalar_bounds_NOT_modified':True,'no_new_coefficient_job_or_values':True}
 (out/'codec-phase-height-and-existing-RAW-cap-proof.json').write_bytes(enc(evidence));owned=25000+scalars(program)+scalars(proof)+scalars(rows)+scalars(evidence)+4096;assert owned<150000;raw=sum(p.stat().st_size for p in inputs)+sum(p.stat().st_size for p in out.iterdir());assert raw+65536<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'all16_codec_qG_qcG_normalization_and_RAW_checks':True,'all700_q_num_max_bits':allnum,'all700_q_den_max_bits':allden,'encoded16_q_num_max_bits':exactEncodedNum,'encoded16_q_den_max_bits':exactEncodedDen,'exact_profile_codec_temporary_bits_upper':codecmax,'broad_codec_temporary_bits_upper':718,'uniform_CODEC_phase_envelope_bits':2048,'DAG174_and_map272_unchanged':True,'all16_existing_RAW_caps_unchanged':True,'new_matrix_coefficients':0,'SOURCE_owned_upper':owned,'actual_READ_plus_evidence_RAW':raw,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
