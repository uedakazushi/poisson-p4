#!/usr/bin/env python3
"""Source-only full G3-linear kernel reduction; no new rank or Taylor2 query."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT
STEM='AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_INPUT_GATE_20261105_v63'
SOURCE='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'source_only_full_G3_zero3_split_no_rank'})
 hs={SOURCE:H.sha(ROOT/SOURCE)}
 peers=['reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
        'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json',
        'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json',
        'reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63_REVIEW/report.json']
 for p in peers:
  q=json.loads((ROOT/p).read_text());assert q['all_passed'];hs[p]=H.sha(ROOT/p)
 active=[0,1,4,5,7,8];zero=[2,3,6]
 columns=list(map(H.parse,H.selected(SOURCE,'actual_full_G3linear102_by418_columns')))
 assert len(columns)==418
 a=[]
 for c in columns[:238]:
  v={}
  for k,x in c.items():
   row,var=divmod(k,9);assert row<42 and var in active
   v[row*6+active.index(var)]=x
  a.append(v)
 pinned=H.selected(SOURCE,'actual_active_G3linear42_by238_columns')
 assert all({(i//9)*6+active.index(i%9):F(x) for i,x in c}==v for c,v in zip(pinned,a))
 # Each nonactive source has one distinct zero-variable target coordinate,
 # absent from every other source. This is a fixed diagonal source block.
 rows=[];signs=[]
 for iz,z in enumerate(zero):
  for j in range(42):rows.append(j*9+z);signs.append(F(1))
 pairs=list(itertools.combinations(range(3),2))
 for ia,ib in pairs:
  for j in range(17):rows.append((42+17*ib+j)*9+zero[ia]);signs.append(F(-1))
 for j in range(3):rows.append((93+3*pairs.index((1,2))+j)*9+zero[0]);signs.append(F(1))
 assert len(rows)==len(set(rows))==180
 for j,c in enumerate(columns):
  got={i:c[k] for i,k in enumerate(rows) if c.get(k)}
  assert got==({} if j<238 else {j-238:signs[j-238]}),('fixed_zero3_diagonal',j,got)
 input_nz=sum(map(len,a));input_bits=max([0]+[max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for c in a for x in c.values()])
 assert input_bits<=384
 # At rank r, keep r normalized252-vectors and 238 full238-source vectors
 # (r basis representations +238-r kernel columns). Do not store a second
 # dense factor matrix or duplicate operation trace.
 bounds={'source_active_columns':input_nz,'normalized_basis252_by238_upper':252*238,
  'basis_representation_plus_kernel_partition238_by238_upper':238*238,
  'pivot_metadata':4*238,'sparse_scratch_and_one_verification_product':4*(252+238),
  'retained_split_rows_signs_and_source_metadata':5000}
 whole=sum(bounds.values());assert whole<150000
 # A 384-bit cap is a fail-closed dynamic guard, not a predicted output bound.
 # Sparse scalar [index,"signed_n/d"] has index<=251; 24 extra bytes
 # cover its index, sign, slash, JSON punctuation, and record separator.
 raw_per_entry=len(str(2**384))*2+24
 raw=(252*238+238*238)*raw_per_entry+2_000_000
 assert raw<32*1024*1024
 proof={'full_source_target_order':'102 targets, each coefficients zA0,...,zA6,rA,ell',
  'active_source_first238':'the238 surviving active F3 source columns of the exact unit-transported resolution',
  'full_kernel_equals_active_kernel':'the displayed180x418 coordinate projection is [0,diag(+/-1)]; all nonactive constants are forced0',
  'full_rank_formula':'rank(full_G3_linear)=180+rank(active252x238)',
  'projected393_kernel_is_not_used':True,
  'quadratic_projection_sufficiency':'For a constant full source c with G3_1*c=0, G3_0=0 and derivative of the99-unit quotient acts on G3_1*c=0; final top quadratic equals PiY0*(G3_2*c), bottomzero. Coherent constant-unit normalization only changes the constant target frame.',
  'not_a_Taylor2_availability_claim':'G3_2 includes true active unit-transport derivatives. A separate source coefficient/budget gate must bind them before any quadratic values.',
  'rank_certificate_plan':'normalized252-basis plus its full238-source representations; every remaining source has one triangular free-index kernel column. Membership and all zero products are replayed; selected source independence follows pivot leading positions. No full factor-matrix duplicate.',
  'actual_T2_HF2_or_length_not_inferred':True}
 data={'active_coordinates6':active,'zero_variable_coordinates3':zero,'active_shape':[252,238],
  'actual_active_G3linear252_by238_columns':[srec(v) for v in a],
  'full_source418_nonactive_diagonal_projection_rows':rows,'full_source180_nonactive_diagonal_signs':list(map(str,signs)),
  'all418_fixed_diagonal_projection_identities':True,'source_proof':proof,
  'next_rank_bounds':bounds,'next_rank_whole_QQ_slot_bound':whole,'next_rank_raw_byte_bound':raw,
  'input_height_bits':input_bits,'rank_or_kernel_searched':False,'Taylor2_coefficients_queried':False,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_bytes':32*1024*1024,'QQ_slots':150000,'dynamic_bits_stop':384}}
 ep=out/'actual-full-G3-linear-zero3-split-and-active-kernel-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'rank_or_kernel_searched':False,'new_coefficients_queried':False,
  'active_shape':[252,238],'full_nonactive_diagonal_rank':180,'next_rank_whole_QQ_slot_bound':whole,
  'next_rank_raw_byte_bound':raw,'input_height_bits':input_bits,'input_sha256':hs,
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','active_shape','full_nonactive_diagonal_rank','next_rank_whole_QQ_slot_bound','next_rank_raw_byte_bound','input_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
