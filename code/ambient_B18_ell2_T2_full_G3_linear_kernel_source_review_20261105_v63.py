#!/usr/bin/env python3
"""Independent full-kernel zero3 source split and bounded-rank input review."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,time
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual-full-G3-linear-zero3-split-and-active-kernel-input.json',
'source':'reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json',
'source_peer':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json',
'transport_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
'minimal_peer':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k!='source'};assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in ['source_peer','transport_peer','minimal_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 S=fields(ROOT/P['source'],['actual_full_G3linear102_by418_columns','actual_active_G3linear42_by238_columns']);cols=list(map(sp,S['actual_full_G3linear102_by418_columns']));E=D['gate'];active=[0,1,4,5,7,8];zero=[2,3,6];ac=[]
 assert len(cols)==418
 for p in cols[:238]:
  q={}
  for k,c in p.items():i,v=divmod(k,9);assert i<42 and v in active;q[i*6+active.index(v)]=c
  ac.append(q)
 assert list(map(enc,ac))==E['actual_active_G3linear252_by238_columns']
 for p,q in zip(S['actual_active_G3linear42_by238_columns'],ac):assert {(i//9)*6+active.index(i%9):Q(c) for i,c in p}==q
 rows=[9*i+z for z in zero for i in range(42)];signs=[Q(1)]*126;pairs=list(itertools.combinations(range(3),2))
 for a,b in pairs:
  rows.extend((42+17*b+i)*9+zero[a] for i in range(17));signs.extend([Q(-1)]*17)
 rows.extend((93+3*pairs.index((1,2))+i)*9+zero[0] for i in range(3));signs.extend([Q(1)]*3)
 assert len(rows)==len(set(rows))==180 and rows==E['full_source418_nonactive_diagonal_projection_rows'] and list(map(str,signs))==E['full_source180_nonactive_diagonal_signs']
 for j,p in enumerate(cols):
  projected={i:p[k] for i,k in enumerate(rows) if p.get(k,Q(0))};assert projected==({} if j<238 else {j-238:signs[j-238]})
 expected={'source_active_columns':sum(map(len,ac)),'normalized_basis252_by238_upper':252*238,'basis_representation_plus_kernel_partition238_by238_upper':238*238,'pivot_metadata':4*238,'sparse_scratch_and_one_verification_product':4*(252+238),'retained_split_rows_signs_and_source_metadata':5000}
 assert expected==E['next_rank_bounds'] and sum(expected.values())==124906==E['next_rank_whole_QQ_slot_bound']<150000
 bits=max(max(abs(c.numerator).bit_length(),c.denominator.bit_length()) for p in ac for c in p.values());assert bits==28==E['input_height_bits']
 raw=(252*238+238*238)*(2*len(str(2**384))+24)+2_000_000;assert raw==31854720==E['next_rank_raw_byte_bound']<32*1024*1024
 assert E['limits']['dynamic_bits_stop']==384 and not E['rank_or_kernel_searched'] and not E['Taylor2_coefficients_queried']
 ep=o/'independent-full418-zero3-diagonal-and-active252-kernel-gate.json';dump(o,ep.name,{'inputs':inputs,'all418_fixed_projection_identities_verified':True,'actual_active_shape':[252,238],'fixed_nonactive_rank':180,
'full_kernel_and_rank_proof':'ProjectfullG3linear onto180displayedcoordinates. Itis[0,diag(+/-1)],soanyconstantfullkernelvectorhasitslast180coefficients0. Thefirst238columnsarethentheactive252matrixembeddedinthetarget; thusfullkernelisexactlyitskernelandfullrankis180+activerank.',
'quadratic_projection_conditional_proof':'TheactualG3constanttermis0. ForaconstantcintheFULLlinearkernel,theorder2partofPiY(R)G3(R)c equalsPiY0G3_2c, sincetheunknownPiY1multipliesfullG3_1c=0. Bottomzeroallordersremovesfirst27cancellationeffects; constantcoherentnormalizationdoesnotintroduceaderivative. No claimis made for thewithdrawnprojected393kernel.',
'no_internal_Taylor2_shortcut':'ActualactiveunittransportG3_2 muststillbesource-boundseparately; rawD3coefficientsorindividualunitderivativesarenotautomaticallysubstitutable.',
'bounded_rank_representation':'Atmost238normalized252vectors plus238total238-sourcevectorspartitionedintobasisrepresentationsandkernelvectors. No secondfactor matrix/tracecopyisbudgeted. Dynamic384bitguardmuststopandpreserveifexceeded.',
'whole_numeric_slot_bound':124906,'raw_byte_bound':31854720,'no_rank_kernel_or_new_coefficient_values_computed':True})
 report={'id':STEM,'all_passed':True,'full_G3_kernel_source_and_active_rank_gate_accepted':True,'actual_rank_kernel_or_Taylor2_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
