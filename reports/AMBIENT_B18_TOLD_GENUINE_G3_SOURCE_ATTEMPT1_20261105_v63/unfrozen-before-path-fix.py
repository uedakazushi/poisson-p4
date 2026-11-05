#!/usr/bin/env python3
"""Source-only gate for the actual Told-to-genuine-connecting target."""
import hashlib,json,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOLD_GENUINE_G3_SOURCE_20261105_v63'
PATHS={
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'naturalpeer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
 'naturalnote':'research/notes/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63.md',
 'head4':'reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63/actual-head4-z5-map-and-F2-kernel-coefficient.json',
 'head4peer':'reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63_REVIEW/report.json',
 'head4code':'code/ambient_B18_Told_Tnil_head4_z5_values_20261105_v63.py',
 'edge':'reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63/actual-five-edge-functionals-and-genuine-Nbar14-cycles.json',
 'edgepeer':'reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'G3peer':'reports/AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63_REVIEW/report.json',
 'finite':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/report.json',
 'finitepeer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'Loewy':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'connectingpeer':'reports/AMBIENT_B18_NIL_CONNECTING16_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k not in ['naturalnote','head4code']}
 assert all(v['all_passed']for k,v in D.items()if k.endswith('peer')or k in ['finite','Loewy'])
 M=D['natural']['actual_origin_data']['actual_origin_map_columns5'];assert len(M)==5 and all(len(c)==5 and c[0]=='0'for c in M)
 assert D['edge']['genuine_Nbar_dimension']==14
 assert D['finite']['actual_ell_regular_whole_E4_and_Ered']and D['finite']['actual_normalized_nil_DVR_free_rank']==9
 assert D['finite']['actual_NilD0_noncanonical_DVR_split']=='K plus D^7'
 key='actual_projected_natural_head4_first_target0_r';assert key in D['head4']
 # It is the SAME actual natural column4/currentminimal head0 derivative,
 # not its z5-product image, nor the unknown canonicallong edge beta.
 code=(ROOT/PATHS['head4code']).read_text();assert "'actual_projected_natural_head4_first_target0_r':Jr"in code
 b=D['edge']['five_actual_cached_functional_values_a_b1_b2_b3_b4'][1:];xi=D['edge']['xi_equals_minus11a_by16']
 data={
  'same_actual_modules':'S=R0/(ell); Told,Tnil from the two accepted same-D0 comparisons; genuineG3=Cprime[ell].',
  'all_five_natural_origin_columns5':M,'all_five_head0_origin_entries_zero':True,
  'actual_genuine_target':'basis(K,Q1), r*K=-(11/16)*Q1, allz andell kill, m²G3=0',
  'rank_at_most1_proof':'The K coefficient of genuineedge is a timesminimalhead0 origin. It vanishes on allfive natural generators, hence onTold. Image is contained k*Q1; its R0rank means finitek-image dimension, notDVRgenericrank.',
  'one_cached_witness':{'old_generator_index':4,'old_actual_free56_label':55,'natural_origin_column5':M[4],'genuine_edge_b1_b2_b3_b4':b,'genuine_edge_rh0_coefficient_xi':xi,'natural_head0_r_first_callback_path':PATHS['head4'],'natural_head0_r_first_callback_key':key},
  'exact_witness_formula':'sum(j=1..4) b_j*M0[j,4] + xi*cachedJr; K coefficient0. Canonical beta of edge(h0) is irrelevant because M0[0,4]=0. Allzfirsts andallhigherterms are killed by genuineG3.',
  'cache_scope':'Jr was already reconstructed and independently reviewed with the actual oldn1/18unit/81unit cleanup in HEAD4_Z5_VALUES. Reuse this saved first coefficient alone, not the quadratic productclass orR35 correction. No new naturalfirst/source/gauge coefficient is required.',
  'source_diagram':{
   'full_complexes':'Lold=[D0→E4], Lnew=[D0→Ered] with identity onD0 and actual E4 quotientEred. Kernel is NE=Nil(E4) in the targetterm.',
   'ell_fiber_exactness':'Both E4,Ered andNE are ellregular, NE Dfree9. Derivedellfiber gives0→Told→Tnil→NE/ellNE→C/ellC→Cprime/ellCprime→0. CommonTor1(D0) is retained/cancels by identity; no H2 inference is made.',
   'genuine_square':'The injection N/ellN→Tnil from the actual nilcomparison followed by Tnil→NE/ellNE is exactly reduction of the actual D0→E4 map onnil. Hence itskernel is genuineNbar intersectimageTold.',
   'conditional_rank1_dimensions':'If the witness is nonzero, δ′|Told rank1 and Told length11 giveintersectionlength10. Nbar14 then hasactualimage inNEbar ofdimension4.',
   'lattice_consequence':'I=image(N→NE)=kernel(NE→C) is the accepted Dfree7. Nbar→Ibar is surjective since N/K=I is Dfree; so Ibar→NEbar hasrank4. For L=coker(I→NE), Dfree resolutions0→D7→D9→L→0 giveDgenericrank2 andTor1^D(L,k)dimension7-4=3. Thus Ltors has3minimalDgenerators, notknownlength/exponents.',
   'notation_warning':'L here is coker(I→NE), not the earlier transfer512 subideal often namedP; no coverage/saturation identity is inferred.'},
  'future_job_profile':{'cached_scalar_products':5,'static_owned_numeric_slots':15000,'hard_owned_numeric_slots':30000,'source_derived_bit_bound':512,'hard_bit_limit':2048,'seconds':60,'memory_bytes':2<<30,'input_plus_output_bytes':2<<20},
  'witness_pairing_rank_and_dimensions_not_evaluated':True,'no_new_source_gauge_jet_or_full_nil_lifts':True}
 ev=out/'actual-old-to-genuine-edge-one-cache-witness-and-diagram.json';ev.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_GATE_actual_Told_edge_rank_atmost1_and_one_existing_cached_witness','actual_edge_image_rank_upper':1,'witness_evaluated':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
