#!/usr/bin/env python3
"""Bound-only dense-method clarification within the unchanged old caps."""
import json,math,time
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_ADJOINT_DENSE_GUARD_20261105_v63'
def main():
 t=time.monotonic();prefix=ROOT/'reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_MIXED_ADJOINT_INPUT_20261105_v63';rp=prefix/'report.json';r=json.loads(rp.read_text());assert r['all_passed']
 p=prefix/'actual136-mixed-output210-slot-adjoint-source-input.json';d=json.loads(p.read_text());assert sha(p)==r['evidence_sha256'][p.name]
 n=len(d['matrix_columns_on210_domain']);m=len(d['full918_domain_slots210']);assert(n,m)==(91,210)
 categories={'decoded_source_and_metadata_plus_parse_overlap':12000,'dense_normalized_basis':n*m,'dense_original_source_representations':n*n,
 'all_pivot_step_records_and_column_metadata':n*(n-1)+3*n,'current_and_uncommitted_dense_vector_and_representation':2*(m+n),
 'two_final_dense_source_adjoint_vectors':2*n,'selected_original_and_pivot_indices':2*n,'streamed_one_row_serialization':2*(m+n),'small_phase_and_frontier_metadata':1500,'scalar_loop_and_verification_registers':640}
 owned=sum(categories.values());assert owned==51562<r['whole_solver_numeric_bound']==58603
 # At384bits, signed ordinary numerator/positive denominator strings use
 # at most235 characters.256 includes quotes/separators. Pending1153-bit
 # scalar, all labels/reports and result records fit the200k margin.
 qchars=2*math.ceil(384*math.log10(2))+5;assert qchars<=256
 raw=2*(n*m+n*n+2*(m+n)+2*n)*256+(n*(n-1)//2)*280+3*n*256+1_000_000+200_000
 assert raw==16842088<r['whole_solver_RAW_bound_at384']==16940016
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 report={'id':STEM,'all_passed':True,'bound_only_no_signed_elimination':True,'method':'dense coefficient vectors; original sparse columns retained by reference; row-streamed single atomic frontier plus one append-only pivot history; no repeated saved fullfrontiers',
 'source_metadata_numeric_cap_before_arithmetic':6000,'source_RAW_input_cap_before_arithmetic':1000000,'numeric_categories':categories,'proven_whole_numeric':owned,'unchanged_numeric_cap':58603,
 'proven_whole_RAW':raw,'unchanged_RAW_cap':16940016,'committed_dynamic_bit_cap':384,'uncommitted_scalar_bit_cap':1153,'no_completed_elimination_height_bound_claimed':True,
 'old_formula_clarification':'Earlier sparse58603 expression omitted some numeric keys/history in an unspecified method. It is not reused as a sparse-method proof. The concrete dense method above fits its unchanged58603 cap; no solve ran before this clarification.',
 'input_sha256':{str(rp.relative_to(ROOT)):sha(rp),str(p.relative_to(ROOT)):sha(p)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'owned':owned,'RAW':raw,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
