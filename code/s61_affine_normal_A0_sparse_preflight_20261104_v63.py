#!/usr/bin/env python3
"""Read-only affine normal-matrix sparsity and constant-block normalization."""
import argparse, collections, gzip, hashlib, json, math, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_AFFINE_NORMAL_A0_SPARSE_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s read-only cap')));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False)
 src=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz';prp=src.parent/'report.json';pr=json.loads(prp.read_text());assert pr['all_passed']
 for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in pr['evidence_sha256'].items():assert sha(src.parent/rel)==h,rel
 data=json.loads(gzip.decompress(src.read_bytes()));A=data['all129_linear_factor_rows'];assert len(A)==129 and all(len(r)==42 for r in A)
 certp=src.parent/'residual42-certificate.json';cert=json.loads(certp.read_text());selected=data['selected_residual_rows_zero_based'];assert selected==cert['selected_residual_rows_zero_based']==list(range(42))
 knownp=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json';known=json.loads(knownp.read_text())['known_five_mod_32003'];assert known[0]==[1,0,0,0]and cert['known_five_determinants_mod_prime'][0]==1680
 priorp=ROOT/'reports/S61_AFFINE_NORMAL_LEFTINVERSE_ASSESSMENT_20261103_v63.json';prior=json.loads(priorp.read_text());assert prior['all_passed']and prior['D1_leftinverse_excluded_in_frozen_F32003_fiber']and not prior['new_D1_solver_launched']
 for rel,h in prior['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 counts=collections.Counter();variables=[0]*45;bits=[0,0];rows=[];cols=[collections.Counter()for _ in range(42)];active=set();a0=[];bad_den=0;max_entry=0
 multipliers=[()]+[(v,)for v in range(45)]
 for i,row in enumerate(A):
  rc=collections.Counter();c0=[]
  for j,p in enumerate(row):
   seen=set();has0='0';counts['entries']+=1;counts['nonzero_entries']+=bool(p);rc['nonzero_entries']+=bool(p);max_entry=max(max_entry,len(p))
   for m,c in p:
    mon=tuple(m);q=Q(c);assert len(mon)<=1 and mon not in seen and q;seen.add(mon)
    if mon:assert 0<=mon[0]<45 and(mon[0]<3 or mon[0]>=j+3)
    else:has0=str(q)
    kind='constant'if not mon else'plane_z'if mon[0]<3 else'normal_n';counts[kind]+=1;rc[kind]+=1;cols[j][kind]+=1
    if mon:variables[mon[0]]+=1
    counts['terms']+=1;rc['terms']+=1
    bits=[max(bits[0],abs(q.numerator).bit_length()),max(bits[1],q.denominator.bit_length())];bad_den+=q.denominator%32003==0
    for mult in multipliers:active.add((j,tuple(sorted(mon+mult))))
   c0.append(has0)
  rows.append(dict(rc));a0.append(c0)
 assert bad_den==0 and counts['terms']==41249
 active_by_degree=[sum(len(mon)==k for j,mon in active)for k in range(3)]
 unknowns=129*46;allconstraints=42*math.comb(47,2);nnz=counts['terms']*46
 assert unknowns==5934 and allconstraints==45402
 dense_raw=allconstraints*unknowns;dense_active=len(active)*unknowns
 block={'selected_original_residual_rows':selected,'coordinate_order':'z0,z1,z2,n0,...,n41; centered deltaq coordinates','constant42x42':[[a0[i][j]for j in range(42)]for i in selected],'inherited_exact_32003_unit_minor':1680,'inherited_minor_source':str(certp.relative_to(ROOT)),'center_plane_coordinates':[0,0,0],'new_inverse_or_determinant_not_computed':True}
 raw=json.dumps(block,separators=(',',':')).encode();(out/'constant42-block.json.gz').write_bytes(gzip.compress(raw,mtime=0))
 result={'source_shape':[129,42],'source_gzip_bytes':src.stat().st_size,'source_json_bytes':len(gzip.decompress(src.read_bytes())),
 'source_counts':dict(counts),'max_entry_terms':max_entry,'row_min_max_terms':[min(r['terms']for r in rows),max(r['terms']for r in rows)],
 'max_coefficient_numerator_denominator_bits':bits,'all_coefficients_32003_integral':True,'first_normal_triangular_support_verified':True,
 'coefficient_counts_by_variable':variables,'row_statistics':rows,'column_statistics':[dict(c)for c in cols],
 'degree1_cofactor_map':{'unknowns_per_target':unknowns,'target_normal_coordinates':42,'all_constraints':allconstraints,'actual_nonzero_constraint_rows':len(active),'actual_active_constraints_by_total_degree':active_by_degree,'exact_initial_nonzeros':nnz,
 'sparse_triplet_uint32_row_col_fmpq16_byte_base_bytes':nnz*24,'dense_all_fmpq16_byte_base_bytes':dense_raw*16,'dense_active_fmpq16_byte_base_bytes':dense_active*16,'dense_all_uint32_bytes':dense_raw*4,'dense_active_uint64_nmod_bytes':dense_active*8,'sparse_elimination_fill_and_QQ_height_not_bounded_by_initial_nnz':True},
 'constant_block_normalization':{'C':'A_selected(0), selected=0..41','inherited_determinant_mod32003':1680,'C_inverse_beta_only_recipe':True,'matrix_recipe':'first42 rows C^-1*A_selected; other87 rows A_i-A_i(0)*C^-1*A_selected','normalized_form':'[I42+L(q);N(q)], L,N homogeneous linear','new_normalized_coefficient_expansion_not_computed':True,
 'degree1_unknowns_after_eliminating_constant_and_linear_constraints':87*46,'remaining_degree2_constraints_upper_bound':42*math.comb(46,2),'degreeD_remaining_unknown_formula':'87*binom(45+D,D)','degree2_remaining_unknowns':87*math.comb(47,2),'finite_polynomial_inverse_not_implied_by_formal_inverse':True},
 'prior_no_repeat_gate':{'assessment_report':str(priorp.relative_to(ROOT)),'D0_generic_obstruction_already_known':True,'D1_same_mod32003_center_obstruction_already_known':True,'QQ_or_generic_beta_D1_obstruction_not_claimed':True,'no_D1_solver_or_new_negative_theorem':True},
 'generic_beta_lifting':{'center_positive_alone_is_not_open':True,'constant_rank_or_free_cokernel_alone_not_sufficient':True,'counterexample':'C(beta)=(1,beta)^T, b=(1,0)^T; solvable at beta=0 and not at generic beta','sufficient_gate':'full row rank of the actual universal coefficient map, or an exact beta-rational solution with every coefficient identity checked; alternatively prove universal compatibility with a smaller surjective equation space','one_nonzero_column_minor_does_not_prove_augmented_generic_membership':True},
 'no_solver_GB_new_field_point_source_regeneration_or_inverse_computation':True}
 (out/'statistics.json').write_text(json.dumps(result,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),src,prp,certp,knownp,priorp,ROOT/'research/notes/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63.md',ROOT/'research/notes/S61_AFFINE_NORMAL_FINITE5_20261103_v63.md']
 report={'id':ID,'all_read_only_source_sparsity_and_normalization_recipe_checks_passed':True,'actually_new_data_only':'source sparsity, active coefficient constraints and constant42 normalization recipe','source_terms':counts['terms'],'degree1_map_shape':[len(active),unknowns],'degree1_map_initial_nonzeros':nnz,'normalized_degree1_unknowns':4002,'no_coefficient_solver_launched':True,'no_new_D1_negative_or_generic_exhaustion_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-started,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);args=p.parse_args();main(args.output.resolve())
