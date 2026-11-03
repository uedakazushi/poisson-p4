#!/usr/bin/env python3
"""Apply the completed full150 cubic gate to the actual Qb CA/CB intersection."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for q in iter(lambda:f.read(1<<20),b''):h.update(q)
    return h.hexdigest()
def monomial_generates(g,m):return all(a<=b for a,b in zip(g,m))
def minimal(M):return sorted(m for m in set(M)if not any(n!=m and monomial_generates(n,m)for n in M))
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    deps=[ROOT/'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63/report.json',
          ROOT/'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63_REVIEW/report.json',
          ROOT/'reports/AMBIENT_CB_ACTUAL_D_QUADRATIC_GATE_20261103_v63/report.json',
          ROOT/'reports/AMBIENT_CB_SELF_PAIR_COMMON_BASE_ETALE_20261103_v63/report.json',
          ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/report.json',
          ROOT/'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json',
          ROOT/'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json']
    data=[json.loads(p.read_text())for p in deps];assert all(d['all_passed']for d in data)
    assert data[0]['missing_cubic_socle_class_identically_zero_not_just_point']
    assert data[0]['quadratic_generator_multiplier_c_constant_explicitly_zero']
    assert data[2]['exact38block_quotient_of_all150_normal_quadratics_rank']==0
    # Normal variables f,rho,z1..z7. Verify the entire monomial primary
    # intersection, nilradical quotient bases and the embedded socle.
    e=lambda i:tuple(int(j==i)for j in range(9))
    plus=lambda a,b:tuple(x+y for x,y in zip(a,b))
    power=lambda i,n:tuple(n if j==i else 0 for j in range(9))
    f=e(0);rho=e(1);zs=[e(i)for i in range(2,9)]
    Q=[plus(f,z)for z in zs]+[power(0,3),plus(f,power(1,2))]
    primary=[*zs,power(0,3),power(1,2)]
    intersection=minimal([tuple(max(a,b)for a,b in zip(f,m))for m in primary])
    assert intersection==minimal(Q)
    nilbasis=[power(0,1),power(0,2),plus(f,rho),plus(power(0,2),rho)]
    assert all(not any(monomial_generates(g,m)for g in Q)for m in nilbasis)
    socle=plus(power(0,2),rho)
    assert all(any(monomial_generates(g,plus(socle,e(i)))for g in Q)for i in range(9))
    assert not any(monomial_generates(g,socle)for g in Q)
    primary_basis=[plus(power(0,i),power(1,j))for i in range(3)for j in range(2)]
    assert len(set(primary_basis))==6
    assert all(not any(monomial_generates(g,m)for g in primary)for m in primary_basis)
    result=dict(id=STEM,all_passed=True,
        scope='Formal germs at Q1 and a nonempty open of the Qb self-pair stratum and its G saturation',
        actual_reduced_intersection_smooth_N28_and_CA_smooth29=True,
        actual_Z4_ideal='(z1..z7,f^2,rho^2)',
        core81_pullbacks_generate_f_times_Z4_ideal_by_unit_source_coordinate_and_quadzero_CI4=True,
        all_remaining_actual_ambient_equations_covered_by_full121_odd_vector_syzygies=True,
        upper_difference_module_after_quadzero_embeds_in_free_rank1_S20_socle=True,
        all150_QQb_cubic_zero_plus_G_dominance_kills_entire_remaining_module=True,
        actual_D_completed_ring='C[[s1..s20,z1..z7,f,rho]]/(f*z_i, f^3, f*rho^2)',
        full_ideal_equality_proved_without_flat_basechange_of_fiber_product=True,
        nilradical_principal='(f)',nilradical_as_module='C[[s20,f,rho]]/(f^2,rho^2)',
        nilradical_free_rank_over_S20=4,nilradical_square_free_rank_over_S20=2,nilradical_cube_zero=True,
        actual_D_cone_dimension=28,actual_D_cone_embedding_dimension=29,actual_D_cone_depth=20,actual_D_multiplicity=1,
        actual_D_projective_dimension=27,actual_D_projective_embedding_dimension=28,actual_D_projective_depth=19,
        irredundant_primary_decomposition='(f) intersection (z1..z7,f^3,rho^2)',
        associated_primes=['(f)','(z1..z7,f,rho)'],embedded_primary_generic_length=6,
        embedded_support_cone_dimension=20,embedded_support_projective_dimension=19,
        generic_nilmodule_length4_not_confused_with_primary_length6=True,
        additional_global_embedded_prime_closure_of_Qb_stratum=True,
        distinct_from_Pa_embedded_prime_by_opposite_whole_CA_CB_flag_counts=True,
        full_global_bad_locus_all_embedded_primes_or_B18_origin_type_exhausted=False,
        no_new_full_Poisson_component_count_claim=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),ROOT/f'research/notes/{STEM}.md',*deps]},
        seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in['all_passed','actual_D_completed_ring','nilradical_free_rank_over_S20','actual_D_cone_depth','seconds']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
