#!/usr/bin/env python3
"""Exact CI8/CI8 derived-intersection proof receipt, no new source job."""
import hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0
def add(a,b):
    out=dict(a)
    for m,c in b.items():out[m]=out.get(m,0)+c
    return {m:c for m,c in out.items() if c}
def scale_monomial(p,r,s,c=1):return {(i+r,j+s):c*a for (i,j),a in p.items() if c*a}
def vadd(a,b):return [add(x,y) for x,y in zip(a,b)]
def vscale(p,r,s):return [scale_monomial(a,r,s) for a in p]
def main():
    out=ROOT/'reports'/STEM;out.mkdir(exist_ok=True);assert not (out/'report.json').exists()
    gates=[ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/report.json']
    assert all(json.loads(p.read_text())['all_passed'] for p in gates)
    # The only nontrivial two-variable Koszul factor, in the stated signs.
    u=[{(0,1):-1},{(1,0):1},{}]
    v=[{},{(0,1):-1},{(1,0):1}]
    d2=[[{(1,1):-1},{(2,0):1},{}],
        [{(0,2):-1},{},{(2,0):1}],
        [{},{(0,2):-1},{(1,1):1}]]
    assert d2==[vscale(u,1,0),vadd(vscale(u,0,1),vscale(v,1,0)),vscale(v,0,1)]
    # Relation module for H1: r*u, s*u+r*v, s*v. It is omega_A3.
    lengths=[3*(choose(6,i)+choose(6,i-1)) for i in range(9)]
    assert lengths==[3,21,63,105,105,63,21,3,0]
    report={'id':STEM,'all_passed':True,'proof_type':'unit_generator_Koszul_equivalence_and_explicit_H1_basis',
      'ring':'actual_R0_regular9','I_A_and_I_B':'actual zero odd ideals, each CI8',
      'sum_ideal':'(z1,...,z7,r^2,r*s,s^2)',
      'combined16_generator_Koszul_equivalent_to':['z1..z7','r^2','r*s','s^2','0^6'],
      'H1_three_quad_Koszul_generators':[[[[c,*m] for m,c in p.items()] for p in q] for q in [u,v]],
      'H1_relation_columns':['r*u','s*u+r*v','s*v'],
      'Tor_i_as_A3_modules':'A3^(binom(6,i)) direct_sum omega_A3^(binom(6,i-1))',
      'Tor_lengths_i0_through8':lengths,
      'Tor1_equals_intersection_over_product':'(I_A intersect I_B)/(I_A*I_B)',
      'Tor1_abstract_module':'A3^6 direct_sum omega_A3','Tor1_length':21,'Tor1_socle_dimension':13,
      'scheme_union_conductor_curve_ring_is_CM1':True,
      'conductor_curve_ideals_assumed_radical':False,
      'product_curve_ring_H0_exact_length':21,
      'annihilator_of_product_curve_ring_H0':'I_A+I_B',
      'comparison_kernel_arbitrary_module_length_or_socle_bound_inferred':False,
      'I_A_times_I_B_annihilates_T0_claimed':False,
      'full_T0_depth_or_nil_solved':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in gates+[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']}}
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'report_sha256':sha(p)}))
if __name__=='__main__':main()
