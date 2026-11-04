#!/usr/bin/env python3
"""Read-only exact source metadata for a complete localized H12 adapter.

No H12 transformation entries are expanded or restricted to the physical
free-four algebra in this preflight.  The distinction is an explicit gate.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_PROPER_SYZYGY_ADAPTER_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    hp=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json'
    kp=ROOT/'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json'
    cp=ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/membership-and-source-cofactors.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json'
    gates=[ROOT/'reports'/f'{x}_20261104_v63_REVIEW/report.json' for x in [
      'AMBIENT_B18_PROPER_POLYNOMIAL_MODEL','AMBIENT_B18_PROPER_RATIONAL_ELIMINATION',
      'AMBIENT_B18_PROPER_CORE_ELIMINATION','AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP',
      'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL']]
    assert all(read(p)['all_passed'] for p in gates)
    H,K,C,E=map(read,[hp,kp,cp,ep])
    order=H['remaining9_variable_order'];assert order==['a','d','v2','v3','v4','v5','v6','v7','v8']
    selected=[0,4,1,2,3];pivots=['v4','v3','v5','v7','v8'];remaining=[5,6,7,8,9,10,11]
    assert K['selected5_source_H_rows']==selected and K['selected5_pivot_order']==pivots
    hs=[q['divided_primitive_polynomial'] for q in H['all12_actual_remaining_and_divided_records']]
    # Check the exact triangular structure without evaluating any root.
    triangular=[]
    units_expected=[(F(2),0,1),(F(8),1,0),(F(80),2,0),(F(16),2,0),(F(192),3,0)]
    for j,(row,pivot,unit) in enumerate(zip(selected,pivots,units_expected)):
        loc=order.index(pivot);future=[order.index(q) for q in pivots[j+1:]]
        assert max(m[loc] for m,n,d in hs[row])==1
        assert all(not any(m[k] for k in future) for m,n,d in hs[row])
        assert all(not any(m[order.index(q)] for q in pivots if q!=pivot)
                   for m,n,d in hs[row] if m[loc])
        record=K['unit_tests'][j];assert record['row']==row and record['pivot']==pivot
        c,ed,epow=unit
        assert record['D_power']==ed and record['P_power']==epow and F(record['constant_factor'])==c
        triangular.append({'source_row':row,'pivot':pivot,'coefficient_constant':str(c),
          'D_power':ed,'P_power':epow,'origin':str(c*5**epow)})
    lambdas=[];lambda_inv_origin=F(1)
    for q in K['remaining7_actual_core_polynomials']:
        c=F(q['QQ_unit_content_numerator'],q['QQ_unit_content_denominator'])
        dp=q['D_denominator_power'];pp=q['P_denominator_power']
        lambdas.append({'QQ_constant':str(c),'D_denominator_power':dp,'P_denominator_power':pp})
        lambda_inv_origin*=5**pp/c
    detA_origin=F(1)
    for q in triangular:detA_origin*=F(q['origin'])
    # Π has sign -1; M3 has determinant -2/(5*D^3) at the origin.
    determinant_origin=-F(1)/detA_origin*lambda_inv_origin*(-F(2,5))
    assert determinant_origin==F(5**10,2359296)
    ht=sum(map(len,hs));rt=sum(len(q['numerator']) for q in K['five_actual_rational_roots'].values())
    ct=sum(len(q) for row in C['cofactor3_by4'] for q in row)
    assert (ht,rt,ct)==(1188,248,4248)
    # Factored scalar-circuit ceiling, not an expanded-polynomial bound.
    # 50 per H monomial: at most 5 telescopes, each at most9 factors+sum.
    # 6 per four-core-variable root/cofactor term, plus shared-power and
    # triangular unit matrix overhead5000. Applying five12x12 factors to
    # 32syzygy vectors adds at most 5*12^2*32 nodes.
    scalar_bound=5000+50*ht+6*(rt+ct)+5*12*12*32
    assert scalar_bound==114416 and scalar_bound<150000
    interface={'actual_R0_rows9':E['actual_R0_rows9'],
      'R0_variable_nodes':E['R0_variable_nodes'],'free4_basis':E['free4_basis'],
      'nu_CA_node':E['nu_CA_node'],'nu_CB_node':E['nu_CB_node'],
      'F12_order':[q['physical_row'] for q in E['all12_zero_physical_difference_columns']],
      'free4_relation_column_rule':'column(4*i+a)=F_i*basis_a',
      'proper_S9_to_physical_C_coefficient_map_serialized':False,
      'proper_H12_to_physical_F12_invertible_row_map_serialized':False,
      'why_two_maps_required':'A source-gauge quotient isomorphism alone is not a literal lift of all coefficients and relations in the regular ambient C.',
      'permitted_finite_interface':'actual proper-source/group factorization, pointed unit Hensel coordinate adapter and unit12 target-difference row adapter',
      'not_permitted':'type-only identification, dividing a conductor/root at B18, assuming a rational global source inverse'}
    recipe={'source_H12_file':str(hp.relative_to(ROOT)),
      'regular_source_ring':'S=QQ[a,d,v2,v6,u0..u4,D^-1,P^-1] completed at origin',
      'D':'1-d^2','P':'d^2+5','permutation_selected_then_rest':selected+remaining,
      'pivots':pivots,'u_definition':'u_i=pivot_i-frozen_root_i(a,d,v2,v6)',
      'telescoping':'H_j=R_j+sum_i A_ji*u_i; A_ji=sum_m coefficient_m*base^m*prod_{l<i}(root_l+u_l)^e_l*sum_{k=0}^{e_i-1}(root_i+u_i)^(e_i-1-k)*root_i^k*prod_{l>i}root_l^e_l',
      'triangular_diagonal':triangular,'residual_Rj_equals_lambda_j_Kj':lambdas,
      'factors_in_action_order':[
        'Π: selected5 followed by remaining7',
        'E1=[[A^-1,0],[-B*A^-1,I7]]: maps ΠH to (u5,R7)',
        'E2=diag(I5,lambda_0^-1,...,lambda_6^-1): gives (u5,K7)',
        'E3=diag(I5,M3,I4), M3=[[1/D,0,0],[0,1/D,0],[-r/(PD),0,-2/(PD)]]: gives (u5,Q0,Q1,s*r²,K3..K6)',
        'E4=unit lower block subtracting Uj*C_ij times the three CI rows from each last graph row: gives (u5,Q0,Q1,s*r²,0^4)'],
      'full_transform':'T=E4*E3*E2*E1*Π','inverse':'Π^-1*E1^-1*E2^-1*E3^-1*E4^-1',
      'determinant':'P^10/(2359296*D^11)','determinant_origin':str(determinant_origin),
      'cofactors_file':str(cp.relative_to(ROOT)),
      'cofactors_key':'cofactor3_by4; Kj=Uj*sum_i C_ij*G_i',
      'whole_E4_preserved':'T*(d*H12)=d*(u5,Q0,Q1,s*r²,0^4)',
      'syzygy_orientation':'For equation columns c=T*h, syz(h^t)=T^t*syz(c^t)',
      'complete_source_syzygy_count':32,
      'complete_source_syzygy_generators':'28 Koszul on regular8 (u5,Q0,Q1,s*r²) plus4 unit coordinate vectors on zero rows, then T^t',
      'why_d_does_not_change_syzygies':'S is regular domain and d nonzero, so d*h^t*v=0 iff h^t*v=0',
      'regular8_proof_gate':str(gates[3].relative_to(ROOT)),
      'source_no_unknown_formal_inverse':'all root rational functions and all inverse diagonals use only D,P,QQ units',
      'factor_circuit_scalar_node_ceiling':scalar_bound,
      'expanded_T12_or128_R0_syzygies_not_constructed':True,
      'physical_free4_interface':interface}
    p=out/'source12-to-regular8-adapter-preflight.json';p.write_text(json.dumps(recipe,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',hp,kp,cp,ep,*gates]
    report={'id':STEM,'all_passed':True,'status':'complete_localized_source12_adapter_formula_and_explicit_physical_restriction_gate',
      'complete_source_syzygy_recipe_count':32,'regular_generators8':True,
      'all5_d_times_u_retained':True,'source_T12_determinant':recipe['determinant'],
      'source_T12_origin':str(determinant_origin),'finite_scalar_circuit_ceiling':scalar_bound,
      'literal_T12_entries_constructed':False,'physical_free4_coefficients_or128_syzygies_claimed':False,
      'D0_K_ell_torsion_or_depth_claimed':False,
      'source_counts':{'H12_terms':ht,'root_numerator_terms':rt,'cofactor_terms':ct},
      'input_sha256':{str(q.relative_to(ROOT)):sha(q) for q in inputs},
      'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','source_T12_origin','finite_scalar_circuit_ceiling','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
