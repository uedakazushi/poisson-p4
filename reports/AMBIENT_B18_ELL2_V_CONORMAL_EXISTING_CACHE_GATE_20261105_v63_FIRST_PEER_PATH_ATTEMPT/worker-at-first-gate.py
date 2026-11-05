#!/usr/bin/env python3
"""Final source/cache/size gate only; no V coefficient evaluation."""
import argparse,hashlib,itertools,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_GATE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hs={}
    def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    profile=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/report.json');assert profile['all_passed']
    adj=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json')
    full=read('reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    joint=read('reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    oldP=read('reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json')
    Y=read('reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json')
    source=read('reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    chain=read('reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json')
    for p in ['reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json']:
        assert read(p)['all_passed']
    assert joint['all105_equations_zero_through_degree3']
    assert all(full['both57_source2'][f]['all57_equations_all55_coefficients_verified'] for f in ('CA','CB'))
    assert set(joint['all105_unknown_original_node_order'])<=set(full['both57_source2']['CA']['unknown57'])|set(full['both57_source2']['CB']['unknown57'])
    assert full['both57_source2']['CA']['F9_nodes'][:7]==[7148,7199,7454,7505,7688,7943,8045]
    assert chain['all48_existing_Q_frame_columns_verified'] and chain['exact_existing_U9_C4'][0]==[1,0,0,0]
    assert chain['physical_target_rows12']==[q['physical_row'] for q in oldP['physical12_coefficient_roots']]
    for q in oldP['physical12_coefficient_roots']:
        r=q['physical_row']
        for f,an,bn in [('CA','aA_node','bA_node'),('CB','aB_node','bB_node')]:
            assert [q[an],q[bn]]==source['families'][f]['complete_normalized_physical_w_coefficients'][0][r]
    basis=[]
    for degree in range(3):
        for inds in itertools.combinations_with_replacement(range(9),degree):
            m=[0]*9
            for i in inds:m[i]+=1
            if m[8]<2:basis.append(m)
    basis+=[[0]*7+[3,0],[0]*7+[2,1]];assert len(basis)==56
    bs={tuple(m) for m in basis}
    for m in basis:
        for i,a in enumerate(m):
            if a:n=list(m);n[i]-=1;assert tuple(n) in bs
    assert all(tuple(v['monomial9']) in bs for v in adj['scalar_monomial_adjoint_library'])
    assert not any(sum(v['monomial9'])==4 for v in adj['scalar_monomial_adjoint_library'])
    cubic=[v for v in adj['scalar_monomial_adjoint_library'] if sum(v['monomial9'])==3]
    assert len(cubic)==4 and all(v['physical_target'] in [0,1] and v['monomial9'] in basis[-2:] for v in cubic)
    monsource=full['source2_monomials'];assert len(monsource)==55
    assert len(joint['monomials'])==10 and len(joint['all105_unknown_jets'])==105
    budgets={'source_closure_metadata':18584,'source_live414_dense56':414*56,
      'cached105_values':105*56,'all105_equation_check':105*56,'all164_roots':164*56,
      'full9_source2_and_bivariate3_cache_inputs':2*57*55+105*10+9*55,
      'nonlinear_YB_copy':9*56,'monomial_multiplication_lookup':56*56,
      'all56_primitiveV_columns_dense':56*3*56,'actual26_by56_matrix':26*56,
      'callback_metadata':1224,'saved105_cache_values':105*56,'source_inverse_checks':6*2*56}
    total=sum(budgets.values());assert total<150000
    evidence={'coefficient_ring_basis56':basis,'all_omitted_monomials_form_an_ideal':True,
      'Pi26_factors_through_B56':True,'all_order4_scalar_callbacks_annihilated':True,
      'new_implicit_source_jet_gauge_or_function_value_to_solve':False,
      'cache_recipe_full9_second':'CA own57; CB own57 composed with FULL9 YB1+YB2 including quadratic correction, then sA=r+ell/5',
      'cache_recipe_plane_third':'COMMON105 bivariateCA(rA,sA) jets, substitute sA=r+ell/5; retainr³,ellr²',
      'required_checks':['all105 original equations in B56','both57 cached equations/full9 nonlinear physical inverse through2',
        'plane coeff0..2 agree betweenbothcaches','all6unitinverseproducts','alloriginalprimitivephysical12 root bindings',
        'exactsame48Qframe source identity','gA505938 prefix0f63=7790 literalbinding','all56conormalcolumnfactorizations/nonzero26minor'],
      'primitive_V_per_physical_row':['(bA,-bB,0)','(delta,0,-bB)','(0,delta,bA)','(-nuB*bB,nuA*bA,delta)'],
      'primitive_U8':['(-zA0,0,0)', '...', '(-zA6,0,0)', '(-gA_plane,0,0)'],
      'no_differentiation_of_expanded_membership_cofactors':'accepted allorder Q48frame equality allows direct original primitive physical scalar roots',
      'gA_actual_source_node':505938,'gA_original_CA_f63':7790,
      'gA_plane_definition':'prefix0[7790], seven CA F9 target coordinates zero; same retainedr/s unit57equations',
      'gA_through3_cache':joint['g_CA_all10_coefficients'],
      'source_scalar_roots_catalog':oldP['physical12_coefficient_roots'],'nu_roots':oldP['nu_roots_CA_CB'],
      'budget':budgets,'total_slots_bound':total,
      'limits':{'seconds':60,'bytes':2_000_000_000,'artifact_bytes':32*1024*1024,'QQ_numerator_or_denominator_bits':384,'total_slots':150000},
      'source_coefficients_or_conormal_values_evaluated':False,
      'still_no_actual_T2_or_H0_length':True}
    ep=out/'whole-existing-cache-V-conormal-source-gate.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    report={'id':STEM,'all_passed':True,'coefficient_ring_dimension':56,'whole_slot_bound':total,
      'new_source_coefficients_to_solve':0,'V_conormal_coefficients_evaluated':False,
      'input_sha256':hs,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)}}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({'all_passed':True,'B56':56,'whole_slot_bound':total,'source_values_evaluated':False}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
