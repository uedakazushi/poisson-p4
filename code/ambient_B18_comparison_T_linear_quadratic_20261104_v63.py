#!/usr/bin/env python3
"""Small exact T/m² lower bound and sufficient m²T=0 certificate.

Only frozen X0/X1 and literal Koszul d3 first coefficients are used. There
is no m1, X2, higher source solve,348 matrix, or unknown quadratic relation.
"""
import argparse,hashlib,itertools,json,resource,signal,sys,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import sparse_rank,read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63'
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
      'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
      'native128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
      'muT':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
      'depth':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json'}
    D={k:read(ROOT/p) for k,p in paths.items()};X0,A,N,M,R,G=[D[k] for k in paths]
    peers={
      'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW/report.json':'baaec4550a00c25ee24b2ab5d03e9eebaf67943caa216a01abe624bf19f3dd17',
      'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63_REVIEW/report.json':'0000a1aa2790891a329f7ff708ef6dcfcaea8b04571de07b79ee5b503caabcb1',
      'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63_REVIEW/report.json':'6b240478705fee2088f5eab02db02573c099ff240b91c82431b2150ce125f587'}
    for p,h in peers.items():assert sha(ROOT/p)==h and read(ROOT/p)['all_passed']
    assert M['mu_Rbar_T_exact']==5 and M['rank_w0']==84 and X0['rank_X0']==27
    assert G['rank_D_K']==4 and G['dim_k_C_ell_torsion']==4 and G['length_k_T_lower_bound']==8
    for p,h in G['input_sha256'].items():assert sha(ROOT/p)==h
    inputs=A['source_columns640_rows56_sparse_dual9']+N['source_columns128_rows48_sparse_dual9']
    X=[{i:rec(p) for i,p in col} for col in inputs]
    constants=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations'] if r.get('zero')]
    assert len(X)==768 and len(constants)==741
    U=X0['free_kernel_rows32'];pivots=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
    pivot_basis={r['pivot']:rec(r['normalized']) for r in pivots};assert len(pivot_basis)==27
    tail=[i for i in range(32) if i not in pivot_basis];assert len(tail)==5
    source_nz=sum(len(p) for col in inputs for i,p in col)
    preflight={'actual_source_QQ_nonzeros':source_nz,'linear_target_dimension':40,'linear_relation_columns':657,
      'normalized_Q1_slots_upper':40*84,'linear_matrix_slots_upper':40*657,
      'quadratic_target_dimension':5*36,'quadratic_columns_if_linear_rank37':37*8+168,
      'quadratic_matrix_slots_upper':180*(37*8+168),'coefficient_cap':150000,
      'no_fresh_scalar_coefficient_query':True}
    assert source_nz+40*84+40*657+180*(37*8+168)<150000
    (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    (out/'phase.json').write_text(json.dumps({'phase':'all_frozen_source_arrays_and_size_guard_bound'})+'\n')
    # Constant row transform sends the27 origin columns to(I27,0); retain its
    # five quotient rows. For Xprime0=0, no derivative of this row unit occurs.
    Q1=[]
    for r in constants:
        combo=rec(r['constant_original_source_combination']);y={}
        for j,c in combo.items():
            for i,row in enumerate(U):
                for k,a in X[j].get(row,{}).items():
                    if k:y[8*i+k-1]=y.get(8*i+k-1,F(0))+c*a
        y={i:a for i,a in y.items() if a}
        for pivot in sorted(pivot_basis):
            p=pivot_basis[pivot]
            for k in range(8):
                coefficient=y.get(8*pivot+k,F(0))
                if coefficient:
                    for row,a in p.items():
                        q=8*row+k;v=y.get(q,F(0))-coefficient*a
                        if v:y[q]=v
                        else:y.pop(q,None)
        assert all(not y.get(8*i+k,F(0)) for i in pivot_basis for k in range(8))
        q={8*j+k:y[8*i+k] for j,i in enumerate(tail) for k in range(8) if y.get(8*i+k,F(0))}
        Q1.append(q)
    # Apply the already frozen W0 elimination to Q1, no new84 solve/search.
    wcert=M['complete_rank_certificate'];qnorm=[];linear=[];linearlabels=[]
    for j,(q,op) in enumerate(zip(Q1,wcert['fixed_elimination_receipts'])):
        v=dict(q)
        for p,c in op['subtract']:addscaled(v,qnorm[p],-F(c))
        if 'pivot_row' in op:qnorm.append({i:a/F(op['divide']) for i,a in v.items()})
        else:linear.append(v);linearlabels.append(constants[j]['source_column'])
    assert len(qnorm)==84 and len(linear)==657
    qsave={'five_remaining_free_source_rows':[U[i] for i in tail],
      'constant_projection_pivot_rows27':sorted(pivot_basis),'constant_projection_tail5':tail,
      'all741_Q1_columns_sparse':list(map(encode,Q1)),
      'normalized84_unitcycle_Q1_columns_sparse':list(map(encode,qnorm)),
      'actual657_minimal_relation_source_labels':linearlabels,
      'actual40_by657_linear_relation_columns':list(map(encode,linear))}
    (out/'phase1-actual-five-target-linear-relations.json').write_text(json.dumps(qsave,separators=(',',':'))+'\n')
    lcert=sparse_rank(linear,40);lr=lcert['rank'];layer=40-lr;lower=5+layer
    (out/'phase.json').write_text(json.dumps({'phase':'complete_five_target_linear_rank','rank':lr,'length_T_mod_m2':lower})+'\n')
    quadcert=None;quadratic=[];quadraticlabels=[];beta_columns=[];full=False
    if lr==37:
        # Products of a complete37-element linear span with all8 variables.
        mons=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(mons)}
        for j in lcert['selected_original_columns']:
            for k in range(8):
                v={}
                for row,c in linear[j].items():
                    target,variable=divmod(row,8);i=36*target+mi[tuple(sorted((variable,k)))];v[i]=v.get(i,F(0))+c
                quadratic.append(v);quadraticlabels.append({'kind':'variable_times_actual_linear_relation','source_relation':linearlabels[j],'variable':k})
        # Decompose d3_linear in the fixed84 normalized wedge unit basis.
        wbasis=list(map(rec,wcert['basis_normalized_columns']));leading={min(v):j for j,v in enumerate(wbasis)};assert len(leading)==84
        for d3j,col in enumerate(R['differentials8_column_sparse'][2]):
            byvar=[{} for _ in range(8)]
            for row,p in col:
                for m,c in p:
                    if sum(m)==1:byvar[m.index(1)][row]=F(c)
            beta={}
            for k,v0 in enumerate(byvar):
                v=dict(v0)
                while v:
                    i=min(v);p=leading[i];c=v[i];beta[8*p+k]=beta.get(8*p+k,F(0))+c;addscaled(v,wbasis[p],-c)
            beta_columns.append(encode(beta));v={}
            for index,c in beta.items():
                p,k=divmod(index,8)
                for row,a in qnorm[p].items():
                    target,var=divmod(row,8);i=36*target+mi[tuple(sorted((var,k)))];v[i]=v.get(i,F(0))-c*a
            v={i:c for i,c in v.items() if c};quadratic.append(v);quadraticlabels.append({'kind':'certified_actual_d3_quadratic','literal_d3_column':d3j})
        assert len(quadratic)==464 and len(beta_columns)==168
        support=source_nz+sum(map(len,Q1))+sum(map(len,qnorm))+sum(map(len,linear))+sum(map(len,quadratic))+sum(map(len,beta_columns))
        assert support<150000
        (out/'phase2-certified-quadratic-input.json').write_text(json.dumps({'quadratic_monomials36':mons,'source_columns464':quadraticlabels,
          'all464_certified_quadratic_columns_sparse':list(map(encode,quadratic)),
          'all168_d3_linear_coordinates_in_unit_wedge_basis':beta_columns,
          'only_first_two_source_levels_used':True,'unknown_cycle_quadratic_coefficients_not_used':True},separators=(',',':'))+'\n')
        quadcert=sparse_rank(quadratic,180);full=quadcert['rank']==180
        assert not full or lower==8
    exact_length=8 if full else None;depth=20 if full else 19 if lower>8 else None
    ev={'actual_mu_Rbar_T':5,'actual_linear_relation_rank':lr,'actual_dim_mT_mod_m2T':layer,
      'actual_length_T_mod_m2T':lower,'complete_linear_rank_certificate':lcert,
      'known_d3_quadratic_formula':'-Qstar_1 * Wstar_0^-1 * d3_linear',
      'first27_unit_elimination_has_no_unknown_linear_correction':'B_p0=0 and Xprime0=0',
      'full84_wedge_elimination_has_no_unknown_linear_correction':'Qstar0=0',
      'd3_has_no_constant_or_linear_minimal_relation_part':True,
      'quadratic_certificate':quadcert,'certified_quadratic_full180':full,
      'm_squared_T_zero_proved':full,'actual_total_length_T_proved':exact_length,
      'actual_origin_cone_depth_proved':depth,
      'proof_if_full':'certified degree2 relations span all5*Sym2(k8), so m²T=m³T; Nakayama gives m²T=0. Linear layer then gives total length5+3=8.',
      'no_lower_length_claim_from_partial_quadratic_span':True,
      'no_X2_m1_new_source_jet_or_348_matrix':True}
    ep=out/'actual-linear-quadratic-T-certificate.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    (out/'phase.json').write_text(json.dumps({'phase':'FINAL','linear_rank':lr,'quadratic_rank':None if quadcert is None else quadcert['rank'],'total_length_T':exact_length,'origin_cone_depth':depth})+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_mu_Rbar_T':5,'actual_linear_relation_rank':lr,'actual_length_T_mod_m2T':lower,
      'certified_quadratic_rank':None if quadcert is None else quadcert['rank'],'m_squared_T_zero_proved':full,
      'actual_total_length_T':exact_length,'actual_origin_cone_depth':depth,
      'no_X2_m1_new_source_solve_jet_point_GB_348_matrix':True,
      'input_sha256':{p:sha(ROOT/p) for p in paths.values()}|peers,
      'code_sha256':sha(Path(__file__)),'rank_helper_sha256':sha(ROOT/'code/ambient_B18_comparison_T_constant_rank_20261104_v63.py'),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','actual_linear_relation_rank','actual_length_T_mod_m2T','certified_quadratic_rank','actual_total_length_T','actual_origin_cone_depth','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
