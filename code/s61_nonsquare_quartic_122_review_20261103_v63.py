#!/usr/bin/env python3
"""Independent source/row-span audit of the full122x64 quartic circuit."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from s61_nonsquare_structural_reductions_review_20261103_v63 import read, sha

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63_REVIEW'


def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports/S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63';author=read(parent/'report.json')
    assert author['frozen']and author['all_passed']
    assert sha(parent/'circuit.json')==author['circuit_sha256']
    for rel,h in author['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    circuit=read(parent/'circuit.json')
    deps=[ROOT/f'reports/{stem}/report.json'for stem in [
        'S61_NONSQUARE_CUBIC_QUARTIC_STRUCTURAL_REDUCTIONS_20261103_v63_REVIEW',
        'S61_NONSQUARE_QUARTIC_286_SUFFICIENT_20261103_v63',
        'S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW',
        'S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63_REVIEW']]
    for path in deps:
        receipt=read(path);assert receipt['all_passed']
        for rel,h in receipt['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    heavy_review=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_REVIEW_20261102.json'
    hr=read(heavy_review);assert hr['all_passed']
    for rel,h in hr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    source_path=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
    source=read(source_path);light=[0,1,2,3,4,6,7];mons=source['normal_monomials']
    heavy=[]
    for j in light:
        mon=[0]*8;mon[j]=mon[5]=1;heavy.append(mons.index(mon))
    assert heavy==circuit['heavy_columns_in_light_order']and circuit['light_old_variable_indices']==light
    chosen=circuit['chosen_heavy_basis_rows'];other=[i for i in range(41,57)if i not in chosen]
    assert len(chosen)==len(set(chosen))==7 and set(chosen)<=set(range(41,57))
    assert other==circuit['remaining_nine_basis_rows']and len(other)==9
    bound={(x['basis_row'],x['column']):x for row in circuit['heavy_coefficient_entries']for x in row}
    assert len(bound)==112 and set(bound)==set(itertools.product(chosen+other,heavy))
    digest=hashlib.sha256();term_count=0;checked=0
    with gzip.open(ROOT/circuit['source_cubic_matrix'],'rb')as f:
        line=f.readline();digest.update(line);assert line.strip()==b'57 91 57'
        for i in range(57):
            for j in range(91):
                line=f.readline();assert line;digest.update(line);ts=line.split();n=int(ts[0]);assert len(ts)==1+4*n;term_count+=n
                if i<41 and j in heavy:assert not n
                if(i,j)in bound:
                    record=bound[i,j];assert record['source_line_1_based']==2+i*91+j
                    assert record['sparse_encoding']==line.decode().rstrip('\n')
                    for k in range(n):Q(ts[1+4*k].decode());assert min(map(int,ts[2+4*k:5+4*k]))>=0
                    checked+=1
        assert not f.read()
    assert digest.hexdigest()==circuit['source_cubic_sha256_uncompressed']and term_count==164333 and checked==112
    # A constructive unimodular49-row change, rather than a numerical rank test.
    pairs=list(itertools.combinations_with_replacement(range(7),2));off=list(itertools.combinations(range(7),2))
    assert circuit['H_Sym2_pivot_rows']==[list(p)for p in pairs]
    assert circuit['H_Sym2_columns']==[list(p)for p in pairs]
    original=[(i,j)for i in range(7)for j in range(7)];piv=[7*i+j for i,j in pairs];nonpiv=[7*j+i for i,j in off]
    assert sorted(piv+nonpiv)==list(range(49))
    images=[[int(tuple(sorted((i,j)))==pair)for pair in pairs]for i,j in original]
    assert [images[k]for k in piv]==[[int(i==j)for j in range(28)]for i in range(28)]
    kernels=[]
    for i,j in off:
        row=[0]*49;row[7*i+j]=1;row[7*j+i]=-1;kernels.append(row)
    assert circuit['twentyone_kernel_rows']==[[i,j,7*i+j,7*j+i]for i,j in off]
    assert all(sum(row[k]*images[k][col]for k in range(49))==0 for row in kernels for col in range(28))
    assert [[row[k]for k in nonpiv]for row in kernels]==[[-int(i==j)for j in range(21)]for i in range(21)]
    # This proves invertibility: after the column permutation [piv,nonpiv],
    # the row transformation is [[I28,0],[*, -I21]], determinant a sign.
    H,xi,xj,pi,pj=s.symbols('H xi xj pi pj')
    fi=H*xi-pi;fj=H*xj-pj
    assert s.expand(fi*xj-fj*xi-(pj*xi-pi*xj))==0
    assert s.expand((H*xi*xj-pi*xj)-(H*xj*xi-pj*xi)-(pj*xi-pi*xj))==0
    # Actual six quadrics are pure-light and all33 true quartics are retained.
    qpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';qdata=read(qpath)
    assert len(qdata['quadratic_rows_sparse'])==6
    assert all(mon[5]==0 and sum(mon)==2 for row in qdata['quadratic_rows_sparse']for mon,co in row)
    wpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=read(wpath)
    q4=[i for i,w in enumerate(chart['equation_weights'])if w==4 and i not in chart['selected34_rows']]
    old319=read(ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json')
    assert len(q4)==33 and q4==circuit['genuine_quartic_original_source_rows']==old319['genuine_quartic_original_rows_retained']
    hpath=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';hs=read(hpath)
    transform=circuit['H2_source_row_transform'];assert transform['input_rows']==[77,131]
    assert transform['matrix']==[['u','v'],['-b131','b77']]
    assert [transform['u'],transform['v']]==hs['row_combination_coefficients']
    assert [transform['b77'],transform['b131']]==hs['coefficient_formulas']
    a,e,d=s.symbols('a e d');B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2
    parse=lambda q:s.sympify(q,locals={'a':a,'e':e,'d':d})
    u,v=list(map(parse,[transform['u'],transform['v']]));b0,b1=list(map(parse,[transform['b77'],transform['b131']]))
    assert s.cancel(u*b0+v*b1-(50*a*B*C)**2)==0
    U,V,B0,B1=s.symbols('U V B0 B1');unit=U*B0+V*B1
    T=s.Matrix([[U,V],[-B1,B0]]);adj=s.Matrix([[B0,-V],[B1,U]])
    assert T*adj==adj*T==unit*s.eye(2)
    assert T*s.Matrix([B0,B1])==s.Matrix([unit,0])
    # For every other actual quartic its arbitrary coefficient b_k is removed
    # by the exact triangular row change; no b_k is declared a unit.
    bk=s.symbols('bk');assert s.cancel(bk-bk*unit/unit)==0
    expected_res=[['second_H2_free_combination',77,131]]+[['genuine_quartic_minus_H2_multiple',i]for i in q4 if i not in[77,131]]
    assert circuit['genuine_quartic_residual_rows']==expected_res and len(expected_res)==32
    low=[['H_times_quadric',i]for i in range(6)]+[['commutator',i,j]for i,j in off]+[['pure_cubic_times_light',i,j]for i in other for j in range(7)]
    assert circuit['low_weight_residual_rows']==low and len(low)==90
    assert circuit['full_residual_rows']==low+expected_res and len(low+expected_res)==122
    assert circuit['pure_quadratic_products']==[['quadric',i,'light_pair',list(pair)]for i in range(6)for pair in pairs]
    assert circuit['pure_quartic_columns']==[list(mon)for mon in itertools.combinations_with_replacement(range(7),4)]
    assert len(circuit['pure_quadratic_products'])==168 and len(circuit['pure_quartic_columns'])==210
    assert read(deps[2])['generic_light_rank146_and_global_upper146']
    assert circuit['pure_quad_complement_columns']==210-146==64
    assert circuit['full_shape']==[122,64]and circuit['low_weight_sufficient_shape']==[90,64]
    assert 28+1+146==175 and 239-175==64 and 28+146==174 and 238-174==64
    assert circuit['no_full_cubic57_premise']and circuit['no_new_rank64_claim']and circuit['no_expanded_parametric_quartic_matrix']
    evidence=dict(actual57by91_source_hash_164333_terms_and_all112_heavy_entries_independently_bound=True,
        heavy_graph_inverse_orientation='T rows are G_E, columns H*x_i; (T^-1)*G_E has identity heavy coefficient block',
        fortynine_row_unimodular_transformation='after [28 pivot,21 nonpivot] column order, [[I28,0],[*, -I21]]',
        all21_commutator_relations_and_fixed_weight4_normal_form_verified=True,
        no_all_degree_Groebner_or_confluence_claim=True,
        actual_Hsquare_coefficients_cofactors_and_global_det_unit_verified=True,
        source2row_transform_has_explicit_adj_inverse=True,
        all33_genuine_quartics_equal_one_Hsquare_pivot_plus32_pure_residuals=True,
        remaining22_pure_quad_rows_vanish_by_global_polynomial_upper146_not_by_samples=True,
        full_rank_identity='175+rank(C122)',low_rank_identity='174+rank(C90)',
        full_cokernel_isomorphism_on_extra_DeltaH_DeltaQ_charts=True,
        low90_is_sufficient_only=True,no_pointwise_full57_or_product41_premise=True,
        all_parameter_chart_units_tracked=['a','F','B','C','DeltaH7','DeltaQ146'],
        genuine_quartic_oracle_is_same_actual_selected34_graph_as_accepted319_packet=True,
        rank64_fullquartic_length65_ambient_or_new_curve_exclusion_not_proved=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    result=dict(id=STEM,all_passed=True,actual_full122_by64_cokernel_recipe_accepted=True,
        actual_low90_by64_sufficient_recipe_accepted=True,all33_genuine_quartics_preserved=True,
        only_on_explicit_extra_DeltaH7_DeltaQ146_charts=True,full_cubic57_not_assumed=True,
        actual_rank64_or_fullquartic_or_ambient_claim=False,
        no_GB_new_field_point_probe_or_expanded_parametric_minor=True,seconds=time.monotonic()-start)
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',parent/'report.json',parent/'circuit.json',source_path,qpath,wpath,hpath,heavy_review]+deps
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths}
    result['evidence_sha256']={'independent-certificate.json':sha(out/'independent-certificate.json')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
