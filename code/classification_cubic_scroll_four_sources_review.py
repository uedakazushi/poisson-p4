#!/usr/bin/env python3
"""Independent Fraction replay of the four whole smooth-scroll zero sources.

No producer, SymPy, CAS, point sample, rank or minor computation.  Only the
previous independent review's elementary sparse polynomial helpers are reused.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

from classification_cubic_scroll_abelian_review import add, scale, mul, diff, parse, monomials

ROOT = Path(__file__).resolve().parents[1]
PAIRS = list(it.combinations(range(5), 2))
CHECKS = 0


def same(a, b, label):
    global CHECKS
    assert a == b, label
    CHECKS += 1


def product(items, n):
    out = {(0,) * n: Q(1)}
    for v in items:
        out = mul(out, v)
    return out


def substitute(poly, images, n):
    out = {}
    for ex, co in poly.items():
        out = add(out, scale(product([images[i] for i, power in enumerate(ex)
                                     for _ in range(power)], n), co))
    return out


def whole_split_identities(gens):
    names = [f'r{i}' for i in range(12)]
    R = lambda text: parse(text, names)
    # Independent displayed primitive identities supporting every unit split.
    first = [(19, '(r11-r2)*(r11-2*r2)'),
             (1, '(r11-2*r2)*r5'), (9, '(r11-2*r2)*r8'),
             (20, '(r11-2*r2)*(r10-r1)'),
             (4, 'r4*(-3*(r10-r1)+4*r5)'),
             (8, 'r4*((r10-r1)-r5)'),
             (5, 'r4*(r4+r9)'), (13, 'r4*r7')]
    for i, text in first:
        same(gens[i-1], R(text), f'primitive forcing g{i}')
    ident = [R(f'r{i}') for i in range(12)]
    ident[11], ident[4] = R('2*r2'), {}
    stage = [substitute(g, ident, 12) for g in gens]
    same(add(stage[1], scale(stage[2], -1)),
         R('(r10-r1+2*r5)*(r1-2*r5)'), 'stage g2-g3')
    same(add(stage[5], scale(stage[2], 6)),
         R('-2*(r1-2*r5)*(r10-r1+r5)'), 'stage g6+6g3')
    same(stage[21], R('2*r2*(r1-2*r5)'), 'stage g22')
    same(stage[6], R('-r9*(r1-2*r5)'), 'stage g7')
    # a=2f,b=2f+delta; use r10 as the independent delta in this substitution.
    ident[1], ident[10] = R('2*r5'), R('2*r5+r10')
    reduced = [substitute(g, ident, 12) for g in gens]
    nine = [R(z) for z in [
        'r2*r9-r5*r10', 'r5*r8', 'r2*r7-r10*r8',
        'r0*r5+r2*(r3+r7)', 'r8*r9-r5*r7',
        'r0*r9+r10*r3+(r10-2*r5)*r7',
        'r0*r8+r2*r6-2*r8**2',
        '-r0*r7-(r10+2*r5)*r6+2*r3*r8+4*r7*r8',
        'r7*(r3+r7)-r6*r9']]
    forms = [R(z) for z in [
        '0','r0','r0','0','0','-6*r0','0','0','0',
        '2*r3-4*r1+r2','-2*r3+4*r1-r2','-2*r5+2*r4','0',
        'r5-r4','r5+r4','r6','r7','r8','0','0','6*r0','0','-4*r0','r2']]
    for i, form in enumerate(forms):
        same(reduced[i], substitute(form, nine + [{}]*3, 12), f'all24 to R9 {i+1}')
    backward = [reduced[1], reduced[23],
                scale(add(reduced[9],scale(reduced[23],-1)),Q(1,2)),
                scale(add(reduced[14], scale(reduced[13],-1)),Q(1,2)),
                scale(add(reduced[14],reduced[13]),Q(1,2)),
                reduced[15],reduced[16],reduced[17]]
    targets=[nine[0],nine[2],add(nine[3],scale(nine[1],-2)),*nine[4:]]
    for i in range(8):
        same(backward[i], targets[i], f'correct eight-system constant inverse {i+1}')
    same(add(reduced[9],reduced[10]),{},'g10+g11 identically zero')
    return dict(all_unit_split_primitive_identities=True,
                displayed_nine_system_equivalence=False,
                correct_eight_system_and_constant_inverse=True,
                exact_defect='g11=-2R4+4R2-R3, so the recovered equation is R4-2R2, not R4 and R2 separately',
                localized_parameter_denominator_used=False)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args(); begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115))
    resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    paths=[Path(__file__).relative_to(ROOT).as_posix(),
           'code/classification_cubic_scroll_abelian_review.py',
           'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json',
           'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_DIMENSION_20261107.json',
           'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_REVIEW_20261107.json',
           'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_FOUR_SOURCES_20261108.json',
           'research/classification/CUBIC_SCROLL_ZERO_FOUR_WHOLE_SOURCES.md',
           'research/classification/CUBIC_SCROLL_ZERO_ACTUAL_FACTOR.md',
           'research/classification/inputs/cubic_scroll_zero_four_sources_manifest_20261108.json',
           'research/classification/inputs/cubic_scroll_four_review_first_failure_20261108.json']
    source=json.loads((ROOT/paths[2]).read_text())['cases'][0]
    dim=json.loads((ROOT/paths[3]).read_text())['cases'][0]
    oldreview=json.loads((ROOT/paths[4]).read_text())
    producer=json.loads((ROOT/paths[5]).read_text())
    assert source['name']=='smooth_scroll_12' and oldreview['all_passed']
    assert oldreview['all350_original_and_reduced_Jacobi_reconstructed']
    assert oldreview['all24_generators_original_Jacobi_constant_multiples']
    rnames=[f'r{i}'for i in range(12)]
    gens=[parse(t,rnames)for t in dim['generators']]
    split=whole_split_identities(gens)
    maps=[('L6',[f'l{i}'for i in range(6)],
           ['l0','l1','l2','l3','l4','0','l5','0','0','-l4','l1','l2']),
          ('F_Segre4',[f't{i}'for i in range(5)],
           ['t1*t4','2*t0*t3','t0*t4','-t1*t3','0','t0*t3','0','0','0',
            't2*t3','2*t0*t3+t2*t4','2*t0*t4']),
          ('U_Segre4',[f't{i}'for i in range(5)],
           ['t1*t4+2*t0*t3','0','t0*t4','-t2*t3','0','0','-t1*t3',
            't2*t3','t0*t3','0','t2*t4','2*t0*t4']),
          ('S_Segre4',[f't{i}'for i in range(5)],
           ['t1*t4','0','0','t1*t3-t0*t2','0','0','t1*t2','t0*t2','0',
            't0*t3','-t0*t4','0']),
          ('COUNTEREXAMPLE_OUTSIDE_FOUR',['q0'],
           ['2','2','0','q0','0','1','q0','0','1','0','2','0'])]
    B=[[Q(v)for v in row]for row in source['complete_original30_by_reduced_source_matrix']]
    assert len(B)==30 and all(len(row)==12 for row in B)
    rmon=source['all_reduced_source_quadratic_monomials']
    assert rmon==[list(v)for v in it.combinations_with_replacement(range(12),2)]
    savedcases={c['source']:c for c in producer['cases']};results=[]
    for name,pars,texts in maps:
        saved=savedcases.get(name);names=pars+[f'x{i}'for i in range(5)];n=len(names);np=len(pars)
        r=[parse(t,names)for t in texts]
        if saved is not None:
            same(pars,saved['all_parameters'],name+' parameter order')
            same(r,[parse(t,names)for t in saved['full_original_r0_to_r11_map']],name+' all maps')
        for i,g in enumerate(gens):
            same(substitute(g,r,n),{},name+f' original g{i+1}')
        quads=[parse(t,names)for t in ['x0*x3-x1*x2','x0*x4-x1*x3','x2*x4-x3**2']]
        P=[[{}for _ in range(5)]for _ in range(5)]
        for pos,(i,j)in enumerate(PAIRS):
            P[i][j]=add(*(mul(quads[g],add(*(scale(r[k],B[10*g+pos][k])
                             for k in range(12))))for g in range(3)))
            P[j][i]=scale(P[i][j],-1)
            if saved is not None:
                same(P[i][j],parse(saved['complete_actual_brackets'][str((i,j))],names),name+' full150')
        deriv=[[[diff(P[i][j],np+l)for l in range(5)]for j in range(5)]for i in range(5)]
        strict={}
        for tri in it.combinations(range(5),3):
            strict[tri]=add(*(mul(P[i][l],deriv[j][k][l])
                for i,j,k in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]
                for l in range(5)))
        for quad in it.combinations(range(5),4):
            value=add(*(scale(mul(parse(f'x{quad[i]}',names),strict[quad[:i]+quad[i+1:]]),(-1)**i)
                        for i in range(4)))
            same(value,{},name+' full350 direct')
        rproducts=[mul(r[i],r[j])for i,j in rmon]
        evaluate=lambda row:add(*(scale(rproducts[int(k)],Q(co))for k,co in row.items()))
        for row in source['all350_reduced_Jacobi_rows']:
            same(evaluate(row),{},name+' all350 original saved rows')
        D={}
        for missing in range(5):
            i,j,k,l=[z for z in range(5)if z!=missing]
            pf=add(mul(P[i][j],P[k][l]),scale(mul(P[i][k],P[j][l]),-1),mul(P[i][l],P[j][k]))
            D=add(D,scale(mul(parse(f'x{missing}',names),pf),(-1)**missing))
        if saved is not None:
            same(D,parse(saved['actual_signed_Euler_Pfaffian'],names),name+' signedD')
        for ex,row in zip(monomials(5,5),source['all126_reduced_actual_D_rows']):
            # Keep parameter monomials intact while comparing each of the126 physical coefficients.
            coeff={e[:np]+(0,)*5:co for e,co in D.items()if e[np:]==ex}
            same(coeff,evaluate(row),name+' all126 savedD')
        if saved is None:
            same(r[5],parse('1',names),'counterexample f1')
            same(r[8],parse('1',names),'counterexample u1')
            results.append(dict(source=name,full_original_r_map=texts,
                all24_original_generators_zero=True,all350_direct_and_original_rows_zero=True,
                all126_original_Pfaffian_coefficients_match=True,
                outside_each_four_sources='L requires f=u=0, F requires u=0, U requires f=0, S requires f=u=0; here f=u=1',
                complete_actual_brackets_sparse={str((i,j)):{','.join(map(str,e)):str(co)for e,co in P[i][j].items()}for i,j in PAIRS},
                full_actual_D_sparse={','.join(map(str,e)):str(co)for e,co in D.items()},
                factorization_not_claimed=True))
            continue
        prod=parse(saved['factor_coefficient'],names);degrees=[];paramfactors=[]
        for item in saved['full_polynomial_factors']:
            factor=parse(item['factor'],names);assert factor
            degree={sum(e[np:])for e in factor}
            same(degree,{item['physical_coordinate_degree']},name+' homogeneous factor degree')
            for _ in range(item['multiplicity']):
                prod=mul(prod,factor)
                if item['physical_coordinate_degree']>0:degrees.append(item['physical_coordinate_degree'])
                else:paramfactors.append(item['factor'])
        same(prod,D,name+' literal full factor product')
        assert D and len(degrees)>=2 and sum(degrees)==5
        results.append(dict(source=name,all_parameters=pars,all12_original_maps=True,
            all24_original_generators_zero=True,whole30by12_tensor_reconstructed=True,
            full150_coefficients_match=True,all350_direct_and_original_rows_zero=True,
            all126_signed_Pfaffian_coefficients_match=True,full_factor_product_identity=True,
            physical_factor_degrees=degrees,parameter_only_factors=paramfactors,
            every_nonzero_specialization_reducible=True,all_zero_parameter_complements_retained=True))
    result=dict(task_id='E_reduced_nonnormal',all_passed=False,all_arithmetic_checks_passed=True,
        four_displayed_sources_positive_checks_passed=True,whole_four_source_cover_rejected=True,
        scope='Four displayed whole sources have proper factor identities, but their cover is false: exact all350 counterexample outside four. Whole smoothscroll irreduciblequintic emptiness not proved.',
        no_producer_or_CAS_import=True,no_new_rank_minor_operator_or_point=True,
        whole_reduced_source_cover_hand_review='REJECTED: incorrect fu0 inference; an additional open fu unit branch is present',
        split_identities=split,cases=results,new_component_count=0,independent_identity_count=CHECKS,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun,
        maxrss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k]for k in ['all_passed','all_arithmetic_checks_passed','whole_four_source_cover_rejected','independent_identity_count','seconds','maxrss_kb']}))


if __name__=='__main__':
    main()
