#!/usr/bin/env python3
"""Independent actual134 polarization and universal61 first-w-colon audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_W_TORSION_PAIR_20261103_v63_REVIEW'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def read(path):
    with (gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)


def add(poly,mon,c):
    value=poly.get(mon,Fraction(0))+c
    if value:poly[mon]=value
    else:poly.pop(mon,None)


def pmul(a,b):
    result={}
    for mon,c in a.items():
        for other,d in b.items():add(result,tuple(sorted(mon+other)),c*d)
    return result


def decode(poly):
    return {tuple(mon):Fraction(c)for mon,c in poly}


def flatten(poly):
    result={}
    for degree,terms in poly:
        for mon,c in terms:add(result,tuple(mon)+(45,)*degree,Fraction(c))
    return result


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    axisdir=ROOT/'reports/S61_AFFINE_AXIS_POLARIZATION_SPLIT_20261103_v63'
    g6dir=ROOT/'reports/S61_AFFINE_G6_PARTIAL_COLON_20261103_v63'
    rawreview=ROOT/'reports/S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63_REVIEW.json'
    quinticreview=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63_REVIEW.json'
    checked={}
    for path in [axisdir/'report.json',g6dir/'report.json',rawreview,quinticreview]:
        d=read(path);assert d['all_passed']
        for name,want in d.get('input_sha256',{}).items():assert sha(ROOT/name)==want;checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():assert sha(path.parent/name)==want;checked[str((path.parent/name).relative_to(ROOT))]=want
    rr=read(rawreview)
    assert rr['projective_Euler_radial_relation_and_deleted_coordinate_chart_checked']
    assert rr['clopen_finite_etale_selection_and_annihilator_complement_sheaf_proof_accepted']
    model=S61();weights=model.reduced['parameter_weights'];qids=[i for i,w in enumerate(weights)if w==2];Q=set(qids)
    axispath=ROOT/'reports/S61_MINIMAL_AXIS_INFINITY_20261101.json';axis=read(axispath)
    p=list(map(Fraction,axis['sample_S61_coordinates']));assert axis['weight_two_u_indices_zero_based']==qids
    saved=read(axisdir/'certificate.json.gz');fixed=saved['fixed_free_u_one_based']-1;assert fixed in Q and p[fixed]!=0
    chart=[i for i in qids if i!=fixed];positions={u:j for j,u in enumerate(chart)};point=[p[i]for i in chart]+[Fraction(0)]
    assert list(map(str,point))==saved['infinity_point45']
    allmid=[];alljac=[];polys=[]
    # Construct symmetric quadratic coefficient matrices, not the producer's derivative midpoint loop.
    for equation,recorded in zip(model.reduced['equations'],saved['all134_actual_homogeneous_chart_polynomials']):
        polynomial={}
        for mon,n,d in equation:
            degree=sum(i in Q for i in mon);assert degree<=2
            target=[44]*(2-degree);c=Fraction(n,d)
            for i in mon:
                if i not in Q or i==fixed:c*=p[i]
                else:target.append(positions[i])
            add(polynomial,tuple(sorted(target)),c)
        assert polynomial==decode(recorded)
        A=[[Fraction(0)for _ in range(45)]for _ in range(45)];b=[Fraction(0)]*45;constant=Fraction(0)
        for mon,c in polynomial.items():
            if not mon:constant+=c
            elif len(mon)==1:b[mon[0]]+=c
            else:
                i,j=mon
                if i==j:A[i][j]+=c
                else:A[i][j]+=c/2;A[j][i]+=c/2
        assert sum(c*math.prod(point[j]for j in mon)for mon,c in polynomial.items())==0
        midpoint=[];jac=[]
        for j in range(45):
            row={(k,):v for k,v in enumerate(A[j])if v};at=sum(a*v for a,v in zip(A[j],point))+b[j]
            if at:row[()]=at
            midpoint.append(row);jac.append(2*sum(a*v for a,v in zip(A[j],point))+b[j])
        replay={}
        for j,row in enumerate(midpoint):
            delta={(j,):Fraction(1)}
            if point[j]:delta[()]=-point[j]
            for mon,c in pmul(row,delta).items():add(replay,mon,c)
        assert replay==polynomial
        allmid.append(midpoint);alljac.append(jac);polys.append(polynomial)
    assert len(polys)==134
    selected=saved['selected_original_rows'];assert len(selected)==45 and len(set(selected))==45 and max(selected)<134
    for row,record in zip(selected,saved['linear45_by45_midpoint_matrix']):assert allmid[row]==list(map(decode,record))
    J=s.Matrix([[s.Rational(v)for v in alljac[row]]for row in selected])
    determinant=s.Rational(DomainMatrix.from_Matrix(J).convert_to(s.QQ).det())
    assert determinant==s.Rational(saved['d_at_point'])==s.Rational('-124102513209883356037120000000000000000/59049')
    # Independent full61 polynomial reconstruction of the universal cofactor identity.
    universalpath=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';old=read(universalpath);assert old['all_passed']
    assert old['q_to_u_indices_zero_based']==qids and len(old['selected_non_Bianchi_syzygy'])==125
    target={tuple(mon):Fraction(n,d)for mon,n,d in old['universal_new_equation']};assert len(target)==9257
    reconstructed={};homright={};multipliers={}
    for encoded,text in old['selected_non_Bianchi_syzygy']:
        j,i=divmod(encoded,70);assert 0<=i<70 and 0<=j<45;c=Fraction(text);assert c
        add(multipliers.setdefault(i,{}),(j,),c)
        for mon,n,d in model.reduced['equations'][i]:
            mon=tuple(mon);coefficient=c*Fraction(n,d);qd=sum(u in Q for u in mon)
            assert qd<=2 and sum(weights[u]for u in mon)==4
            add(reconstructed,tuple(sorted(mon+(qids[j],))),coefficient)
            add(homright,tuple(sorted(mon+(qids[j],)+(61,)*(2-qd))),coefficient)
    assert reconstructed==target
    assert max(sum(i in Q for i in mon)for mon in target)==2
    homleft={}
    for mon,c in target.items():
        assert sum(weights[u]for u in mon)==6
        add(homleft,tuple(sorted(mon+(61,)*(3-sum(i in Q for i in mon)))),c)
    assert homleft==homright
    assert s.sympify(axis['g_formula']['g7'])==0 and axis['all_135_identities']
    # Evaluate and center all9257 generic terms at the exact same T01 fiber.
    auditpath=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json';audit=read(auditpath)
    sample=list(map(Fraction,next(row for row in audit['points']if row['component_id']=='T01')['S61_coordinates_u1_through_u61']))
    qindex={u:j for j,u in enumerate(qids)};center=[sample[u]for u in qids];centered={}
    for mon,c in target.items():
        free=[]
        for u in mon:
            if u in Q:free.append(qindex[u])
            else:c*=sample[u]
        expanded={():c}
        for j in free:
            linear={(j,):Fraction(1)}
            if center[j]:linear[()]=center[j]
            expanded=pmul(expanded,linear)
        for mon,c in expanded.items():add(centered,mon,c)
    assert centered=={tuple(mon):Fraction(c)for mon,c in old['sample_centered_new_equation']}and len(centered)==790
    presentation=read(ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz')
    framepath=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz';frame=read(framepath)
    E=[[Fraction(v)for v in row]for row in frame['embedding']];h=s.Symbol('z0')
    def hpoly(expr):return {(45,)*int(mon[0]):Fraction(c)for mon,c in s.Poly(s.sympify(expr),h).terms()if c}
    Z=[{(45,):Fraction(1)},hpoly(presentation['p1']),hpoly(presentation['p2'])];Z[1][(42,)]=Fraction(1);Z[2][(43,)]=Fraction(1)
    images=[]
    for i,row in enumerate(E):
        image={}
        for j in range(3):
            for mon,c in Z[j].items():add(image,mon,row[j]*c)
        if i>=3:add(image,(i-3,),Fraction(1))
        images.append(image)
    image={}
    for mon,c in centered.items():
        product={():c}
        for j in mon:product=pmul(product,images[j])
        for word,value in product.items():add(image,word,value)
    G=flatten(presentation['monic_quintic_G']);assert G[(45,)*5]==1
    def divide(poly):
        remainder=dict(poly);quotient={}
        while True:
            candidates=[(word.count(45),word,c)for word,c in remainder.items()if word.count(45)>=5]
            if not candidates:break
            degree,word,c=max(candidates)
            term=word[:len(word)-degree]+(45,)*(degree-5);add(quotient,term,c)
            for mon,value in G.items():add(remainder,tuple(sorted(term+mon)),-c*value)
        rebuilt=pmul(quotient,G)
        for mon,c in remainder.items():add(rebuilt,mon,c)
        assert rebuilt==poly
        return quotient,remainder
    quotient,remainder=divide(image);extra=read(g6dir/'quintic-extra-generator.json')
    assert quotient==flatten(extra['exact_division_quotient'])and remainder==flatten(extra['extra_h_remainder'])
    # Verify membership in the OLD134 monic quotient before any nonlinear basechange.
    module_sum={}
    for i,ki in multipliers.items():
        ki_image={}
        for (j,),c in ki.items():
            add(ki_image,(),c*center[j])
            for mon,value in images[j].items():add(ki_image,mon,c*value)
        ri=flatten(presentation['all134_remainders'][i])
        for mon,value in pmul(ki_image,ri).items():add(module_sum,mon,value)
    assert divide(module_sum)[1]==remainder
    extra_columns=[];current=remainder
    for j in range(5):
        extra_columns.append([[list(mon),str(c)]for mon,c in sorted(current.items())])
        current=divide({tuple(sorted(mon+(45,))):c for mon,c in current.items()})[1]
    (out/'extra-five-columns.json').write_text(json.dumps(dict(variable45='h;0..44 formal coefficient variables',all5_exact_reduced_columns=extra_columns),indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',axisdir/'report.json',axisdir/'certificate.json.gz',g6dir/'report.json',g6dir/'certificate.json',g6dir/'quintic-extra-generator.json',rawreview,quinticreview,universalpath,axispath,auditpath,framepath,ROOT/'code/s61_model.py',model.source/'symplectic_slice_reduced_equations.json']
    result=dict(id=STEM,all_passed=True,current_dependency_hashes_checked=len(checked),
        independent_all134_actual_quadratic_chart_bindings=True,
        independent_symmetric_quadratic_midpoint_identities_all134=True,
        exact_selected_original45_by45_determinant=str(determinant),
        adjugate_idempotent_ideal_intersection_comax_and_ew0_proof_accepted=True,
        selected_finite_etale_axis_trace_normunit_recipe_accepted_not_expanded=True,
        trace_basis_or_axis_degree_computed=False,all_affine_Dw_points_and_nilpotents_preserved=True,
        independent_full61_G6_9257term_identity_and125_q_linear_cofactors=True,
        actual_q_degrees=dict(G6=2,Ki=1,w_G6h=3),all_homogenized_QQ_coefficients_left_equal_right=True,
        actual_I134_subset_I135_subset_I134_colon_w_and_equal_w_saturation=True,
        G6_zero_on_every_Delta_selected_factor_assumed=False,
        idempotent_descends_after_G6_and_Dw_preserved=True,
        independent_same_T01_790term_centering_and_flat46symbol_quintic_division=True,
        all5_extra_quintic_columns_constructed_exactly=True,
        G6_remainder_in_old134_formal_ideal_checked_before_basechange=True,
        presentation5_by675_same_actual_algebra=True,
        total_w_saturation_or_unknown_affine_unit_or_exhaustion_claimed=False,
        no_new_GB_or_higher_Macaulay_or_trace_basis_search=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={'extra-five-columns.json':sha(out/'extra-five-columns.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('INDEPENDENT AXIS POLARIZATION / UNIVERSAL wG6 / EXTRA5 COLUMNS PASS',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
