#!/usr/bin/env python3
"""Independent full-source/full-ambient review of the nil31 half component.

No producer, producer helper, CAS, NumPy or legacy arithmetic imports.
All elementary sparse and Fraction matrices come from prior independent
review libraries, and the source, permutation and derivatives are rebuilt.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import re
import time

from classification_minimal_quartic_review import LaurentRing
import classification_translated_minimal_quartic_review as independent

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
PARAMETERS=tuple(f'A{i}' for i in range(6))+('u','v','t','p')


def saved_parse(ring,value):
    value=re.sub(r'\bq\b','p',value)
    for old,new in [('x1','a'),('x2','b'),('x3','c'),('x4','d')]:
        value=value.replace(old,new)
    return ring.parse(value)


def derivative(ring,polynomial,position):
    result={}
    for powers,value in polynomial.items():
        if powers[position]:
            target=list(powers)
            target[position]-=1
            result[tuple(target)]=value*powers[position]
    return result


def matrix(ring,entries):
    P=[[{} for j in range(5)] for i in range(5)]
    for (i,j),value in zip(itertools.combinations(range(1,5),2),entries):
        P[i][j],P[j][i]=value,ring.scale(value,-1)
    return P


def evaluate(ring,polynomial,point,scale=Q(1)):
    result={}
    for powers,value in polynomial.items():
        value*=scale
        for name,power in zip(ring.names[:ring.n],powers[:ring.n]):
            value*=point[name]**power
        monomial=powers[ring.n:]+(0,0)
        result[monomial]=result.get(monomial,Q(0))+value
    return {e:value for e,value in result.items() if value}


def point_matrix(ring,P,point,scale=Q(1)):
    return independent.arithmetic.antisymmetric({ij:value for ij in PAIRS
                 if (value:=evaluate(ring,P[ij[0]][ij[1]],point,scale))})


def source_review(saved,old_source,third=False):
    ring=LaurentRing(PARAMETERS,('h','a','b','c','d'))
    r=Q(-1,3)if third else Q(-1,2)
    A=ring.parse('A0*h**2+A1*a*h+A2*d*h+A3*a**2+A4*a*d+A5*d**2')
    K=ring.parse('t')if third else ring.parse('t-1/2')
    S=ring.parse('u*a+v*h+(t+1/3)*d')if third else ring.parse('u*a+v*h+(t+1/2)*d')
    g=ring.parse('p*h*a')if third else ring.parse('p*a**2')
    entries=[ring.scale(ring.parse('a*c'),r),ring.scale(ring.parse('h*a'),r),{},
             ring.add(A,ring.mul(K,ring.parse('h*b-c**2/2'))),
             ring.add(g,ring.mul(ring.variable('c'),S)),
             ring.mul(ring.variable('h'),S)]
    P=matrix(ring,entries)
    for i,j in PAIRS:
        assert P[i][j]==saved_parse(ring,saved['whole_homogeneous_h_Casimir_brackets'].get(f'{i},{j}','0'))
    strict=[ring.add(*(ring.mul(P[a][k],ring.derivative(P[b][d],k))
              for a,b,d in [(i,j,l),(j,l,i),(l,i,j)] for k in range(5)))
            for i,j,l in itertools.combinations(range(5),3)]
    assert strict==[{}]*10
    pf=ring.add(ring.mul(entries[0],entries[5]),ring.scale(ring.mul(entries[1],entries[4]),-1),
                 ring.mul(entries[2],entries[3]))
    expected_pf='p*h**2*a**2/3'if third else'p*h*a**3/2'
    expected_D='p*h**3*a**2/3'if third else'p*h**2*a**3/2'
    assert pf==ring.parse(expected_pf)==saved_parse(ring,saved['actual_affine_homogeneous_pfaffian'])
    assert ring.mul(ring.variable('h'),pf)==ring.parse(expected_D)

    # Reconstruct the entire old15-source, homogenize independently, and
    # explicitly apply the projective permutation before Euler subtraction.
    oldring=LaurentRing(tuple(f'n{i}' for i in range(15)),('H','X','Y','Z','W'))
    binding=[ring.variable('A3'),ring.variable('A0'),ring.variable('A2'),ring.variable('A5'),ring.variable('A1'),
             {}if third else ring.variable('p'),ring.variable('p')if third else {},{},ring.variable('A4'),{},ring.constant(r),{},
             ring.variable('u'),ring.variable('v'),ring.variable('t')]
    assert len(saved['old_full_source_parameter_substitution'])==15
    for i,value in enumerate(binding):
        assert value==saved_parse(ring,saved['old_full_source_parameter_substitution'][f'n{i}'])
    target_coordinates=[ring.variable(name) for name in ('a','b','c','h','d')]
    oldP=[[{} for _ in range(5)]for _ in range(5)]
    for ij,expression in zip(itertools.combinations(range(1,5),2),old_source['full_six_brackets']):
        for old,new in [('x1','X'),('x2','Y'),('x3','Z'),('x4','W')]:
            expression=expression.replace(old,new)
        affine=oldring.parse(expression)
        polynomial={}
        for powers,value in affine.items():
            coordinate_powers=list(powers[15:])
            assert coordinate_powers[0]==0 and sum(coordinate_powers)<=2
            coordinate_powers[0]=2-sum(coordinate_powers)
            term=ring.constant(value)
            for factor,power in zip(binding,powers[:15]):
                for _ in range(power):term=ring.mul(term,factor)
            for factor,power in zip(target_coordinates,coordinate_powers):
                for _ in range(power):term=ring.mul(term,factor)
            polynomial=ring.add(polynomial,term)
        oldP[ij[0]][ij[1]]=polynomial
        oldP[ij[1]][ij[0]]=ring.scale(polynomial,-1)
    order=[3,0,1,2,4]
    assert (-1)**sum(order[i]>order[j] for i,j in itertools.combinations(range(5),2))==-1
    gauge=[{},{},ring.scale(ring.variable('c'),-r),ring.scale(ring.variable('h'),-r),{}]
    assert gauge==[saved_parse(ring,value)for value in saved['Euler_gauge_vector']]
    E=[ring.variable(name)for name in ('h','a','b','c','d')]
    for i,j in PAIRS:
        transport=oldP[order[i]][order[j]]
        actual=ring.add(transport,ring.scale(ring.mul(E[i],gauge[j]),-1),ring.mul(gauge[i],E[j]))
        assert actual==P[i][j]
    # Divergence before and after gauge for every source parameter. The
    # intrinsic triple hyperplane is a in the half branch and h in third.
    old_coordinate_positions=[1,2,3,0,4]
    old_modular=[ring.add(*(ring.derivative(oldP[i][j],old_coordinate_positions[i])for i in range(5)))for j in range(5)]
    assert old_modular==[{},ring.variable('c'),ring.variable('h'),{},{}]
    new_modular=[ring.add(*(ring.derivative(P[i][j],i)for i in range(5)))for j in range(5)]
    N=[{},{},ring.variable('c'),ring.variable('h'),{}]
    assert new_modular==[ring.scale(value,5*r+1)for value in N]

    point=dict(zip(PARAMETERS,map(Q,[1,2,3,5,7,11,2,3,5,13])))
    assert point=={('p'if key=='q'else key):Q(value)for key,value in saved['witness_parameter_values'].items()}
    witness=point_matrix(ring,P,point,Q(12))
    expected=independent.arithmetic.antisymmetric({ij:value for ij in PAIRS
                   if(value:=evaluate(ring,saved_parse(ring,saved['witness_brackets'].get(f'{ij[0]},{ij[1]}','0')),point))})
    assert witness==expected
    actual_derivatives=[]
    for position,name in enumerate(PARAMETERS):
        data={}
        for i,j in PAIRS:
            value=evaluate(ring,derivative(ring,P[i][j],position),point,Q(12))
            expected=evaluate(ring,saved_parse(ring,saved['all_ten_actual_parameter_derivatives'][position].get(f'{i},{j}','0')),point)
            assert value==expected
            if value:data[i,j]=value
        actual_derivatives.append(independent.arithmetic.antisymmetric(data))
    return witness,actual_derivatives,dict(full_parameter_count=10,source_irreducible=True,
        all_six_full_source_brackets_reconstructed=True,all_ten_strict_Jacobi_zero=True,
        universal_pfaffian=expected_pf,universal_projective_divisor=expected_D,
        all15_old_parameters_and_all150_PGL_Euler_coefficients_checked=True,
        projective_permutation_determinant='-1',gauge=['0','0',str(-r)+'*c',str(-r)+'*h','0'],
        whole_old_triple_a_modular_field=['0','0','c','h','0'],
        whole_new_h_modular_field_scalar=str(5*r+1),whole_new_h_modular_field_direction=['0','0','c','h','0'],
        all_ten_actual_derivatives_reconstructed=True,auxiliary_kappa_direction_used=False,
        p_nonzero_rank_four_complement_retained=True,p_zero_boundary_in_same_polynomial_source=True)


def ambient_review(point,directions,saved,third=False):
    D,G,sparse_rows=independent.matrices(point)
    V=[row[:25]for row in G]
    columns=[independent.vector(value)for value in directions]
    Fmatrix=[row+[col[i]for col in columns]for i,row in enumerate(G)]
    assert all(sum(value*col[j]for j,value in row)==0 for row in sparse_rows for col in columns)
    certificates=saved['certificates']
    ranks=dict(D=98 if third else 100,vertical=24,orbit_vertical_scale=48 if third else 47,whole_source=50)
    matrices=dict(D=D,vertical=V,orbit_vertical_scale=G,whole_source=Fmatrix)
    checked={}
    for name,M in matrices.items():
        certificate=certificates[name]
        assert certificate['rank_QQ']==ranks[name]
        determinant=independent.minor(M,certificate['rows'],certificate['columns'])
        assert determinant==Q(certificate['determinant']) and determinant
        encoded=json.dumps(M,default=str,separators=(',',':')).encode()
        assert hashlib.sha256(encoded).hexdigest()==saved['full_matrix_sha256'][name]
        if name in ['vertical','orbit_vertical_scale']:
            assert independent.arithmetic.rank(M)[0]==ranks[name]
        checked[name]=dict(shape=[len(M),len(M[0])],rank_QQ=ranks[name],minor_determinant=str(determinant),
                           all_matrix_coefficients_match=True)
    if third:
        assert independent.arithmetic.rank(Fmatrix)[0]==50
        assert independent.arithmetic.rank(D)[0]==98
    assert len(Fmatrix)==150 and len(Fmatrix[0])==61
    result=dict(all150_variables_and350_equations=True,all61_actual_columns_killed_by_D=True,
                independent_matrices=checked,full_D_rank=ranks['D'],whole_source_rank=50,
                projective_tangent_dimension=27 if third else 25,actual_source_saturation_dimension_lower=25,
                smooth_unique_full_component_dimension=None if third else 25,
                full_source_closure_is_the_entire_component=not third)
    return (D,G,Fmatrix,result)if third else result


def catalogue_review(snapshot):
    rows=snapshot['components']
    assert len(rows)==92 and len({row['id']for row in rows})==92
    assert all(type(row['projective_dimension'])is int for row in rows)
    same={row['id']:row for row in rows if row['projective_dimension']==25}
    expected={'T01','T39','T40','R4_POSITIVE_WEIGHT6_25','R4_MODULAR_NIL31_DOUBLE_DOUBLE_H_25'}
    assert set(same)==expected
    for key in ['T01','T39','T40']:
        assert same[key]['meets_SNC_log_symplectic_locus'] is True
    positive=same['R4_POSITIVE_WEIGHT6_25']['generic_actual_degeneracy_divisor']
    assert positive['irreducible_factor_degrees']==[1,2] and positive['corresponding_multiplicities']==[3,1]
    assert 'rank5' in positive['description']
    cn=same['R4_MODULAR_NIL31_DOUBLE_DOUBLE_H_25']['generic_actual_degeneracy_divisor']
    assert cn['irreducible_factor_degrees']==[1,1,1] and cn['corresponding_multiplicities']==[2,2,1]
    assert cn['hyperplane_span_dimension']==3
    exclusions={row['id']:('different full-component dimension'if row['projective_dimension']!=25
                else 'closed H³L² divisor locus contains no SNC point'if row['id']in ['T01','T39','T40']
                else 'generic residual smooth rank5 quadric is not a square'if row['id']=='R4_POSITIVE_WEIGHT6_25'
                else 'generic three distinct 2,2,1 factors contain no triple factor')for row in rows}
    return dict(prior_component_count=92,all_prior_dimensions_present=True,same_dimension_components=sorted(same),
                exclusions_for_every_prior_component=exclusions,
                closed_factor_image='P(V*)×P(V*) -> P(Sym⁵V*), (H,L) -> H³L²',
                includes_zero_degeneracy_by_affine_cone=True,
                nonsquare_free_divisor_condition_preserved_on_entire_saturated_source_closure=True,
                proper_closed_incidence_has_no_SNC_divisor=True,
                distinct_from_every_prior92_component=True,
                mathematical_scope='Dimensions and adopted generic divisor/SNC properties come from the audited prior92catalogue. Properness of the multiplication map and UFD factorization give the stated distinctions after the independent full-component proof. No claim of publication novelty.')


def third_quadratic_review(point,directions,saved,component):
    D,G,Fmatrix,ambient=ambient_review(point,directions,saved,third=True)
    ring=LaurentRing(PARAMETERS,('h','a','b','c','d'))
    values={name:Q(0)for name in PARAMETERS}
    normals=[]
    for encoding in component['normals']:
        data={}
        for i,j in PAIRS:
            value=evaluate(ring,saved_parse(ring,encoding.get(f'{i},{j}','0')),values)
            if value:data[i,j]=value
        normals.append(independent.arithmetic.antisymmetric(data))
    assert len(normals)==2
    vectors=[independent.vector(N)for N in normals]
    for vector in vectors:
        assert all(sum(a*b for a,b in zip(row,vector))==0 for row in D)
    FN=[row+[column[i]for column in vectors]for i,row in enumerate(Fmatrix)]
    assert independent.arithmetic.rank(FN)[0]==52
    coordinates=component['literal_slice_coordinate_columns']
    rows=component['implicit_D_rows']
    assert len(coordinates)==len(rows)==98 and len(set(coordinates))==98
    determinant=independent.minor(D,rows,coordinates)
    assert determinant==Q(component['implicit_D_minor'])and determinant
    slice_columns=[[Q(i==column)for i in range(150)]for column in coordinates]+vectors
    full=[row+[column[i]for column in slice_columns]for i,row in enumerate(Fmatrix)]
    assert independent.arithmetic.rank(full)[0]==150
    assert independent.arithmetic.rank([list(row)for row in zip(*slice_columns)])[0]==100
    q0=[value.get((0,0),Q(0))for value in independent.jacobi_vector(normals[0])]
    q2=[value.get((0,0),Q(0))for value in independent.jacobi_vector(normals[1])]
    added=independent.arithmetic.antisymmetric({ij:independent.arithmetic.add(normals[0].get(ij,{}),normals[1].get(ij,{}))for ij in PAIRS})
    both=[value.get((0,0),Q(0))for value in independent.jacobi_vector(added)]
    q1=[a-b-c for a,b,c in zip(both,q0,q2)]
    Qcolumns=[q0,q1,q2]
    pairings=[]
    for encoding in component['normalized_left_functionals']:
        functional={int(i):Q(value)for i,value in encoding.items()}
        assert all(sum(value*D[i][column]for i,value in functional.items())==0 for column in range(150))
        pairings.append([sum(value*column[i]for i,value in functional.items())for column in Qcolumns])
    A=Q(113435,3825851)
    B=Q(214061,7651702)
    assert pairings==[[Q(1),Q(0),-A],[Q(0),Q(1),B]]
    coefficient=B*B-A
    assert coefficient==Q(component['radical_nonzero_coefficient'])==Q(-1690119521019,58548543496804)
    assert coefficient
    expression=f'v*(u**2-({A})*v**2)-(u-({B})*v)*(u*v+({B})*v**2)'
    assert independent.parse(expression)==independent.parse(f'({coefficient})*v**3')
    ambient.update(full_kernel_dimension=52,genuine_normal_count=2,F_plus_normals_rank=52,
        literal_slice_dimension=100,actual_independent_linear_cut_count=50,
        F_plus_literal_slice_rank=150,implicit_minor_determinant=str(determinant),
        both_normalized_left_functionals_kill_all150_D_columns=True,
        all_three_projected_quadratic_coefficients=[[str(value)for value in row]for row in pairings],
        radical_certificate_nonzero_coefficient=str(coefficient),
        quadratic_initial_ideal_has_maximal_radical=True,true_formal_cut_has_dimension_zero=True,
        full_local_projective_dimension=25,full_source_closure_is_the_entire_component=True,
        point_smoothness_or_local_unique_component_claimed=False)
    return ambient


def half_third_distinction():
    # Both projective cofactor directions have this same nilpotent N, but the
    # intrinsic triple hyperplane has a different position in its dual flag.
    N=[[Q(0)for _ in range(5)]for _ in range(5)]
    N[2][3]=N[3][0]=Q(1)
    square=[[sum(N[i][k]*N[k][j]for k in range(5))for j in range(5)]for i in range(5)]
    cube=[[sum(square[i][k]*N[k][j]for k in range(5))for j in range(5)]for i in range(5)]
    assert all(not value for row in cube for value in row)
    assert independent.arithmetic.rank(N)[0]==2 and independent.arithmetic.rank(square)[0]==1
    covector_a=[Q(i==1)for i in range(5)]
    covector_h=[Q(i==0)for i in range(5)]
    transpose=[list(row)for row in zip(*N)]
    transpose_square=[list(row)for row in zip(*square)]
    assert independent.arithmetic.rank([row+[value]for row,value in zip(transpose,covector_a)])[0]==3
    assert independent.arithmetic.rank([row+[value]for row,value in zip(transpose,covector_h)])[0]==2
    assert independent.arithmetic.rank([row+[value]for row,value in zip(transpose_square,covector_h)])[0]==1
    assert all(sum(a*b for a,b in zip(row,covector_a))==0 for row in N)
    return dict(whole_dense_source_triple_hyperplane_flags_verified=True,
        half_intrinsic_triple_hyperplane='a=0, equivalently old_h=0',
        third_intrinsic_triple_hyperplane='h=0, equivalently old_z=0',
        common_projective_cofactor_direction=[[str(v)for v in row]for row in N],
        half_triple_covector_not_in_image_N_transpose=True,
        third_triple_covector_in_image_N_transpose_squared=True,
        half_affine_modular_field_has_a_zero=True,third_affine_modular_field_has_no_zero=True,
        half_affine_type='nil31,constant0',third_affine_type='nil211,long-top nonzero translation',
        distinction_argument='Each entire nonzero rank4 source has its stated cofactor flag at the unique triple hyperplane. PGL/scalar equivalence preserves this intrinsic pair. If the two same-dimensional component closures coincided, their dense constructible saturated images would contain intersecting nonempty opens, contradicting the two flag values. No global closedness of the fixed-rank image flag is assumed.',
        half_and_third_components_distinct=True)


def review_third(output):
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic()
    producer_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json'
    component_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_THIRD_COMPONENT_20261107.json'
    old_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
    catalogue_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_PRIOR92_20261107.json'
    packet=json.loads(producer_path.read_text())
    component=json.loads(component_path.read_text())
    old=json.loads(old_path.read_text())['cases'][0]
    # Recheck both universal modular fields for every source parameter before
    # using the cofactor flag to distinguish the two components.
    _,_,half_source=source_review(packet['branches'][0],old)
    point,directions,source=source_review(packet['branches'][1],old,third=True)
    ambient=third_quadratic_review(point,directions,packet['branches'][1],component)
    distinction=half_third_distinction()
    catalogue=catalogue_review(json.loads(catalogue_path.read_text()))
    inputs=[Path(__file__),Path(independent.__file__),ROOT/'code/classification_minimal_quartic_review.py',
        ROOT/'code/classification_nil22_weighted_slice_review.py',producer_path,component_path,old_path,catalogue_path,
        ROOT/'research/classification/KNOWN_COMPONENT_AUDIT.md',ROOT/'research/classification/MODULAR_NIL31_COMPONENT.md']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_or_helper_import=True,no_CAS_or_NumPy_import=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},
        full_source=source,full_ambient_and_true_slice=ambient,whole_half_source_modular_check=half_source,
        prior92_distinctness=catalogue,half_third_intrinsic_distinction=distinction,
        componenthood_mathematical_review='The actual fullsource rank50 supplies a25-dimensional projective image lower bound. A literal100-dimensional complement to that actual50-space is checked in the entireC150. A nonzero98derivative minor solves98directions; the two functionals kill all150Dcolumns, so their pairings are genuine quadratic initial forms of this actualslice. The exact nonzero B²-A identity makes the normal initialsupport only0, hence the cut dimension0. Krull bounds the full localdimension by50 before quotient. Equality with the actualsource lower bound makes the irreducible whole-source saturation closure a25fullcomponent. This does not imply smoothness, local uniqueness, or source-normal product decomposition.',
        accepted_new_component_relative_to_prior92_and_half=True,projective_component_dimension=25,
        full_projective_tangent_dimension=27,point_smoothness_or_uniqueness_claimed=False,
        scope='Independent full-source and true-slice proof of the third25component, distinct from half and all prior92. Global classification and exact scheme thickness remain unproved.',
        reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun,
        maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x')as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({key:result[key]for key in ['all_passed','projective_component_dimension','seconds','maximum_RSS_KiB']}))


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic()
    producer_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json'
    old_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
    catalogue_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_PRIOR92_20261107.json'
    saved=json.loads(producer_path.read_text())['branches'][0]
    assert saved['fixed_r']=='-1/2'
    point,directions,source=source_review(saved,json.loads(old_path.read_text())['cases'][0])
    ambient=ambient_review(point,directions,saved)
    catalogue=catalogue_review(json.loads(catalogue_path.read_text()))
    inputs=[Path(__file__),Path(independent.__file__),ROOT/'code/classification_minimal_quartic_review.py',
            ROOT/'code/classification_nil22_weighted_slice_review.py',producer_path,old_path,catalogue_path,
            ROOT/'research/classification/KNOWN_COMPONENT_AUDIT.md',ROOT/'research/classification/MODULAR_NIL31_COMPONENT.md']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_import=True,no_producer_helper_import=True,
                no_CAS_or_NumPy_import=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},
                full_source=source,full_ambient=ambient,prior92_distinctness=catalogue,
                componenthood_mathematical_review='The ten-parameter polynomial source and connected PGL/vertical/scalar saturation are irreducible. Its actual tangent rank50 and the entire350x150 derivative rank100 give the same dimension50 before quotient, hence a regular point and a unique projective25component. The irreducible whole-source image closure has that dimension and equals this component, including p=0 boundaries. No containment in CA/CB is inferred from a singular special-point intersection.',
                accepted_new_component_relative_to_prior92=True,projective_dimension=25,generic_rank=4,
                generic_actual_divisor='H³L², two distinct hyperplanes',
                third_r_minus_one_third_branch_componenthood_claimed=False,
                scope='Independent acceptance of the entire half source as a new smooth25fullcomponent distinct from prior92. Global classification, scheme structure and r=-1/3branch are not proved.',
                reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,
                maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x')as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({key:result[key]for key in ['all_passed','projective_dimension','seconds','maximum_RSS_KiB']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--third',action='store_true',help='Review the complete third source, genuine two-normal support, and distinction from half')
    args=parser.parse_args()
    (review_third if args.third else main)(args.output)
