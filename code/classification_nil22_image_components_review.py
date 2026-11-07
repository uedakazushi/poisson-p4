#!/usr/bin/env python3
"""Independent all-source and genuine-cut review of two nil22 components.

Standard-library Fraction sparse arithmetic only. No producer, producer
helper, CAS, NumPy or legacy imports. This does not classify the other nine
nil22 image sources.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

from classification_minimal_quartic_review import LaurentRing
import classification_translated_minimal_quartic_review as independent

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2))
PARAMETERS=tuple(f'm{i}'for i in range(6))


def parse(ring,text):
    for old,new in [('x1','q1'),('x2','q2'),('x3','x'),('x4','y')]:
        text=re.sub(r'\b'+old+r'\b',new,text)
    return ring.parse(text)


def evaluate(ring,value,point,scale=Q(1)):
    answer={}
    for powers,coefficient in value.items():
        coefficient*=scale
        for name,power in zip(ring.names[:ring.n],powers[:ring.n]):coefficient*=point[name]**power
        exponent=powers[ring.n:]+(0,0)
        answer[exponent]=answer.get(exponent,Q(0))+coefficient
    return {e:c for e,c in answer.items()if c}


def at_point(ring,entries,point,scale=Q(1)):
    return independent.arithmetic.antisymmetric({ij:value for ij in PAIRS if(value:=evaluate(ring,entries.get(ij,{}),point,scale))})


def derivative_parameter(value,position):
    answer={}
    for powers,c in value.items():
        if powers[position]:
            new=list(powers);new[position]-=1
            answer[tuple(new)]=c*powers[position]
    return answer


def full_source(saved,index):
    ring=LaurentRing(PARAMETERS,('h','q1','q2','x','y'))
    h,q1,q2,x,y=(ring.variable(a)for a in ('h','q1','q2','x','y'))
    A=index==4
    # Reconstruct M=-adj(S+Hess(x^3)/2)+p tensor V, then the entire F0
    # polynomial and the two q coefficients from the actual divergences.
    texts=['6*h**2+h*x-x**2/3','-x*y/2','h*y-x*y/3','18*h**2-3*h*x-y**2/2']if A else ['-x**2/3','h**2-2*x*y/3','h**2-x*y/3','-3*h*x-2*y**2/3']
    M=[ring.parse(text)for text in texts]
    F0=ring.parse('m0*h**2+m1*h*x+m2*h*y+m3*x**2+m4*x*y+m5*y**2')
    div_second=ring.add(ring.derivative(M[2],3),ring.derivative(M[3],4))
    div_first=ring.add(ring.derivative(M[0],3),ring.derivative(M[1],4))
    F=ring.add(F0,ring.mul(q1,ring.add(y,div_second)),ring.scale(ring.mul(q2,ring.add(x,div_first)),-1))
    entries={(1,2):F,(1,3):M[0],(1,4):M[1],(2,3):M[2],(2,4):M[3]}
    for ij in PAIRS:assert entries.get(ij,{})==parse(ring,saved['full_six_magnetic_source'].get(f'{ij[0]},{ij[1]}','0'))
    scale=Q(6 if A else 3);assert scale==Q(saved['source_integer_scale'])
    for ij in PAIRS:assert ring.scale(entries.get(ij,{}),scale)==parse(ring,saved['universal_brackets'].get(f'{ij[0]},{ij[1]}','0'))
    P=[[{}for j in range(5)]for i in range(5)]
    for i,j in PAIRS:P[i][j],P[j][i]=entries.get((i,j),{}),ring.scale(entries.get((i,j),{}),-1)
    strict=[ring.add(*(ring.mul(P[a][k],ring.derivative(P[b][c],k))for a,b,c in [(i,j,l),(j,l,i),(l,i,j)]for k in range(5)))for i,j,l in it.combinations(range(5),3)]
    assert strict==[{}]*10
    modular=[ring.add(*(ring.derivative(P[i][j],i)for i in range(1,5)))for j in range(1,5)]
    assert modular==[x,y,{},{}]
    pf=ring.add(ring.mul(M[1],M[2]),ring.scale(ring.mul(M[0],M[3]),-1))
    cubic=ring.parse('x**3-9*h*x**2-3*h*y**2+108*h**3'if A else'x**3+h*x*y-h**3')
    assert pf==ring.scale(ring.mul(h,cubic),-1)
    D=ring.mul(h,pf)
    assert D==parse(ring,saved.get('whole_source_actual_quintic','-h**2*(x**3-9*h*x**2-3*h*y**2+108*h**3)'if A else'-h**2*(x**3+h*x*y-h**3)'))
    point=dict(zip(PARAMETERS,map(Q,[1,2,3,5,7,11])))
    assert point=={name:Q(v)for name,v in saved['full_source_specification']['witness'].items()}
    witness=at_point(ring,entries,point,scale)
    expected=at_point(ring,{ij:parse(ring,saved['witness_brackets'].get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS},point)
    assert witness==expected
    directions=[at_point(ring,{ij:derivative_parameter(value,i)for ij,value in entries.items()},point,scale)for i in range(6)]
    assert independent.arithmetic.rank([independent.vector(p)for p in directions])[0]==6
    return ring,entries,witness,directions,dict(full_six_magnetic_source_reconstructed=True,
        all_ten_strict_Jacobi_zero=True,whole_modular_field=['x','y','0','0'],
        actual_quintic='-h^2*(x^3-9hx^2-3hy^2+108h^3)'if A else'-h^2*(x^3+hxy-h^3)',
        all_six_actual_parameter_derivatives=True,parameter_domain='A^6, no magnetic localization',
        full_source_irreducible=True,source_polynomial_parameter_rank=6,
        generic_global_Poisson_rank=4,all_parameter_boundaries_retained=True)


def certify(M):
    rank,columns=independent.arithmetic.rank(M)
    rows=independent.arithmetic.rank([[row[j]for row in M]for j in columns])[1]
    determinant=independent.minor(M,rows,columns)
    assert determinant
    return dict(shape=[len(M),len(M[0])],rank=rank,rows=rows,columns=columns,determinant=str(determinant),
        full_matrix_sha256=hashlib.sha256(json.dumps(M,default=str,separators=(',',':')).encode()).hexdigest())


def true_cut(ring,point,directions,source,component,index):
    A=index==4;d_rank=101 if A else 100;normal_count=1 if A else 2
    D,G,sparse_D=independent.matrices(point)
    V=[row[:25]for row in G]
    vectors=[independent.vector(p)for p in directions]
    F=[row+[col[i]for col in vectors]for i,row in enumerate(G)]
    certificates={name:certify(M)for name,M in [('D',D),('vertical',V),('GL5_vertical_scalar',G),('whole_source',F)]}
    assert [certificates[a]['rank']for a in certificates]==[d_rank,24,48,48]
    assert all(sum(v*column[j]for j,v in row)==0 for row in sparse_D for column in vectors)
    values={name:Q(0)for name in PARAMETERS}
    normals=[at_point(ring,{ij:parse(ring,encoding.get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS},values)for encoding in component['normals']]
    assert len(normals)==normal_count
    normal_vectors=[independent.vector(value)for value in normals]
    assert all(sum(v*column[j]for j,v in row)==0 for row in sparse_D for column in normal_vectors)
    FN=[row+[column[i]for column in normal_vectors]for i,row in enumerate(F)]
    assert independent.arithmetic.rank(FN)[0]==150-d_rank
    coordinates=component['literal_slice_coordinate_columns'];rows=component['implicit_D_rows']
    assert len(coordinates)==len(rows)==d_rank
    determinant=independent.minor(D,rows,coordinates)
    assert determinant==Q(component['implicit_D_minor'])and determinant
    cut_columns=[[Q(i==j)for i in range(150)]for j in coordinates]+normal_vectors
    cut=[list(row)for row in zip(*cut_columns)]
    assert len(cut[0])==102 and independent.arithmetic.rank(cut)[0]==102
    assert independent.arithmetic.rank([a+b for a,b in zip(F,cut)])[0]==150
    q0=[p.get((0,0),Q(0))for p in independent.jacobi_vector(normals[0])]
    quadratic_columns=[q0]
    if not A:
        q2=[p.get((0,0),Q(0))for p in independent.jacobi_vector(normals[1])]
        added=independent.arithmetic.antisymmetric({ij:independent.arithmetic.add(normals[0].get(ij,{}),normals[1].get(ij,{}))for ij in PAIRS})
        both=[p.get((0,0),Q(0))for p in independent.jacobi_vector(added)]
        quadratic_columns.extend([[a-b-c for a,b,c in zip(both,q0,q2)],q2])
    pairings=[]
    for encoding in component['normalized_left_functionals']:
        functional={int(i):Q(v)for i,v in encoding.items()}
        assert all(sum(v*D[i][j]for i,v in functional.items())==0 for j in range(150))
        pairings.append([sum(v*column[i]for i,v in functional.items())for column in quadratic_columns])
    expected=[[Q(1)]]if A else[[Q(1),Q(0),-Q(3,124)],[Q(0),Q(1),Q(5,31)]]
    assert pairings==expected
    if not A:
        coefficient=Q(5,31)**2-Q(3,124)
        assert coefficient==Q(7,3844)==Q(component['radical_nonzero_coefficient'])
        identity=independent.parse('v*(u**2-3*v**2/124)-(u-5*v/31)*(u*v+5*v**2/31)-7*v**3/3844')
        assert not identity
    return dict(independent_full_matrices=certificates,all150_variables_and350_equations=True,
        full_D_times_actual_source_and_normals_zero=True,full_kernel_dimension=150-d_rank,
        F_plus_normal_rank=150-d_rank,genuine_normal_count=normal_count,
        literal_true_slice_dimension=102,independent_actual_cut_count=48,F_plus_literal_slice_rank=150,
        implicit_D_minor=str(determinant),all_functionals_kill_all150_D_columns=True,
        all_normal_quadratic_pairings=[[str(v)for v in row]for row in pairings],
        radical_nonzero_coefficient='1'if A else'7/3844',true_normal_initial_support_is_origin=True,
        full_ambient_local_dimension=48,whole_irreducible_saturated_source_is_a_full_component=True,
        projective_component_dimension=23,projective_tangent_dimension=125-d_rank,
        smoothness_or_local_component_uniqueness_claimed=False)


def geometry():
    ring=LaurentRing((),('h','x','y'))
    f=ring.parse('x**3-9*h*x**2-3*h*y**2+108*h**3');g=ring.parse('x**3+h*x*y-h**3')
    assert f==ring.parse('(x-6*h)**2*(x+3*h)-3*h*y**2')
    gradients=[[ring.derivative(p,i)for i in range(3)]for p in [f,g]]
    def restrict(p,values):
        answer={}
        for powers,c in p.items():
            remaining=list(powers)
            for i,a in values.items():c*=a**powers[i];remaining[i]=0
            e=tuple(remaining);answer[e]=answer.get(e,Q(0))+c
        return {e:c for e,c in answer.items()if c}
    at_infinity=[restrict(a,{0:0,1:0,2:1})for a in gradients[0]]
    assert at_infinity==[ring.constant(-3),{},{}]
    assert all(not restrict(a,{0:1,1:6,2:0})for a in gradients[0])
    assert all(not restrict(a,{0:0,1:0,2:1})for a in gradients[1])
    tangent={e:c for e,c in restrict(g,{2:1}).items()if e[0]+e[1]==2}
    assert tangent==ring.parse('h*x')
    hessian=ring.add(ring.mul(ring.derivative(gradients[0][1],1),ring.derivative(gradients[0][2],2)),
        ring.scale(ring.mul(ring.derivative(gradients[0][1],2),ring.derivative(gradients[0][2],1)),-1))
    assert restrict(hessian,{0:1,1:6,2:0})==ring.constant(-108)
    # Independence of the three partial derivatives proves the essential
    # variable count3 and makes the cubic cone vertex intrinsic.
    for grad in gradients:
        monomials=sorted(set().union(*(set(a)for a in grad)))
        assert independent.arithmetic.rank([[a.get(e,Q(0))for e in monomials]for a in grad])[0]==3
    return dict(A_cubic_irreducible_by_quadratic_discriminant='(x-6)^2*(x+3)/3 is not a square in C(x)',
        B_cubic_irreducible_by_Gauss='gcd(h*x,x^3-h^3)=1; primitive linear polynomial in y',
        A_unique_plane_cubic_node='[1:6:0]',B_unique_plane_cubic_node='[0:0:1]',
        A_node_Hessian_determinant='-108',B_node_tangent_cone='h*x',
        essential_variable_count_for_both=3,whole_source_double_H_unique=True,
        A_cubic_cone_singular_plane_contained_in_double_H=False,
        B_cubic_cone_singular_plane_contained_in_double_H=True,
        all_six_magnetic_parameters_leave_actual_quintic_unchanged=True,
        pairwise_distinct_full_components=True,
        proof='The singular-plane/double-H incidence is intrinsic for every point of each saturated constructible dense source image. Equality of the two component closures would force these two dense images to intersect on nonempty opens, contradicting the opposite incidence values.')


def catalogue(snapshot):
    rows=snapshot['components']
    assert len(rows)==94 and len({row['id']for row in rows})==94
    assert all(type(row['projective_dimension'])is int for row in rows)
    exclusions={};same=[]
    for row in rows:
        name=row['id']
        if row['projective_dimension']!=23:
            exclusions[name]='different full-component dimension after componenthood';continue
        same.append(name)
        if row['meets_SNC_log_symplectic_locus']:
            exclusions[name]='closed squared-hyperplane divisor locus contains no reduced SNC divisor';continue
        divisor=row['generic_actual_degeneracy_divisor']
        assert divisor
        if name=='R4_ZERO_LOG_5H':
            assert divisor['irreducible_factor_degrees']==[1]and divisor['corresponding_multiplicities']==[5]
            exclusions[name]='prior whole closure satisfies closed D=H^5; new source has irreducible cubic factor'
        else:
            assert divisor['reduced']is True
            assert all(a==1 for a in divisor['corresponding_multiplicities'])
            exclusions[name]='prior general divisor is reduced, violating the new entire-closure closed D=H^2*F3 condition'
    assert len(exclusions)==94
    return dict(prior_count=94,all_prior_dimensions_present=True,same23_components=sorted(same),
        exclusions_for_every_prior_component=exclusions,relative_distinctness_from_all_prior94=True,
        new_entire_closure_closed_factor_image='P(V*) x P(Sym^3 V*) -> P(Sym^5 V*), (H,C) -> H^2 C',
        prior_zero_log_closed_factor_image='P(V*) -> P(Sym^5 V*), H -> H^5',
        zero_quintic_allowed_by_affine_cone=True,
        prior_componenthood_and_generic_divisor_data='Adopted from the audited immutable prior94 snapshot; not independently reproved here.',
        publication_novelty_claimed=False)


def substitute_sparse(oldring,polynomial,ring,parameter_values,coordinate_values,homogenize=False):
    answer={}
    for powers,coefficient in polynomial.items():
        spatial=list(powers[oldring.n:])
        if homogenize:
            assert not spatial[0] and sum(spatial)<=2
            spatial[0]=2-sum(spatial)
        term=ring.constant(coefficient)
        for value,power in zip(parameter_values,powers[:oldring.n]):
            for _ in range(power):term=ring.mul(term,value)
        for value,power in zip(coordinate_values,spatial):
            for _ in range(power):term=ring.mul(term,value)
        answer=ring.add(answer,term)
    return answer


def matrix_entries(ring,upper):
    P=[[{}for _ in range(5)]for _ in range(5)]
    for i,j in PAIRS:P[i][j],P[j][i]=upper.get((i,j),{}),ring.scale(upper.get((i,j),{}),-1)
    return P


def strict_and_projective(ring,P):
    J={triple:ring.add(*(ring.mul(P[a][k],ring.derivative(P[b][c],k))for a,b,c in [triple,triple[1:]+triple[:1],triple[2:]+triple[:2]]for k in range(5)))for triple in it.combinations(range(5),3)}
    assert all(not value for value in J.values())
    for quad in it.combinations(range(5),4):
        wedge=ring.add(*(ring.scale(ring.mul(ring.variable(ring.names[ring.n+quad[i]]),J[quad[:i]+quad[i+1:]]),(-1)**i)for i in range(4)))
        assert not wedge


def limits_review(output,double_d=False):
    begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
    limit_path=ROOT/('reports/CLASSIFICATION_NIL22_IMAGE_DOUBLE_D_LIMIT_20261107.json'if double_d else'reports/CLASSIFICATION_NIL22_IMAGE_COMPONENT_LIMITS_20261107.json')
    nil31_path=ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json'
    old_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json'
    sources=json.loads(source_path.read_text())['cases'];packet=json.loads(limit_path.read_text());limits=[packet]if double_d else packet['cases']
    branches=json.loads(nil31_path.read_text())['branches'];old=json.loads(old_path.read_text())['cases'][0]
    parameters=PARAMETERS+(('a','u','e')if double_d else('n','c','d','e'))
    ring=LaurentRing(parameters,('h','q1','q2','x','y'),invert=('e',))
    spatial=[ring.variable(a)for a in ('h','q1','q2','x','y')]
    records=[]
    eposition=parameters.index('e')
    for index,item in zip([3]if double_d else[8,6],limits):
        half=index==8
        # Entire original9-source reconstructed from its canonical S,V and
        # unrestricted F0, before changing the infinity hyperplane.
        Mp=['a*x**2','-(2*a+1)*x*(y+u*h/2)','h*x+a*x*y','-u*h**2-2*h*y-(2*a+1)*y*(y+u*h/2)']if double_d else ['-x**2/2','x*(n*h+c*x+d*y)','-x*y/2','-h**2+y*(n*h+c*x+d*y)']if half else['-x**2/3','x*(n*h+c*x+d*y)','-x*y/3','-3*h*x+y*(n*h+c*x+d*y)']
        M=[ring.parse(text)for text in Mp]
        F0=ring.parse('m0*h**2+m1*h*x+m2*h*y+m3*x**2+m4*x*y+m5*y**2')
        F=ring.add(F0,ring.mul(spatial[1],ring.add(spatial[4],ring.derivative(M[2],3),ring.derivative(M[3],4))),ring.scale(ring.mul(spatial[2],ring.add(spatial[3],ring.derivative(M[0],3),ring.derivative(M[1],4))),-1))
        target=matrix_entries(ring,{(1,2):F,(1,3):M[0],(1,4):M[1],(2,3):M[2],(2,4):M[3]})
        for i,j in PAIRS:assert target[i][j]==parse(ring,sources[index]['full_six_magnetic_source'].get(f'{i},{j}','0'))
        strict_and_projective(ring,target)
        if index==6:
            # Exact PGL permutation followed by a polynomial Euler gauge.
            perm=[3,1,2,0,4]
            old_coordinate_values=[spatial[3],spatial[1],spatial[2],spatial[0],spatial[4]]
            raw=[[substitute_sparse(ring,target[perm[i]][perm[j]],ring,[ring.variable(a)for a in parameters],old_coordinate_values)for j in range(5)]for i in range(5)]
            gauge=[{},ring.parse('h/3'),ring.parse('y/3'),{},{}]
            assert item['changed_infinity']['permutation']==perm
            assert gauge==[ring.parse(v)for v in item['changed_infinity']['Euler_gauge_vector']]
            target=[[ring.add(raw[i][j],ring.scale(ring.mul(spatial[i],gauge[j]),-1),ring.mul(gauge[i],spatial[j]))for j in range(5)]for i in range(5)]
            assert all(not value for value in target[0])
            strict_and_projective(ring,target)
        if half or double_d:
            source_ring=LaurentRing(tuple(f'n{i}'for i in range(15)),('H','X','Y','Z','W'))
            source_parameters={f'n{i}':'0'for i in range(15)}
            if double_d:
                source_parameters.update(n0='-m0/e-u/e**2',n1='-m3/e',n2='-m4/e',n3='-m5/e-(5*a/2+1)/e**2',n4='-m1/e',n5='-u/e',n7='-2/e',n8='-m2/e-(3+(a+1/2)*u)/e**2',n10='a',n12='-(2*a+1)*u/2',n14='-1-2*a')
                for relation in ['(2*n10+n14+1)*n7+2*n11*n6','n11*n7','n11*n5+n7*n9','(6*n10+2)*n6+3*n13*n7','(2*n10+1)*n5+n12*n7+2*n6*n9','n5*n9']:
                    assert not substitute_sparse(source_ring,source_ring.parse(relation),ring,[ring.parse(source_parameters[name])for name in source_ring.names[:15]],[{}]*5)
            else:
                source_parameters.update(n0='-(e*m0+1)/e**2',n1='-m3/e',n2='(c-e*m4)/e**2',n3='(6*d-4*e*m5+1)/(4*e**2)',n4='-m1/e',n5='-1/e',n8='(n-e*m2)/e**2',n10='-1/2',n12='n',n13='c',n14='d')
            values=[ring.parse(source_parameters[name])for name in source_ring.names[:15]]
            old_coordinates=[ring.variable('h'),ring.parse('q2/e'),ring.parse('q1+y/e'),ring.variable('x'),ring.variable('y')]
            source_upper={}
            for ij,text in zip(it.combinations(range(1,5),2),old['full_six_brackets']):
                for a,b in [('x1','X'),('x2','Y'),('x3','Z'),('x4','W')]:text=re.sub(r'\b'+a+r'\b',b,text)
                source_upper[ij]=substitute_sparse(source_ring,source_ring.parse(text),ring,values,old_coordinates,True)
            Btexts=[['0','1','0','-1/e'],['e','0','0','0'],['0','0','1','0'],['0','0','0','1']]
            expected_mod=['x','y+e*q1','0','0'];expected_pf='-h*x**2*(y+u*h/2)'if double_d else'-h**2*x**2/2'
        else:
            source_ring=LaurentRing(tuple(f'A{i}'for i in range(6))+('u','v','t','q'),('H','a','b','C','D'))
            source_parameters=dict(A0='-m3/e',A1='-(e*m1+3)/e**2',A2='(c-e*m4)/e**2',A3='-m0/e',A4='(n-e*m2)/e**2',A5='(9*d-6*e*m5+2)/(6*e**2)',u='n',v='c',t='d',q='-3/e')
            values=[ring.parse(source_parameters[name])for name in source_ring.names[:10]]
            old_coordinates=[ring.variable('h'),ring.variable('x'),ring.parse('q2/e'),ring.parse('q1+y/e'),ring.variable('y')]
            A='A0*H**2+A1*a*H+A2*D*H+A3*a**2+A4*a*D+A5*D**2'
            S='(u*a+v*H+(t+1/3)*D)'
            formula={(1,2):'-a*C/3',(1,3):'-H*a/3',(2,3):A+'+t*(H*b-C**2/2)',(2,4):'q*H*a+C*'+S,(3,4):'H*'+S}
            source_upper={}
            for ij in PAIRS:
                text=branches[1]['whole_homogeneous_h_Casimir_brackets'].get(f'{ij[0]},{ij[1]}','0')
                for a,b in [('h','H'),('x1','a'),('x2','b'),('x3','C'),('x4','D')]:text=re.sub(r'\b'+a+r'\b',b,text)
                original=source_ring.parse(text)
                assert original==source_ring.parse(formula.get(ij,'0'))
                source_upper[ij]=substitute_sparse(source_ring,original,ring,values,old_coordinates)
            Btexts=[['0','0','1','-1/e'],['0','e','0','0'],['1','0','0','0'],['0','0','0','1']]
            expected_mod=['-2*h/3','-2*(y+e*q1)/3','0','0'];expected_pf='h**2*x**2'
        for name,text in item['whole_old_source_Laurent_parameters'if double_d else'source_Laurent_parameters'].items():assert ring.parse(text)==ring.parse(source_parameters[name])
        assert Btexts==item['actual_GL_matrix']
        B=[[ring.constant(i==j)for j in range(5)]for i in range(5)]
        for i,j in it.product(range(4),repeat=2):B[i+1][j+1]=ring.parse(Btexts[i][j])
        from classification_minimal_quartic_review import determinant
        assert determinant(ring,B)==ring.parse('-e')
        source=matrix_entries(ring,source_upper)
        family=[[ring.add(*(ring.mul(ring.mul(B[i][a],source[a][b]),B[j][b])for a,b in it.product(range(5),repeat=2)))for j in range(5)]for i in range(5)]
        assert all(powers[eposition]>=0 for row in family for value in row for powers in value)
        assert all(sum(powers[len(parameters):])==2 for row in family for value in row for powers in value)
        for i,j in PAIRS:
            assert family[i][j]==ring.parse(item['whole_regular_homogeneous_family'].get(f'{i},{j}','0'))
            endpoint={powers:coefficient for powers,coefficient in family[i][j].items()if powers[eposition]==0}
            assert endpoint==target[i][j]
        strict_and_projective(ring,family)
        pf=ring.add(ring.mul(family[1][2],family[3][4]),ring.scale(ring.mul(family[1][3],family[2][4]),-1),ring.mul(family[1][4],family[2][3]))
        assert pf==ring.parse(expected_pf)
        modular=[ring.add(*(ring.derivative(family[i][j],i)for i in range(1,5)))for j in range(1,5)]
        assert modular==[ring.parse(v)for v in expected_mod]
        records.append(dict(name='double_D_to_CN31_25'if double_d else item['name'],entire_original_parameters=list(parameters[:-1]),all_zero_boundaries_retained=True,
            sole_parameter_inverted='e',nonzero_source_period='n7=-2/e'if double_d else'n5=-1/e'if half else'q=-3/e',
            exact_GL_determinant='-e',all150_punctured_GL_source_coefficients_checked=True,
            all_Laurent_negative_e_powers_cancel=True,whole_regular_homogeneous_family_endpoint_checked=True,
            full_parameter_strict10_and_projective350_Jacobi_zero=True,universal_homogeneous_Pfaffian=expected_pf,
            universal_modular=expected_mod,changed_infinity_PGL_and_Euler_identity_checked=index==6,
            e_zero_complement='The exact whole original parameter endpoint; no parameterexception remains.',
            entire_original_source_covered_by_accepted_component=True,
            full_n7_open_six_compatibilities_checked=double_d,
            source_component='accepted CN31_25'if double_d else'accepted half25'if half else'accepted third25'))
    inputs=[Path(__file__),ROOT/'code/classification_minimal_quartic_review.py',ROOT/'code/classification_quadratic_cofactor_review.py',source_path,limit_path,nil31_path,old_path,
        ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_REVIEW_20261107.json',ROOT/'reports/CLASSIFICATION_NIL31_IMAGE_THIRD_COMPONENT_REVIEW_20261107.json']
    if double_d:inputs.append(ROOT/'reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_REVIEW_20261107.json')
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_CAS_or_NumPy_import=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},cases=records,
        new_component_count=0,proof='For each arbitrary original parameter choice, the punctured family is an actual GL transform of the accepted entire source with nonzero period. Every negative e power cancels in all150coefficients, giving a regular polynomial Poisson family with the exact original endpoint (up to the checked PGL/Euler change, if any). Closedness of that fixed accepted component gives endpoint containment. No smoothpoint uniqueness, divisor-only recognition, or finite jet inference is used.',
        reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({key:result[key]for key in ['all_passed','new_component_count','seconds','maximum_RSS_KiB']}))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--catalogue',type=Path)
    ap.add_argument('--limits',action='store_true',help='Review whole9parameter quadratic_B/triple_C regular limits to half/third25')
    ap.add_argument('--double-d',action='store_true',help='Review whole8parameter double_D regular limit to accepted CN31_25')
    args=ap.parse_args()
    if args.double_d:return limits_review(args.output,True)
    if args.limits:return limits_review(args.output)
    assert args.catalogue,'--catalogue is required for component review'
    args.catalogue=args.catalogue.resolve()
    begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
    component_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_COMPONENTS_20261107.json'
    sources=json.loads(source_path.read_text())['cases'];components=json.loads(component_path.read_text())['cases']
    records=[]
    for index,component in zip([4,5],components):
        ring,entries,point,directions,source=full_source(sources[index],index)
        ambient=true_cut(ring,point,directions,sources[index],component,index)
        records.append(dict(name=sources[index]['name'],full_source=source,true_full_ambient_support=ambient))
    inputs=[Path(__file__),Path(independent.__file__),ROOT/'code/classification_minimal_quartic_review.py',ROOT/'code/classification_quadratic_cofactor_review.py',ROOT/'code/classification_nil22_weighted_slice_review.py',source_path,component_path,args.catalogue]
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_or_producer_helper_import=True,no_CAS_or_NumPy_import=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},cases=records,
        intrinsic_pair_geometry=geometry(),all_prior94_distinctness=catalogue(json.loads(args.catalogue.read_text())),
        accepted_two_new_full_components_relative_to_prior94=True,projective_dimensions=[23,23],generic_ranks=[4,4],
        full_nine_remaining_sources_membership_or_exhaustion_claimed=False,
        support_proof='For each irreducible whole6source, actualF48 supplies the lower bound. The actual102dim linear complement and full nonzero D minor give a genuine formal cut. FullD-annihilating functionals have maximal-radical quadratic initial forms in one or two variables, making this actualcut zero dimensional. Krull gives the upper48, so the whole irreducible saturated source closure is a full component; Euler24 plus scalar1 leaves23. No smoothness or local uniqueness follows.',
        reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({key:result[key]for key in ['all_passed','projective_dimensions','seconds','maximum_RSS_KiB']}))


if __name__=='__main__':main()
