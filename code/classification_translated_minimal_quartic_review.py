#!/usr/bin/env python3
"""Fraction-only replay of the complete mixed slice and SS20 degeneration.

No producer, producer sparse engine, SymPy, NumPy or legacy arithmetic imports.
Uses the separately authored nil22 reviewer's elementary sparse operations.
"""
import argparse
import ast
from fractions import Fraction as F
import hashlib
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from pathlib import Path
import resource
import time

import classification_nil22_weighted_slice_review as arithmetic

arithmetic.N = 7
ROOT = Path(__file__).resolve().parents[1]
NAMES = ('h','x','y','z','w','u','v')
PAIRS = list(combinations(range(5),2))
QUADS = list(combinations(range(5),4))


def parse(expression):
    def walk(node):
        if isinstance(node,ast.Constant) and type(node.value) is int:
            return arithmetic.atom(node.value)
        if isinstance(node,ast.Name):
            return arithmetic.atom(1,5 if node.id=='t' else NAMES.index(node.id))
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
            return arithmetic.scale(walk(node.operand),-1)
        if isinstance(node,ast.BinOp):
            if isinstance(node.op,ast.Add):
                return arithmetic.add(walk(node.left),walk(node.right))
            if isinstance(node.op,ast.Sub):
                return arithmetic.add(walk(node.left),arithmetic.scale(walk(node.right),-1))
            if isinstance(node.op,ast.Mult):
                return arithmetic.mul(walk(node.left),walk(node.right))
            if isinstance(node.op,ast.Div):
                denominator = walk(node.right)
                assert set(denominator)=={(0,)*7}
                return arithmetic.scale(walk(node.left),1/denominator[(0,)*7])
            if isinstance(node.op,ast.Pow):
                assert isinstance(node.right,ast.Constant) and type(node.right.value) is int and node.right.value>=0
                result = arithmetic.atom(1)
                for _ in range(node.right.value):
                    result = arithmetic.mul(result,walk(node.left))
                return result
        raise ValueError(ast.dump(node))
    return walk(ast.parse(expression,mode='eval').body)


def decode(data):
    return arithmetic.antisymmetric({tuple(map(int,key.split(','))):parse(value) for key,value in data.items()})


def monomials(degree):
    result = []
    for term in combinations_with_replacement(range(5),degree):
        powers = [0]*5
        for i in term:
            powers[i] += 1
        result.append(tuple(powers))
    return sorted(result,reverse=True)


COLUMNS = [(ij,e) for ij in PAIRS for e in monomials(2)]
ROWS = [(quad,e) for quad in QUADS for e in monomials(4)]


def vector(point):
    return [point.get(ij,{}).get(e+(0,0),F(0)) for ij,e in COLUMNS]


def vector_polynomials(point):
    return [{powers[-2:]:value for powers,value in point.get(ij,{}).items() if powers[:5]==e}
            for ij,e in COLUMNS]


def jacobi_vector(point):
    J = arithmetic.jacobi(point)
    return [{powers[-2:]:value for powers,value in J[quad].items() if powers[:5]==e}
            for quad,e in ROWS]


def parameter_polynomial(expression):
    result = parse(expression)
    assert all(powers[:5]==(0,)*5 for powers in result)
    return {powers[-2:]:value for powers,value in result.items()}


def add_parameters(*values):
    result = {}
    for polynomial in values:
        for powers,value in polynomial.items():
            result[powers] = result.get(powers,F(0))+value
    return {e:value for e,value in result.items() if value}


def matrices(point):
    assert all(not value for value in arithmetic.jacobi(point).values())
    columns = []
    for ij,e in COLUMNS:
        changed = {key:dict(value) for key,value in point.items()}
        changed[ij] = arithmetic.add(changed.get(ij,{}),{e+(0,0):F(1)})
        changed[ij[::-1]] = arithmetic.scale(changed[ij],-1)
        J = jacobi_vector(changed)
        assert all(not value or set(value)=={(0,0)} for value in J)
        columns.append([value.get((0,0),F(0)) for value in J])
    D = [list(row) for row in zip(*columns)]
    actual = []
    for a,b in product(range(5),repeat=2):
        vertical = {}
        for i,j in PAIRS:
            vertical[i,j] = arithmetic.add(arithmetic.atom(1,i,b) if j==a else {},
                                             arithmetic.atom(-1,j,b) if i==a else {})
        actual.append(vector(vertical))
    for a,b in product(range(5),repeat=2):
        action = {}
        for i,j in PAIRS:
            action[i,j] = arithmetic.add(point.get((b,j),{}) if i==a else {},
                                         point.get((i,b),{}) if j==a else {},
                                         arithmetic.scale(arithmetic.mul(arithmetic.atom(1,b),
                                                          arithmetic.der(point.get((i,j),{}),a)),-1))
        actual.append(vector(action))
    actual.append(vector(point))
    G = [list(row) for row in zip(*actual)]
    sparse_rows = [[(j,value) for j,value in enumerate(row) if value] for row in D]
    assert all(sum(value*column[j] for j,value in row)==0 for row in sparse_rows for column in actual)
    return D,G,sparse_rows


def minor(matrix,rows,columns):
    return arithmetic.determinant([[matrix[i][j] for j in columns] for i in rows])


def verify_weighted_model(saved,point):
    D,G,sparse_D = matrices(point)
    weights = saved['coordinate_weights']
    character = saved['P_character']
    relative = [weights[i]+weights[j]-sum(a*b for a,b in zip(weights,e))-character for (i,j),e in COLUMNS]
    row_weights = [sum(weights[i] for i in quad)-sum(a*b for a,b in zip(weights,e))-2*character for quad,e in ROWS]
    assert relative==saved['all150_coordinate_relative_weights'] and row_weights==saved['all350_equation_relative_weights']
    assert (min(relative),max(relative))==(-4,3)
    assert all(not value or relative[i]==0 for i,value in enumerate(vector(point)))
    assert all(len({relative[i] for i,value in enumerate(column) if value})<=1 for column in zip(*G))
    assert all(row_weights[i]==relative[j] for i,row in enumerate(sparse_D) for j,value in row)
    group_rows,group_columns = saved['actual_group_cut_rows'],saved['actual_group_independent_columns']
    assert len(group_rows)==len(group_columns)==47
    assert arithmetic.rank(G)[0]==47
    group_det = minor(G,group_rows,group_columns)
    assert group_det==F(saved['actual_group_cut_minor_determinant']) and group_det
    sl = saved['weight_stable_slice_coordinate_indices']
    implicit,rows = saved['implicit_coordinate_indices'],saved['selected_equation_rows']
    free = saved['free_normal_coordinate_indices']
    assert sorted(group_rows+sl)==list(range(150)) and len(sl)==103
    assert len(implicit)==len(rows)==101 and sorted(implicit+free)==sorted(sl)
    det = minor(D,rows,implicit)
    assert det==F(saved['implicit_minor_determinant']) and det
    normal = [[F(encoding.get(str(i),'0')) for i in range(150)] for encoding in saved['normal_vectors']]
    assert len(normal)==2
    assert all(sum(value*column[j] for j,value in row)==0 for row in sparse_D for column in normal)
    assert all(column[i]==0 for column in normal for i in group_rows)
    assert [[column[i] for i in free] for column in normal]==[[F(1),F(0)],[F(0),F(1)]]
    assert [{relative[i] for i,value in enumerate(column) if value} for column in normal]==[{-1},{-2}]
    assert max([0]+[-relative[i] for i in implicit])==saved['finite_weight_bound_proved_before_solving']==4
    curve = decode(saved['full_polynomial_slice_brackets'])
    assert all(sum(powers[:5])==2 for ij in PAIRS for powers in curve.get(ij,{}))
    values = vector_polynomials(curve)
    base = vector(point)
    differences = [add_parameters(value,{(0,0):-b} if b else {}) for value,b in zip(values,base)]
    assert all(not differences[i] for i in group_rows)
    assert [differences[i] for i in free]==[{(1,0):F(1)},{(0,1):F(1)}]
    nonlinear = []
    for i,value in enumerate(differences):
        for powers,coefficient in value.items():
            assert powers[0]+2*powers[1]==-relative[i]
        correction = add_parameters(value,{(1,0):-normal[0][i]} if normal[0][i] else {},
                                      {(0,1):-normal[1][i]} if normal[1][i] else {})
        assert not correction or i in implicit
        assert all(sum(powers)>=2 and powers[0]+2*powers[1]<=4 for powers in correction)
        nonlinear.append(correction)
    recovered = [{} for _ in range(150)]
    assert [item['weighted_degree'] for item in saved['finite_corrections']]==[2,3,4]
    for item in saved['finite_corrections']:
        for key,value in item['coefficients'].items():
            polynomial = parameter_polynomial(value)
            assert int(key) in implicit
            assert all(a+2*b==item['weighted_degree'] for a,b in polynomial)
            recovered[int(key)] = add_parameters(recovered[int(key)],polynomial)
    assert recovered==nonlinear
    residual = jacobi_vector(curve)
    expected = {int(i):parameter_polynomial(value) for i,value in saved['all350_literal_residuals'].items()}
    assert {i:value for i,value in enumerate(residual) if value}==expected
    assert all(not residual[i] for i in rows)
    assert all(a+2*b==-row_weights[i] for i,value in enumerate(residual) for a,b in value)
    uv_rows = [i for i,value in enumerate(residual) if value and set(value)=={(1,1)}]
    v2_rows = [i for i,value in enumerate(residual) if value and set(value)=={(0,2)}]
    assert uv_rows and v2_rows
    assert all(b>=1 and (a>=1 or b>=2) for value in residual for a,b in value)
    axis_upper = {ij:{powers:coefficient for powers,coefficient in curve.get(ij,{}).items()
                          if powers[-1]==0} for ij in PAIRS}
    axis = arithmetic.antisymmetric({ij:value for ij,value in axis_upper.items() if value})
    assert axis==decode(saved['actual_v_zero_family_brackets'])
    assert all(not value for value in arithmetic.jacobi(axis).values())
    return axis,dict(full_D_exact_rank=101,actual_group_rank=47,whole_kernel_dimension=49,
                     all_D_times_group_and_normals_zero=True,true_C150_weight_stable_direct_sum=True,
                     true_group_cut_minor=str(group_det),implicit_minor=str(det),normal_relative_weights=[-1,-2],
                     a_priori_polynomial_bound=4,complete_model_and_all350_residuals_replayed=True,
                     complete_formal_slice_ideal=['u*v','v^2'],uv_generator_witness_rows=uv_rows,v2_generator_witness_rows=v2_rows,
                     every_residual_in_uv_v2_ideal=True,v_zero_full_actual_Poisson_curve=True,
                     nonzero_u_one_actual_orbit_by_verified_torus_equivariance=True)


def univariate_add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for power,value in polynomial.items():
            out[power] = out.get(power,F(0))+value
    return {power:value for power,value in out.items() if value}


def univariate_mul(*polynomials):
    out = {0:F(1)}
    for polynomial in polynomials:
        result = {}
        for p,a in out.items():
            for q,b in polynomial.items():
                result[p+q] = result.get(p+q,F(0))+a*b
        out = {power:value for power,value in result.items() if value}
    return out


def modular_zero_value(polynomial,x_constant=-16):
    result = {}
    for powers,value in polynomial.items():
        if powers[2] or powers[3] or powers[4] or powers[6]:
            continue
        result[powers[5]-2*powers[1]] = result.get(powers[5]-2*powers[1],F(0))+value*x_constant**powers[1]
    return {power:value for power,value in result.items() if value}


def ss20_membership(axis):
    div = [arithmetic.add(*(arithmetic.der(axis.get((i,j),{}),i) for i in range(5))) for j in range(5)]
    assert div==[{},parse('4*y'),parse('4*h+u**2*x/4'),parse('-u*z'),parse('u*w')]
    assert all(not modular_zero_value(value) for value in div)
    C = [[{} for _ in range(4)] for _ in range(4)]
    for i,j in combinations(range(4),2):
        C[i][j] = modular_zero_value(axis.get((i+1,j+1),{}))
        C[j][i] = {power:-value for power,value in C[i][j].items()}
    assert C[0][1]=={-2:F(-16)} and C[2][3]=={-2:F(-32)}
    assert all(not C[i][j] for i,j in [(0,2),(0,3),(1,2),(1,3)])
    pf = univariate_add(univariate_mul(C[0][1],C[2][3]),
                         {p:-v for p,v in univariate_mul(C[0][2],C[1][3]).items()},
                         univariate_mul(C[0][3],C[1][2]))
    assert pf=={-4:F(512)}
    A = [[{}, {0:F(4)}, {}, {}], [{2:F(1,4)}, {}, {}, {}],
         [{},{},{1:F(-1)},{}], [{},{},{},{1:F(1)}]]
    square = [[univariate_add(*(univariate_mul(A[i][k],A[k][j]) for k in range(4))) for j in range(4)] for i in range(4)]
    assert square==[[{2:F(1)} if i==j else {} for j in range(4)] for i in range(4)]
    T = [[{1:F(1,4)},{0:F(-1)},{},{}], [{},{},{0:F(1)},{}],
         [{1:F(1,4)},{0:F(1)},{},{}], [{},{},{},{0:F(1)}]]
    eigenvalues = [-1,-1,1,1]
    assert all(univariate_add(*(univariate_mul(T[i][k],A[k][j]) for k in range(4)))==
               univariate_mul({1:F(eigenvalues[i])},T[i][j]) for i,j in product(range(4),repeat=2))
    determinant = {}
    for permutation in permutations(range(4)):
        inversions = sum(permutation[i]>permutation[j] for i,j in combinations(range(4),2))
        determinant = univariate_add(determinant,univariate_mul({0:F((-1)**inversions)},
                                                               *(T[i][permutation[i]] for i in range(4))))
    assert determinant=={1:F(-1,2)}
    transformed = [[univariate_add(*(univariate_mul(T[i][a],C[a][b],T[j][b])
                                       for a,b in product(range(4),repeat=2))) for j in range(4)] for i in range(4)]
    assert transformed[0][2]=={-1:F(-8)} and transformed[1][3]=={-2:F(-32)}
    assert all(not transformed[i][j] for i,j in [(0,1),(0,3),(1,2),(2,3)])
    return dict(full_family_modular_field=['0','4y','4h+u²x/4','-uz','uw'],
                modular_zero='(x,y,z,w)=(-16/u²,0,0,0)',
                constant_P12='-16/u²',constant_P34='-32/u²',constant_pfaffian='512/u⁴',
                centered_A_squared_equals_u_squared_identity=True,
                exact_eigenspace_coordinates=['u(x+16/u²)/4-y','z','u(x+16/u²)/4+y','w'],
                eigenspace_frame_determinant='-u/2',negative_positive_crossing_matrix=['-8/u','-32/u²'],
                crossing_rank_two_for_every_nonzero_u=True,
                whole_family_to_SS20_rank_two_branch_verified=True)


def generic_ranks(axis,report):
    saved = report['results'][0]
    scale = F(saved['integer_scale'])
    generic = {}
    for ij,value in axis.items():
        out = {}
        for powers,coefficient in value.items():
            target = powers[:5]+(0,0)
            out[target] = out.get(target,F(0))+coefficient*scale
        generic[ij] = {powers:value for powers,value in out.items() if value}
    assert generic==decode(saved['full_homogeneous_witness'])
    D,G,_ = matrices(generic)
    assert arithmetic.rank(G)[0]==48
    certificates = saved['certificates']
    assert minor(D,certificates['D']['rows'],certificates['D']['columns'])==F(certificates['D']['determinant'])
    assert minor(G,certificates['group']['rows'],certificates['group']['columns'])==F(certificates['group']['determinant'])
    assert certificates['D']['rank_QQ']==102 and certificates['group']['rank_QQ']==48
    return dict(full350_D_rank=102,actual_orbit_vertical_scale_rank=48,
                projective_tangent_dimension=23,actual_projective_orbit_dimension=23,
                smooth_unique_component_dimension=23,known_MPS_membership_uses_SS20_theorem=True)


def review_cofactor_curves(output):
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started = time.monotonic()
    normal_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS_20261107.json'
    normal = json.loads(normal_path.read_text())
    records = []
    paths = [Path(__file__),Path(arithmetic.__file__),normal_path]
    for name,ratio,expected_source in [('cubic',F(1,3),7),('linear',F(1),20)]:
        path = ROOT/f'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_{name.upper()}_COFACTOR_CURVE_20261107.json'
        saved = json.loads(path.read_text())
        paths.append(path)
        assert saved['actual_curve']['kappa']==str(ratio)
        curve = decode(saved['actual_curve']['full_polynomial_brackets'])
        # Strict Jacobi, not only the Euler/projective image.
        strict = {triple:arithmetic.add(*(arithmetic.mul(curve.get((a,k),{}),
                        arithmetic.der(curve.get((b,c),{}),k))
                        for a,b,c in [triple,(triple[1],triple[2],triple[0]),
                                     (triple[2],triple[0],triple[1])] for k in range(5)))
                  for triple in combinations(range(5),3)}
        assert all(not value for value in strict.values())
        assert all(not value for value in arithmetic.jacobi(curve).values())
        assert all(sum(powers[:5])==2 for ij in PAIRS for powers in curve.get(ij,{}))
        target_model = next(model for model in normal['three_actual_Poisson_models'] if model['name']==name)
        target = {}
        for ij,value in zip(combinations(range(1,5),2),target_model['brackets']):
            polynomial = {}
            for powers,coefficient in parse(value).items():
                assert powers[0]==0 and powers[-2:]==(0,0)
                transformed = list(powers)
                transformed[0] = 2-sum(powers[1:5])
                assert transformed[0]>=0
                polynomial[tuple(transformed)] = coefficient
            if polynomial:
                target[ij] = polynomial
        target = arithmetic.antisymmetric(target)
        zero = arithmetic.antisymmetric({ij:{powers:coefficient for powers,coefficient in curve.get(ij,{}).items()
                                               if not powers[5] and not powers[6]} for ij in PAIRS
                                         if any(not powers[5] and not powers[6] for powers in curve.get(ij,{}))})
        assert zero==target
        divergence = [arithmetic.add(*(arithmetic.der(curve.get((i,j),{}),i) for i in range(5))) for j in range(5)]
        expected_divergence = [{},parse('y'),parse('h+t**2*x'),arithmetic.scale(parse('t*z'),ratio),
                              arithmetic.scale(parse('t*w'),-ratio)]
        assert divergence==expected_divergence
        weights = saved['coordinate_weights']
        character = saved['P_character']
        for ij in PAIRS:
            for powers,coefficient in curve.get(ij,{}).items():
                assert powers[6]==0
                assert weights[ij[0]]+weights[ij[1]]-sum(a*b for a,b in zip(weights,powers[:5]))-character==-character*powers[5]
        pf = arithmetic.add(arithmetic.mul(curve.get((1,2),{}),curve.get((3,4),{})),
                             arithmetic.scale(arithmetic.mul(curve.get((1,3),{}),curve.get((2,4),{})),-1),
                             arithmetic.mul(curve.get((1,4),{}),curve.get((2,3),{})))
        assert pf==parse(saved['actual_curve']['full_actual_pfaffian'])
        assert all(not modular_zero_value(value,x_constant=-1) for value in divergence)
        C = [[{} for _ in range(4)] for _ in range(4)]
        for i,j in combinations(range(4),2):
            C[i][j] = modular_zero_value(curve.get((i+1,j+1),{}),x_constant=-1)
            C[j][i] = {power:-value for power,value in C[i][j].items()}
        PfC = univariate_add(univariate_mul(C[0][1],C[2][3]),
                              {p:-v for p,v in univariate_mul(C[0][2],C[1][3]).items()},
                              univariate_mul(C[0][3],C[1][2]))
        assert PfC=={-4:F(1,8)}
        A = [[{}, {0:F(1)}, {}, {}], [{2:F(1)}, {}, {}, {}],
             [{},{},{1:ratio},{}], [{},{},{},{1:-ratio}]]
        assert all(not univariate_add(*(univariate_add(univariate_mul(A[i][k],C[k][j]),
                              univariate_mul(C[i][k],A[j][k])) for k in range(4)))
                   for i,j in product(range(4),repeat=2))
        frame = [[{1:F(1)},{0:F(-1)},{},{}],[{1:F(1)},{0:F(1)},{},{}],
                 [{},{},{0:F(1)},{}],[{},{},{},{0:F(1)}]]
        eig = [-F(1),F(1),ratio,-ratio]
        assert all(univariate_add(*(univariate_mul(frame[i][k],A[k][j]) for k in range(4)))==
                   univariate_mul({1:eig[i]},frame[i][j]) for i,j in product(range(4),repeat=2))
        frame_det = {}
        for permutation in permutations(range(4)):
            inversions = sum(permutation[i]>permutation[j] for i,j in combinations(range(4),2))
            frame_det = univariate_add(frame_det,univariate_mul({0:F((-1)**inversions)},
                                                         *(frame[i][permutation[i]] for i in range(4))))
        assert frame_det=={1:F(2)}
        assert saved['actual_curve']['known_semisimple_source_type']==expected_source
        records.append(dict(normal_form=name,kappa=str(ratio),all_t_strict_Jacobi_and_projective_Jacobi_zero=True,
                            all_homogeneous_quadratic_coefficients_checked=True,t_zero_entire_representative_identity=True,
                            full_modular_field_checked=True,whole_curve_torus_equivariance_checked=True,
                            exact_actual_pfaffian_checked=True,modular_zero='x=-1/t²,y=z=w=0',
                            constant_Pfaffian='1/(8t⁴)',constant_matrix_nondegenerate_for_every_nonzero_t=True,
                            exact_semisimple_frame=['t(x+1/t²)-y','t(x+1/t²)+y','z','w'],
                            frame_determinant='2t',semisimple_eigenvalues=list(map(str,eig)),
                            known_semisimple_whole_source_type=expected_source,whole_original_orbit_MPS_coverage=True))
    paths.extend([ROOT/'research/classification/HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT.md',
                  ROOT/'research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md'])
    report = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_producer_engine_import=True,
                  no_CAS_import=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  records=records,reviewer_source_bytes=Path(__file__).read_text(),
                  mathematical_review='For each fixed normal form, the actual all-t polynomial curve contains the entire original representative at0. The irreducible quartic at0 belongs to the open complement of the projective closed reducible-quartic locus. Exact torus equivariance makes all nonzerot members one actual PGL/scalar orbit, so all are rank4 with irreducible quartic. The explicit modularzero/eigenframe identifies the entire nonzerot orbit with an accepted covered semisimple source (type7 or20). The finite MPS closed union includes t=0; PGL invariance then covers the whole original normal-form orbit. No graph-negative inference or new component count is used.',
                  external_dependencies=['Complete type10 three-orbit reduction','Whole semisimple type7 and20 support coverage','Adopted MPS SNC component theorem'],
                  scope='Both remaining linear and cubic type10 branches covered by knownMPS. Together with the earlier mixed proof this closes the entire type10 simpleH/irreduciblequartic source. Other nilpotent types and global classification remain open.',
                  seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x') as handle:
        json.dump(report,handle,indent=2)
        handle.write('\n')
    print(json.dumps({key:report[key] for key in ('all_passed','seconds','maximum_RSS_KiB')}))


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started = time.monotonic()
    weighted_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_WEIGHTED_MIXED_20261107.json'
    ambient_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_AMBIENT_20261107.json'
    generic_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_GENERIC_MIXED_20261107.json'
    weighted = json.loads(weighted_path.read_text())
    point = next(value for value in json.loads(ambient_path.read_text())['results'] if value['name']=='mixed')
    axis,replay = verify_weighted_model(weighted,decode(point['full_homogeneous_witness']))
    ss20 = ss20_membership(axis)
    generic = generic_ranks(axis,json.loads(generic_path.read_text()))
    paths = [Path(__file__),Path(arithmetic.__file__),weighted_path,ambient_path,generic_path,
             ROOT/'research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md']
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_producer_engine_import=True,
                  no_CAS_import=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  full_weighted_model=replay,whole_nonzero_family_SS20_membership=ss20,
                  independent_generic_ambient_ranks=generic,
                  mathematical_review='The entire C150 action, true actual-group coordinate cut and nonzero101minor make the formal IFT model equivariant. All normal weights are negative and all coordinate weights bounded below by-4, so its unique implicit series is a polynomial of weighted degree<=4. The independently expanded full model solves every selected equation, and its full350residual ideal is exactly(uv,v²). This proves the actual complete formal slice and exhibits the whole polynomial u-axis family. For u!=0 the exact modular-zero/eigenframe/crossing-rank calculation applies the accepted wholeSS20 rank2 MPS theorem. Closedness of the finite MPS union and PGL invariance include u=0 and its entire orbit, covering the mixed type10 branch.',
                  external_dependency='WholeSS20 rank2 source coverage theorem in HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md, and its adopted MPS SNC component theorem; these are not reproved here.',
                  known_component_count_increment=0,
                  scope='Mixed type10 branch completely covered by known MPS; linear/cubic type10 branches remain separate. Exact nil index and global classification not claimed.',
                  reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,
                  maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({key:result[key] for key in ('all_passed','seconds','maximum_RSS_KiB')}))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cofactor-curves',action='store_true',help='Independent all-t review of the cubic and linear semisimple deformations')
    args=parser.parse_args()
    (review_cofactor_curves if args.cofactor_curves else main)(args.output)
