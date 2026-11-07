#!/usr/bin/env python3
"""Independent fullambient/S61 recognition of the second regular-nil source.

No producer, S61 arithmetic engine or CAS import.  The two imports are older
independent standard-library sparse/Fraction arithmetic implementations.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import itertools as it
import json
import math
from pathlib import Path
import re
import resource
import time

import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_quadratic_cofactor_review import Ring, PAIRS as AFFINE_PAIRS, MONS


def simultaneous(R,poly,values):
    out = {}
    for e,c in poly.items():
        term = {((0,)*R.n)+e[-4:]:c}
        for i,power in enumerate(e[:-4]):
            for _ in range(power):
                term = R.mul(term,values[i])
        out = R.add(out,term)
    return out


def convert_family(R,entries):
    """Homogenize full family; n10=a,n0=b,n4=c become independent tau0:3."""
    out = {}
    for (i,j),p in zip(AFFINE_PAIRS,entries):
        poly = {}
        for e,c in p.items():
            assert all(not value for k,value in enumerate(e[:-4]) if k not in (0,4,10))
            exponent = (2-sum(e[-4:]),)+e[-4:]+(e[10],e[0],e[4])
            assert exponent[0]>=0
            poly[exponent] = poly.get(exponent,F(0))+c
        out[i+1,j+1] = poly
    return q.antisymmetric(out)


def family_checks(source):
    row = source['full_nilpotent_sources'][0]
    R = Ring(13)
    raw = [R.parse(v) for v in row['brackets']]
    n = [R.variable(f'n{i}') for i in range(13)]
    normalized_values = [R.constant(0) for _ in range(13)]
    for i in (0,4,10):
        normalized_values[i] = n[i]
    normalized_values[1] = R.scale(R.mul(n[10],n[4]),-6)
    normalized_values[6] = R.scale(R.mul(n[10],n[10]),3)
    normalized_values[8] = R.constant(F(-1,4))
    normalized_values[12] = R.constant(F(-1,2))
    family = [simultaneous(R,p,normalized_values) for p in raw]
    full = convert_family(R,family)
    assert all(not j for j in q.jacobi(full).values())
    # Also check all ten raw homogeneous Jacobi, not just their E-wedge.
    rawJ = [q.add(*(q.mul(full.get((a,l),{}),q.der(full.get((b,c),{}),l))
                  for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(5)))
            for i,j,k in it.combinations(range(5),3)]
    assert all(not j for j in rawJ)
    # Specialize a=1,b=c=0 in the literal entire family.
    point = {}
    for ij,p in full.items():
        point[ij] = {}
        for ex,c in p.items():
            if ex[6]==ex[7]==0:
                point[ij][ex[:5]+(0,0,0)] = point[ij].get(ex[:5]+(0,0,0),F(0))+c

    # Complete N^2 normalization on the second scalar branch, before imposing
    # its four Jacobi equations. No a-dependent division is made.
    initial = list(n)
    for idx,value in {8:F(-1,4),12:F(-1,2),11:0,7:0}.items():
        initial[idx] = R.constant(value)
    initial[9] = R.scale(n[5],-3)
    initial[2] = R.scale(R.mul(n[5],n[5]),-12)
    before = [simultaneous(R,p,initial) for p in raw]
    k = R.scale(n[5],6)
    transformed = [R.add(before[0],R.mul(k,before[2]),R.scale(R.mul(k,before[3]),-1),R.mul(k,k,before[5])),
                   before[1],R.add(before[2],R.mul(k,before[5])),
                   R.add(before[3],R.scale(R.mul(k,before[5]),-1)),before[4],before[5]]
    # Old x=x'-k*z', old y=y'-k*w', with N^4=0.
    def coord_substitute(poly,pos,value):
        out = {}
        for ex,c in poly.items():
            e = list(ex)
            power = e[pos]
            e[pos] = 0
            term = {tuple(e):c}
            for _ in range(power):
                term = R.mul(term,value)
            out = R.add(out,term)
        return out
    transformed = [coord_substitute(coord_substitute(p,13,R.add(R.variable('x'),R.scale(R.mul(k,R.variable('z')),-1))),
                                    14,R.add(R.variable('y'),R.scale(R.mul(k,R.variable('w')),-1))) for p in transformed]
    coeff = []
    for p in transformed:
        for mon in MONS:
            coeff.append({ex[:-4]+(0,)*4:c for ex,c in p.items() if ex[-4:]==mon})
    new = [coeff[j] for j in row['free_coefficient_indices']]
    expected = list(n)
    expected[0] = R.add(n[0],R.scale(R.mul(n[5],n[6]),-12))
    expected[3] = R.add(n[3],R.scale(R.mul(n[10],n[5]),-12))
    for idx in (2,5,7,9,11):
        expected[idx] = {}
    expected[8],expected[12] = R.constant(F(-1,4)),R.constant(F(-1,2))
    assert new==expected
    assert transformed==[simultaneous(R,p,new) for p in raw]
    return point,dict(full_Gm_a_times_A2_source_raw_Jacobi_zero=True,
                      parameter_domain='Gm(a) x A2(n0,n4), a!=0',
                      normalized_parameters=[R.encode(v) for v in normalized_values],
                      canonical_point='a=1,n0=n4=0',
                      N_squared_matrix='G=I+6*n5*N^2, G^-1=I-6*n5*N^2',
                      all13_transformed_parameters=[R.encode(v) for v in new],
                      full90_coefficient_centralizer_reconstruction=True,
                      centralizer_denominator='4*n12+1=-1, so no excluded locus on this source')


def vector_to_poly(v):
    out = {}
    for value,(ij,e) in zip(v,q.COLS):
        if value:
            out.setdefault(ij,{})[e+(0,0,0)] = value
    return q.antisymmetric(out)


def parse_saved(text):
    return q.parse(re.sub(r'\b([xyzw])\b',lambda m:dict(x='x1',y='x2',z='x3',w='x4')[m[1]],text))


def ambient_check(point,saved):
    scale = F(saved['integer_scale'])
    assert scale==12
    P = {ij:q.scale(p,scale) for ij,p in point.items()}
    for ij in q.PAIRS:
        assert P.get(ij,{})==parse_saved(saved['witness_brackets'].get(f'{ij[0]},{ij[1]}','0'))
    Dcols = []
    for ij,e in q.COLS:
        variation = dict(P)
        variation[ij] = q.add(P.get(ij,{}),{e+(0,0,0):F(1)})
        variation[ij[::-1]] = q.scale(variation[ij],-1)
        Dcols.append(q.jvec(q.jacobi(variation)))
    D = q.transpose(Dcols)
    vertical,orbit = [],[]
    for a,b in it.product(range(5),repeat=2):
        V,O = {},{}
        for i,j in q.PAIRS:
            V[i,j] = q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
            O[i,j] = q.scale(q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},
                                   q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1)),-1)
        vertical.append(q.vec(V))
        orbit.append(q.vec(O))
    Mcols = vertical+orbit+[q.vec(P)]
    normals = [[F(v.get(str(i),'0')) for i in range(150)] for v in saved['independent_normal_sparse_vectors']]
    Kcols = Mcols+normals
    assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in Kcols)
    Dcert = saved['derivative']
    det = q.determinant([[D[i][j] for j in Dcert['columns']] for i in Dcert['rows']])
    assert det==F(Dcert['determinant']) and det
    K = q.transpose(Kcols)
    Kcert = saved['whole_kernel']
    kernel_det = q.determinant([[K[i][j] for j in Kcert['columns']] for i in Kcert['rows']])
    assert kernel_det and q.rank(K)[0]==50
    assert q.rank(q.transpose(vertical))[0]==24
    assert q.rank(q.transpose(Mcols))[0]==47
    return dict(full_derivative_shape=[350,150],exact_derivative_rank=100,entire_kernel_dimension=50,
                independent_100_minor=dict(rows=Dcert['rows'],columns=Dcert['columns'],determinant=str(det)),
                independent_50_kernel_minor=dict(rows=Kcert['rows'],columns=Kcert['columns'],determinant=str(kernel_det)),
                all150_times50_kernel_products_zero=True,vertical_rank=24,orbit_vertical_scale_rank=47,
                projective_tangent_dimension=25)


def transform(P,M):
    Mi = q.inverse(M)
    lin = [q.add(*(q.atom(v,j) for j,v in enumerate(row))) for row in Mi]+[q.atom(1,i) for i in range(5,8)]
    out = {ij:q.substitute(q.add(*(q.scale(P.get((a,b),{}),M[ij[0]][a]*M[ij[1]][b])
                                  for a,b in it.product(range(5),repeat=2))),lin) for ij in q.PAIRS}
    return q.antisymmetric(out)


def same_projective(A,B):
    delta = q.antisymmetric({ij:q.add(A.get(ij,{}),q.scale(B.get(ij,{}),-1)) for ij in q.PAIRS})
    for i,j,k in it.combinations(range(5),3):
        assert not q.add(q.mul(q.atom(1,i),delta[j,k]),q.scale(q.mul(q.atom(1,j),delta[i,k]),-1),q.mul(q.atom(1,k),delta[i,j]))


def graph_check(point,graph,cert,model,reduced,T39,T40):
    K = [[point.get((i+1,j+1),{}).get((2,0,0,0,0,0,0,0),F(0)) for j in range(4)] for i in range(4)]
    assert K==[[F(v) for v in row] for row in graph['nondegenerate_constant_matrix']]
    B = [[F(v) for v in row] for row in graph['Darboux_frame']]
    B4 = [row[1:] for row in B[1:]]
    J = [[F((i,j) in [(0,1),(2,3)])-F((i,j) in [(1,0),(3,2)]) for j in range(4)] for i in range(4)]
    assert q.matmul(q.matmul(B4,K),q.transpose(B4))==J
    A = [[F(i==j) for j in range(5)] for i in range(5)]
    A[0][1:] = list(map(F,graph['perspective']))
    M = q.matmul(A,B)
    Mi = q.inverse(M)
    assert q.matmul(M,Mi)==[[F(i==j) for j in range(5)] for i in range(5)]
    pushed = transform(point,M)
    params = list(map(F,graph['all61_parameters']))
    assert len(params)==61
    def evaluate(data):
        out = F(0)
        for mon,num,den in data:
            v = F(num,den)
            for k in mon:
                v *= params[k]
            out += v
        return out
    vector = q.vec({(1,2):q.atom(1,0,0),(3,4):q.atom(1,0,0)})
    for at,v in zip(cert['raw_free_indices'],params):
        vector[at] = v
    for entry in model['eliminated_raw_variables']:
        vector[entry['raw_index']] = evaluate(entry['polynomial'])
    represented = vector_to_poly(vector)
    same_projective(pushed,represented)
    same_projective(transform(represented,Mi),point)
    assert len(reduced['equations'])==134 and all(evaluate(e)==0 for e in reduced['equations'])
    values = {f'u{i+1}':v for i,v in enumerate(params)}
    residuals = {}
    for label,section in [('T39',T39),('T40',T40)]:
        residual = [params[k]-q.evaluate(v,values) for k,v in
                    zip(section['weight2_u_indices_zero_based'],section['quadratic_section_q'])]
        assert len(residual)==45
        assert list(map(str,residual))==graph['graph_membership'][label]['all45_residuals']
        residuals[label] = list(map(str,residual))
    assert all(F(v)==0 for v in residuals['T39'])
    assert any(F(v)!=0 for v in residuals['T40'])
    return dict(Darboux_value_exactly_J=True,finite_perspective_invertible=True,projective_scalar=1,
                complete150_coefficient_quotient_forward_and_inverse=True,
                all_ten_E_wedge_relations_each_direction=True,
                all134_original_reduced_equations_zero=True,all61_parameters=list(map(str,params)),
                all45_graph_residuals=residuals,T39_whole_graph_membership=True,T40_whole_graph_membership=False,
                exact_matrix=[[str(v) for v in row] for row in M],exact_inverse=[[str(v) for v in row] for row in Mi])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    started = time.monotonic()
    paths = ['reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json',
             'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_SECOND_INITIAL_20261107.json',
             'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_GRAPH_20261107.json',
             'legacy/v02/data/symplectic_slice_certificate.json','legacy/v02/data/symplectic_slice_equations.json',
             'legacy/v02/data/symplectic_slice_reduced_equations.json',
             'reports/Q51_UNIVERSAL_SECTION_20261101.json','reports/Q52_UNIVERSAL_SECTION_20261101.json']
    source,second,graph,cert,model,reduced,T39,T40 = [json.loads(Path(p).read_text()) for p in paths]
    point,family = family_checks(source)
    ambient = ambient_check(point,second['nil_regular_second_component'])
    recognition = graph_check(point,graph['nil_regular_graph'],cert,model,reduced,T39,T40)
    catalogue_path = Path('reports/CLASSIFICATION_BASELINE89_COMPONENT_CATALOGUE_20261107.json.gz')
    catalogue = json.loads(gzip.decompress(catalogue_path.read_bytes()))
    T39row = next(r for r in catalogue['components'] if r['id']=='T39')
    assert T39row['projective_dimension']==25 and T39row['generic_rank']==4
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
                  complete_source=family,fullambient=ambient,whole_graph=recognition,
                  existing_T39_dimension=25,canonical_point_local_ring_regular=True,
                  unique_component_at_point='T39',entire_normalized_source_contained_in_T39=True,
                  after_inverse_N_squared_and_connected_PGL_entire_source_contained_in_T39=True,
                  geometric_bridge='T39 through the point gives local dimension>=25; independently verified fullambient tangent dimension<=25 gives regularity and uniqueness. The irreducible full Gm(a)xA2 source contains this point, so its image closure is in T39. Connected group actions preserve every component.',
                  universal_Q51_graph_equals_T39_accepted_as_existing_dependency=True,
                  global_exhaustion_claimed=False,
                  input_sha256={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in
                                paths+[str(catalogue_path),str(Path(q.__file__).relative_to(Path.cwd())),'code/classification_quadratic_cofactor_review.py']},
                  reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  reviewer_source_bytes=Path(__file__).read_text(),
                  seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','maximum_RSS_KiB','unique_component_at_point']}))


if __name__=='__main__':
    main()
