#!/usr/bin/env python3
"""Independent full second nil31-long source -> T40 review.

No producer, S61 engine, CAS, NumPy, or legacy arithmetic imports.
Only previously independent Fraction/Laurent arithmetic is reused.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

from classification_minimal_quartic_review import LaurentRing
from classification_quartic_nil31_long_review import substitute
from classification_quadratic_cofactor_review import verify_source
import classification_quartic_second_graph_review as independent

q=independent.q
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2))


def mproduct(R,A,B):
    return [[R.add(*(R.mul(a,b)for a,b in zip(row,col)))for col in zip(*B)]for row in A]


def family_review(source,support,normal):
    original=source['full_nilpotent_sources'][4]
    original_review=verify_source(original,5)
    assert len(original_review['free_indices'])==15
    R=LaurentRing([f'n{i}'for i in range(15)],['h','x','y','z','w'],invert=['n14'])
    p=R.parse;a=R.variable('n14')
    replacement={i:R.variable(f'n{i}')for i in range(15)}
    replacement.update({2:p('1/(216*n14**2)'),8:p('-n12/(6*n14)'),13:p('-1/6'),
        1:p('-27*n14*n4**2-n6/(12*n14)-n12/(36*n14**2)'),
        0:p('54*n14*n4**2*n12-18*n14*n3*n4+n6*n12/(2*n14)'),
        5:p('n12**2/n14-6*n6*n12'),9:p('1/(36*n14)'),
        11:p('-9*n14*n4'),7:p('-6*n14*n6+2*n12'),10:p('6*n14*(3*n4*n12-n3)')})
    packet=support['nil31_long_support']['all_two_complete_families'][1]
    import re
    parse_saved=lambda t:p(re.sub(r'\ba\b','n14',t))
    assert replacement=={int(k[1:]):parse_saved(v)for k,v in packet['complete_source_substitutions'].items()}
    aff=[substitute(R,p(t),replacement)for t in original['brackets']]
    assert aff==[parse_saved(t)for t in packet['universal_brackets']]
    def homogeneous(poly):
        out={}
        for e,c in poly.items():
            f=list(e);assert f[15]==0;f[15]=2-sum(f[16:]);assert f[15]>=0;out[tuple(f)]=c
        return out
    P={(i+1,j+1):homogeneous(v)for (i,j),v in zip(it.combinations(range(4),2),aff)}
    for i,j in PAIRS:P[j,i]=R.scale(P.get((i,j),{}),-1)
    canonical={}
    at={i:R.constant(1 if i==14 else 0)for i in [14,3,4,6,12]}
    for ij in PAIRS:canonical[ij]=substitute(R,P.get(ij,{}),at)
    for i,j in PAIRS:canonical[j,i]=R.scale(canonical[i,j],-1)
    expected=normal['nil31_long_normalform']['all_two_actual_orbit_normalforms'][1]
    assert all(canonical[ij]==p(expected['canonical_homogeneous_brackets'].get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS)
    d=p('-3*n6/2+n12/(2*n14)');s=p('18*n14*n6');f=p('-3*n14*n4');e=p('-2*n3-54*n14*n4*n6+6*n4*n12')
    z,o=R.constant(0),R.constant(1)
    M=[[o,z,z,z,z],[e,o,z,d,f],[d,z,o,z,z],[z,z,z,o,z],[s,z,z,z,a]]
    Mi=[[o,z,z,z,z],[R.add(R.mul(f,s,p('1/n14')),R.scale(e,-1)),o,z,R.scale(d,-1),R.scale(R.mul(f,p('1/n14')),-1)],
        [R.scale(d,-1),z,o,z,z],[z,z,z,o,z],[R.scale(R.mul(s,p('1/n14')),-1),z,z,z,p('1/n14')]]
    identity=[[R.constant(i==j)for j in range(5)]for i in range(5)]
    assert mproduct(R,M,Mi)==mproduct(R,Mi,M)==identity
    N=[[R.constant((i,j)in[(1,2),(2,3),(3,0)])for j in range(5)]for i in range(5)]
    assert mproduct(R,M,N)==mproduct(R,N,M)
    assert M==[[parse_saved(t)for t in row]for row in expected['canonical_to_universal_affine_matrix']]
    assert Mi==[[parse_saved(t)for t in row]for row in expected['universal_to_canonical_affine_matrix']]
    coordinates=[R.variable(t)for t in ['h','x','y','z','w']]
    newcoords=[R.add(*(R.mul(c,v)for c,v in zip(row,coordinates)))for row in M]
    for i,j in PAIRS:
        lhs=substitute(R,P.get((i,j),{}),{15+k:v for k,v in enumerate(newcoords)})
        rhs=R.add(*(R.mul(M[i][a],M[j][b],canonical.get((a,b),{}))for a,b in it.product(range(5),repeat=2)))
        assert lhs==rhs
    strict=[R.add(*(R.mul(P.get((a,l),{}),R.derivative(P.get((b,c),{}),l))for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in it.combinations(range(5),3)]
    assert all(not v for v in strict)
    point=q.antisymmetric({ij:{e[-5:]+(0,0,0):c for e,c in canonical[ij].items()}for ij in PAIRS})
    assert all(not value for value in q.jacobi(point).values())
    return point,dict(original_full90_source_independently_reconstructed=original_review,
        entire15_parameter_substitution_reconstructed=True,parameter_domain='Gm(n14) x A4(n3,n4,n6,n12)',
        universal_full_six_brackets_exactly_match=True,all_ten_strict_Jacobi_zero=True,
        whole_source_actual_affine_equivalence=True,tensor_scalar='1',matrix_determinant='n14',
        exact_canonical_to_source_matrix=[[R.encode(v)for v in row]for row in M],exact_inverse=[[R.encode(v)for v in row]for row in Mi],
        inverse_both_directions=True,commutes_with_full_homogeneous_cofactor=True,
        excluded_locus='n14=0 only, already separately excluded for the whole H+irreducible quartic support')


def certify(matrix):
    rank,columns=q.rank(matrix);rows=q.rank([[row[j]for row in matrix]for j in columns])[1]
    minor=q.determinant([[matrix[i][j]for j in columns]for i in rows]);assert minor
    return dict(shape=[len(matrix),len(matrix[0])],rank=rank,rows=rows,columns=columns,determinant=str(minor),sha256=hashlib.sha256(json.dumps(matrix,default=str,separators=(',',':')).encode()).hexdigest())


def ambient_review(point,saved):
    scale=Q(saved['integer_scale']);assert scale==216
    P=q.antisymmetric({ij:q.scale(point.get(ij,{}),scale)for ij in PAIRS})
    assert all(P.get(ij,{})==independent.parse_saved(saved['witness_brackets'].get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS)
    columns=[]
    for ij,e in q.COLS:
        variation=dict(P);variation[ij]=q.add(P.get(ij,{}),{e+(0,0,0):Q(1)});variation[ij[::-1]]=q.scale(variation[ij],-1)
        columns.append(q.jvec(q.jacobi(variation)))
    D=q.transpose(columns);Dcert=certify(D);assert Dcert['rank']==100
    old=saved['derivative'];minor=q.determinant([[D[i][j]for j in old['columns']]for i in old['rows']]);assert minor==Q(old['determinant'])and minor
    Kcols=[[Q(v.get(str(i),'0'))for i in range(150)]for v in saved['full_exact_kernel_sparse_columns']]
    assert len(Kcols)==50 and all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in Kcols)
    Kcert=certify(q.transpose(Kcols));assert Kcert['rank']==50
    E=[];G=[]
    for a,b in it.product(range(5),repeat=2):
        vertical={};action={}
        for i,j in PAIRS:
            vertical[i,j]=q.add(q.atom(1,i,b)if j==a else {},q.atom(-1,j,b)if i==a else {})
            action[i,j]=q.add(P.get((b,j),{})if i==a else {},P.get((i,b),{})if j==a else {},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
        E.append(q.vec(vertical));G.append(q.vec(action))
    actual=E+G+[q.vec(P)];assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in actual)
    Ecert=certify(q.transpose(E));Gcert=certify(q.transpose(actual));assert Ecert['rank']==24 and Gcert['rank']==47
    return dict(derivative=Dcert,saved_full100_minor_independently_reconstructed=True,whole_kernel=Kcert,
        all150_times50_kernel_products_zero=True,Euler=Ecert,actual_GL5_Euler_scalar=Gcert,
        projective_tangent_dimension=25,projective_actual_orbit_dimension=22,
        orbit_alone_does_not_prove_componenthood=True)


def gauge(A,B):
    delta=q.antisymmetric({ij:q.add(A.get(ij,{}),q.scale(B.get(ij,{}),-1))for ij in PAIRS})
    independent.same_projective(A,B)
    divergence=[q.add(*(q.der(delta.get((i,j),{}),i)for i in range(5)))for j in range(5)]
    for i,j in PAIRS:
        assert delta.get((i,j),{})==q.scale(q.add(q.mul(q.atom(1,i),divergence[j]),q.scale(q.mul(q.atom(1,j),divergence[i]),-1)),Q(1,5))
    return dict(all150_coefficients_exact_Euler_gauge=True,linear_gauge_field=[{','.join(map(str,e)):str(c/5)for e,c in p.items()}for p in divergence])


def graph_review(point,saved,cert,model,reduced,T39,T40):
    translation=[[Q(v)for v in row]for row in saved['translation']]
    assert translation==[[Q(i==j)-Q(i==2 and j==0)for j in range(5)]for i in range(5)]
    translated=independent.transform(point,translation)
    K=[[translated.get((i+1,j+1),{}).get((2,0,0,0,0,0,0,0),Q(0))for j in range(4)]for i in range(4)]
    assert K==[[Q(v)for v in row]for row in saved['nondegenerate_constant_matrix']]
    pf=K[0][1]*K[2][3]-K[0][2]*K[1][3]+K[0][3]*K[1][2];assert pf==Q(4,3)
    B=[[Q(v)for v in row]for row in saved['Darboux_frame']];B4=[row[1:]for row in B[1:]]
    J=[[Q((i,j)in[(0,1),(2,3)])-Q((i,j)in[(1,0),(3,2)])for j in range(4)]for i in range(4)]
    assert q.matmul(q.matmul(B4,K),q.transpose(B4))==J
    A=[[Q(i==j)for j in range(5)]for i in range(5)];A[0][1:]=list(map(Q,saved['perspective']))
    M=q.matmul(A,q.matmul(B,translation));Mi=q.inverse(M);identity=[[Q(i==j)for j in range(5)]for i in range(5)]
    assert q.matmul(M,Mi)==q.matmul(Mi,M)==identity
    pushed=independent.transform(point,M);params=list(map(Q,saved['all61_parameters']));assert len(params)==61
    def evaluate(data):
        out=Q(0)
        for mon,num,den in data:
            v=Q(num,den)
            for k in mon:v*=params[k]
            out+=v
        return out
    vector=q.vec({(1,2):q.atom(1,0,0),(3,4):q.atom(1,0,0)})
    for at,v in zip(cert['raw_free_indices'],params):vector[at]=v
    for entry in model['eliminated_raw_variables']:vector[entry['raw_index']]=evaluate(entry['polynomial'])
    represented=independent.vector_to_poly(vector)
    forward=gauge(pushed,represented);backward=gauge(independent.transform(represented,Mi),point)
    assert len(reduced['equations'])==134 and all(evaluate(v)==0 for v in reduced['equations'])
    values={f'u{i+1}':v for i,v in enumerate(params)};residuals={}
    for label,section in [('T39',T39),('T40',T40)]:
        residual=[params[k]-q.evaluate(v,values)for k,v in zip(section['weight2_u_indices_zero_based'],section['quadratic_section_q'])]
        assert len(residual)==45 and list(map(str,residual))==saved['graph_membership'][label]['all45_residuals']
        residuals[label]=list(map(str,residual))
    assert all(Q(v)==0 for v in residuals['T40'])and any(Q(v)!=0 for v in residuals['T39'])
    return dict(actual_evaluation=['0','1','0','0'],actual_affine_Pfaffian='4/3',Darboux_value_exactly_J=True,
        translation= saved['translation'],full_linear_PGL_matrix=[[str(v)for v in row]for row in M],full_inverse=[[str(v)for v in row]for row in Mi],
        finite_perspective_invertible=True,complete_150_forward=forward,complete_150_inverse=backward,
        all134_original_equations_zero=True,all61_parameters=list(map(str,params)),all45_graph_residuals=residuals,
        whole_T40_graph_membership=True,whole_T39_graph_membership=False,
        external_input='The accepted universal Q52 graph equals the reduced S61 closure of the known25-dimensional T40 component; this theorem is used, not reproved by this point check.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_SECOND_REVIEW_20261107.json');args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(179,180));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);start=time.monotonic()
    names=['reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json','reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_SUPPORT_20261107.json','reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_NORMALFORM_20261107.json','reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_ORBIT2_AMBIENT_20261107.json','reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_ORBIT2_GRAPH_20261107.json','legacy/v02/data/symplectic_slice_certificate.json','legacy/v02/data/symplectic_slice_equations.json','legacy/v02/data/symplectic_slice_reduced_equations.json','reports/Q51_UNIVERSAL_SECTION_20261101.json','reports/Q52_UNIVERSAL_SECTION_20261101.json']
    source,support,normal,ambient,graph,cert,model,reduced,T39,T40=[json.loads((ROOT/p).read_text())for p in names]
    for packet in [support,normal,ambient,graph]:assert hashlib.sha256(packet['source_bytes'].encode()).hexdigest()==packet['source_sha256']
    point,family=family_review(source,support,normal)
    ambient_result=ambient_review(point,ambient['source_point'])
    recognized=graph_review(point,graph['source_graph'],cert,model,reduced,T39,T40)
    imports=['code/classification_minimal_quartic_review.py','code/classification_quartic_nil31_long_review.py','code/classification_quadratic_cofactor_review.py','code/classification_quartic_second_graph_review.py','code/classification_degenerate_quadric_q4_nil31_I_review.py',str(Path(__file__).relative_to(ROOT))]
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_or_CAS_import=True,whole_source=family,full_ambient=ambient_result,faithful_full_graph=recognized,
        geometric_implication=dict(known_T40_dimension=25,ambient_tangent_dimension=25,canonical_point_is_regular=True,unique_ambient_component_at_point='T40',entire_second_source_and_its_PGL_saturation_in_T40=True,
            proof='The faithful full graph places the canonical point in the known25-dimensional T40 component. Full tangent25 forces local regularity and uniqueness. Every point of the complete a!=0 source is an explicitly verified invertible affine image of this point; PGL preserves T40. No new component is counted.',additional_component_count=0,first_source_not_reviewed=True),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in names+imports},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest=ROOT/args.output
    if dest.exists():raise FileExistsError(str(dest))
    dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':main()
