#!/usr/bin/env python3
"""Independent Fraction addendum for the two-normal nodal plane-cubic source.
Only arithmetic from our independent five-source checker is reused; no producer
engine, SymPy, NumPy, or claimed matrix is imported.
"""
import hashlib,json,resource,time
from pathlib import Path
from itertools import product
import classification_two_lines_plane_cubic_component_review as q
ROOT=Path(__file__).resolve().parents[1]
def main():
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();paths=[ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json',ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_TWO_NORMALS_20261107.json',ROOT/'code/classification_two_lines_plane_cubic_component_review.py']
    source=next(r for r in json.loads(paths[0].read_text())['cases'] if r['case']=='node_tangent_flex')
    saved=next(r for r in json.loads(paths[1].read_text())['cases'] if r['case']=='node_tangent_flex')
    P,directions,source_check=q.full_source_check(source);R=[q.decode(v) for v in saved['normals']];assert len(R)==2
    Dcols=[]
    for ij,e in q.COLS:
        v={key:dict(value) for key,value in P.items()};inc={e+(0,):q.F(1)}
        v[ij]=q.add(v.get(ij,{}),inc);v[ij[::-1]]=q.scale(v[ij],-1);Dcols.append(q.jvec(q.jacobi(v)))
    D=q.transpose(Dcols);VV=[];OO=[]
    for a,b in product(range(5),repeat=2):
        V={};O={}
        for i,j in q.PAIRS:
            V[i,j]=q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
            O[i,j]=q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
        VV.append(q.vec(V));OO.append(q.vec(O))
    G=q.transpose(VV+OO+[q.vec(P)]);F=q.transpose(VV+OO+[q.vec(P)]+[q.vec(v) for v in directions])
    rv=[q.vec(v) for v in R];K=[row+[v[i] for v in rv] for i,row in enumerate(F)]
    assert q.rank(q.transpose(VV))[0]==24 and q.rank(G)[0]==47 and q.rank(F)[0]==48 and q.rank(K)[0]==50
    assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in q.transpose(K))
    cert=source['certificates']['D'];rows=cert['rows'];cols=cert['columns'];assert len(rows)==len(cols)==100
    det=q.determinant([[D[i][j] for j in cols] for i in rows]);assert det and det==q.F(cert['determinant'])
    Q0=q.jvec(q.jacobi(R[0]));Q1=q.jvec(q.jacobi(R[1]));total=q.jvec(q.jacobi({ij:q.add(R[0].get(ij,{}),R[1].get(ij,{})) for ij in set(R[0])|set(R[1])}))
    cross=[v-a-b for v,a,b in zip(total,Q0,Q1)];Qcols=[Q0,cross,Q1];pairings=[]
    for saved_lambda in saved['quadratic_cokernel_functionals']:
        left={int(i):q.F(v) for i,v in saved_lambda['nonzero_coefficients'].items()}
        assert all(sum(v*D[i][j] for i,v in left.items())==0 for j in range(150))
        vals=[sum(v*Q[i] for i,v in left.items()) for Q in Qcols]
        assert vals==[q.F(v) for v in saved_lambda['quadratic_pairing']];pairings.append(vals)
    pdet=q.determinant(pairings);assert pdet==q.F(-9011301580800000000000,7)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'case':'node_tangent_flex',**source_check,
      'method':'Independent standard-library Fraction reconstruction, reusing only the independent five-source review arithmetic. No producer or symbolic engine imported.',
      'full_D_shape':[350,150],'full_D_exact_rank':100,'Euler_vertical_exact_rank':24,'actual_group_vertical_scale_rank':47,'actual_full_source_rank':48,
      'full_kernel_actual_source_plus_two_normals_rank':50,'full_D_times_all57_actual_directions_and_two_normals_zero':True,
      'derivative_minor_rows':rows,'derivative_minor_columns':cols,'derivative_minor_determinant':str(det),
      'literal_two_normals':saved['normals'],'all_three_functionals_annihilate_all150_D_columns':True,
      'quadratic_monomial_order':['n0^2','n0*n1','n1^2'],'quadratic_pairing_matrix':[[str(v) for v in row] for row in pairings],
      'quadratic_pairing_determinant':str(pdet),'full_quadratic_initial_space':'Sym^2(C^2)*','reduced_normal_support':'origin',
      'true_slice_dimension':102,'true_slice_basis':'The two literal normals and the100 coordinate directions in the checked derivative minor complement the actual48-dimensional source tangent.',
      'projective_full_source_closure_component_dimension':23,
      'component_implication':'After48 true linear cuts and100 IFT eliminations, the initial ideal contains all binary quadratics. Hence the cut support is zero-dimensional. Krull gives ambient dimension<=48, equal to the actual source lower bound; its irreducible closure is a component. No dense orbit or local uniqueness claimed.',
      'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
      'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    out=ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_TWO_NORMAL_REVIEW_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'pairing_determinant':str(pdet)}))
if __name__=='__main__':main()
