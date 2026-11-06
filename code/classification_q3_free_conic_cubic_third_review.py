#!/usr/bin/env python3
"""Independent finite cubic obstruction of the remaining three complete sources.
Uses only our own Fraction review arithmetic; no producer or symbolic engine.
"""
import hashlib,json,resource,time
from pathlib import Path
from itertools import product,combinations,combinations_with_replacement
import classification_q3_free_conic_cubic_component_review as q
ROOT=Path(__file__).resolve().parents[1]
def mixed(A,B):
    cubic={}
    for ij in combinations(range(5),3):
        cubic[ij]=q.add(*(term for a,b,c in [ij,(ij[1],ij[2],ij[0]),(ij[2],ij[0],ij[1])] for l in range(5) for term in [q.mul(A.get((a,l),{}),q.der(B.get((b,c),{}),l)),q.mul(B.get((a,l),{}),q.der(A.get((b,c),{}),l))]))
    return q.jvec({inds:q.add(*(q.scale(q.mul(q.atom(1,inds[k]),cubic[inds[:k]+inds[k+1:]]),(-1)**k) for k in range(4))) for inds in q.QUADS})
def main():
    resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();kinds=['SOURCES','CANONICAL','NORMALS','CURVES','THIRD'];paths=[ROOT/('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_'+k+'_20261107.json') for k in kinds]
    all_data={k:{row['case']:row for row in json.loads(p.read_text())['cases']} for k,p in zip(kinds,paths)};records=[]
    for name in ['cusp_outside_3','cusp_at_cusp_2','cusp_at_cusp_3']:
        begin=time.monotonic();src=all_data['SOURCES'][name];canon=all_data['CANONICAL'][name];norm=all_data['NORMALS'][name];third=all_data['THIRD'][name];curve=all_data['CURVES'][name]
        source_check=q.full_source_check(src,canon);P=q.decode(canon['canonical_point_brackets']);R=[q.decode(v) for v in norm['literal_normals']];axis=third['free_axis_normal_index'];assert axis==curve['normal_index'];T=q.decode(curve['second_order_correction'])
        Dcols=[]
        for ij,e in q.COLS:
            v={key:dict(value) for key,value in P.items()};inc={e+(0,):q.F(1)}
            v[ij]=q.add(v.get(ij,{}),inc);v[ij[::-1]]=q.scale(v[ij],-1);Dcols.append(q.jvec(q.jacobi(v)))
        D=q.transpose(Dcols);VV=[];OO=[];actual=[]
        for a,b in product(range(5),repeat=2):
            V={};O={}
            for i,j in q.PAIRS:
                V[i,j]=q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
                O[i,j]=q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
            V=q.antisymmetric(V);O=q.antisymmetric(O);VV.append(q.vec(V));OO.append(q.vec(O));actual.extend([V,O])
        actual.append(P);G=q.transpose(VV+OO+[q.vec(P)]);rv=[q.vec(v) for v in R];K=[row+[v[i] for v in rv] for i,row in enumerate(G)]
        assert q.rank(q.transpose(VV))[0]==24 and q.rank(G)[0]==48 and q.rank(K)[0]==48+len(R)
        assert all(sum(a*b for a,b in zip(row,col))==0 for row in D for col in q.transpose(K))
        rows=norm['implicit_D_minor_rows'];cols=norm['implicit_D_minor_columns'];drank=150-48-len(R);assert len(rows)==len(cols)==drank
        det=q.determinant([[D[i][j] for j in cols] for i in rows]);assert det and det==q.F(norm['implicit_D_minor_determinant'])
        monpairs=list(combinations_with_replacement(range(len(R)),2));Qcols=[q.jvec(q.jacobi(R[i])) if i==j else mixed(R[i],R[j]) for i,j in monpairs]
        second_pairs=[]
        for witness in norm['left_null_witnesses']:
            left={int(i):q.F(v) for i,v in witness['functional'].items()}
            assert all(sum(v*D[i][j] for i,v in left.items())==0 for j in range(150))
            pairing=[sum(v*column[i] for i,v in left.items()) for column in Qcols]
            assert pairing==[q.F(v) for v in witness['pairings']];second_pairs.append(pairing)
        if len(R)==2:
            assert q.rank(second_pairs)[0]==2
            axis_square=monpairs.index((axis,axis));assert all(row[axis_square]==0 for row in second_pairs)
        else:assert not second_pairs
        Qa=q.jvec(q.jacobi(R[axis]));tv=q.vec(T)
        assert all(sum(a*b for a,b in zip(row,tv))+v==0 for row,v in zip(D,Qa))
        cubic=mixed(R[axis],T)
        assert {str(i):str(v) for i,v in enumerate(cubic) if v}==third['cubic_vector_sparse']
        lam={int(i):q.F(v) for i,v in third['left_null_functional'].items()}
        assert all(sum(v*D[i][j] for i,v in lam.items())==0 for j in range(150))
        assert all(sum(v*column[i] for i,v in lam.items())==0 for column in Qcols)
        for A in actual:
            cross=mixed(R[axis],A)
            assert sum(v*cross[i] for i,v in lam.items())==0
        pairing=sum(v*cubic[i] for i,v in lam.items());assert pairing and pairing==q.F(third['nonzero_cubic_pairing'])
        rec={'case':name,'all_passed':True,**source_check,'full_D_exact_rank':drank,'group_vertical_scale_rank':48,
          'full_kernel_group_plus_normals_rank':48+len(R),'normal_count':len(R),'free_axis_normal_index':axis,
          'D_times_all_actual_directions_and_normals_zero':True,'derivative_minor_rows':rows,'derivative_minor_columns':cols,'derivative_minor_determinant':str(det),
          'literal_normals':norm['literal_normals'],'quadratic_monomial_pairs':monpairs,'independently_reconstructed_quadratic_pairings':[[str(v) for v in row] for row in second_pairs],
          'quadratic_reduced_support':'the selected free normal axis' if len(R)==2 else 'the sole normal line',
          'literal_second_order_correction':curve['second_order_correction'],'full_D_T_plus_Jacobi_R_zero':True,
          'all350_cubic_vector_entries_reconstructed':True,'lambda_annihilates_all150_D_columns':True,'lambda_annihilates_all_normal_quadratics':True,
          'lambda_annihilates_all51_group_cross_terms':True,'arbitrary_full_kernel_correction_ambiguity_removed':True,
          'cubic_left_functional':third['left_null_functional'],'nonzero_cubic_pairing':str(pairing),
          'genuine_quadratic_plus_cubic_initial_support':'origin','projective_orbit_closure_component_dimension':23,
          'logical_scope':'The cubic functional kills every lower normal quadratic and the entire derivative, so its residual has no term below degree3. The nonzero cubic restriction to the quadratic-support axis is independent of the second-order lift because all group and normal kernel cross terms are killed. This proves reduced support at the origin, not existence of an actual curve.',
          'seconds':time.monotonic()-begin}
        records.append(rec);print(json.dumps({'case':name,'all_passed':True,'pairing':str(pairing),'seconds':rec['seconds']}),flush=True)
    helper=ROOT/'code/classification_q3_free_conic_cubic_component_review.py';paths.append(helper)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'method':'Independent standard-library Fraction reconstruction of all full-source identities, derivatives, genuine quadratics, second-order lift, cubic vector and all kernel ambiguities. No producer or symbolic engine arithmetic imported.',
      'cases':records,'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    out=ROOT/'reports/CLASSIFICATION_Q3_FREE_CONIC_CUBIC_THIRD_REVIEW_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':main()
