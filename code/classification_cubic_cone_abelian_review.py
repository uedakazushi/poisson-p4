#!/usr/bin/env python3
"""Independent corrected whole-cone base/orbit/linear-fiber replay.

Uses only accepted reviewer Fraction/sparse arithmetic, never producer/CAS.
The invalid first orbit report is rejected explicitly; V2 is mandatory.
"""
import argparse,hashlib,itertools as it,json,resource,time
from pathlib import Path
from classification_cubic_scroll_abelian_review import Q,PAIRS,add,scale,mul,diff,parse,rref,determinant,bivectors,jacobi_and_divisor,decode

ROOT=Path(__file__).resolve().parents[1]


def mm(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]


def inverse(A):
    n=len(A);R=rref([row+[Q(i==j)for j in range(n)]for i,row in enumerate(A)])
    assert sorted(R)==list(range(n));return [[R[i].get(n+j,Q(0))for j in range(n)]for i in range(n)]


def linear_substitution(poly,matrix):
    n=len(matrix[0]);zero=(0,)*n;variables=[]
    for row in matrix:
        v={}
        for i,c in enumerate(row):
            if c:e=list(zero);e[i]=1;v[tuple(e)]=c
        variables.append(v)
    out={}
    for e,c in poly.items():
        term={zero:c}
        for i,k in enumerate(e):
            for _ in range(k):term=mul(term,variables[i])
        out=add(out,term)
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    paths={k:ROOT/('reports/'+v)for k,v in dict(old='CLASSIFICATION_CUBIC_SCROLL_ZERO_JACOBI_20261107.json',source='CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json',geometry='CLASSIFICATION_CUBIC_CONE_BASE_GEOMETRY_20261107.json',orbits='CLASSIFICATION_CUBIC_CONE_BASE_ORBITS_V2_20261107.json',erratum='CLASSIFICATION_CUBIC_CONE_BASE_ORBITS_ERRATUM_20261107.json').items()}
    data={k:json.loads(p.read_text())for k,p in paths.items()};case=next(v for v in data['source']['cases']if v['name']=='twisted_cubic_cone_03');old=next(v for v in data['old']['cases']if v['name']=='twisted_cubic_cone_03');geometry=data['geometry'];orbits=data['orbits'];erratum=data['erratum']
    assert geometry['input_sha256']==hashlib.sha256(paths['source'].read_bytes()).hexdigest();assert orbits['input_sha256']==hashlib.sha256(paths['geometry'].read_bytes()).hexdigest()
    assert erratum['corrected_report']==str(paths['orbits'].relative_to(ROOT))
    invalid=ROOT/erratum['invalid_universal_scope_report'];assert hashlib.sha256(invalid.read_bytes()).hexdigest()==erratum['sha256']
    quadrics=[parse(t,[f'x{i}'for i in range(5)])for t in ['x0*x2-x1**2','x0*x3-x1*x2','x1*x3-x2**2']]
    mapping=[(3,0,0),(2,1,0),(1,2,0),(0,3,0),(0,0,1)]
    def pull(poly):
        out={}
        for e,c in poly.items():
            f=tuple(sum(e[i]*mapping[i][j]for i in range(5))for j in range(3));out[f]=out.get(f,Q(0))+c
        return {e:c for e,c in out.items()if c}
    assert all(not pull(q)for q in quadrics)
    derivatives=[[diff(q,i)for i in range(5)]for q in quadrics];rows=[];labels=[]
    for a,b in it.combinations(range(3),2):
        for k in range(5):
            columns=[]
            for g in range(3):
                for i,j in PAIRS:columns.append(pull(mul(derivatives[g][k],add(mul(derivatives[a][i],derivatives[b][j]),scale(mul(derivatives[a][j],derivatives[b][i]),-1)))))
            for e in sorted(set().union(*(set(v)for v in columns)),reverse=True):
                rows.append([v.get(e,Q(0))for v in columns]);labels.append(dict(generator_pair=[a,b],ambient_partial=k,parametrization_exponents=list(e)))
    assert len(rows)==78 and rows==[[Q(v)for v in row]for row in case['full_linear_matrix']]and labels==case['all_generator_pair_all_firstpartial_equations']
    R=rref(rows);assert len(R)==12;free=[i for i in range(30)if i not in R];B=[[Q(0)for _ in free]for i in range(30)]
    for k,j in enumerate(free):
        B[j][k]=1
        for p in R:B[p][k]=-R[p].get(j,Q(0))
    assert B==[[Q(v)for v in row]for row in case['complete_original30_by_reduced_source_matrix']]
    cert=case['rank_minor'];assert determinant([[rows[i][j]for j in cert['columns']]for i in cert['rows']])==Q(cert['determinant'])!=0
    identity=[[Q(i==j)for j in range(30)]for i in range(30)];J30,D30=jacobi_and_divisor(bivectors(quadrics,identity));assert J30==decode(old['all350_Jacobi_rows'])and D30==decode(old['all126_actual_D_rows'])
    bases=bivectors(quadrics,B);J,D=jacobi_and_divisor(bases);assert J==decode(case['all350_reduced_Jacobi_rows'])and D==decode(case['all126_reduced_actual_D_rows'])
    vertex=[k for k in range(18)if all(not P[i][j]for i,j in it.combinations(range(4),2)for P in [bases[k]])];base=[k for k in range(18)if k not in vertex]
    assert base==[1,6,10,12,13,15]==geometry['base_six_original18_indices'];assert vertex==geometry['vertex_twelve_original18_indices']and len(vertex)==12
    assert len(rref([[B[i][j]for i in range(30)if 4 in PAIRS[i%10]]for j in vertex]))==12
    for i,j in it.combinations(range(4),2):
        expected={}
        for k,column in enumerate(base):
            for exponent,c in bases[column][i][j].items():
                assert exponent[4]==0
                parameter=tuple(int(t==k)for t in range(6));expected[parameter+exponent[:4]]=c
        assert expected==parse(geometry['complete_base_bivector'][f'{i},{j}'],[f'b{k}'for k in range(6)]+[f'x{k}'for k in range(4)])
    cols18=list(it.combinations_with_replacement(range(18),2));cols6=list(it.combinations_with_replacement(range(6),2));ci6={v:i for i,v in enumerate(cols6)};base_rows=[];fiber=[]
    for row in J:
        a={};m=[{}for _ in vertex]
        for index,c in row.items():
            i,j=cols18[index]
            if i in base and j in base:a[ci6[tuple(sorted((base.index(i),base.index(j))))]]=c
            else:
                assert(i in base)!=(j in base);b=i if i in base else j;v=j if i in base else i;e=[0]*6;e[base.index(b)]=1;m[vertex.index(v)][tuple(e)]=c
        assert not a or not any(m);base_rows.append(a);fiber.append(m)
    assert base_rows==decode(geometry['all350_base_quadratic_rows']);assert all(fiber[i][j]==parse(geometry['all350_linear_vertex_fiber_matrix'][i][j],[f'b{k}'for k in range(6)])for i in range(350)for j in range(12))
    assert len(rref(base_rows))==9
    equations=['(2*b0-b5)*(3*b0-b5)','b1*(3*b0-b5)','2*b0*b2-2*b1*b4+b3*b5','6*b0*b3+2*b1*b4-3*b3*b5','b1*b2','b2*(b2+b3)','(2*b0-b5)*b4+2*b1*b5','12*b1**2-4*b1*b4+(b2+b3)*b5','6*b1*b3-(b2+b3)*b4'];b_names=[f'b{k}'for k in range(6)];eq=[parse(t,b_names)for t in equations]
    def quadratic_row(poly):return {ci6[tuple(i for i,n in enumerate(e)for _ in range(n))]:c for e,c in poly.items()}
    assert rref([quadratic_row(p)for p in eq])==rref(base_rows)
    for text in geometry['base_rowbasis_polynomials']:
        v=quadratic_row(parse(text,b_names));p=min(v);weights=geometry['base_rowbasis_original350_functionals'][str(p)];assert v==add(*(scale(base_rows[int(i)],Q(c))for i,c in weights.items()))
    L=[[Q(v)for v in row]for row in [[1,0,0],[0,0,0],[0,1,0],[0,-1,0],[0,0,1],[2,0,0]]];W=[[Q(v)for v in row]for row in [[1,0,0],[0,1,0],[0,0,0],[0,0,1],[0,6,0],[3,0,0]]]
    assert L==[[Q(v)for v in row]for row in orbits['L3_basis']]and W==[[Q(v)for v in row]for row in orbits['Q2_linear_span_basis']];assert all(not linear_substitution(p,L)for p in eq)
    nilquad=parse('z0*z2-4*z1**2',['z0','z1','z2'])
    for p in eq:
        value=linear_substitution(p,W)
        if value:
            ratio=value[next(iter(value))]/nilquad[next(iter(value))];assert value==scale(nilquad,ratio)
    change=[a+b for a,b in zip(L,W)];inv=inverse(change);assert determinant(change)
    fields={'h':[[Q(v)for v in row]for row in [[3,0,0,0,0],[0,1,0,0,0],[0,0,-1,0,0],[0,0,0,-3,0],[0,0,0,0,0]]],
        'e':[[Q(0)for j in range(5)]for i in range(5)],'f':[[Q(0)for _ in range(5)]for _ in range(5)]}
    for i in range(1,4):fields['e'][i][i-1]=Q(i);fields['f'][i-1][i]=Q(4-i)
    actions={};blocks={}
    pivot_exponents=[(1,0,1,0,0),(1,0,0,1,0),(0,1,0,1,0)]
    for name,T in fields.items():
        X=[{tuple(int(k==j)for k in range(5)):T[i][j]for j in range(5)if T[i][j]}for i in range(5)];rho=[[Q(0)for _ in base]for _ in base]
        for col,k in enumerate(base):
            P=bases[k];image=[[{}for _ in range(5)]for _ in range(5)];coefficients=[Q(0)]*30
            for t,(i,j)in enumerate(PAIRS):
                image[i][j]=add(*(mul(X[l],diff(P[i][j],l))for l in range(5)),*(scale(P[l][j],-T[i][l])for l in range(5)),*(scale(P[i][l],-T[j][l])for l in range(5)))
                for g,e in enumerate(pivot_exponents):coefficients[10*g+t]=image[i][j].get(e,Q(0))
                assert image[i][j]==add(*(scale(quadrics[g],coefficients[10*g+t])for g in range(3)))
            coordinates=[coefficients[j]for j in free];assert all(not coordinates[j]for j in vertex)
            for i,j in enumerate(base):rho[i][col]=coordinates[j]
        assert rho==[[Q(v)for v in row]for row in orbits['SL2_action_on_full_base6'][name]];actions[name]=rho;block=mm(mm(inv,rho),change);blocks[name]=block
        assert block==[[Q(v)for v in row]for row in orbits['SL2_block_actions'][name]];assert all(block[i][j]==0 for i in range(6)for j in range(6)if(i<3)!=(j<3))
    invariantL=parse('4*z0*z1+z2**2',['z0','z1','z2'])
    for block in blocks.values():
        for offset,invariant in [(0,invariantL),(3,nilquad)]:
            vector=[]
            for i in range(3):vector.append({tuple(int(k==j)for k in range(3)):block[offset+i][offset+j]for j in range(3)if block[offset+i][offset+j]})
            assert not add(*(mul(diff(invariant,i),vector[i])for i in range(3)))
    assert [blocks['h'][i][i]for i in range(3)]==[2,-2,0]and[blocks['h'][i][i]for i in range(3,6)]==[2,0,-2]
    for offset in (0,3):
        assert len(rref([[block[offset+i][offset+j]for i in range(3)for j in range(3)]for block in blocks.values()]))==3
    records=[]
    for saved,point,base_dim,rank in zip(orbits['cases'],[[0,0,0,0,1,0],[0,0,1,-1,0,0],[1,0,0,0,0,3]],[3,2,2],[6,6,9]):
        assert list(map(Q,saved['base_representative']))==point
        M=[[sum(c*point[e.index(1)]for e,c in poly.items())for poly in row]for row in fiber]
        assert M==[[Q(v)for v in row]for row in saved['full350_by12_vertex_matrix']]and len(rref(M))==rank
        K=[[Q(v)for v in row]for row in saved['full_vertex_kernel_basis']];dimension=12-rank;assert len(K)==12 and len(K[0])==dimension and len(rref(list(map(list,zip(*K)))))==dimension
        assert mm(M,K)==[[0]*dimension for _ in range(350)]
        cert=saved['rank_minor'];det=determinant([[M[i][j]for j in cert['columns']]for i in cert['rows']]);assert det==Q(cert['determinant'])and det
        total=base_dim+dimension;bound=24+total-1-8;assert total==saved['whole_fixed_affine_source_dimension_upper']and bound==saved['whole_projective_PGL_image_dimension_upper']
        records.append(dict(name=saved['name'],correct_representative=point,base_stratum_dimension=base_dim,full_fiber_rank=rank,complete_fiber_dimension=dimension,rank_minor=str(det),fixed_source_dimension_upper=total,whole_projective_PGL_image_dimension_upper=bound))
    inputs=[Path(__file__),ROOT/'code/classification_cubic_scroll_abelian_review.py',*paths.values(),invalid]
    result=dict(task_id='G',all_passed=True,scope='Whole twisted-cubic-cone rankzero irreducible-reduced-quintic G branch only; no plane/rank2 restriction/E or emptiness assertion.',original30_and_all350_and126_reconstructed=True,linear_shape=[78,30],linear_rank=12,complete_kernel_dimension=18,base6_vertex12_split_verified=True,full350_purebase_and_linearvertex_verified=True,whole_nine_quadratic_support_L3_union_Q2_verified=True,Azero_complement_included=True,actual_SL2_actions_and_two_Sym2_blocks_verified=True,correct_L_invariant='4AC+E^2',correct_Q_invariant='AD-4B^2',invalid_first_semisimple_report_explicitly_rejected=True,cases=records,whole_projective_PGL_image_dimension_upper=24,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,whole_projective_PGL_image_dimension_upper=24,seconds=result['seconds'])))


if __name__=='__main__':main()
