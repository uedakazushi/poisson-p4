#!/usr/bin/env python3
"""One finite complete three-model bound for special-j squarefree m4 triangles.

The preceding geometry classifies the three whole leading pencils. No point
scan, GB, jet iteration or producer report is imported. Necessary quadratic
compatibilities are kept, including their original-row lift functionals.
"""
import argparse
from fractions import Fraction as F
import hashlib,itertools as it,json,resource,time
from pathlib import Path
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]


def rref(M,n):
    A=[[F(v)for v in row]for row in M];r=0;piv=[]
    for j in range(n):
        k=next((k for k in range(r,len(A))if A[k][j]),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r];d=A[r][j];A[r]=[v/d for v in A[r]]
        for k in range(len(A)):
            if k!=r and A[k][j]:
                d=A[k][j];A[k]=[a-d*b for a,b in zip(A[k],A[r])]
        piv.append(j);r+=1
    return A[:r],piv


def kernel(M,n):
    A,piv=rref(M,n);out=[]
    for j in range(n):
        if j not in piv:
            v=[F()for _ in range(n)];v[j]=F(1)
            for row,p in zip(A,piv):v[p]=-row[j]
            assert all(not sum(a*b for a,b in zip(row,v))for row in M)
            out.append(v)
    return out


def matrix_cert(M):
    rank,cols=q.rank(M);_,rows=q.rank(q.transpose(M))
    det=q.determinant([[M[i][j]for j in cols]for i in rows]);assert det
    return dict(rank=rank,rows=rows,columns=cols,minor_determinant=str(det))


def plus(*Ps):
    return q.antisymmetric({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in q.PAIRS})


def actual(P):
    out={}
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
        out=q.add(out,q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    X=[q.atom(1,i)for i in range(5)];f=q.mul(q.mul(X[0],X[1]),X[2])
    Pf=q.antisymmetric({(0,1):q.der(f,2),(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)})
    first=[]
    for normal in [3,4]:
        for ij in it.combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,mon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in it.combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    assert len(first)==30 and len(second)==9
    d3mon=[e for e in q.monomials(5)if e[3]+e[4]==3];assert len(d3mon)==24
    d4mon=[e for e in q.monomials(5)if e[3]+e[4]==4];assert len(d4mon)==15
    models=[dict(name='J0_secant',beta_u=[F(1,8),F(3,8),F(1,2)],beta_v=[F(3,8),F(1,8),F(1,2)],leading=[(4,0),(0,1)]),
            dict(name='J0_compression',beta_u=[F(),F(1,2),F(1,2)],beta_v=[F(1,2),F(),F(1,2)],leading=[(3,0),(1,1)]),
            dict(name='I0_mixed',beta_u=[F(1,6),F(1,3),F(1,2)],beta_v=[F(1,2),F(),F(1,2)],leading=[(4,0),(1,1)])]
    records=[]
    for model in models:
        bu,bv=model['beta_u'],model['beta_v'];assert sum(bu)==sum(bv)==1
        parallel=[]
        for k in range(5):
            ex=[1-(k-1)*bu[i]-(3-k)*bv[i]for i in range(3)]
            allowed=next((i for i in range(3)if ex==[F(j==i)for j in range(3)]),None)
            parallel.append(dict(normal_u_degree=k,required_base_monomial_exponents=list(map(str,ex)),linear_coefficient_coordinate=allowed))
        assert [(r['normal_u_degree'],r['linear_coefficient_coordinate'])for r in parallel if r['linear_coefficient_coordinate']is not None]==sorted(model['leading'])
        vectors=[]
        for beta in [bu,bv]:
            tau=[v-F(1,3)for v in beta]
            vectors.append([tau[2]-tau[1],tau[0]-tau[2],tau[1]-tau[0]])
        raw=plus(Pf,q.antisymmetric({(i,normal):q.scale(q.atom(1,i,normal),-vectors[j][i])for j,normal in enumerate([3,4])for i in range(3)}))
        rawD=actual(raw);target=q.atom(1,0,1,2,3,4);assert set(rawD)<=set(target)
        kuv=-rawD.get(next(iter(target)),F())/3
        P0=plus(raw,q.antisymmetric({(3,4):q.atom(kuv,3,4)}))
        assert not actual(P0)and all(not v for v in q.jacobi(P0).values())
        J1=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in first])
        D3=[[actual(plus(P0,P)).get(e+(0,0,0),F())for P in first]for e in d3mon]
        necessary=J1+D3;K1=kernel(necessary,30);k1=len(K1)
        allowed=[]
        for k,coordinate in model['leading']:
            e=[0]*5;e[coordinate]=1;e[3]=k;e[4]=4-k;allowed.append(tuple(e))
        shape_mon=[e for e in d4mon if e not in allowed];assert len(shape_mon)==13
        shape=lambda P:[actual(P).get(e+(0,0,0),F())for e in shape_mon]
        J2=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in second])
        shapeL=q.transpose([shape(plus(P0,P))for P in second]);L=J2+shapeL
        lcert=matrix_cert(L);rank=lcert['rank']
        family=[q.antisymmetric({ij:q.add(*(q.scale(P.get(ij,{}),c)for P,c in zip(first,row)if c))for ij in q.PAIRS})for row in K1]
        pairs=list(it.combinations_with_replacement(range(k1),2));quadratic_columns=[]
        own=[q.jvec(q.jacobi(P))+shape(P)for P in family]
        for i,j in pairs:
            if i==j:quadratic_columns.append(own[i])
            else:
                v=q.jvec(q.jacobi(plus(family[i],family[j])))+shape(plus(family[i],family[j]))
                quadratic_columns.append([a-b-c for a,b,c in zip(v,own[i],own[j])])
        Q=q.transpose(quadratic_columns)
        rows,cols=lcert['rows'],lcert['columns'];inv=q.inverse([[L[i][j]for j in cols]for i in rows])
        compatibility=[];lifts=[]
        for i,row in enumerate(L):
            coeff=[sum(row[c]*inv[k][j]for k,c in enumerate(cols))for j in range(rank)]
            assert all(row[c]==sum(v*L[j][c]for v,j in zip(coeff,rows))for c in range(9))
            residual=[v-sum(c*Q[j][at]for c,j in zip(coeff,rows))for at,v in enumerate(Q[i])]
            if any(residual):
                functional={i:F(1)}
                for j,c in zip(rows,coeff):functional[j]=functional.get(j,F())-c
                compatibility.append(residual);lifts.append({str(j):str(c)for j,c in functional.items()if c})
        # A finite rowspace of true ideal members supplies an initial-ideal
        # dimension upper bound. No Buchberger completion is performed.
        exponents=[tuple(int(k==i)+int(k==j)for k in range(k1))for i,j in pairs]
        order=sorted(range(len(pairs)),key=lambda j:tuple(-e for e in reversed(exponents[j])),reverse=True)
        ordered=[[row[j]for j in order]for row in compatibility]
        RR,piv=rref(ordered,len(pairs));leading=[exponents[order[j]]for j in piv]
        supports=[{j for j,e in enumerate(mon)if e}for mon in leading]
        allowed_variables=set(range(k1))-{next(iter(v))for v in supports if len(v)==1}
        independent=None
        for d in range(len(allowed_variables),-1,-1):
            independent=next((v for v in it.combinations(sorted(allowed_variables),d)if all(not e<=set(v)for e in supports)),None)
            if independent is not None:break
        assert independent is not None
        p1_upper=len(independent);fixed_upper=p1_upper+9-rank;image_upper=24+fixed_upper-10
        records.append(dict(type=model['name'],beta_u_logarithmic_coefficients=list(map(str,bu)),beta_v_logarithmic_coefficients=list(map(str,bv)),all_five_parallel_monomial_tests=parallel,
            raw_mixed_diagonal_base_fields=[list(map(str,v))for v in vectors],forced_raw_normal_uv_coefficient=str(kuv),full150_P0=list(map(str,q.vec(P0))),P0_full350_Jacobi_zero=True,P0_actual_Pf_zero=True,
            complete350_weight1_matrix=[[str(v)for v in row]for row in J1],all24_D3_rows=[[str(v)for v in row]for row in D3],weight1_rank_certificate=matrix_cert(necessary),full30_weight1_kernel=[list(map(str,row))for row in K1],weight1_kernel_dimension=k1,
            allowed_two_leading_degree5_monomials=allowed,all13_excluded_leading_monomials=shape_mon,full350_weight2_matrix=[[str(v)for v in row]for row in J2],all13_leading_shape_weight2_rows=[[str(v)for v in row]for row in shapeL],full_weight2_rank_certificate=lcert,weight2_fiber_dimension_upper=9-rank,
            all_quadratic_parameter_pairs=pairs,full363_by_quadratic_monomials_matrix=[[str(v)for v in row]for row in Q],all_nonzero_compatibility_quadrics=[[str(v)for v in row]for row in compatibility],compatibility_original363_row_functionals=lifts,grevlex_quadratic_column_order=order,full_quadratic_rowspace_RREF=[[str(v)for v in row]for row in RR],all_certified_leading_monomials=leading,max_independent_variable_set=list(independent),P1_necessary_dimension_upper=p1_upper,
            fixed_P0_full_source_dimension_upper=fixed_upper,actual_connected_group_fiber_dimension=10,whole_projective_image_dimension_upper=image_upper,G24_closed=image_upper<=24))
    paths=[Path(__file__),ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_M4_SCALAR_REVIEW.md',ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md']
    result=dict(task_id='G',all_arithmetic_passed=True,all_three_whole_leading_pencils_included=True,cases=records,all_special_nonconcurrent_G24_closed=all(r['G24_closed']for r in records),
        coordinates=['x','y','z','u','v'],full30_weight1_basis=[list(map(str,q.vec(P)))for P in first],full9_weight2_basis=[list(map(str,q.vec(P)))for P in second],full350_row_labels=[dict(quad=ij,monomial=e)for ij,e in q.ROWS],
        scope='Only three complete squarefree special-I0/J0 rank2 leading pencils on nonconcurrent triangle. Rank1 scalar and conic+line geometric complements are separate. Concurrent, repeated roots, m2/m3, nonreduced R and whole E are not closed.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_special_nonconcurrent_G24_closed'],cases=[{k:r[k]for k in ['type','weight1_kernel_dimension','weight2_fiber_dimension_upper','P1_necessary_dimension_upper','whole_projective_image_dimension_upper']}for r in records],seconds=result['seconds'])))


if __name__=='__main__':main()
