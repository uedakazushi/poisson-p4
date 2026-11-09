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


def concurrent_main():
    """One whole four-pencil concurrent bound; all old functions stay immutable."""
    ap=argparse.ArgumentParser(description=concurrent_main.__doc__)
    ap.add_argument('--concurrent',action='store_true');ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    X=[q.atom(1,i)for i in range(5)]
    f=q.mul(q.mul(X[0],X[1]),q.add(X[0],q.scale(X[1],-1)))
    base=lambda G:q.add(q.mul(q.der(f,1),q.der(G,0)),q.scale(q.mul(q.der(f,0),q.der(G,1)),-1))
    raw_basis=[(a,b,k)for a in [3,4]for b in [3,4]for k in [0,1]]
    action=lambda G,a,b,k:q.add(q.mul(q.atom(1,k,b),q.der(G,a)),q.scale(q.mul(X[k],G),-int(a==b)))
    mon6=[e+(0,0,0)for e in q.monomials(6)]
    leading_monomials=[e for e in q.monomials(5)if e[3]+e[4]==4]
    H=[{e+(0,0,0):F(1)}for e in leading_monomials];assert len(H)==15
    first=[]
    for normal in [3,4]:
        for ij in it.combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,mon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in it.combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    assert len(first)==30 and len(second)==9
    A=q.atom(4,3,3,3,4)
    B=q.add(q.atom(F(3,4),3,3,3,3),q.atom(6,3,3,4,4),q.atom(-4,4,4,4,4))
    three=q.add(q.mul(q.add(q.scale(X[0],2),q.scale(X[1],-1)),A),q.mul(X[1],B))
    models=[('J0_secant',[q.atom(1,0,3,3,3,3),q.atom(1,1,4,4,4,4)]),
            ('J0_compression',[q.atom(1,0,3,3,3,4),q.atom(1,1,3,4,4,4)]),
            ('I0_mixed',[q.atom(1,0,3,3,3,3),q.atom(1,1,3,4,4,4)]),
            ('I0_three_distinct',[three])]
    encode=lambda G:[dict(exponent=list(e[:5]),coefficient=str(c))for e,c in sorted(G.items())]
    records=[]
    for name,Gs in models:
        # A full linear system, not a chosen residue/diagonal-connection ansatz.
        connection=[];rhs=[]
        for G in Gs:
            cols=[action(G,*key)for key in raw_basis];target=q.scale(base(G),-1)
            connection.extend([[col.get(e,F())for col in cols]for e in mon6])
            rhs.extend([target.get(e,F())for e in mon6])
        augmented,piv=rref([row+[b]for row,b in zip(connection,rhs)],9)
        assert piv==list(range(8));Mcoef=[row[-1]for row in augmented]
        assert all(sum(a*b for a,b in zip(row,Mcoef))==v for row,v in zip(connection,rhs))
        M={(a,b):q.add(*(q.scale(X[k],c)for (aa,bb,k),c in zip(raw_basis,Mcoef)if aa==a and bb==b))for a in [3,4]for b in [3,4]}
        trace=q.add(M[3,3],M[4,4])
        parallel=lambda G:q.add(base(G),*(q.mul(q.mul(M[a,b],X[b]),q.der(G,a))for a in [3,4]for b in [3,4]),q.scale(q.mul(trace,G),-1))
        assert all(not parallel(G)for G in Gs)
        parallel_matrix=[[parallel(G).get(e,F())for G in H]for e in mon6]
        KL=kernel(parallel_matrix,15);assert len(KL)==len(Gs)
        annihilator=kernel(KL,15);assert len(annihilator)==15-len(KL)
        for G in Gs:
            gv=[G.get(e+(0,0,0),F())for e in leading_monomials]
            assert all(not sum(a*b for a,b in zip(row,gv))for row in annihilator)
        P0=q.antisymmetric({(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0),
                (2,3):q.add(*(q.mul(M[3,b],X[b])for b in [3,4])),
                (2,4):q.add(*(q.mul(M[4,b],X[b])for b in [3,4]))})
        assert not actual(P0)and all(not v for v in q.jacobi(P0).values())
        J1=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in first])
        d3mon=[e for e in q.monomials(5)if e[3]+e[4]==3]
        D3=[[actual(plus(P0,P)).get(e+(0,0,0),F())for P in first]for e in d3mon]
        necessary=J1+D3;K1=kernel(necessary,30);k1=len(K1)
        shape=lambda P:[sum(a*actual(P).get(e+(0,0,0),F())for a,e in zip(row,leading_monomials))for row in annihilator]
        J2=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in second])
        shapeL=q.transpose([shape(plus(P0,P))for P in second]);L=J2+shapeL
        lcert=matrix_cert(L);rank=lcert['rank']
        family=[q.antisymmetric({ij:q.add(*(q.scale(P.get(ij,{}),c)for P,c in zip(first,row)if c))for ij in q.PAIRS})for row in K1]
        pairs=list(it.combinations_with_replacement(range(k1),2));quadratic_columns=[]
        own=[q.jvec(q.jacobi(P))+shape(P)for P in family]
        for i,j in pairs:
            if i==j:quadratic_columns.append(own[i])
            else:
                vv=q.jvec(q.jacobi(plus(family[i],family[j])))+shape(plus(family[i],family[j]))
                quadratic_columns.append([a-b-c for a,b,c in zip(vv,own[i],own[j])])
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
        exponents=[tuple(int(k==i)+int(k==j)for k in range(k1))for i,j in pairs]
        order=sorted(range(len(pairs)),key=lambda j:tuple(-e for e in reversed(exponents[j])),reverse=True)
        RR,piv=rref([[row[j]for j in order]for row in compatibility],len(pairs));leading=[exponents[order[j]]for j in piv]
        supports=[{j for j,e in enumerate(mon)if e}for mon in leading]
        independent=None
        for d in range(k1,-1,-1):
            independent=next((v for v in it.combinations(range(k1),d)if all(not e<=set(v)for e in supports)),None)
            if independent is not None:break
        assert independent is not None
        # Actual connected normal centralizer, not all GL2 for non-scalar M.
        centralizer=[]
        for k in [0,1]:
            Mk=[[next((c for (aa,bb,kk),c in zip(raw_basis,Mcoef)if (aa,bb,kk)==(a,b,k)),F())for b in [3,4]]for a in [3,4]]
            for a,b in it.product(range(2),repeat=2):
                centralizer.append([int(a==c)*Mk[d][b]-Mk[a][c]*int(d==b)for c,d in it.product(range(2),repeat=2)])
        CK=kernel(centralizer,4);group_dimension=3+6+len(CK)
        fixed_upper=len(independent)+9-rank;image_upper=24+fixed_upper-group_dimension
        records.append(dict(type=name,whole_leading_generators=[encode(G)for G in Gs],
            full_connection_basis_labels=raw_basis,full_connection_equations=[[str(v)for v in row]for row in connection],full_connection_rhs=list(map(str,rhs)),unique_raw_normal_matrix_linear_coefficients=list(map(str,Mcoef)),connection_rank_certificate=matrix_cert(connection),
            all15_leading_monomials=leading_monomials,full_parallel_leading_matrix=[[str(v)for v in row]for row in parallel_matrix],full_parallel_leading_kernel=[list(map(str,row))for row in KL],parallel_leading_rank_certificate=matrix_cert(parallel_matrix),leading_annihilator=[list(map(str,row))for row in annihilator],
            full150_P0=list(map(str,q.vec(P0))),raw_normal_normal_bracket_zero_forced_by_actual_D2=True,P0_full350_Jacobi_zero=True,P0_actual_Pf_zero=True,
            complete350_weight1_matrix=[[str(v)for v in row]for row in J1],all24_D3_rows=[[str(v)for v in row]for row in D3],weight1_rank_certificate=matrix_cert(necessary),full30_weight1_kernel=[list(map(str,row))for row in K1],weight1_kernel_dimension=k1,
            full350_weight2_matrix=[[str(v)for v in row]for row in J2],all_leading_shape_weight2_rows=[[str(v)for v in row]for row in shapeL],full_weight2_rank_certificate=lcert,weight2_fiber_dimension_upper=9-rank,
            all_quadratic_parameter_pairs=pairs,full_original_by_quadratic_monomials_matrix=[[str(v)for v in row]for row in Q],all_nonzero_compatibility_quadrics=[[str(v)for v in row]for row in compatibility],compatibility_original_row_functionals=lifts,
            grevlex_quadratic_column_order=order,full_quadratic_rowspace_RREF=[[str(v)for v in row]for row in RR],all_certified_leading_monomials=leading,max_independent_variable_set=list(independent),P1_necessary_dimension_upper=len(independent),
            normal_centralizer_matrix=[[str(v)for v in row]for row in centralizer],full_normal_centralizer_kernel=[list(map(str,row))for row in CK],normal_centralizer_rank_certificate=matrix_cert(centralizer),actual_connected_group_fiber_dimension=group_dimension,
            fixed_P0_full_source_dimension_upper=fixed_upper,whole_projective_image_dimension_upper=image_upper,G24_closed=image_upper<=24))
        print(json.dumps(dict(type=name,P1=k1,compatibility_rank=len(RR),P1_bound=len(independent),P2_fiber=9-rank,group=group_dimension,image_bound=image_upper)),flush=True)
    paths=[Path(__file__),ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_PENCILS.md',ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md',ROOT/'research/classification/PLANE_CONCURRENT_M4_SPECIAL_PENCILS.md']
    result=dict(task_id='G',all_arithmetic_passed=True,all_four_whole_leading_pencils_included=True,cases=records,
        coordinates=['x','y','z','u','v'],base_cubic=encode(f),full30_weight1_basis=[list(map(str,q.vec(P)))for P in first],full9_weight2_basis=[list(map(str,q.vec(P)))for P in second],full350_row_labels=[dict(quad=ij,monomial=e)for ij,e in q.ROWS],
        all_concurrent_special_rank2_G24_closed=all(r['G24_closed']for r in records),scope='Only concurrent reduced cubic rank2 restriction, m4 generic-squarefree special I0/J0 rank2 leading pencil. No generic-j/rank1 reopening, repeated roots, m2/m3 or nonreduced R claims.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_concurrent_special_rank2_G24_closed'],seconds=result['seconds'])))


def _exact_polynomial_quotient(p,d):
    out={};rem=dict(p);de=max(d);dc=d[de]
    while rem:
        e=max(rem);assert all(a>=b for a,b in zip(e,de)),(p,d,rem)
        term={tuple(a-b for a,b in zip(e,de)):rem[e]/dc}
        out=q.add(out,term);rem=q.add(rem,q.scale(q.mul(d,term),-1))
    assert q.mul(out,d)==p
    return out


def _whole_filtered_plane_certificate(P0,annihilator,group_dimension):
    """Complete necessary source, constant fibers, original-row quadrics only."""
    first=[]
    for normal in [3,4]:
        for ij in it.combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,mon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in it.combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    d3mon=[e for e in q.monomials(5)if e[3]+e[4]==3]
    d4mon=[e for e in q.monomials(5)if e[3]+e[4]==4]
    shape=lambda P:[sum(c*actual(P).get(e+(0,0,0),F())for c,e in zip(row,d4mon))for row in annihilator]
    J1=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in first])
    D3=[[actual(plus(P0,P)).get(e+(0,0,0),F())for P in first]for e in d3mon]
    K1=kernel(J1+D3,30);k1=len(K1)
    J2=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in second])
    shapeL=q.transpose([shape(plus(P0,P))for P in second]);L=J2+shapeL
    lcert=matrix_cert(L);rank=lcert['rank'];K2=kernel(L,9)
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
    exponents=[tuple(int(k==i)+int(k==j)for k in range(k1))for i,j in pairs]
    order=sorted(range(len(pairs)),key=lambda j:tuple(-e for e in reversed(exponents[j])),reverse=True)
    RR,piv=rref([[row[j]for j in order]for row in compatibility],len(pairs));leading=[exponents[order[j]]for j in piv]
    supports=[{j for j,e in enumerate(mon)if e}for mon in leading]
    independent=None
    for d in range(k1,-1,-1):
        independent=next((v for v in it.combinations(range(k1),d)if all(not e<=set(v)for e in supports)),None)
        if independent is not None:break
    assert independent is not None
    K2family=[q.antisymmetric({ij:q.add(*(q.scale(P.get(ij,{}),c)for P,c in zip(second,row)if c))for ij in q.PAIRS})for row in K2]
    qzero=not any(v for row in Q for v in row)
    common_z=qzero and all(not v or(2 in ij and not e[2])for P in [P0]+family+K2family for (ij,e),v in zip(q.COLS,q.vec(P)))
    fixed=len(independent)+9-rank;image=24+fixed-group_dimension
    return dict(full30_weight1_basis=[list(map(str,q.vec(P)))for P in first],full9_weight2_basis=[list(map(str,q.vec(P)))for P in second],
        full350_weight1_matrix=[[str(v)for v in row]for row in J1],all24_D3_rows=[[str(v)for v in row]for row in D3],weight1_rank_certificate=matrix_cert(J1+D3),full30_weight1_kernel=[list(map(str,row))for row in K1],weight1_kernel_dimension=k1,
        full350_weight2_matrix=[[str(v)for v in row]for row in J2],all_leading_shape_weight2_rows=[[str(v)for v in row]for row in shapeL],weight2_rank_certificate=lcert,full9_weight2_kernel=[list(map(str,row))for row in K2],weight2_fiber_dimension_upper=9-rank,
        all_quadratic_parameter_pairs=pairs,full_original_by_quadratic_monomials_matrix=[[str(v)for v in row]for row in Q],all_nonzero_compatibility_quadrics=[[str(v)for v in row]for row in compatibility],compatibility_original_row_functionals=lifts,
        grevlex_quadratic_column_order=order,full_quadratic_rowspace_RREF=[[str(v)for v in row]for row in RR],all_certified_leading_monomials=leading,max_independent_variable_set=list(independent),P1_necessary_dimension_upper=len(independent),
        all_original_quadratic_coefficients_zero=qzero,entire_necessary_source_common_partial_z_factor=common_z,rank4_source_empty=common_z,
        actual_connected_group_fiber_dimension=group_dimension,fixed_P0_full_source_dimension_upper=fixed,whole_projective_image_dimension_upper=image,G24_closed=common_z or image<=24)


def repeated_211_main():
    """One whole six-model generic-211 repeated-quartic source test."""
    ap=argparse.ArgumentParser(description=repeated_211_main.__doc__)
    ap.add_argument('--repeated-211',action='store_true');ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    X=[q.atom(1,i)for i in range(5)];x,y,z=X[:3]
    tri=q.mul(q.mul(x,y),z);con=q.mul(q.mul(x,y),q.add(x,q.scale(y,-1)))
    tan=q.mul(z,q.add(q.mul(x,z),q.scale(q.mul(y,y),-1)))
    models=[('triangle_diagonal',tri,x,{},y,10,2),('triangle_Jordan',tri,{},q.scale(y,F(1,2)),x,10,2),
            ('concurrent_diagonal',con,x,{},y,11,2),('concurrent_Jordan',con,{},q.scale(y,F(1,2)),x,11,2),
            ('concurrent_three_used',con,x,q.scale(x,-1),q.add(x,q.scale(y,-1)),10,1),
            ('tangent_conic_net',tan,x,y,z,9,1)]
    mon2=[e for e in q.monomials(2)if not e[3]and not e[4]]
    VF=[]
    for a,b in it.permutations(range(3),2):VF.append([X[b]if i==a else{}for i in range(3)])
    VF.extend([[x,q.scale(y,-1),{}],[x,{},q.scale(z,-1)]])
    contraction=lambda V:[q.add(q.mul(z,V[1]),q.scale(q.mul(y,V[2]),-1)),q.add(q.mul(x,V[2]),q.scale(q.mul(z,V[0]),-1)),q.add(q.mul(y,V[0]),q.scale(q.mul(x,V[1]),-1))]
    form_labels=[(i,e+(0,0,0))for i in range(3)for e in mon2]
    field_map=[[contraction(V)[i].get(e,F())for V in VF]for i,e in form_labels]
    records=[];encode=lambda G:[dict(exponent=list(e[:5]),coefficient=str(c))for e,c in sorted(G.items())]
    for name,f,A,B,C,group,leading_dim in models:
        T=q.add(q.mul(A,C),q.scale(q.mul(B,B),-1));assert T
        dt=[_exact_polynomial_quotient(q.mul(f,q.der(T,i)),T)for i in range(3)]
        fc=_exact_polynomial_quotient(f,C)
        fbt=_exact_polynomial_quotient(q.mul(f,B),q.mul(C,T))if B else{}
        Ba=[q.add(q.scale(q.der(f,i),F(1,2)),q.scale(dt[i],F(-1,4)))for i in range(3)]
        Bd=[q.add(q.scale(q.der(f,i),F(1,2)),q.scale(q.mul(fc,q.der(C,i)),-1),q.scale(dt[i],F(1,4)))for i in range(3)]
        Bb=[q.add(q.scale(q.mul(fbt,q.der(T,i)),F(1,2)),q.scale(q.mul(fc,q.der(B,i)),-1))for i in range(3)]
        beta={(3,3):Ba,(3,4):[{}, {}, {}],(4,3):Bb,(4,4):Bd}
        # Full rational curvature numerator: d beta - beta wedge beta.
        curvature=[]
        for a,b in it.product([3,4],repeat=2):
            for i,j in it.combinations(range(3),2):
                v=q.add(q.mul(f,q.add(q.der(beta[a,b][j],i),q.scale(q.der(beta[a,b][i],j),-1))),
                    q.scale(q.mul(q.der(f,i),beta[a,b][j]),-1),q.mul(q.der(f,j),beta[a,b][i]),
                    *(q.scale(q.add(q.mul(beta[a,c][i],beta[c,b][j]),q.scale(q.mul(beta[a,c][j],beta[c,b][i]),-1)),-1)for c in [3,4]))
                assert not v;curvature.append(v)
        fields={};field_coeff=[]
        for a,b in it.product([3,4],repeat=2):
            target=[q.scale(q.add(beta[a,b][i],q.scale(q.der(f,i),-F(a==b,3))),3)for i in range(3)]
            rhs=[target[i].get(e,F())for i,e in form_labels]
            rr,piv=rref([row+[v]for row,v in zip(field_map,rhs)],9);assert piv==list(range(8))
            coef=[row[-1]for row in rr];field_coeff.append(list(map(str,coef)))
            fields[a,b]=[q.add(*(q.scale(V[i],c)for V,c in zip(VF,coef)))for i in range(3)]
            assert contraction(fields[a,b])==target
        Pf=q.antisymmetric({(0,1):q.der(f,2),(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)})
        raw=plus(Pf,q.antisymmetric({(i,a):q.scale(q.add(*(q.mul(X[b],fields[a,b][i])for b in [3,4])),-1)for i in range(3)for a in [3,4]}))
        normal=[q.atom(1,3,3),q.atom(1,3,4),q.atom(1,4,4)]
        d0=actual(raw);dcols=[q.add(actual(plus(raw,q.antisymmetric({(3,4):m}))),q.scale(d0,-1))for m in normal]
        dmon=[e+(0,0,0)for e in q.monomials(5)]
        DK=[[p.get(e,F())for p in dcols]for e in dmon];rhs=[-d0.get(e,F())for e in dmon]
        rr,piv=rref([row+[v]for row,v in zip(DK,rhs)],4);assert piv==list(range(3))
        k=[row[-1]for row in rr];P0=plus(raw,q.antisymmetric({(3,4):q.add(*(q.scale(m,c)for m,c in zip(normal,k)))}))
        assert not actual(P0)and all(not v for v in q.jacobi(P0).values())
        leading=[e for e in q.monomials(5)if e[3]+e[4]==4]
        H=[{e+(0,0,0):F(1)}for e in leading]
        operator=[]
        for G in H:
            operator.append([q.add(q.mul(f,q.der(G,i)),q.scale(q.mul(q.der(f,i),G),-1),
                *(q.mul(q.mul(beta[a,b][i],X[b]),q.der(G,a))for a,b in it.product([3,4],repeat=2)),
                q.scale(q.mul(q.add(Ba[i],Bd[i]),G),-1))for i in range(3)])
        pmon=[e+(0,0,0)for e in q.monomials(7)]
        PM=[[col[i].get(e,F())for col in operator]for i in range(3)for e in pmon]
        KL=kernel(PM,15);assert len(KL)==leading_dim
        ann=kernel(KL,15)
        source=_whole_filtered_plane_certificate(P0,ann,group)
        records.append(dict(type=name,base_cubic=encode(f),Q_A=encode(A),Q_B=encode(B),Q_C=encode(C),Q_determinant=encode(T),full_f_beta={str(ij):[encode(p)for p in form]for ij,form in beta.items()},all12_curvature_numerators_zero=True,
            full18x8_form_field_map=[[str(v)for v in row]for row in field_map],full4x8_tracefree_field_coefficients=field_coeff,forced_raw_uv_quadratic_coefficients=list(map(str,k)),raw_uv_forcing_rank_certificate=matrix_cert(DK),full150_P0=list(map(str,q.vec(P0))),P0_full350_Jacobi_zero=True,P0_actual_Pf_zero=True,
            all15_leading_monomials=leading,full_parallel_leading_matrix=[[str(v)for v in row]for row in PM],parallel_leading_rank_certificate=matrix_cert(PM),full_parallel_leading_kernel=[list(map(str,row))for row in KL],leading_annihilator=[list(map(str,row))for row in ann],**source))
        print(json.dumps(dict(type=name,P1=source['weight1_kernel_dimension'],P2_fiber=source['weight2_fiber_dimension_upper'],P1_bound=source['P1_necessary_dimension_upper'],image_bound=source['whole_projective_image_dimension_upper'],common_z_empty=source['rank4_source_empty'],G24_closed=source['G24_closed'])),flush=True)
    paths=[Path(__file__),ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_M4_REPEATED_LEADING_GEOMETRY.md',ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md']
    result=dict(task_id='G',all_arithmetic_passed=True,all_six_whole_generic211_sources_included=True,cases=records,
        all_generic211_high_rank_G24_closed=all(c['G24_closed']for c in records),coordinates=['x','y','z','u','v'],full350_row_labels=[dict(quad=ij,monomial=e)for ij,e in q.ROWS],
        scope='Only six complete reduced-reducible-cubic generic211 image-rank2/3 m4 leading normal forms. All image-rank1/22/4 reduced separately to accepted essential3 source theorem. Generic31 rank2, m2/m3, nonreducedR remain explicit.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_generic211_high_rank_G24_closed'],seconds=result['seconds'])))


def triple_root_quintics_main():
    """One full150 linear strict-Casimir source batch for both actual F31 forms."""
    ap=argparse.ArgumentParser(description=triple_root_quintics_main.__doc__)
    ap.add_argument('--triple-root-quintics',action='store_true');ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    X=[q.atom(1,i)for i in range(5)]
    mon6=[e+(0,0,0)for e in q.monomials(6)]
    mon1=[e+(0,0,0)for e in q.monomials(1)]
    labels=[dict(type='Casimir',component=i,monomial=e[:5])for i in range(5)for e in mon6]
    labels +=[dict(type='divergence',component=i,monomial=e[:5])for i in range(5)for e in mon1]
    records=[]
    for eps in [0,1]:
        Fp=q.add(q.atom(1,0,2,2,2,2),q.atom(1,1,2,2,2,3),q.atom(eps,2,2,3,3,3),q.atom(1,3,3,3,3,3))
        gradient=[q.der(Fp,i)for i in range(5)]
        assert not q.add(q.mul(X[3],gradient[0]),q.scale(q.mul(X[2],gradient[1]),-1))
        cols=[]
        for (i,j),ex in q.COLS:
            m={ex+(0,0,0):F(1)};cas=[{}for _ in range(5)];dv=[{}for _ in range(5)]
            cas[i]=q.mul(m,gradient[j]);cas[j]=q.scale(q.mul(m,gradient[i]),-1)
            dv[i]=q.scale(q.der(m,j),-1);dv[j]=q.der(m,i)
            cols.append([cas[a].get(e,F())for a in range(5)for e in mon6]+[dv[a].get(e,F())for a in range(5)for e in mon1])
        M=q.transpose(cols);K=kernel(M,150);cert=matrix_cert(M)
        assert len(K)==150-cert['rank']
        degrees=[]
        for d in [0,1,2]:
            mons=[e for e in q.monomials(d)if not e[4]]
            outmons=[e+(0,0,0)for e in q.monomials(d+4)if not e[4]]
            C=[q.mul({e+(0,0,0):F(1)},gradient[i])for i in range(4)for e in mons]
            GM=[[p.get(e,F())for p in C]for e in outmons]
            GK=kernel(GM,len(C))
            degrees.append(dict(degree=d,coefficient_labels=[dict(component=i,monomial=e)for i in range(4)for e in mons],full_gradient_map=[[str(v)for v in row]for row in GM],full_kernel=[list(map(str,row))for row in GK],rank_certificate=matrix_cert(GM),kernel_dimension=len(GK)))
        if not eps:
            assert [v['kernel_dimension']for v in degrees]==[0,2,8]
            Z=[X[3],q.scale(X[2],-1),{},{}]
            Y=[q.scale(X[0],-4),q.scale(X[1],-3),X[2],{}]
            mons=[v['monomial']for v in degrees[2]['coefficient_labels'][:10]]
            expected=[]
            for V in [Z,Y]:
                for k in range(4):
                    field=[q.mul(X[k],p)for p in V]
                    expected.append([field[i].get(e+(0,0,0),F())for i in range(4)for e in mons])
            assert q.rank(expected)[0]==8
            assert all(not sum(a*b for a,b in zip(row,col))for row in [[F(v)for v in row]for row in degrees[2]['full_gradient_map']]for col in expected)
        dim=len(K);image=18+dim-1
        records.append(dict(epsilon=eps,F=[dict(exponent=e[:5],coefficient=str(c))for e,c in sorted(Fp.items())],full1075x150_strict_Casimir_matrix=[[str(v)for v in row]for row in M],rank_certificate=cert,full150_kernel=[list(map(str,row))for row in K],strict_Casimir_source_dimension=dim,all_kernel_original_rows_zero=True,complete_gradient_syzygies_degrees_0_1_2=degrees,
            actualD_Ga_symmetry_Z_checked=True,actualD_family_dimension_upper=18,necessary_projective_Poisson_image_upper=image,F0_all_base_gradient_rows_in_pointwise_two_distribution=(not eps),F0_vertex_adds_at_most_one_then_skew_rank_at_most2=(not eps),rank4_source_empty=(not eps),G24_closed=(not eps)or image<=24))
        print(json.dumps(dict(epsilon=eps,strict_Casimir_dimension=dim,gradient_syzygy_dimensions=[v['kernel_dimension']for v in degrees],projective_image_upper=image,rank4_empty=(not eps),G24_closed=(not eps)or image<=24)),flush=True)
    paths=[Path(__file__),ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_M4_TRIPLE_ROOT_QUINTICS.md',ROOT/'research/classification/ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md']
    result=dict(task_id='G',all_arithmetic_passed=True,all_two_whole_actualD_normalforms_included=True,cases=records,all_m4_generic31_rank2_G24_closed=all(c['G24_closed']for c in records),coordinates=['A','B','u','v','w'],full1075_row_labels=labels,full150_column_labels=[dict(pair=ij,monomial=e)for ij,e in q.COLS],
        scope='Whole actual irreducible m4 generic31 leading-image-rank2 quintics, arbitrary rank2-plane zero cubic R including nonreducedR. Reduced image-rank1 is separate accepted essential3 source. No m2/m3 or other leading/support claim.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_m4_generic31_rank2_G24_closed'],seconds=result['seconds'])))


def m2_double_tangent_main():
    """One bounded full-source test for three fixed and one whole tangent model."""
    from classification_unimodular_logarithmic_review import RF,ua,um,ud,ug
    import sympy as sy
    ap=argparse.ArgumentParser(description=m2_double_tangent_main.__doc__)
    ap.add_argument('--m2-double-tangent',action='store_true');ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jobpath=ROOT/'research/classification/jobs/plane_m2_double_tangent.json'
    job=json.loads(jobpath.read_text())
    for row in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    assert not args.output.exists()
    def compact(v):
        return (v.p[0]if v.p else F())if isinstance(v,RF)and v.q==(F(1),)and len(v.p)<=1 else v
    def encode(v):
        if isinstance(v,RF):return dict(numerator=list(map(str,v.p)),denominator=list(map(str,v.q)))
        if isinstance(v,AF):return dict(residue=list(map(str,v.p)),modulus=list(map(str,v.mod)))
        return str(v)
    def polyencode(p):return [dict(exponent=e[:5],coefficient=encode(v))for e,v in sorted(p.items())]
    def menc(M):return [[encode(v)for v in row]for row in M]
    class AF:
        """QQ[t]/irreducible wall, including every complex embedding at once."""
        mod=()
        def __init__(self,v=0):
            if isinstance(v,AF):self.p=v.p;return
            if isinstance(v,RF):
                a=AF(v.p);b=AF(v.q);self.p=(a/b).p;return
            self.p=ud(tuple(v)if isinstance(v,(list,tuple))else(F(v),),self.mod)[1]
        def __add__(self,v):return AF(ua(self.p,AF(v).p))
        __radd__=__add__
        def __neg__(self):return AF(tuple(-v for v in self.p))
        def __sub__(self,v):return self+-AF(v)
        def __rsub__(self,v):return AF(v)+-self
        def __mul__(self,v):return AF(um(self.p,AF(v).p))
        __rmul__=__mul__
        def inv(self):
            assert self.p;a,b=self.mod,self.p;s,t=(),(F(1),)
            while b:
                quotient,rem=ud(a,b);a,b=b,rem;s,t=t,ua(s,tuple(-v for v in um(quotient,t)))
            assert len(a)==1;return AF(tuple(v/a[0]for v in s))
        def __truediv__(self,v):return self*AF(v).inv()
        def __rtruediv__(self,v):return AF(v)*self.inv()
        def __bool__(self):return bool(self.p)
        def __eq__(self,v):return self.p==AF(v).p
    def cert(M,n):
        B=[[compact(v)for v in row]for row in M];ids=list(range(len(B)));piv=[];rows=[];det=F(1)
        for j in range(n):
            r=len(piv);k=next((k for k in range(r,len(B))if B[k][j]),None)
            if k is None:continue
            B[r],B[k]=B[k],B[r];ids[r],ids[k]=ids[k],ids[r];d=B[r][j]
            rows.append(ids[r]);piv.append(j);det=compact(det*d)
            B[r]=[compact(v/d)if v else F()for v in B[r]]
            active=[l for l in range(j,n)if B[r][l]]
            for k in range(r+1,len(B)):
                c=B[k][j]
                if c:
                    for l in active:B[k][l]=compact(B[k][l]-c*B[r][l])
        rank=len(piv)
        for r in range(rank-1,-1,-1):
            j=piv[r];active=[l for l in range(j+1,n)if B[r][l]]
            for k in range(r):
                c=B[k][j]
                if c:
                    for l in active:B[k][l]=compact(B[k][l]-c*B[r][l])
                    B[k][j]=F()
        K=[]
        for j in range(n):
            if j not in piv:
                v=[F()for _ in range(n)];v[j]=F(1)
                for row,i in zip(B,piv):v[i]=-row[j]
                K.append(v)
        for row in M:
            nz=[(i,v)for i,v in enumerate(row)if v]
            assert all(not sum((c*v[i]for i,c in nz if v[i]),F())for v in K)
        return dict(rank=rank,rows=rows,columns=piv,minor_determinant=encode(det),nullity=n-rank,
                    complete_kernel=menc(K),all_original_rows_annihilate_kernel=True),K,det
    X=[q.atom(1,i)for i in range(5)];x,y,z,u,v=X
    first=[]
    for normal in [3,4]:
        for ij in it.combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,mon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in it.combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    assert len(first)==30 and len(second)==9
    universal=[];basis=first+second
    for i,j in it.combinations_with_replacement(range(39),2):
        P=plus(basis[i],basis[j]);J=q.jvec(q.jacobi(P));D=actual(P)
        if i==j:J=[a/4 for a in J];D=q.scale(D,F(1,4))
        if any(J)or D:universal.append(dict(pair=[i,j],all_nonzero_Jacobi_rows={str(k):encode(a)for k,a in enumerate(J)if a},actualD=polyencode(D)))
    f=q.mul(z,q.add(q.mul(x,z),q.scale(q.mul(y,y),-1)))
    Pf=q.antisymmetric({(0,1):q.der(f,2),(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)})
    N=[q.scale(y,2),z,{}]
    def fixed(name,cube,c):
        XX=[q.scale(a,-8 if cube else-2)for a in N];YY=[q.scale(a,-5 if cube else-2)for a in N]
        ZZ=[z,{},{}]if cube else[q.scale(a,c)for a in N]
        P=plus(Pf,q.antisymmetric({(i,3):q.mul(u,XX[i])for i in range(3)}),
               q.antisymmetric({(i,4):q.add(q.mul(u,ZZ[i]),q.mul(v,YY[i]))for i in range(3)}),
               q.antisymmetric({(3,4):{}if cube else q.mul(u,u)}))
        target=q.atom(-9,2,2,2,3,3)if cube else q.scale(q.mul(f,q.mul(u,u)),3)
        return name,P,target
    T=RF((F(),F(1)))
    moving=q.antisymmetric({(0,1):q.scale(q.add(q.mul(z,x),q.scale(q.mul(y,y),-1)),3),
        (0,3):q.add(q.scale(q.mul(y,u),2*T+3),q.scale(q.mul(z,v),-3)),
        (0,4):q.add(q.scale(q.mul(x,u),-3),q.scale(q.mul(y,v),2*T+3)),
        (1,3):q.scale(q.mul(z,u),T+3),(1,4):q.add(q.scale(q.mul(y,u),-3),q.scale(q.mul(z,v),T)),
        (3,4):q.scale(q.mul(u,u),T+3)})
    target_moving=q.scale(q.mul(z,q.mul(q.add(q.mul(y,u),q.mul(z,v)),q.add(q.mul(y,u),q.mul(z,v)))),3*T)
    models=[fixed('constant_hf_Z0',False,0),fixed('constant_hf_ZN',False,1),fixed('constant_hcube',True,0),('moving_all_t_nonzero',moving,target_moving)]
    eulers=[]
    for a,b in it.product(range(5),repeat=2):
        eulers.append(q.vec(q.antisymmetric({(i,j):q.add(q.atom(1,i,b)if j==a else{},q.atom(-1,j,b)if i==a else{})for i,j in q.PAIRS})))
    assert cert(q.transpose(eulers),25)[0]['rank']==24
    levi=[(a,b)for a,b in it.product(range(5),repeat=2)if(a<3)==(b<3)]
    assert len(levi)==13
    def operators(P):
        J1=q.transpose([q.jvec(q.jacobi(plus(P,Q)))for Q in first]);J2=q.transpose([q.jvec(q.jacobi(plus(P,Q)))for Q in second])
        orbit=[]
        for a,b in levi:
            O=q.antisymmetric({(i,j):q.add(P.get((b,j),{})if i==a else{},P.get((i,b),{})if j==a else{},q.scale(q.mul(X[b],q.der(P.get((i,j),{}),a)),-1))for i,j in q.PAIRS})
            orbit.append(q.vec(O))
        G=q.transpose(eulers+orbit+[q.vec(P)])
        return J1,J2,G
    def process(name,P,target,pdim,field='QQ',prepared=None):
        assert all(not a for a in q.jacobi(P).values())and actual(P)==target,(name,'P0 Jacobi/actualPf mismatch')
        J1,J2,G=prepared if prepared is not None else operators(P)
        c1,K1,d1=cert(J1,30);c2,K2,d2=cert(J2,9);cg,_,dg=cert(G,39)
        # Euler rank24, and projective Levi orbit rank=rank(G)-25.
        sL=37-cg['rank'];group_fiber=sL+6
        upper=24-group_fiber+pdim+len(K1)+len(K2)
        u_divisible=all(e[3]for i in range(5)for e in P.get((min(i,3),max(i,3)),{})if i!=3)
        bad_at=[at for at,Q in enumerate(first)if any(e[3]==0 for i in range(5)for e in Q.get((min(i,3),max(i,3)),{})if i!=3)]
        u_factor=u_divisible and all(not row[j]for row in K1 for j in bad_at)
        linD=[polyencode(q.add(actual(plus(P,Q)),q.scale(target,-1)))for Q in basis]
        out=dict(type=name,field=field,source_parameter_dimension=pdim,full150_P0=list(map(encode,q.vec(P))),P0_full350_Jacobi_zero=True,
            P0_actualG2=polyencode(target),full350by30_J1=menc(J1),full350by9_J2=menc(J2),weight1_certificate=c1,weight2_certificate=c2,
            full150_Euler_Levi_scale_matrix=menc(G),normalizer_certificate=cg,graded_Levi_projective_stabilizer_dimension=sL,
            all_six_filtered_base_shears_in_group_fiber=True,actual_filtered_group_fiber_dimension=group_fiber,
            fixed_P0_higher_necessary_dimension_upper=len(K1)+len(K2),whole_projective_image_dimension_upper=upper,
            all_full39_actualD_linear_terms=linD,full_actualD_quadratic_terms_in_universal_table=True,
            all_full39_Jacobi_quadratic_terms_in_universal_table=True,all_bad_urow_coefficients_vanish=u_factor,
            actual_u_hyperplane_factor_for_every_full_solution=u_factor,rank4_irreducible_source_empty=u_factor,G24_closed=u_factor or upper<=24)
        print(json.dumps(dict(type=name,field=field,k1=len(K1),k2=len(K2),Levi_stabilizer=sL,group_fiber=group_fiber,pdim=pdim,image_upper=upper,actual_u_factor=u_factor,G24_closed=out['G24_closed'])),flush=True)
        return out,(d1,d2,dg),K1+K2
    records=[];generic_data=None
    for name,P,target in models:
        out,dets,ker=process(name,P,target,int(name.startswith('moving')),field='QQ(t)'if name.startswith('moving')else'QQ')
        records.append(out)
        if name.startswith('moving'):generic_data=(out,dets,ker,operators(P))
    symbol=sy.Symbol('t');factors={}
    for a in list(generic_data[1])+[v for row in generic_data[2]for v in row]:
        r=RF(a)
        for polynomial in (r.p,r.q):
            if len(polynomial)>1:
                expr=sum(sy.Rational(v.numerator,v.denominator)*symbol**i for i,v in enumerate(polynomial))
                for fac,mult in sy.factor_list(expr,symbol)[1]:
                    coeff=tuple(F(v)for v in reversed(sy.Poly(fac,symbol).monic().all_coeffs()))
                    if coeff!=(F(),F(1)):factors[coeff]=str(fac)
    walls=[]
    for modulus,text in sorted(factors.items()):
        if len(modulus)==2:
            value=-modulus[0]/modulus[1];convert=lambda v:sum((c*value**i for i,c in enumerate(RF(v).p)),F())/sum((c*value**i for i,c in enumerate(RF(v).q)),F())
            label='QQ';field_data=dict(t=str(value),minimal_polynomial=list(map(str,modulus)))
        else:
            AF.mod=modulus;convert=AF;label='QQ[t]/('+text+')';field_data=dict(minimal_polynomial=list(map(str,modulus)),all_complex_embeddings_covered=True)
        P=q.antisymmetric({ij:{e:convert(c)for e,c in p.items()if convert(c)}for ij,p in moving.items()if ij[0]<ij[1]})
        target={e:convert(c)for e,c in target_moving.items()if convert(c)}
        prepared=tuple([[convert(v)for v in row]for row in M]for M in generic_data[3])
        out,_,_=process('moving_wall_'+text,P,target,0,field=label,prepared=prepared);out.update(field_data);walls.append(out)
    records[-1]['whole_nonzero_t_walls']=walls;records[-1]['generic_open_excluded_irreducible_factors']=[list(map(str,v))for v in factors]
    records[-1]['whole_nonzero_t_G24_closed']=records[-1]['G24_closed']and all(v['G24_closed']for v in walls)
    closed=all(c['G24_closed']for c in records[:3])and records[-1]['whole_nonzero_t_G24_closed']
    result=dict(task_id='G',all_exact_arithmetic_passed=True,whole_tangent_double_G24_closed=closed,cases=records,
        full30_weight1_basis=[list(map(str,q.vec(P)))for P in first],full9_weight2_basis=[list(map(str,q.vec(P)))for P in second],
        full350_row_labels=[dict(quad=ij,monomial=e)for ij,e in q.ROWS],universal_full39_quadratic_Jacobi_and_actualD=universal,
        source_geometry_status='Producer pureproof reduction; separate independent review required for the three fixed models. Arithmetic does not certify classification.',
        scope='Only reduced tangent conic+line m2 normal-double; three constant-root fixed whole sources plus one whole moving-root t!=0 family. Triangle/concurrent/nonreducedR/m3 remain separate.',
        next_if_not_closed='No further profile or automatic jet/GB. Preserve all equations and turn to structural Ore/source geometry.',
        inputs=job['exact_inputs_and_hashes'],producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_source_bytes=Path(__file__).read_text(),
        seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_exact_arithmetic_passed=True,whole_tangent_double_G24_closed=closed,seconds=result['seconds'])))


if __name__=='__main__':
    import sys
    if '--m2-double-tangent'in sys.argv:m2_double_tangent_main()
    elif '--triple-root-quintics'in sys.argv:triple_root_quintics_main()
    elif '--repeated-211'in sys.argv:repeated_211_main()
    elif '--concurrent'in sys.argv:concurrent_main()
    else:main()
