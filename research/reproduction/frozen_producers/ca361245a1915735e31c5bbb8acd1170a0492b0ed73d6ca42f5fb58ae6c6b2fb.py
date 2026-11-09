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


if __name__=='__main__':
    import sys
    if '--repeated-211'in sys.argv:repeated_211_main()
    elif '--concurrent'in sys.argv:concurrent_main()
    else:main()
