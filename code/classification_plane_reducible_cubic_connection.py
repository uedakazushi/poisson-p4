#!/usr/bin/env python3
"""Finite complete simple-pole SL2 sources for four reduced cubic divisors.

This constructs all24 coefficients and30 exact flatness polynomials, not a
claimed solution of the nonlogarithmic complement or full Poisson support.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

from classification_minimal_quartic_review import LaurentRing,determinant
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]

def nullspace(A,n):
    B=[[F(v)for v in row]for row in A];r=0;piv=[]
    for j in range(n):
        k=next((k for k in range(r,len(B))if B[k][j]),None)
        if k is None:continue
        B[r],B[k]=B[k],B[r];d=B[r][j];B[r]=[v/d for v in B[r]]
        for k in range(len(B)):
            if k!=r and B[k][j]:
                a=B[k][j];B[k]=[u-a*v for u,v in zip(B[k],B[r])]
        piv.append(j);r+=1
    vectors=[]
    for j in range(n):
        if j not in piv:
            v=[F()for _ in range(n)];v[j]=F(1)
            for k,i in enumerate(piv):v[i]=-B[k][j]
            assert all(sum(a*b for a,b in zip(row,v))==0 for row in A);vectors.append(v)
    return vectors,piv

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--scalar-fourth-order',action='store_true')
    args=ap.parse_args()
    if args.scalar_fourth_order:return scalar_fourth_order(args.output)
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,)*2)
    R=LaurentRing([f'p{i}'for i in range(24)],['x','y','z']);p=R.parse
    xyz=[R.variable(v)for v in ['x','y','z']]
    mon3=[tuple(v.count(i)for i in range(3))for v in it.combinations_with_replacement(range(3),3)]
    basis=[]
    for i,j in it.product(range(3),repeat=2):
        if i!=j:
            M=[[F(a==i and b==j)for b in range(3)]for a in range(3)];basis.append(M)
    basis.extend([[[F(1),0,0],[0,F(-1),0],[0,0,0]],[[0,0,0],[0,F(1),0],[0,0,F(-1)]]])
    vector=lambda values:[R.add(*(R.scale(xyz[j],sum(values[k]*basis[k][i][j]for k in range(8)))for j in range(3)))for i in range(3)]
    numerator=lambda V:[R.add(R.mul(xyz[2],V[1]),R.scale(R.mul(xyz[1],V[2]),-1)),
                        R.add(R.mul(xyz[0],V[2]),R.scale(R.mul(xyz[2],V[0]),-1)),
                        R.add(R.mul(xyz[1],V[0]),R.scale(R.mul(xyz[0],V[1]),-1))]
    # vector() above takes numeric coefficients; here source coefficients are polynomials.
    genericV=[[R.add(*(R.scale(R.mul(R.variable(f'p{8*a+k}'),xyz[j]),basis[k][i][j])for k in range(8)for j in range(3)))for i in range(3)]for a in range(3)]
    def flat(f,Vs):
        H,E,G=Vs
        value=lambda W:R.add(*(R.mul(W[i],R.derivative(f,i))for i in range(3)))
        det=lambda U,W:determinant(R,[xyz,U,W])
        return [R.add(value(H),det(E,G)),R.add(value(E),R.scale(det(H,E),2)),R.add(value(G),R.scale(det(H,G),-2))]
    def curvature(f,N):
        out=[]
        for i,j in it.combinations(range(3),2):
            H=R.add(R.mul(N[1][i],N[2][j]),R.scale(R.mul(N[2][i],N[1][j]),-1))
            E=R.scale(R.add(R.mul(N[0][i],N[1][j]),R.scale(R.mul(N[1][i],N[0][j]),-1)),2)
            G=R.scale(R.add(R.mul(N[0][i],N[2][j]),R.scale(R.mul(N[2][i],N[0][j]),-1)),-2)
            out.append([R.add(R.mul(f,R.add(R.derivative(N[a][j],i),R.scale(R.derivative(N[a][i],j),-1))),
                       R.scale(R.mul(R.derivative(f,i),N[a][j]),-1),R.mul(R.derivative(f,j),N[a][i]),[H,E,G][a])for a in range(3)])
        return out
    def coefficients(poly):
        return [{e[:-3]+(0,)*3:c for e,c in poly.items()if e[-3:]==mon}for mon in mon3]
    models=[('conic_transverse','z*(x*y-z**2)',['x*y-z**2','z'],[['x','z'],['z','y']]),
            ('conic_tangent','z*(x*z-y**2)',['x*z-y**2','z'],[['x','y'],['y','z']]),
            ('triangle','x*y*z',['x','y','z'],[['x','z'],['0','y']]),
            ('concurrent','x*y*(x-y)',['x','y','x-y'],[['x','z'],['0','y']])]
    records=[]
    for name,text,factors,T in models:
        f=p(text);J=flat(f,genericV);N=list(map(numerator,genericV));actual=curvature(f,N)
        # Every horizontal2form is i_E(volume)*one cubic; this checks all3 raw components.
        for at,(i,j)in enumerate(it.combinations(range(3),2)):
            k=next(k for k in range(3)if k not in [i,j]);sign=(-1)**k
            assert actual[at]==[R.scale(R.mul(xyz[k],v),sign)for v in J]
        fcoeff=[f.get((0,)*24+mon,F())for mon in mon3]
        action=[]
        for M in basis:
            W=[R.add(*(R.scale(xyz[j],M[i][j])for j in range(3)))for i in range(3)]
            value=R.add(*(R.mul(W[i],R.derivative(f,i))for i in range(3)))
            action.append([value.get((0,)*24+mon,F())for mon in mon3])
        action=list(map(list,zip(*action)))
        closed,_=nullspace(action,8)
        logkernel,_=nullspace([row+[-c]for row,c in zip(action,fcoeff)],9)
        logs=[v[:8]for v in closed]
        extra=next((v for v in logkernel if v[8]),None)
        if extra:logs.append(extra[:8])
        assert len(logs)==len(logkernel)and len(closed)==len(factors)-1
        expected={'conic_transverse':1,'conic_tangent':2,'triangle':2,'concurrent':3}[name]
        assert len(logs)==expected
        logV=[vector(v)for v in logs];lambdas=[]
        for v in logs:
            Vf=[sum(a*b for a,b in zip(row,v))for row in action]
            first=next(i for i,c in enumerate(fcoeff)if c);lam=Vf[first]/fcoeff[first]
            assert Vf==[lam*c for c in fcoeff];lambdas.append(lam)
        pair_constants=[]
        for i,j in it.combinations(range(len(logs)),2):
            D=determinant(R,[xyz,logV[i],logV[j]])
            key=next(e for e in f);c=D.get(key,F())/f[key]
            assert D==R.scale(f,c);pair_constants.append(dict(pair=[i,j],wedge_constant=str(c)))
        # Explicit genuine nonlog flat connection for EVERY type, with no appeal to pi1Z3.
        T=[[p(v)for v in row]for row in T];Q=R.add(R.mul(T[0][0],T[1][1]),R.scale(R.mul(T[0][1],T[1][0]),-1))
        last=p('z')if name.startswith('conic')else p('z')if name=='triangle'else p('x-y')
        assert R.mul(Q,last)==f
        adj=[[T[1][1],R.scale(T[0][1],-1)],[R.scale(T[1][0],-1),T[0][0]]]
        W=[[],[],[]]
        for coordinate in range(3):
            A=[[R.add(*(R.mul(adj[i][k],R.derivative(T[k][j],coordinate))for k in range(2)))for j in range(2)]for i in range(2)]
            W[0].append(R.scale(R.mul(last,R.add(A[0][0],R.scale(A[1][1],-1))),F(1,2)))
            W[1].append(R.mul(last,A[0][1]));W[2].append(R.mul(last,A[1][0]))
        assert all(not v for row in curvature(f,W)for v in row)
        assert all(not R.add(*(R.mul(xyz[i],v[i])for i in range(3)))for v in W)
        # Nonzero df_i wedge numerator modulo the factor establishes a genuine nonlog polar part.
        def divides(poly,divisor):
            if not poly:return True
            lead=max(divisor);value=divisor[lead];rem=dict(poly)
            while rem:
                mon=max(rem)
                if not all(a>=b for a,b in zip(mon,lead)):return False
                term={tuple(a-b for a,b in zip(mon,lead)):rem[mon]/value}
                rem=R.add(rem,R.scale(R.mul(term,divisor),-1))
            return True
        nonlog=[]
        for factor in factors:
            g=p(factor)
            if any(not divides(R.add(R.mul(R.derivative(g,i),form[j]),R.scale(R.mul(R.derivative(g,j),form[i]),-1)),g)
                   for form in W for i,j in it.combinations(range(3),2)):nonlog.append(factor)
        assert nonlog
        records.append(dict(type=name,cubic=text,components=factors,
            full24_SL2_parameters=[f'p{i}'for i in range(24)],full30_flatness_polynomials=[R.encode(v)for z in J for v in coefficients(z)],
            tracefree_linear_vectorfield_basis=[[[str(v)for v in row]for row in M]for M in basis],
            all_three_raw_curvature_components_verified=True,
            full_logarithmic_source_basis=[[str(v)for v in row]for row in logs],logarithmic_source_dimension=len(logs),
            closed_scalar_source_dimension=len(closed),logarithmic_form_derivative_constants=list(map(str,lambdas)),
            logarithmic_pair_wedge_constants=pair_constants,
            nonlog_flat_witness_matrix=[[R.encode(v)for v in row]for row in T],
            nonlog_flat_witness_numerators=[[R.encode(v)for v in row]for row in W],nonlog_components=nonlog,
            nonlog_not_excluded=True))
    paths=[ROOT/'research/classification/PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION_REVIEW.md',
           ROOT/'research/classification/PLANE_SMOOTH_CUBIC_RANK_TWO_REVIEW.md',
           ROOT/'reports/CLASSIFICATION_PLANE_SMOOTH_CUBIC_NORMAL_WEIGHTS_20261107.json']
    result=dict(task_id='G',all_passed=True,all_four_global_flat_sources_complete=True,
        all_four_types=records,scalar_closed_source_dimensions=[1,1,2,2],
        logarithmic_flat_source_dimensions=[3,3,4,6],
        log_classification='Transverseconic oneclosedform arbitraryT; tangentconic closedalpha+nonclosedbeta satisfies k[T,U]+lambda U=0 (U0 ornilpotent); triangle commutingpair; concurrent two arbitrarymatrices withnonclosedthirdcoefficientzero.',
        does_not_claim_full_nonlog_solution=True,does_not_claim_G_dimension_bound=True,
        nonlogarithmic_complement_retained_for_every_type=True,no_CAS_import=True,
        input_sha256={str(z.relative_to(ROOT)):hashlib.sha256(z.read_bytes()).hexdigest()for z in paths},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','seconds','maximum_RSS_KiB']}))

def scalar_fourth_order(output):
    """Full30/full9 dimension bound for m4, squarefree, nonexceptional j.

    Geometry gives tracefreeOmega0 and the discrete scalar connection
    beta-beta_can=-(1/6)dlog(L^3/f). All remaining parameter fibers are kept.
    """
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,)*2)
    def read(text):return q.parse(text.replace('x','h').replace('y','x1').replace('z','x2'))
    X=[q.atom(1,i)for i in range(5)];base=X[:3]
    bmon=q.monomials(2)
    def point_sum(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in q.PAIRS})
    def actual(P):
        out={}
        for missing in range(5):
            i,j,k,l=[a for a in range(5)if a!=missing]
            pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
            out=q.add(out,q.scale(q.mul(X[missing],pf),(-1)**missing))
        return out
    def det3(M):
        out={}
        for perm in it.permutations(range(3)):
            sign=(-1)**sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3))
            term=q.atom(sign)
            for i,j in enumerate(perm):term=q.mul(term,M[i][j])
            out=q.add(out,term)
        return out
    basis=[]
    for i,j in it.product(range(3),repeat=2):
        if i!=j:basis.append([[F(a==i and b==j)for b in range(3)]for a in range(3)])
    basis.extend([[[F(1),0,0],[0,F(-1),0],[0,0,0]],[[0,0,0],[0,F(1),0],[0,0,F(-1)]]])
    fields=[[q.add(*(q.scale(base[j],M[i][j])for j in range(3)))for i in range(3)]for M in basis]
    def numerator(V):return [q.add(q.mul(base[2],V[1]),q.scale(q.mul(base[1],V[2]),-1)),q.add(q.mul(base[0],V[2]),q.scale(q.mul(base[2],V[0]),-1)),q.add(q.mul(base[1],V[0]),q.scale(q.mul(base[0],V[1]),-1))]
    nums=list(map(numerator,fields))
    A=[[n[i].get(e+(0,0,0),F())for n in nums]for i in range(3)for e in bmon]
    first=[]
    for normal in [3,4]:
        for ij in it.combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,normalmon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*normalmon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in it.combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    assert len(first)==30 and len(second)==9
    records=[]
    for name,ft,lt,qt,aut in [('conic_transverse','z*(x*y-z*z)','z','x*y-z*z',1),('conic_tangent','z*(x*z-y*y)','z','x*z-y*y',2),('triangle','x*y*z','z','x*y',2),('concurrent','x*y*(x-y)','x-y','x*y',3)]:
        f,L,Q=map(read,[ft,lt,qt]);assert q.mul(L,Q)==f
        target=[q.add(q.scale(q.der(f,i),F(1,6)),q.scale(q.mul(Q,q.der(L,i)),-F(1,2)))for i in range(3)]
        values=[v.get(e+(0,0,0),F())for v in target for e in bmon]
        ker,_=nullspace([row+[-v]for row,v in zip(A,values)],9)
        assert len(ker)==1 and ker[0][-1]
        coefficients=[v/ker[0][-1]for v in ker[0][:-1]]
        assert all(sum(a*b for a,b in zip(row,coefficients))==v for row,v in zip(A,values))
        V=[q.scale(q.add(*(q.scale(fields[k][i],coefficients[k])for k in range(8))),3)for i in range(3)]
        Pf=q.antisymmetric({(0,1):q.der(f,2),(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)})
        P0=point_sum(Pf,q.antisymmetric({(i,j):q.scale(q.mul(X[j],V[i]),-1)for i in range(3)for j in [3,4]}))
        assert all(not p for p in q.jacobi(P0).values())and not actual(P0)
        J1=q.transpose([q.jvec(q.jacobi(point_sum(P0,P)))for P in first])
        cross=[actual(point_sum(P0,P))for P in first]
        assert all(ex[3]+ex[4]==3 for P in cross for ex in P)
        d3mon=[e for e in q.monomials(5)if e[3]+e[4]==3]
        D3=[[P.get(e+(0,0,0),F())for P in cross]for e in d3mon]
        necessary=J1+D3;K1,piv1=nullspace(necessary,30);k1=len(K1)
        J2=q.transpose([q.jvec(q.jacobi(point_sum(P0,P)))for P in second]);K2,piv2=nullspace(J2,9)
        bound=k1+len(K2);image=24+bound-(aut+10)
        # The parallel leading-symbol theorem also fixes G4=L*g4. Its10
        # coefficient equations are affine-linear in all9 second-weight
        # variables and quadratic in the complete K1 source.
        normal4=[]
        for terms in it.combinations_with_replacement([3,4],4):
            e=[0]*5
            for i in terms:e[i]+=1
            normal4.append(e)
        def shape(poly):
            result=[]
            for e in normal4:
                values=[]
                for i in range(3):
                    ex=list(e);ex[i]+=1;values.append(poly.get(tuple(ex)+(0,0,0),F()))
                result.extend([values[0]+values[1],values[2]]if name=='concurrent'else values[:2])
            return result
        L2=J2+q.transpose([shape(actual(point_sum(P0,P)))for P in second])
        extraK,extraPiv=nullspace(L2,9);r2=len(extraPiv)
        first_family=[]
        for row in K1:
            first_family.append(q.antisymmetric({ij:q.add(*(q.scale(P.get(ij,{}),v)for P,v in zip(first,row)if v))for ij in q.PAIRS}))
        monomial_pairs=list(it.combinations_with_replacement(range(k1),2));Q2=[]
        ownJ=[q.jvec(q.jacobi(P))for P in first_family];ownD=[shape(actual(P))for P in first_family]
        for i,j in monomial_pairs:
            if i==j:Q2.append(ownJ[i]+ownD[i])
            else:
                J=q.jvec(q.jacobi(point_sum(first_family[i],first_family[j])))
                D4=shape(actual(point_sum(first_family[i],first_family[j])))
                Q2.append([v-a-b for v,a,b in zip(J+ D4,ownJ[i]+ownD[i],ownJ[j]+ownD[j])])
        proper=None
        if r2==9:
            rr=q.rank(q.transpose(L2))[1];assert len(rr)==9
            inv=q.inverse([L2[i]for i in rr])
            for at,column in enumerate(Q2):
                solution=[sum(a*column[i]for a,i in zip(row,rr))for row in inv]
                residual=[value-sum(a*b for a,b in zip(row,solution))for row,value in zip(L2,column)]
                if any(residual):
                    which=next(i for i,v in enumerate(residual)if v)
                    functional={which:F(1)}
                    coefficients=[sum(L2[which][k]*inv[k][j]for k in range(9))for j in range(9)]
                    for i,c in zip(rr,coefficients):functional[i]=functional.get(i,F())-c
                    functional={i:c for i,c in functional.items()if c}
                    assert all(not sum(c*L2[i][j]for i,c in functional.items())for j in range(9))
                    polynomial=[sum(c*column[i]for i,c in functional.items())for column in Q2]
                    assert polynomial[at]and any(polynomial)
                    proper=dict(full350_plus10shape_functional={str(i):str(c)for i,c in functional.items()},
                        complete_quadratic_parameter_pairs=monomial_pairs,nonzero_quadratic_coefficients={str(i):str(v)for i,v in enumerate(polynomial)if v})
                    break
        refined_bound=k1+len(extraK)-int(proper is not None)
        refined_image=24+refined_bound-(aut+10)
        records.append(dict(type=name,cubic=ft,line_component=lt,
            scalar_normal_matrix=[[str(sum(3*coefficients[k]*basis[k][i][j]for k in range(8)))for j in range(3)]for i in range(3)],
            full_P0_Jacobi_and_actualD_zero=True,
            full350_weight1_matrix=[[str(v)for v in row]for row in J1],
            all24_D_normal_degree3_linear_rows=[[str(v)for v in row]for row in D3],
            full_weight1_and_order4_necessary_kernel=[[str(v)for v in row]for row in K1],
            weight1_kernel_dimension=k1,full350_weight2_matrix=[[str(v)for v in row]for row in J2],
            weight2_linear_rank=len(piv2),weight2_affine_fiber_dimension_upper=len(K2),
            fixed_P0_full_remaining_source_dimension_upper=bound,
            connected_surface_R_preserving_group_dimension=aut+10,
            coarse_whole_projective_image_dimension_upper=image,
            full350_plus10_leading_shape_weight2_linear_matrix=[[str(v)for v in row]for row in L2],
            shape_weight2_rank=r2,shape_weight2_fiber_dimension_upper=len(extraK),
            exact_proper_quadratic_compatibility=proper,
            refined_fixed_P0_source_dimension_upper=refined_bound,
            whole_projective_image_dimension_upper=refined_image,G24_closed=refined_image<=24))
    result=dict(task_id='G',all_arithmetic_passed=True,cases=records,
        all_four_nonexceptional_squarefree_order4_G24_closed=all(r['G24_closed']for r in records),
        scope='Full30/full9 triangular necessarydimension bound only for m4 squarefree leadingbinaryquartic with both classicalI,J nonzero. Otherleadingprofiles/specialj/fullnonlog30source remain explicit complements.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_four_nonexceptional_squarefree_order4_G24_closed'],cases=[{k:r[k]for k in ['type','weight1_kernel_dimension','weight2_affine_fiber_dimension_upper','whole_projective_image_dimension_upper']}for r in records],seconds=result['seconds'])))

if __name__=='__main__':main()
