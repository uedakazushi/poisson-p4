#!/usr/bin/env python3
"""ONE N/R bivariate source4, exact actual-generator conductor kernel and CI4."""
import argparse
from fractions import Fraction as F
import gzip
import itertools
import json
from math import factorial
from pathlib import Path
import resource
import signal
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own
from ambient_B18_actual_D_tangent_cone_support_review_20261103_v63 import inverse,matrix
from ambient_W1_S61_intersection_recognition_review_20261103_v63 import rref
from ambient_B18_Noether_fiber_even_quadratic_review_20261103_v63 import exact_row_rank

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_DRED_CONDUCTOR_20261103_v63'
MON=[(a,n-a) for n in range(5) for a in range(n,-1,-1)]
MI={e:i for i,e in enumerate(MON)}


def plus(A,B,c=F(1)):
    out={e:v[:] for e,v in A.items()}
    for e,v in B.items():out[e]=[a+c*b for a,b in zip(out.get(e,[F(0)]*150),v)]
    return {e:v for e,v in out.items() if any(v)}


def ma(A,B,c=F(1)):
    out={e:[r[:] for r in M] for e,M in A.items()}
    for e,M in B.items():
        out[e]=[[a+c*b for a,b in zip(ar,br)] for ar,br in zip(out.get(e,[[F(0)]*5 for _ in range(5)]),M)]
    return {e:M for e,M in out.items() if any(any(r) for r in M)}


def mm(A,B):
    out={}
    for e,M in A.items():
        for f,N in B.items():
            g=(e[0]+f[0],e[1]+f[1])
            if sum(g)<=4:out=ma(out,{g:own.matrix_product(M,N)})
    return out


def act(A,P):
    out={}
    for e,M in A.items():
        for f,V in P.items():
            g=(e[0]+f[0],e[1]+f[1])
            if sum(g)<=4:out=plus(out,{g:own.lie(M,V)})
    return out


def exponential(A,P):
    out,term=P,P
    for j in range(1,5):term=act(A,term);out=plus(out,term,F(1,factorial(j)))
    return out


def poly_mul(A,B):
    out={}
    for e,a in A.items():
        for f,b in B.items():
            g=(e[0]+f[0],e[1]+f[1])
            if sum(g)<=4:out[g]=out.get(g,F(0))+a*b
    return {e:a for e,a in out.items() if a}


def c4mul(a,b):
    return [sum((a[j]*b[i-j] for j in range(i+1)),F(0)) for i in range(4)]


def Taction(v):return [F(0),*v[:3],F(0),*v[4:7]]


def tauaction(v,nu):return v[4:]+c4mul(nu,v[:4])


def dump(path,data):
    def serial(x):
        if isinstance(x,F):return str(x)
        if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [serial(v) for v in x]
        return x
    with path.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as gz:
            gz.write((json.dumps(serial(data),separators=(',',':'))+'\n').encode())


def check_physics(q):
    div=[[F(0)]*15 for _ in range(25)]
    for ip,(i,j) in enumerate(own.PAIRS):
        for im,(a,b) in enumerate(own.MONS):
            V=q[15*ip+im]
            for k,rem in [(a,b),(b,a)]:
                if k==i:div[5*j+rem]=[x+y for x,y in zip(div[5*j+rem],V)]
                if k==j:div[5*i+rem]=[x-y for x,y in zip(div[5*i+rem],V)]
    assert all(not any(V) for V in div)
    cubic=list(itertools.combinations_with_replacement(range(5),3));ci={m:i for i,m in enumerate(cubic)}
    jac=[[F(0)]*15 for _ in range(350)];cache={}
    def entry(i,j):return own.PI[tuple(sorted((i,j)))],1 if i<j else -1
    for ti,(a,b,c) in enumerate(itertools.combinations(range(5),3)):
        for i,j,k in [(a,b,c),(b,c,a),(c,a,b)]:
            jb,sb=entry(j,k)
            for l in range(5):
                if i==l:continue
                ia,sa=entry(i,l)
                for am,mon in enumerate(own.MONS):
                    qa=15*ia+am
                    if not any(q[qa]):continue
                    for bm,(u,v) in enumerate(own.MONS):
                        count=int(u==l)+int(v==l);qb=15*jb+bm
                        if not count or not any(q[qb]):continue
                        key=tuple(sorted((qa,qb)))
                        if key not in cache:
                            p=poly_mul({e:x for e,x in zip(MON,q[qa]) if x},{e:x for e,x in zip(MON,q[qb]) if x})
                            cache[key]=[p.get(e,F(0)) for e in MON]
                        rem=v if u==l else u
                        row=35*ti+ci[tuple(sorted((*mon,rem)))];factor=sa*sb*count
                        jac[row]=[x+factor*y for x,y in zip(jac[row],cache[key])]
    assert all(not any(V) for V in jac)


def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
    assert not out.exists();out.mkdir(parents=True)
    gates=[ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json']
    contact,N,firstgate,orbit=[own.accepted(p) for p in gates];assert contact['k_exact']==4
    oldpath=ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63/certificate.json.gz'
    redpath=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    curvepath=ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    firstpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    old,red,curve,first=[own.read(p) for p in [oldpath,redpath,curvepath,firstpath]]
    ca=red['families']['CA'];J=[[F(a) for j,a in enumerate(row) if j!=20] for row in ca['linear_J']]
    oi=red['orbit_GL_generator_indices'];K=matrix(ca['source_kernel_Lie_matrix'])
    assert F(ca['source_velocity'][2])==0
    t=own.symbols(['t'])['t'];rp=[own.Poly(own.evaluate(q,{'t':t})) for q in curve['physical_curve150']]
    target=[[p.terms.get((a,),F(0)) if b==0 else F(0) for a,b in MON] for p in rp]
    P0=[row[0] for row in target];LK0=own.lie(K,P0)
    fixed=[[row[j] for j in range(19,27)] for row in J]
    regular=old['regular27_physical_rows'];fold=old['fold_physical_row'];L=matrix(old['regular27_inverse'])
    assert own.matrix_product(L,[J[i] for i in regular])==[[F(i==j) for j in range(27)] for i in range(27)]
    assert own.matrix_product([J[i] for i in regular],L)==[[F(i==j) for j in range(27)] for i in range(27)]
    def physical(source):
        P={(0,0):P0[:],(0,1):[-a for a in LK0]};M={}
        for n,e in enumerate(MON):
            fv=[row[n] for row in source[19:]]
            if any(fv):P=plus(P,{e:own.matvec(fixed,fv)})
            A=[[F(0)]*5 for _ in range(5)]
            for row,ij in zip(source[:19],oi):A[ij//5][ij%5]+=row[n]
            if any(any(row) for row in A):M[e]=A
        assert (0,0) not in M
        log,power={},M
        for j in range(1,5):log=ma(log,power,F((-1)**(j+1),j));power=mm(power,M)
        Q=exponential(log,exponential({(0,1):K},P))
        return [[Q.get(e,[F(0)]*150)[i] for e in MON] for i in range(150)]
    source=[[F(0)]*15 for _ in range(27)]
    for degree in range(1,5):
        before=physical(source)
        for n,e in enumerate(MON):
            if sum(e)==degree:
                correction=own.matvec(L,[target[i][n]-before[i][n] for i in regular])
                for row,a in zip(source,correction):row[n]+=a
        dump(out/f'stage{degree}.json.gz',dict(degree=degree,source27=source))
    q=physical(source)
    assert all(q[i]==target[i] for i in regular)
    assert [row[:10] for row in source]==matrix(old['source27_IFT_coefficients'])
    assert [row[:10] for row in q]==matrix(old['full150_physical_source_coefficients_afterIFT'])
    check_physics(q)
    f={e:q[fold][n]-target[fold][n] for n,e in enumerate(MON) if q[fold][n]!=target[fold][n]}
    assert {e:a for e,a in f.items() if sum(e)<4}=={(1,1):F(48),(0,2):F(-3)}
    u02=f.get((0,4),F(0));u11=f.get((1,3),F(0))+16*u02;u20=f.get((2,2),F(0))+16*u11
    w31=-(f.get((3,1),F(0))+16*u20)/3;w40=-f.get((4,0),F(0))/3
    W={(0,2):F(1),(1,1):F(-16)}
    if w31:W[3,1]=w31
    if w40:W[4,0]=w40
    U={(0,0):F(-3)}
    for e,a in [((0,2),u02),((1,1),u11),((2,0),u20)]:
        if a:U[e]=a
    assert poly_mul(U,W)==f
    def reduce_term(a,b,v,result):
        if a+b>4 or not v:return
        if b<2:result[a,b]=result.get((a,b),F(0))+v
        else:
            for (i,j),c in W.items():
                if (i,j)!=(0,2):reduce_term(a+i,b-2+j,-v*c,result)
    aa,bb=[],[]
    for row in q:
        rem={}
        for e,v in zip(MON,row):reduce_term(*e,v,rem)
        aa.append([rem.get((i,0),F(0)) for i in range(4)])
        bb.append([rem.get((i,1),F(0)) for i in range(4)])
    center=[F(0),F(8),F(0),-w31/2]
    nu=[F(0),F(0),F(64),F(0)]
    traceaa=[[a[i]+c4mul(center,b)[i] for i in range(4)] for a,b in zip(aa,bb)]
    d=[[target[j][MI[i,0]]-traceaa[j][i] for i in range(4)] for j in range(150)]
    assert all(b[0]==0 for b in bb)
    assert [b[:3] for b in bb]==matrix(old['full150_remainders_odd_b'])
    assert traceaa==matrix(old['full150_trace_zero_even_a'])
    nrows=curve['N37_coordinate_rows'];odd=nrows[28:];m=odd.index(61)
    assert set(nrows[:28])==set(regular+[fold])
    sigma=[bb[i] for i in odd];eta=[d[i] for i in odd]
    assert sigma[m][1]==46 and eta[m][2]==232
    for i in range(9):
        for j in range(9):
            assert c4mul(sigma[i],eta[j])==c4mul(sigma[j],eta[i])
            assert c4mul(eta[i],eta[j])==c4mul(nu,c4mul(sigma[i],sigma[j]))
    # Save the complete source computation BEFORE any conductor postprocessing.
    dump(out/'source4.json.gz',dict(monomials=MON,regular27=regular,fold=fold,source27=source,
        all150_source_coefficients=q,all150_R_target_coefficients=target,
        fold_equation=[[e,a] for e,a in f.items()],Weierstrass_W=[[e,a] for e,a in W.items()],Weierstrass_U=[[e,a] for e,a in U.items()],
        all150_odd_b_throughT3=bb,all150_trace_even_a_throughT3=traceaa,
        all150_phi_odd_throughT3=d,actual_odd9_physical_rows=odd,phi_even9=sigma,phi_odd9=eta,phi_nu=nu,
        full_source_recipe=dict(point150=P0,fixed8_matrix150_by8=fixed,GL19_indices=oi,kernel_GL=K,compensation150=[-a for a in LK0],
            regular27_rows=regular,fold_row=fold,inverse27=L,
            definition='Full rational GL(I+M)*exp(pK) source; unique27implicit solve and distinguished monic quadratic; tau=p-trace(p)/2; bi=tau-coefficients of physical odd9.')))
    gens=[sigma[i]+eta[i] for i in range(9)]+[eta[i]+c4mul(nu,sigma[i]) for i in range(9)]
    columns=[]
    for g in gens:
        v=g
        for _ in range(4):columns.append(v);v=Taction(v)
    M=list(map(list,zip(*columns)));RR,pivots=rref(M);assert len(pivots)==4
    free=[j for j in range(72) if j not in pivots];kernel=[]
    for j in free:
        v=[F(0)]*72;v[j]=1
        for i,p in enumerate(pivots):v[p]=-RR[i][j]
        assert own.matvec(M,v)==[0]*8
        kernel.append(v)
    assert len(kernel)==68
    gm=gens[m];Tg=Taction(gm);taug=tauaction(gm,nu);T2g=Taction(Tg)
    cyclic=list(map(list,zip(gm,Tg,taug,T2g)))
    _,unitrows=rref(list(map(list,zip(*cyclic))));assert len(unitrows)==4
    ic,cdet=inverse([cyclic[i] for i in unitrows])
    h=[]
    for g in gens[:9]:
        v=own.matvec(ic,[g[i] for i in unitrows]);assert own.matvec(cyclic,v)==g;h.append(v)
    assert h[m]==[1,0,0,0]
    mu=F(116,23);assert mu*mu-64==-F(20400,529)
    assert [a-mu*b for a,b in zip(Taction(taug),T2g)]==[0]*8
    assert [a-64*b for a,b in zip(tauaction(taug,nu),T2g)]==[0]*8
    assert Taction(T2g)==[0]*8
    for i,g in enumerate(gens[:9]):
        hv=h[i];z=[sum((hv[j]*col[k] for j,col in enumerate([taug,Taction(taug),tauaction(taug,nu),Taction(Taction(taug))])),F(0)) for k in range(8)]
        assert z==gens[9+i]
    ranks=[];current=M
    for j in range(4):
        ranks.append(exact_row_rank(current));current=[row[:] for row in list(map(list,zip(*(Taction(v) for v in zip(*current)))))]
    assert ranks==[4,2,1,0]
    # Actual even-base coordinate/gauge binding; group equations contain the full curve.
    Rrows=red['common28_linear_target_rows'];OR=[row[:19] for i,row in enumerate(J) if i in Rrows]
    OR=[J[i][:19] for i in Rrows]
    _,gpiv=rref(list(map(list,zip(*OR))));assert len(gpiv)==19
    comp=[i for i in range(28) if i not in gpiv]
    Q=[[F(i==j) for j in comp] for i in range(28)]
    ib,baseunit=inverse([a+b for a,b in zip(OR,Q)]);projection=ib[19:]
    vcurve=[target[i][MI[1,0]] for i in Rrows]
    projected=own.matvec(projection,vcurve);pick=next(i for i,a in enumerate(projected) if a)
    ell=[a/projected[pick] for a in projection[pick]]
    assert own.matvec([ell],vcurve)==[1] and own.matrix_product([ell],OR)==[[0]*19]
    rc=[[target[i][MI[j,0]] for j in range(1,4)] for i in Rrows]
    well=[sum((a*row[j] for a,row in zip(ell,rc)),F(0)) for j in range(3)];assert well[0]==1
    invcoef=[F(1),-well[1],2*well[1]**2-well[2]]
    aug=[row+[F(first['N_full_firstjet_data']['W'][i][0])] for i,row in enumerate(J)]
    ainv,_=inverse([aug[i] for i in Rrows]);nugrad=ainv[-1]
    assert own.matvec([nugrad],vcurve)==[0] and own.matrix_product([nugrad],OR)==[[0]*19]
    coordinate_rows=[nugrad,ell];coords=[]
    for j in gpiv+comp:
        row=[F(i==j) for i in range(28)]
        if exact_row_rank(coordinate_rows+[row])>len(coordinate_rows):coordinate_rows.append(row);coords.append(j)
    assert len(coords)==26 and coords[:19]==gpiv
    _,coordinateunit=inverse(coordinate_rows)
    # Normalizer ann of the cyclic image: 26 actual even curve constraints and2 quadrics.
    # Quotienting19 exact orbit gauges leaves7+2=9 CI generators in the slice.
    dump(out/'conductor-kernel.json.gz',dict(
        actual_conductor_R_generators=[dict(kind=kind,physical_row=i) for kind in ['bi','tau_bi'] for i in odd],
        phi_even9=sigma,phi_odd9=eta,phi_nu=nu,C4_linear_map8_by72=M,
        QQ_kernel_basis68_by72=kernel,QQ_pivot_columns=pivots,
        cyclic_generator_physical_row=61,cyclic_basis_order=['1','u','tau','u^2'],cyclic_basis8_by4=cyclic,cyclic_unit_rows=unitrows,cyclic_minor=cdet,
        all9_coefficients_hi_in_cyclic_basis=h,
        annihilator_quadrics=[dict(monomials=['u*tau','u^2'],coefficients=[1,-mu]),dict(monomials=['tau^2','u^2'],coefficients=[1,-64])],
        base_R28_physical_rows=Rrows,base_Rcurve_coefficients_T_T2_T3=rc,
        uniformizer_linear_ell=ell,ell_of_Rcurve_coefficients=well,uniformizer_inverse_through3=invcoef,
        uniformizer_full_definition='u = formal compositional inverse of exact ell(Rcurve(T)) applied to ell(R); saved exact Rcurve is immutable.',
        source_nu_full_definition='nu=c0+c1^2/4 from the unique full Weierstrass quadratic p^2-c1(R)p-c0(R).',
        source_nu_gradient=nugrad,
        group19_gauge_physical_R_indices=gpiv,group19_gauge_physical_rows=[Rrows[i] for i in gpiv],
        group19_gauge_definition='R_j - exact Rcurve_j(u), no truncation; Jacobian on orbit is OR[gpiv].',
        residual7_even_normal_physical_indices=coords[19:],residual7_even_normal_physical_rows=[Rrows[i] for i in coords[19:]],
        residual7_annihilator_constraints=[dict(physical_R_row=Rrows[i],coefficients_u_u2=rc[i][:2]) for i in coords[19:]],
        actual_base_coordinate_chart_unit=coordinateunit,
        auxiliary_full_R28_curve_kernel_definition='For any j with ell_j nonzero: (u_3^4, R_i-Rcurve_i(u_3) mod u_3^4 for i!=j), with u_3=ellR - b(ellR)^2+(2b^2-c)(ellR)^3.',
        full_R_module_kernel_definition='kerPsi=(kernel_rho)*IB_N + sum_R 68 lifted QQ kernel vectors, T^k lifted by u_3^k; then quotient exact19 group gauges and restore G19.',
        exact_slice_BN_ideal_J_generators='8 differences bi-hi(u,tau)*b61 (i!=61), plus9 annihilator generators timesb61; actual bi are bound by full source recipe.',
        original_reduced_conductor_in_normalization='J_slice direct_sum (T^4), then formal G19 product; NOT I intersect Ann(g61).',
        image_R_module_Tpower_ranks=ranks,image_R_module_Hilbert=[2,1,1],image_BN_cyclic_Hilbert=[1,2,1]))
    inputs=[Path(__file__),Path(own.__file__),ROOT/'code/ambient_B18_actual_D_tangent_cone_support_review_20261103_v63.py',
        ROOT/'code/ambient_W1_S61_intersection_recognition_review_20261103_v63.py',ROOT/'code/ambient_B18_Noether_fiber_even_quadratic_review_20261103_v63.py',
        ROOT/f'research/notes/{STEM}.md',oldpath,redpath,curvepath,firstpath,*gates]
    files=[out/'source4.json.gz',out/'conductor-kernel.json.gz',*[out/f'stage{i}.json.gz' for i in range(1,5)]]
    result=dict(id=STEM,all_passed=True,
        ONE_targeted_N_R_bivariate_degree4_source_with_no_larger_retry=True,
        genuine_fixed_N19_full_rational_source_recipe_bound=True,
        all150_source4_and_original_divergence25_Jacobi350_verified=True,
        full27_IFT_with_two_inverse_products_and_all_old_degree3_coefficients_verified=True,
        exact_Weierstrass_unit_through4_and_full_physical_odd_coefficients_T3=True,
        actual18_conductor_generators_complete_from_actual_N37_and_full_source_trace_recipe=True,
        psi_map_shape=[8,72],psi_QQ_rank=4,psi_QQ_kernel_dimension=68,
        cyclic_BN_image_generator_physical_b_row=61,cyclic_image_length=4,
        image_R_module_Hilbert_function=[2,1,1],image_BN_annihilator_quotient_Hilbert_function=[1,2,1],
        explicit_actual_BN_slice_annihilator_CI_generator_count=9,
        explicit_J_slice_BN_ideal_generator_upper_bound=17,
        original_Dred_conductor_in_normalization='(J_slice, T^4), restored by actual G19 product',
        actual_even_G19_slice_contains_full_Rcurve_and_has_unit_orbit_Jacobian=True,
        no_abstract_normalform_to_physical_guess_no_image_ideal_basechange_or_tau_value_assumed=True,
        full_original_conductor_ideal_given_by_actual_generators_finite_CI_recipe=True,
        full_nonreduced_D_depth_nilpotency_or_higher_source_jets_claimed=False,
        scope='Formal generic B18 reduced intersection; actual even-base transverse slice plus G19 restoration. Full expanded28-variable coefficients are not serialized.',
        input_sha256={str(p.relative_to(ROOT)):own.sha(p) for p in inputs},
        evidence_sha256={p.name:own.sha(p) for p in files},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','psi_QQ_rank','cyclic_image_length','explicit_J_slice_BN_ideal_generator_upper_bound','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    main(parser.parse_args().output)
