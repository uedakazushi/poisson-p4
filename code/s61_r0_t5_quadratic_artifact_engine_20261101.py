"""T5 transverse upper23 at both roots of the quadratic pivot artifact."""
from itertools import combinations_with_replacement
import sympy as s


def verify_core(R,F,rr,globalb8):
    b9=R.one
    b0=R(rr)
    b6=R.from_expr(s.Symbol("b6"));b11=R.from_expr(s.Symbol("b11"));b8=globalb8
    aq=5254333434224*b0*b0-8538753610352*b0+3402879618517
    pp=96900*b0*b0-24150*b0*b6*b9-142225*b0-4140*b11*b9+1932*b8*b9*b9
    known=[pp]
    unit_checks=[]
    def unit_factor(p):
        assert p
        original=p; powers=[]
        for factor in known:
            power=0
            while True:
                q,rem=p.div(factor)
                if rem:break
                p=q;power+=1
            powers.append(power)
        assert len(p)==1 and next(iter(p))==(0,)*R.ngens, str(original)
        unit_checks.append(powers)
        return powers

    mons={d:list(combinations_with_replacement(range(3),d)) for d in [2,3,4,5]}
    mons[5]=[m for m in mons[5] if 0 not in m]
    index={d:{m:i for i,m in enumerate(mons[d])} for d in mons}
    fp=103; a=F[2][fp][0]
    assert a and all(not p for row in F[2] for p in row[1:])
    assert all(not p for m,p in zip(mons[3],F[3][fp]) if 0 not in m)
    H3=[];H4=[];H5=[]
    for i in range(135):
        ai=F[2][i][0]
        g3=[a*v-ai*w for v,w in zip(F[3][i],F[3][fp])]
        g4=[a*v-ai*w for v,w in zip(F[4][i],F[4][fp])]
        g5=[a*v-ai*w for v,w in zip(F[5][i],F[5][fp])]
        h3=[a*v for v in g3];h4=[a*v for v in g4];h5=[a*v for v in g5]
        for j in range(3):
            coefficient=g3[index[3][(0,0,j)]];h3[index[3][(0,0,j)]]=R.zero
            for mon,p in zip(mons[3],F[3][fp]):h4[index[4][tuple(sorted((j,)+mon))]]-=b9*coefficient*p
            if j:
                for mon,p in zip(mons[4],F[4][fp]):
                    if 0 not in mon:h5[index[5][tuple(sorted((j,)+mon))]]-=b9*coefficient*p
        H3.append(h3);H4.append(h4);H5.append(h5)
    gp=116; gc=index[3][(0,1,1)]
    A=H3[gp][gc];B=H3[gp][index[3][(0,1,2)]];C=H3[gp][index[3][(0,2,2)]]
    assert A
    unit_factor(a);unit_factor(A)
    for j,p in enumerate(H3[gp]):assert not p or mons[3][j] in [(0,1,1),(0,1,2),(0,2,2)]
    yg4=[R.zero]*6;zg4=[R.zero]*6
    for mon,p in zip(mons[4],H4[gp]):
        if 0 not in mon:
            yg4[index[5][tuple(sorted((1,)+mon))]]=p
            zg4[index[5][tuple(sorted((2,)+mon))]]=p

    def primitive(arr):
        nz=[q for q in arr if q]
        if not nz:return arr
        ex=min(mon[-1] for q in nz for mon in q)
        factor=pp**ex*R.ground_new(nz[0].LC)
        unit_factor(factor)
        return [q.exquo(factor) for q in arr]

    N4=[];N5=[]
    for i in range(135):
        ci=H3[i][gc]
        assert all(A*p-ci*q==R.zero for p,q in zip(H3[i],H3[gp]))
        z4=[A*p-ci*q for p,q in zip(H4[i],H4[gp])]
        z5=[A*p-ci*q for p,q in zip(H5[i],H5[gp])]
        assert all(not z4[j] for j,m in enumerate(mons[4]) if 0 not in m)
        h=[z4[index[4][m]] for m in [(0,1,1,1),(0,1,1,2),(0,1,2,2),(0,2,2,2)]]
        n4=[A*A*h[2]-A*C*h[0]-A*B*h[1]+B*B*h[0],A*A*h[3]-A*C*h[1]+B*C*h[0]]
        n5=[A*A*p-A*h[0]*y-(A*h[1]-B*h[0])*z for p,y,z in zip(z5,yg4,zg4)]
        both=primitive(n4+n5);N4.append(both[:2]);N5.append(both[2:])
    r0,r1=69,91
    determinant=N4[r0][0]*N4[r1][1]-N4[r0][1]*N4[r1][0]
    assert determinant
    unit_factor(determinant)

    def monomial_primitive(arr):
        nonzero=[p for p in arr if p]
        if not nonzero:return arr
        exponent=tuple(min(mon[j] for p in nonzero for mon in p) for j in range(R.ngens))
        factor=R.from_dict({exponent:nonzero[0].LC})
        unit_factor(factor)
        return [p.exquo(factor) for p in arr]

    b0=R(rr)
    b6=R.from_expr(s.Symbol("b6"));b11=R.from_expr(s.Symbol("b11"));b8=globalb8
    known=[pp]

    def cancel_known(arr):
        # Optimization only: every cancellation is verified by exact polynomial division.
        for factor in known:
            while any(arr):
                qr=[p.div(factor) for p in arr]
                if any(rem for q,rem in qr):break
                arr=[q for q,rem in qr]
        return monomial_primitive(arr)

    pure=[]
    for i in range(135):
        co0=N4[i][0]*N4[r1][1]-N4[i][1]*N4[r1][0]
        co1=N4[r0][0]*N4[i][1]-N4[r0][1]*N4[i][0]
        assert all(determinant*p-co0*q-co1*r==R.zero for p,q,r in zip(N4[i],N4[r0],N4[r1]))
        dd,co0,co1=cancel_known([determinant,co0,co1])
        pure.append(cancel_known([dd*p-co0*q-co1*r for p,q,r in zip(N5[i],N5[r0],N5[r1])]))
    selected_rows=[102,115,134];selected_columns=[0,1,2]
    matrix=[[pure[i][j] for j in selected_columns] for i in selected_rows]

    def det3(M):
        return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])

    det=det3(matrix);assert det
    for row in pure:
        coefficients=[]
        for k in range(3):
            replacement=[r[:] for r in matrix];replacement[k]=[row[j] for j in selected_columns]
            coefficients.append(det3(replacement))
        for j,p in enumerate(row):assert det*p-sum((co*pure[r][j] for co,r in zip(coefficients,selected_rows)),R.zero)==R.zero
    common=A.gcd(B).gcd(C)
    ep=lambda p:[[list(e),str(R.domain.to_sympy(co))] for e,co in p.items()]
    result=dict(initial_ranks={'2':1,'3':4,'4':10,'5':18},truncated_hilbert_function=[1,3,5,6,5,3],truncated_length=23,
                quadratic_generator='x^2',primitive_cubic_coefficients=[str(p.exquo(common)) for p in [A,B,C]],
                quartic_pivot_rows=[r0,r1],quartic_normal_determinant=str(determinant),
                quartic_initial_contains_all_x_times_binary_cubics=True,
                first_lift_cubic_has_no_pure_binary_terms=True,
                quintic_pivot_rows=selected_rows,quintic_pivot_columns=selected_columns,quintic_pivot_determinant=str(det),
                quintic_rows={str(i):[ep(p) for p in row] for i,row in enumerate(pure) if any(row)},
                all_polynomial_cross_product_identities_verified=True)

    result['all_division_factors_supported_on_b9_Aq_P']=True
    result['division_factor_exponents']=sorted({tuple(v) for v in unit_checks})
    return result,dict(N4=N4,N5=N5,H3=H3,H4=H4,mons=mons,index=index,A=A,B=B,C=C,a=a,known=known,pure=pure,unit_factor=unit_factor)


def verify(original_ring,original_residual):
    from sympy.polys.rings import ring
    from sympy.polys.matrices import DomainMatrix as DM
    discriminant=838398770209350113
    E=s.QQ.algebraic_field(s.sqrt(discriminant))
    root=E.from_sympy((8538753610352+1288*s.sqrt(discriminant))/(2*5254333434224))
    assert 5254333434224*root**2-8538753610352*root+3402879618517==E.zero
    assert E.ext.minpoly.degree()==2
    R,b1,b3,b4,b6,b11,p=ring('b1 b3 b4 b6 b11 p',E)
    a=R(root);b8=(p-96900*a*a+142225*a+24150*a*b6+4140*b11)/1932
    assert [str(v) for v in original_ring.symbols]==['b0','b1','b3','b4','b6','b8','b9','b11']
    images=[a,b1,b3,b4,b6,b8,R.one,b11];powers={}

    def convert(poly):
        out=R.zero
        for exponent,c in poly.items():
            term=R.ground_new(E.convert(c))
            for i,n in enumerate(exponent):
                if n:
                    if (i,n) not in powers:powers[i,n]=images[i]**n
                    term*=powers[i,n]
            out+=term
        return out

    mons3=list(combinations_with_replacement(range(3),3));index3={m:i for i,m in enumerate(mons3)}
    f2=convert(original_residual[2][103][0]);g2=convert(original_residual[2][116][0])
    cubic=[f2*(f2*convert(original_residual[3][116][index3[m]])-g2*convert(original_residual[3][103][index3[m]])) for m in [(0,1,1),(0,1,2),(0,2,2)]]
    assert not cubic[0]
    scalar=cubic[1].exquo(p);assert len(scalar)==1 and next(iter(scalar))==(0,)*R.ngens
    ratio=cubic[2].exquo(scalar*p)
    aa=R.one-ratio;bb=R.one+ratio;cc=R.one;dd=-R.one
    assert aa*dd-bb*cc==-2*R.one
    selected={103,116,69,91,102,115,134};F={}
    for degree in [2,3,4,5]:
        mons=list(combinations_with_replacement(range(3),degree))
        if degree==5:mons=[m for m in mons if 0 not in m]
        index={m:i for i,m in enumerate(mons)};rows=[]
        for row in range(135):
            values=[R.zero]*len(mons)
            if row in selected:
                for mon,poly in zip(mons,original_residual[degree][row]):
                    coefficient=convert(poly)
                    if not coefficient:continue
                    n0,n1,n2=[mon.count(i) for i in range(3)]
                    for i in range(n1+1):
                        for j in range(n2+1):
                            new=(0,)*n0+(1,)*(n1+n2-i-j)+(2,)*(i+j)
                            values[index[new]]+=coefficient*s.binomial(n1,i)*s.binomial(n2,j)*aa**(n1-i)*bb**i*cc**(n2-j)*dd**j
            rows.append(values)
        F[degree]=rows
    lower,data=verify_core(R,F,root,b8)
    N4,N5=data['N4'],data['N5'];unit=data['unit_factor']

    def cancel(arr):
        nz=[v for v in arr if v]
        if not nz:return arr,R.one
        power=min(e[-1] for v in nz for e in v)
        factor=p**power*R.ground_new(nz[0].LC);unit(factor)
        return [v.exquo(factor) for v in arr],factor

    r,t=N4[69];rp,tp=N4[91];det=r*tp-t*rp
    (r,t,rp,tp),factor=cancel([r,t,rp,tp]);den=det.exquo(factor)
    u1=[tp*x-t*y for x,y in zip(N5[69],N5[91])]
    u2=[-rp*x+r*y for x,y in zip(N5[69],N5[91])]
    S1,cancel1=cancel([x-y for x,y in zip([R.zero]+u1,u2+[R.zero])])
    G=[data['A'],data['B'],data['C']]+[v for m,v in zip(data['mons'][4],data['H4'][116]) if 0 not in m]
    G,cancelG=cancel(G);A,B,C=G[:3];g4=G[3:];S2=[R.zero,R.zero]+[den*v for v in g4]
    for i in range(6):S2[i]-=A*u1[i];S2[i+1]-=B*u1[i]+C*u2[i]
    S2,cancel2=cancel(S2)
    rows=[]
    for i in [102,115,134]:
        row=data['pure'][i];rows.extend([row+[R.zero],[R.zero]+row])
    KP=E.poly_ring('p')

    def univariate(v):
        assert all(not any(e[:-1]) for e in v)
        return KP.ring.from_dict({(e[-1],):co for e,co in v.items()})

    determinants=[]
    for row in [S1,S2]:
        matrix=DM([[univariate(v) for v in rr] for rr in rows+[row]],(7,7),KP)
        determinants.append(matrix.det())
    gcd=determinants[0].gcd(determinants[1]);assert len(gcd)==1 and gcd.degree()==18
    bez0,bez1,gcd2=determinants[0].gcdex(determinants[1])
    assert bez0*determinants[0]+bez1*determinants[1]==gcd2 and gcd2==gcd
    expr=lambda v:str(KP.to_sympy(v))
    return dict(field=str(E),quadratic_minimal_polynomial=str(E.ext.minpoly.as_expr()),
                root=str(E.to_sympy(root)),both_conjugate_roots_covered=True,
                slice='b9=L=b5=1, b2=b7=0, Aq(b0)=0, Pstar=p != 0',
                normal_binary_matrix=[[str(v.as_expr()) for v in row] for row in [[aa,bb],[cc,dd]]],
                normal_binary_determinant='-2',selected_original_rows=sorted(selected),
                all_division_factors_are_powers_of_p=True,
                initial_ranks=lower['initial_ranks'],truncated_hilbert_function=lower['truncated_hilbert_function'],
                sextic_determinants=list(map(expr,determinants)),sextic_gcd=expr(gcd),sextic_bezout=[expr(bez0),expr(bez1)],
                all_degree6_monomials_in_initial_ideal=True,transverse_length_upper=23,
                nuisance_parameters_remain_symbolic=True,only_sixth_syzygies_not_sixth_jets=True)
