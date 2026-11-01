"""T5 transverse upper23 on b9*Aq*P != 0, using sextic syzygies of five-jets."""
from itertools import combinations_with_replacement
import sympy as s


def verify_core(R,F):
    b9=R.from_expr(s.Symbol('b9'))
    b0,b6,b8,b11=[R.from_expr(s.Symbol(v)) for v in ['b0','b6','b8','b11']]
    aq=5254333434224*b0*b0-8538753610352*b0+3402879618517
    pp=96900*b0*b0-24150*b0*b6*b9-142225*b0-4140*b11*b9+1932*b8*b9*b9
    known=[b9,aq,pp]
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
        nonzero=[p for p in arr if p]
        if not nonzero:return arr
        factor=nonzero[0]
        for p in nonzero[1:]:factor=factor.gcd(p)
        unit_factor(factor)
        return [p.exquo(factor) for p in arr]

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

    b0,b6,b8,b11=[R.from_expr(s.Symbol(v)) for v in ['b0','b6','b8','b11']]
    known=[b9,5254333434224*b0*b0-8538753610352*b0+3402879618517,
           96900*b0*b0-24150*b0*b6*b9-142225*b0-4140*b11*b9+1932*b8*b9*b9]

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
    ep=lambda p:[[list(e),str(co)] for e,co in p.items()]
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


def verify(R,F):
    lower,data=verify_core(R,F)
    N4,N5=data['N4'],data['N5'];known=data['known'];unit=data['unit_factor']
    ep=lambda f:[[list(e),str(c)] for e,c in f.items()]

    def cancel(arr):
        factor=R.one
        for f in known:
            while any(arr):
                qr=[v.div(f) for v in arr]
                if any(rem for q,rem in qr):break
                arr=[q for q,rem in qr];factor*=f
        nz=[v for v in arr if v]
        if nz:
            exponent=tuple(min(e[j] for v in nz for e in v) for j in range(R.ngens))
            f=R.from_dict({exponent:nz[0].LC})
            arr=[v.exquo(f) for v in arr];factor*=f
        unit(factor)
        return arr,factor

    r0,r1=69,91
    r,t=N4[r0];rp,tp=N4[r1];det=r*tp-t*rp
    (r,t,rp,tp),cofactor_cancel=cancel([r,t,rp,tp])
    denominator=det.exquo(cofactor_cancel)
    u1=[tp*x-t*y for x,y in zip(N5[r0],N5[r1])]
    u2=[-rp*x+r*y for x,y in zip(N5[r0],N5[r1])]
    # Leading forms: h1=nu*xi*eta^2, h2=nu*eta^3.
    # S1=eta*h1-xi*h2. Only the pure binary fifth parts contribute in degree6.
    sextic1,cancel1=cancel([x-y for x,y in zip([R.zero]+u1,u2+[R.zero])])
    G=[data['A'],data['B'],data['C']]+[v for m,v in zip(data['mons'][4],data['H4'][116]) if 0 not in m]
    G,gcancel=cancel(G);A,B,C=G[:3];g4=G[3:]
    # S2=eta^2*g-A*xi*h1-B*eta*h1-C*eta*h2, g3=nu*(A xi^2+B xi eta+C eta^2).
    sextic2=[R.zero,R.zero]+[denominator*v for v in g4]
    for i in range(6):
        sextic2[i]-=A*u1[i];sextic2[i+1]-=B*u1[i]+C*u2[i]
    sextic2,cancel2=cancel(sextic2)

    a,b6,b8,b9,b11,p=s.symbols('b0 b6 b8 b9 b11 p')
    P=known[2].as_expr();shear=-5*a*(P-1)/16
    ref={b6:0,b11:0,b9:1,b8:(1-96900*a*a+142225*a)/1932}
    references=[]
    for i in [102,115,134]:
        row=[v.as_expr() for v in data['pure'][i]]
        reference=[s.cancel(c.subs(ref,simultaneous=True)) for c in row]
        for k in range(6):
            actual=sum(reference[j]*s.binomial(5-j,k-j)*b9**(5-k)*shear**(k-j)*P**j for j in range(k+1))
            assert s.cancel(actual-row[k])==0,(i,k)
        references.append(reference)
    one=[v.as_expr() for v in sextic1]
    one_ref=[s.cancel(v.subs(ref,simultaneous=True)) for v in one]
    for k in range(7):
        actual=sum(one_ref[j]*s.binomial(6-j,k-j)*b9**(6-k)*shear**(k-j)*P**j for j in range(k+1))
        assert s.cancel(actual-one[k])==0,k
    rows=[]
    for row in references:rows.extend([row+[s.S.Zero],[s.S.Zero]+row])
    matrix=s.Matrix(rows)
    dm=matrix.col_join(s.Matrix([one_ref])).to_DM()
    first_det=s.factor(dm.domain.to_sympy(dm.det()))
    D1=s.Poly(s.fraction(first_det)[0],a).monic().as_expr()
    assert s.degree(D1,a)==2
    # A primitive dual has no common zero. We certify this explicitly by a polynomial Bezout chain.
    kernel=matrix.to_DM().nullspace().to_Matrix();assert kernel.rows==1
    dual=list(kernel[0,:]);g=dual[0]
    for c in dual[1:]:g=s.gcd(g,c)
    dual=[s.cancel(c/g) for c in dual]
    assert (matrix*s.Matrix(dual)).applyfunc(s.cancel)==s.zeros(6,1)
    gcd=s.Poly(0,a);bez=[s.S.Zero]*7
    for j,c in enumerate(dual):
        h,k,new=s.gcdex(gcd,s.Poly(c,a))
        bez=[s.expand(h.as_expr()*v) for v in bez];bez[j]+=k.as_expr();gcd=new
    assert gcd.as_expr()==1 and s.expand(sum(v*c for v,c in zip(bez,dual)))==1
    # The second relation is not individually GL2-covariant: compare its class modulo old sextics.
    canonical_sub={b6:0,b11:0,b9:1,b8:(p-96900*a*a+142225*a)/1932}
    two=[v.as_expr() for v in sextic2]
    canonical=[s.cancel(c.subs(canonical_sub,simultaneous=True)) for c in two]
    for i,(original,c) in enumerate(zip(two,canonical)):
        assert s.cancel(b9**(6-i)*c.subs(p,P)-original)==0,i
    sh=-5*a*(p-1)/16
    transformed_dual=[sum(dual[k]*s.binomial(6-j,k-j)*(-sh)**(k-j)*p**(6-k) for k in range(j,7)) for j in range(7)]
    value=s.factor(sum(w*c for w,c in zip(transformed_dual,canonical)))
    numerator=s.fraction(value)[0]
    remainder=s.factor(s.rem(numerator,D1,a))
    linear=s.Poly(s.cancel(remainder/p**7),a).monic().as_expr()
    assert s.degree(linear,a)==1 and s.cancel(remainder/(p**7*linear)).is_Rational
    ba,bb,gg=s.gcdex(D1,linear,a);assert gg==1 and s.expand(ba*D1+bb*linear)==1
    # The old six rows have rank6 at both roots of D1. This removes any ambiguity in using their dual.
    cofactor_gcd=s.Poly(0,a)
    for column in range(7):
        minor=matrix[:,[j for j in range(7) if j!=column]].to_DM()
        detminor=minor.domain.to_sympy(minor.det())
        cofactor_gcd=s.gcd(cofactor_gcd,s.Poly(detminor,a))
    assert s.gcd(cofactor_gcd,s.Poly(D1,a)).degree()==0
    lower['sextic_upper_certificate']=dict(
        open_factors=[str(v.as_expr()) for v in known],
        normalized_binary_change=['X=b9*xi-5*b0*(P-1)*eta/16','Y=P*eta'],
        quintic_reference_rows=[[str(v) for v in row] for row in references],
        sextic1=[ep(v) for v in sextic1],sextic2=[ep(v) for v in sextic2],
        cancellation_factors=[str(v.as_expr()) for v in [cofactor_cancel,gcancel,cancel1,cancel2]],
        first_reference_determinant=str(first_det),first_factor=str(D1),
        primitive_dual=list(map(str,dual)),primitive_dual_bezout=list(map(str,bez)),
        second_canonical_coefficients=list(map(str,canonical)),second_dual_class=str(value),
        second_class_remainder=str(remainder),second_linear_factor=str(linear),
        final_bezout=[str(ba),str(bb)],all_degree6_monomials_in_initial_ideal=True,
        transverse_length_upper=23,sixth_jet_not_needed=True,
        second_sextic_individual_covariance_not_assumed=True)
    return lower
