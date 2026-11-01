"""Polynomial certificates for the T5 normal Hilbert function through degree five."""
from itertools import combinations_with_replacement
import sympy as s


def verify(R,F):
    b9=R.from_expr(s.Symbol('b9'))
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

    def monomial_primitive(arr):
        nonzero=[p for p in arr if p]
        if not nonzero:return arr
        exponent=tuple(min(mon[j] for p in nonzero for mon in p) for j in range(R.ngens))
        factor=R.from_dict({exponent:nonzero[0].LC})
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
    return dict(initial_ranks={'2':1,'3':4,'4':10,'5':18},truncated_hilbert_function=[1,3,5,6,5,3],truncated_length=23,
                quadratic_generator='x^2',primitive_cubic_coefficients=[str(p.exquo(common)) for p in [A,B,C]],
                quartic_pivot_rows=[r0,r1],quartic_normal_determinant=str(determinant),
                quartic_initial_contains_all_x_times_binary_cubics=True,
                first_lift_cubic_has_no_pure_binary_terms=True,
                quintic_pivot_rows=selected_rows,quintic_pivot_columns=selected_columns,quintic_pivot_determinant=str(det),
                quintic_rows={str(i):[ep(p) for p in row] for i,row in enumerate(pure) if any(row)},
                all_polynomial_cross_product_identities_verified=True)
