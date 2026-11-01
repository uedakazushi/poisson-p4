"""Exact polynomial arithmetic for the T5 L=0 chart."""
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring


def compute(model, orbit, frame, g6, progress=print):
    b=s.symbols('b0:12'); substitutions={s.Symbol(k):s.sympify(v) for k,v in frame['slice'].items()}
    variables=list(map(s.Symbol,frame['variables'])); R,*_=ring(variables,s.QQ); b9=R.one
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    qi={i:j for j,i in enumerate(two)}; ui={i:j for j,i in enumerate(one)}
    H=s.Matrix(orbit['beta_kernel']).applyfunc(s.Rational)
    beta=[R.from_expr(v) for v in (H*s.Matrix(b)).xreplace(substitutions)]
    q=[s.Rational(orbit['point61'][i]) for i in two]
    fixed=frame['fixed_q']; chart=frame['chart_columns']; rows=frame['pivot_rows']; columns=frame['pivot_columns']
    J=[[R.zero for j in range(45)] for i in range(135)]; terms=[]; cache={}
    for row,equation in enumerate(model.reduced['equations']+[g6]):
        high=R.zero; local={}
        for mon,nu,de in equation:
            qm=tuple(qi[i] for i in mon if i in qi); bm=tuple(ui[i] for i in mon if i in ui)
            if bm not in cache:
                value=R.one
                for i in bm:value*=beta[i]
                cache[bm]=value
            coefficient=R.ground_new(s.Rational(nu,de))*cache[bm]
            if len(qm)==2:
                i,j=qm; high+=coefficient*q[i]*q[j]
                if i!=fixed:J[row][chart.index(i)]+=coefficient*q[j]
                if j!=fixed:J[row][chart.index(j)]+=coefficient*q[i]
            elif len(qm)==1:J[row][44]+=coefficient*q[qm[0]]
            if fixed not in qm:
                ij=tuple([chart.index(i) for i in qm]+[44]*(2-len(qm)))
                local[ij]=local.get(ij,R.zero)+coefficient
        assert not high,(row,'highest equation')
        terms.append([(i,j,v) for (i,j),v in local.items() if v])

    def clear(expression,denominator):
        n,d=s.fraction(s.sympify(expression))
        return (R.from_expr(n)*denominator).exquo(R.from_expr(d))

    fullT=[[clear(v,b9) for v in row] for row in frame['tangent']]
    I=[[clear(v,b9**6) for v in row] for row in frame['inverse']]
    for i,row in enumerate(J):
        for j in range(5):assert not sum((a*fullT[k][j] for k,a in enumerate(row) if a),R.zero)
    for i,row in enumerate(frame['free_columns']):
        assert fullT[row]==[b9 if i==j else R.zero for j in range(5)]
    for i,row in enumerate(rows):
        for j in range(40):
            value=sum((J[row][c]*I[k][j] for k,c in enumerate(columns) if J[row][c] and I[k][j]),R.zero)
            assert value==(b9**6 if i==j else R.zero)
    progress('raw135 linear frame and inverse verified',flush=True)
    W={1:[[row[j] for j in frame['normal_columns']] for row in fullT]}; den={1:b9}
    residual={}; rden={}; monomials={d:list(combinations_with_replacement(range(3),d)) for d in range(1,6)}
    statistics=[]

    def reduce_den(matrix,denominator):
        assert denominator==R.one
        return matrix,denominator

    for degree in [2,3,4]:
        index={m:i for i,m in enumerate(monomials[degree])}; n=len(index)
        raw=[[R.zero]*n for _ in range(135)]; common=b9**(2*degree-2)
        for row,ts in enumerate(terms):
            for i,j,coefficient in ts:
                for a in range(1,degree):
                    scale=common.exquo(den[a]*den[degree-a])
                    for k,ma in enumerate(monomials[a]):
                        if not W[a][i][k]:continue
                        for l,mb in enumerate(monomials[degree-a]):
                            if W[degree-a][j][l]:raw[row][index[tuple(sorted(ma+mb))]]+=coefficient*scale*W[a][i][k]*W[degree-a][j][l]
        wd=[[R.zero]*n for _ in range(45)]
        for i,c in enumerate(columns):
            for j in range(n):wd[c][j]=-sum((I[i][k]*raw[r][j] for k,r in enumerate(rows) if I[i][k] and raw[r][j]),R.zero)
        W[degree],den[degree]=reduce_den(wd,b9**6*common)
        rd=common.lcm(den[degree])
        rp=[[raw[i][j]*rd.exquo(common)+sum((J[i][k]*W[degree][k][j] for k in range(45) if J[i][k] and W[degree][k][j]),R.zero)*rd.exquo(den[degree]) for j in range(n)] for i in range(135)]
        residual[degree],rden[degree]=reduce_den(rp,rd)
        statistics.append(dict(degree=degree,lift_monomials=sum(len(p) for row in W[degree] for p in row),residual_monomials=sum(len(p) for row in residual[degree] for p in row),lift_denominator=str(den[degree]),residual_denominator=str(rden[degree])))
        progress('normal jet '+str(degree)+' verified',flush=True)
    # In degree >=4 all monomials containing the first normal variable belong
    # to the initial ideal. Only the six binary degree-five coefficients are needed.
    binary={d:[i for i,m in enumerate(monomials[d]) if 0 not in m] for d in range(1,6)}
    binary_monomials={d:[monomials[d][i] for i in binary[d]] for d in binary}
    index={m:i for i,m in enumerate(binary_monomials[5])}; raw=[[R.zero]*6 for _ in range(135)]
    for row,ts in enumerate(terms):
        for i,j,coefficient in ts:
            for a in range(1,5):
                scale=(b9**8).exquo(den[a]*den[5-a])
                for k in binary[a]:
                    if not W[a][i][k]:continue
                    for l in binary[5-a]:
                        if W[5-a][j][l]:raw[row][index[tuple(sorted(monomials[a][k]+monomials[5-a][l]))]]+=coefficient*scale*W[a][i][k]*W[5-a][j][l]
    projected=[[R.zero]*6 for _ in range(135)]
    for row in range(135):
        if row in rows:continue
        left=[sum((J[row][c]*I[i][k] for i,c in enumerate(columns) if J[row][c] and I[i][k]),R.zero) for k in range(len(rows))]
        for j in range(6):projected[row][j]=b9**6*raw[row][j]-sum((lc*raw[r][j] for lc,r in zip(left,rows) if lc and raw[r][j]),R.zero)
    residual[5],rden[5]=reduce_den(projected,b9**14)
    assert rden=={2:R.one,3:b9,4:b9**3,5:b9**5}
    statistics.append(dict(degree=5,binary_only=True,residual_monomials=sum(len(p) for row in residual[5] for p in row),residual_denominator=str(rden[5])))
    progress('binary normal jet 5 verified',flush=True)
    return R,residual,statistics
