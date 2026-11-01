#!/usr/bin/env python3
"""Constant quadratic obstruction ideal on the five-variable free-M0 flag chart."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
enc=lambda M:[[str(v)for v in row]for row in M.tolist()]
def polynomial_inverse(Q,a):
    R,*_=ring(','.join(map(str,a)),s.QQ);n=Q.rows;A=[[R.from_expr(Q[i,j])for j in range(n)]+[R.one if i==j else R.zero for j in range(n)]for i in range(n)];perm=list(range(n));pivots=[]
    for col in range(n-1):
        choices=[(i,j)for i in range(col,n)for j in range(col,n)if A[i][j] and A[i][j].is_ground];assert choices;i,j=choices[0];A[col],A[i]=A[i],A[col]
        for row in A:row[col],row[j]=row[j],row[col]
        perm[col],perm[j]=perm[j],perm[col];p=A[col][col].LC;pivots.append(str(p));A[col]=[v.quo_ground(p)for v in A[col]]
        for i in range(n):
            if i==col or not A[i][col]:continue
            c=A[i][col];A[i]=[x-c*y for x,y in zip(A[i],A[col])]
    p=A[-1][n-1];N=s.Matrix([[s.cancel((A[i][j]*p-A[i][n-1]*A[n-1][j]).as_expr()/p.as_expr())for j in range(n,2*n)]for i in range(n-1)]+[[s.cancel(A[n-1][j].as_expr()/p.as_expr())for j in range(n,2*n)]]);out=s.zeros(n)
    for i,j in enumerate(perm):out[j,:]=N[i,:]
    assert(Q.to_DM()*out.to_DM()).to_Matrix()==s.eye(n)
    return out,dict(constant_pivots=pivots,last_pivot=str(p),column_permutation=perm)
def build():
    m=S61();a=s.symbols('a0:5');old=json.loads((ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text());nf=json.loads((ROOT/'reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json').read_text());point=json.loads((ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json').read_text())['strata']['free_M0']['local_data'];beta=s.Matrix(nf['normalized_beta16']).applyfunc(s.sympify);b=s.symbols('b0:15');J=s.Matrix(old['linear']['jacobian']).applyfunc(s.sympify).subs(dict(zip(b,beta[1:,:])));one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];two,F=specialize(m,dict(zip(one,beta)),True);q=s.Matrix(point['point61']).applyfunc(s.Rational)[two,:];fixed=two.index(point['chart_fixed_u']-1);chart=point['chart_columns'];direct=s.zeros(135,45)
    for row,f in enumerate(F):
        value=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q[i]*q[j]
                if i!=fixed:direct[row,chart.index(i)]+=c*q[j]
                if j!=fixed:direct[row,chart.index(j)]+=c*q[i]
            elif len(mon)==1:direct[row,44]+=c*q[mon[0]]
        assert s.expand(value)==0
    assert direct==J;rr,columns=J.to_DM().rref();rr=rr.to_Matrix();columns=list(columns);assert len(columns)==34;free=[i for i in range(45)if i not in columns];T=s.zeros(45,11)
    for j,i in enumerate(free):
        T[i,j]=1
        for k,c in enumerate(columns):T[c,j]=-rr[k,i]
    assert(J.to_DM()*T.to_DM()).is_zero_matrix and T[-1,:]==s.zeros(1,11);assert not s.lcm([s.denom(s.cancel(v))for v in T]).free_symbols;T[:,10]+=s.Rational(125,252)*a[1]*T[:,3];rows=list(J[:,columns].T.to_DM().rref()[1]);iv,gauss=polynomial_inverse(J[rows,columns],a);U=s.zeros(45,34)
    for j,i in enumerate(columns):U[i,:]=iv[j,:]
    return a,J,T,U,rows,columns,free,F,chart,fixed,gauss
def run():
    st=time.time();a,J,T,U,rows,columns,free,F,chart,fixed,gauss=build();mons=list(combinations_with_replacement(range(11),2));R,*_=ring(','.join(map(str,a)),s.QQ);TT=[[R.from_expr(v)for v in row]for row in T.tolist()];HH=[[R.zero for _ in mons]for _ in range(135)]
    for row,f in enumerate(F):
        for mon,c in f.items():
            if fixed in mon:continue
            i,j=[chart.index(k)for k in mon]+[44]*(2-len(mon));co=R.from_expr(c)
            for k,(x,y)in enumerate(mons):HH[row][k]+=co*(TT[i][x]*TT[j][y]+(TT[i][y]*TT[j][x]if x!=y else R.zero))
    Q=s.Matrix([[v.as_expr()for v in row]for row in HH]);JU=J.to_DM()*U.to_DM();A=(Q.to_DM()-JU*Q[rows,:].to_DM()).to_Matrix();sample=dict(zip(a,range(1,6)));selected=list(A.subs(sample).T.to_DM().rref()[1]);assert len(selected)==12;cr,piv=A[selected,:].to_DM().rref();C=cr.to_Matrix().applyfunc(s.factor);assert len(piv)==12 and not any(v.free_symbols for v in C);assert(A.to_DM()-A[:,list(piv)].to_DM()*C.to_DM()).is_zero_matrix
    hessian=s.zeros(12*11,11)
    for row in range(12):
        for k,(i,j)in enumerate(mons):hessian[11*row+i,j]+=C[row,k];hessian[11*row+j,i]+=C[row,k]
    assert hessian[:,:3]==s.zeros(132,3)and hessian.to_DM().rank()==8;det=s.factor(J[rows,columns].det(method='domain-ge'));assert s.factor(det/a[1]).is_Rational
    paths=[Path(__file__)]+[ROOT/f for f in ['code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json','reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json']]
    return dict(id='S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101',all_passed=True,field='QQ(a0,...,a4)',jacobian=enc(J),jacobian_rank=34,jacobian_minor=str(det),pivot_rows=rows,pivot_columns=columns,free_columns=free,tangent_basis_after_shift=enc(T),tangent_shift='v3 -> v3+(125/252)*a1*v10',right_inverse=enc(U),constant_pivot_elimination=gauss,quadratic_monomials=[list(m)for m in mons],selected_quadratic_rows=selected,quadratic_row_coefficients=enc(A[:,list(piv)]),constant_quadratic_rref=enc(C),quadratic_rank=12,common_hessian_rank=8,common_radical_coordinates=[0,1,2],scope=['Exact universal first and second obstruction data on the generic five-variable normal chart.','Proper closed conditions in this five-variable slice have beta images of dimension at most15 after Sp4 and weighted scaling.','No universal cubic/higher normal form, t-nilpotence, or saturation-removal claim.'],proof_note='research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.time()-st)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();result=run()
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
