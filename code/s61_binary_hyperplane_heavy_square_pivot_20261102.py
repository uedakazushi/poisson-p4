#!/usr/bin/env python3
"""A global heavy-square unit pivot on the binary hyperplane chart."""
import argparse,hashlib,json,time,os,subprocess,tempfile
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
root=Path(__file__).resolve().parents[1]
NOTE=root/'research/notes/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.md'
def verify():
    begin=time.time()
    data=json.loads((root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json').read_text())
    old=json.loads((root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();pol=lambda x:R.from_expr(s.sympify(x));frac=lambda x:K.from_expr(s.sympify(x))
    T=[[frac(row[j]) for j in data['normal34']] for row in data['T34']]
    U=[[frac(v) for v in row] for row in data['U34']]
    td=[]
    for j in range(8):
        den=R.one
        for row in T:den=den.lcm(row[j].denom)
        td.append(den)
    ud=R.one
    for row in U:
        for v in row:ud=ud.lcm(v.denom)
    Tn=[[v.numer*td[j].exquo(v.denom) for j,v in enumerate(row)] for row in T]
    Un=[[v.numer*ud.exquo(v.denom) for v in row] for row in U]
    m=S61();beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()};two,F=specialize(m,beta,True)
    fixed=old['fixed'];cols=old['cols'];J=[[R.zero for j in range(45)] for i in range(135)];H=[]
    for i,f in enumerate(F):
        row=[]
        for mon,c in f.items():
            cc=pol(c)
            if len(mon)==2:
                x,y=mon
                if x!=fixed and y==fixed:J[i][cols.index(x)]+=cc
                if y!=fixed and x==fixed:J[i][cols.index(y)]+=cc
            elif len(mon)==1 and mon[0]==fixed:J[i][44]+=cc
            if fixed not in mon:
                ij=[cols.index(j) for j in mon]+[44]*(2-len(mon))
                if cc:row.append((ij,cc))
        H.append(row)


    zero=R.zero
    selected=data['selected34_rows']
    assert all(sum((v*Tn[j][k] for j,v in enumerate(row) if v),zero)==0
               for row in J for k in range(8))
    assert all(sum((v*Un[j][k] for j,v in enumerate(J[i]) if v),zero)==(ud if r==k else zero)
               for r,i in enumerate(selected) for k in range(34))
    Q=[sum((c*Tn[i][5]*Tn[j][5] for (i,j),c in terms),zero) for terms in H]
    W=[-sum((v*Q[selected[j]] for j,v in enumerate(row) if v),zero) for row in Un]
    A=[ud*Q[i]+sum((v*W[j] for j,v in enumerate(row) if v),zero) for i,row in enumerate(J)]
    Fcore=500*a+90*e**2+3*e-24
    B=50*a*d+18*d*e**2-9*d*e-2*d+60*e-40
    C=100*a+18*e**2-9*e-2
    P=125*a**2*d-450*a*e-200*a-54*e**3-27*e**2+33*e+6
    Rpol=15*a*d*e-10*a*d+200*a+18*e**2+6*e-12
    assert td[5]==B/50 and ud==a**2*Fcore*B/25000
    denominator=ud*td[5]**2
    coefficients=[K(v)/K(denominator) for v in A]
    assert coefficients[77]==-3*K(P**2)/(500*K(B**2))
    assert coefficients[131]==3*K(Rpol)/(10*K(B))
    assert data['equation_weights'][77]==data['equation_weights'][131]==4
    assert all(not c or data['equation_weights'][i]==4 for i,c in enumerate(coefficients))
    S=3*e-2;D0=50*a+9*e**2+3*e-6;E=B-S*(d*(6*e+1)+40)
    alpha=-S*E;gamma=25*a*E+10*S*C;V=50*a*B*C
    assert S*P-25*a*Rpol==-D0*C
    assert 50*a*B==D0*E+10*S*Rpol
    assert alpha*P+gamma*Rpol==V
    left=-s.QQ(500,3)*B**2*alpha**2
    right=s.QQ(10,3)*B*gamma*(2*alpha*P+gamma*Rpol)
    assert K(left)*coefficients[77]+K(right)*coefficients[131]==K(V**2)
    nw=data['normal_weights34'];assert nw==[1,1,1,1,1,2,1,1]
    def count(weight,index=0):
        if index==8:return int(weight==0)
        return sum(count(weight-k*nw[index],index+1) for k in range(weight//nw[index]+1))
    assert count(4)==239
    quadpath=root/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
    quad=json.loads(quadpath.read_text());assert quad['all_passed']
    # Every saved quadratic relation is a pure-light quadratic.
    for row in quad['quadratic_rows_sparse']:
        for mon,poly in row:assert mon[5]==0 and sum(mon)==2
    paths=[Path(__file__),NOTE,
           root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json',
           root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',quadpath,
           root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py']
    return dict(id='S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102',all_passed=True,
        core_variables=['a','e','d'],raw_four_jet_cache_used=False,
        original_equations=135,heavy_variable=5,selected_source_rows_zero_based=[77,131],
        linear_section_checked=True,heavy_column_denominator=str(td[5].as_expr()),
        right_inverse_denominator=str(ud.as_expr()),
        coefficient_formulas=[str(coefficients[i].as_expr()) for i in [77,131]],
        P=str(P.as_expr()),R=str(Rpol.as_expr()),alpha=str(alpha.as_expr()),gamma=str(gamma.as_expr()),
        bezout_unit=str(V.as_expr()),row_combination_coefficients=[str(left.as_expr()),str(right.as_expr())],
        combined_heavy_square_coefficient=str((V**2).as_expr()),
        identities_checked=['S*P-25*a*R=-D0*C','50*a*B=D0*E+10*S*R',
                            'alpha*P+gamma*R=50*a*B*C','left*c77+right*c131=(50*a*B*C)^2'],
        weighted_four_columns=239,columns_after_unit_pivot=238,
        heavy_square_absent_in_lower_weight_products=True,
        limitations=['No uniform rank238 assertion for the remaining columns.',
                     'No new length65 or ambient saturation-removal assertion.'],
        proof_note=str(NOTE.relative_to(root)),
        input_sha256={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-begin)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    result=verify()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','weighted_four_columns','columns_after_unit_pivot','seconds']}))
