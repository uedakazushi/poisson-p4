#!/usr/bin/env python3
"""Independent sparse canonical/Jacobi and obstruction audit at a=2/9."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MON2=list(itertools.combinations_with_replacement(range(5),2))
MON3=list(itertools.combinations_with_replacement(range(5),3))
X=s.symbols('x0:5')


def tensor(R,data):
    out=[[R.zero]*5 for _ in range(5)]
    for key,value in data.items():
        i,j=map(int,key.split(','));out[i][j]=R.from_expr(s.sympify(value));out[j][i]=-out[i][j]
    return out


def schouten(P,Q,x):
    R=P[0][1].ring;out=[]
    for i,j,k in itertools.combinations(range(5),3):
        value=R.zero
        for l in range(5):
            for u,v,w in [(i,j,k),(j,k,i),(k,i,j)]:
                value+=P[u][l]*Q[v][w].diff(x[l])
                if P is not Q:value+=Q[u][l]*P[v][w].diff(x[l])
        out.append(value)
    return out


def coefficients(P,Q,x):
    R=P[0][1].ring;div=[sum((Q[i][j].diff(x[i])for i in range(5)),R.zero)for j in range(5)]
    mon3=[tuple(ids.count(i)for i in range(5))for ids in MON3]
    return [v.get(tuple(int(i==k)for i in range(5)),R.domain.zero)for v in div for k in range(5)]+[
        v.get(m,R.domain.zero)for v in schouten(P,Q,x)for m in mon3]


def verify(output):
    path=ROOT/'reports/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62.json';record=json.loads(path.read_text())
    a,eps=s.symbols('a eps');field=s.QQ.frac_field(a);R,*x=ring('x0:5',field)
    qdata={'1,4':'a*x0**2','2,3':'-2*x0**2/9','2,4':'4*x0*x1/3','3,4':'2*x0*x2/3+x1**2'}
    Q=tensor(R,qdata);basis=[];indices=[]
    for index,(pair,mon)in enumerate(itertools.product(PAIRS,MON2)):
        if sum(pair)-sum(mon)==-1:
            basis.append(tensor(R,{f'{pair[0]},{pair[1]}':s.prod(X[i]for i in mon)}));indices.append(index)
    assert indices==record['weight_block_columns']and len(basis)==20
    M=s.Matrix.hstack(*(s.Matrix([field.to_sympy(v)for v in coefficients(Q,B,x)])for B in basis))
    active=[i for i in range(375)if any(M.row(i))];assert active==record['active_rows']
    A=M[active,:];A0=A.subs(a,s.Rational(2,9));assert A.to_DM().rank()==17 and A0.rank()==16
    Kt=[tensor(R,t)for t in record['four_kernel_tensors']]
    constant=[]
    for T in Kt:
        column=[]
        for index in indices:
            pair,mon=list(itertools.product(PAIRS,MON2))[index]
            exponent=tuple(mon.count(i)for i in range(5));column.append(field.to_sympy(T[pair[0]][pair[1]].get(exponent,field.zero)))
        constant.append(s.Matrix(column))
    kernel=s.Matrix.hstack(*constant);assert kernel.rank()==4 and A0*kernel==s.zeros(len(active),4)
    r=kernel*s.Matrix(record['unique_projective_second_order_extendible_weight_block_direction'])
    combination=record['unique_projective_second_order_extendible_weight_block_direction']
    reported_R=tensor(R,record['integrable_R'])
    assert all(reported_R[i][j]==sum((combination[l]*Kt[l][i][j]for l in range(4)),R.zero)
               for i in range(5)for j in range(5))
    W=s.Matrix(list(map(s.Rational,record['derivative_cokernel_left_witness'])))
    assert W.T*A0==s.zeros(1,20)
    value=(W.T*A.diff(a)*r)[0];assert value==-15==s.Rational(record['derivative_cokernel_nonzero_value'])
    assert list(W.T*A.diff(a)*kernel)==list(map(s.Rational,record['late_extension_functional']))
    # Rebuild full parameter curve without the producer's derivative routines.
    F=s.QQ.frac_field(eps);S,*y=ring('x0:5',F)
    Q0=tensor(S,{k:s.sympify(v).subs(a,s.Rational(2,9))for k,v in qdata.items()})
    extra=tensor(S,record['integrable_R']);T=tensor(S,{'0,1':'x4**2'})
    P=[[Q0[i][j]+F.from_sympy(eps)*extra[i][j]-75*F.from_sympy(eps**2)*T[i][j]for j in range(5)]for i in range(5)]
    assert not any(schouten(P,P,y))
    assert all(not sum((P[i][j].diff(y[i])for i in range(5)),S.zero)for j in range(5))
    degeneration=S.zero
    for i in range(5):
        j,k,l,m=[t for t in range(5)if t!=i]
        degeneration+=(-1)**i*y[i]*(P[j][k]*P[l][m]-P[j][l]*P[k][m]+P[j][m]*P[k][l])
    assert degeneration.get((5,0,0,0,0))==F.from_sympy(-s.Rational(4,81))
    assert not(degeneration-S.from_expr(s.sympify(record['degeneracy_polynomial'])))
    # Weight separates all possible second corrections; the needed weight-7
    # space of all150 raw tensors is exactly the displayed one-dimensional T.
    negative_seven=[(pair,mon)for pair,mon in itertools.product(PAIRS,MON2)if sum(pair)-sum(mon)==-7]
    assert negative_seven==[((0,1),(4,4))]
    z=s.symbols('z0:4');G=s.QQ.frac_field(*z);U,*v=ring('x0:5',G)
    Qz=tensor(U,{k:s.sympify(text).subs(a,s.Rational(2,9))for k,text in qdata.items()})
    kernel_tensors=[tensor(U,b)for b in record['four_kernel_tensors']]
    Rz=[[sum((G.gens[l]*kernel_tensors[l][i][j]for l in range(4)),U.zero)for j in range(5)]for i in range(5)]
    image=[s.sympify(G.to_sympy(c))for c in coefficients(Qz,tensor(U,{'0,1':'x4**2'}),v)]
    quadratic=[s.S.Zero]*25+[G.to_sympy(c)for poly in schouten(Rz,Rz,v)for m in MON3
                                          for c in [poly.get(tuple(m.count(i)for i in range(5)),G.zero)]]
    pivot=next(i for i,c in enumerate(image)if c)
    w0,w1,w2,t=s.symbols('w0 w1 w2 t')
    changed=[s.Poly(s.expand((quadratic[i]-image[i]*quadratic[pivot]/image[pivot]).subs(dict(zip(z,[w0+4*t,w1+6*t,w2-4*t,t])))),w0,w1,w2,t)for i in range(375)]
    mons=[w0*w0,w0*w1,w0*w2,w1*w1,w1*w2,w2*w2,w0*t,w1*t,w2*t]
    assert all(not p.coeff_monomial(t*t)for p in changed)
    assert s.Matrix([[p.coeff_monomial(m)for m in mons]for p in changed]).rank()==9
    for name,digest in record['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    paths=[Path(__file__),path,ROOT/'research/notes/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62.md']
    report=dict(id='Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62_REVIEW',all_passed=True,
                independent_raw375_by20_linearization_generic17_special16=True,
                independent_late_cokernel_pairing=-15,
                independent_all350_exact_curve_Jacobi_and_divergence_zero=True,
                independent_degree5_polynomial_and_constant_x0_power_checked=True,
                full_second_order_weight_block_obstruction_ideal_checked=True,
                unique_projective_second_order_direction=[4,6,-4,1],
                full35_variable_formal_model_claimed=False,new_component_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
    with output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS INDEPENDENT EXACT CURVE + RESONANT KERNEL/OBSTRUCTION')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();verify(args.output)
