#!/usr/bin/env python3
"""Independent algebraic audit of the published T5 open sextic certificate."""
import argparse, hashlib, json, time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.time()
    source=ROOT/'reports/S61_R0_T5_OPEN_RANK23_20261101.json'
    r=json.loads(source.read_text());assert r['all_passed']
    for name,digest in r['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    frame_path=ROOT/'reports/S61_R0_T5_NORMAL_FRAME_INPUT_20261101.json'
    variables=list(map(s.Symbol,json.loads(frame_path.read_text())['variables']))
    R,*_=ring(variables,s.QQ)
    a,b6,b8,b9,b11,p,xi,eta=s.symbols('b0 b6 b8 b9 b11 p xi eta')
    c=r['normal_initial_certificate'];e=c['sextic_upper_certificate']
    def decode(poly):return R.from_dict({tuple(m):s.Rational(v) for m,v in poly}).as_expr()
    quintics=[[s.sympify(v) for v in row] for row in e['quintic_reference_rows']]
    P=s.sympify(e['open_factors'][2]);X=b9*xi-5*a*(P-1)*eta/16;Y=P*eta
    ref={b6:0,b11:0,b9:1,b8:(1-96900*a*a+142225*a)/1932}
    for i,row in zip([102,115,134],quintics):
        actual=[decode(v) for v in c['quintic_rows'][str(i)]]
        transformed=s.Poly(sum(row[j]*X**(5-j)*Y**j for j in range(6)),xi,eta)
        assert all(s.expand(transformed.coeff_monomial(xi**(5-j)*eta**j)-actual[j])==0 for j in range(6))
    one=[decode(v) for v in e['sextic1']]
    one_ref=[s.cancel(v.subs(ref,simultaneous=True)) for v in one]
    transformed=s.Poly(sum(one_ref[j]*X**(6-j)*Y**j for j in range(7)),xi,eta)
    assert all(s.expand(transformed.coeff_monomial(xi**(6-j)*eta**j)-one[j])==0 for j in range(7))
    rows=[]
    for row in quintics:rows.extend([row+[s.S.Zero],[s.S.Zero]+row])
    M=s.Matrix(rows)
    dm=M.col_join(s.Matrix([one_ref])).to_DM();det=dm.domain.to_sympy(dm.det())
    assert s.cancel(det-s.sympify(e['first_reference_determinant']))==0
    D=s.sympify(e['first_factor']);assert s.Poly(det,a).monic().as_expr()==D
    dual=list(map(s.sympify,e['primitive_dual']));bez=list(map(s.sympify,e['primitive_dual_bezout']))
    assert (M*s.Matrix(dual)).applyfunc(s.expand)==s.zeros(6,1)
    assert s.expand(sum(x*y for x,y in zip(dual,bez)))==1
    cofactor_gcd=s.Poly(0,a)
    for j in range(7):
        dm=M[:,[k for k in range(7) if k!=j]].to_DM()
        cofactor_gcd=s.gcd(cofactor_gcd,s.Poly(dm.domain.to_sympy(dm.det()),a))
    assert s.gcd(cofactor_gcd,s.Poly(D,a)).degree()==0
    two=[decode(v) for v in e['sextic2']]
    canonical=list(map(s.sympify,e['second_canonical_coefficients']))
    assert all(s.expand(b9**(6-j)*v.subs(p,P)-two[j])==0 for j,v in enumerate(canonical))
    # Independently form the inverse binary transformation as a coefficient matrix.
    back=s.Poly((xi+5*a*(p-1)*eta/16),xi,eta).as_expr()
    transform=s.Matrix([[s.Poly(back**(6-j)*(p*eta)**j,xi,eta).coeff_monomial(xi**(6-k)*eta**k)
                         for k in range(7)] for j in range(7)])
    # The contragredient formula avoids division and is p^6 times the inverse action.
    w=[sum(dual[k]*s.binomial(6-j,k-j)*(5*a*(p-1)/16)**(k-j)*p**(6-k)
           for k in range(j,7)) for j in range(7)]
    value=s.factor(sum(x*y for x,y in zip(w,canonical)))
    assert s.cancel(value-s.sympify(e['second_dual_class']))==0
    remainder=s.rem(s.fraction(value)[0],D,a)
    assert s.expand(remainder-s.sympify(e['second_class_remainder']))==0
    L=s.sympify(e['second_linear_factor'])
    assert s.cancel(remainder/(p**7*L)).is_Rational
    aa,bb=map(s.sympify,e['final_bezout']);assert s.expand(aa*D+bb*L)==1
    known=[R.from_expr(s.sympify(v)) for v in e['open_factors']]
    for value in e['cancellation_factors']:
        f=R.from_expr(s.sympify(value))
        for g in known:
            while True:
                q,rem=f.div(g)
                if rem:break
                f=q
        assert len(f)==1 and next(iter(f))==(0,)*R.ngens
    paths=[Path(__file__),source,frame_path,ROOT/'research/notes/S61_R0_T5_OPEN_RANK23_REVIEW_20261101.md']
    return dict(id='S61_R0_T5_OPEN_RANK23_REVIEW_20261101',all_passed=True,accepted=True,
                independent_checks=['All three universal quintic GL2 coefficient identities.',
                                    'Universal first sextic GL2 identity and exact reference determinant.',
                                    'Primitive dual, its Bezout identity, and rank6 at both determinant roots.',
                                    'Second sextic canonical identity, transformed dual class, remainder and final Bezout.',
                                    'All four recorded cancellation factors supported on b9,Aq,P.'],
                proof_review=['The sextic syzygies cancel their degree5 leading forms exactly.',
                              'Degree6 terms containing nu are already in the known initial ideal.',
                              'Same-surface lower23 and upper23 imply finite freeness; t^23 descends by torsion-freeness.',
                              'The stated smooth equivariant source is required for ambient descent.'],
                limitations=['This audit does not independently regenerate the raw135 universal jets; the primary producer does.',
                             'L=0 and Aq roots are outside this theorem. No assertion t=0 or full exhaustion.'],
                input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},
                seconds=time.time()-started)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
