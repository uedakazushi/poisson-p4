#!/usr/bin/env python3
"""Exact six-light-coordinate chart of the nonsquare quotient by ell.

The API preserves the actual unscaled polynomial images, including v5 of
weight2. The chart is C*(5*a*d-6*e-4) !=0; no all-core rank24/41 is asserted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]


def load_chart():
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    P,*v=ring('v0,v1,v2,v3,v4,v5,v6,v7',R)
    T,*w=ring('w0,w1,w2,w3,w4,w5,H',K)
    qp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
    pp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_PAIR_20261102.json'
    data,pair=[json.loads(p.read_text())for p in [qp,pp]]
    for report in [data,pair]:
        assert report['all_passed']
        for name,want in report['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want,name
    quadrics=[P.from_dict({tuple(m):R.from_dict({tuple(ex):s.QQ(c)for ex,c in terms})
                           for m,terms in row})for row in data['quadratic_rows_sparse']]
    ell=(-(d*e+8)*quadrics[1]+8*quadrics[2]+d*quadrics[3]).exquo(v[3])
    assert s.expand(ell.as_expr()-s.sympify(pair['linear_form']))==0
    coefficients=[ell.get(tuple(int(i==j)for j in range(8)),R.zero)for i in range(8)]
    assert coefficients[5]==0
    F=500*a+90*e**2+3*e-24
    B=(50*a+18*e**2-9*e-2)*d+60*e-40
    C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4
    g=coefficients[6]+coefficients[1]
    assert g==-s.Rational(16,243)*C*J
    # Coordinates are (w0,...,w5,H)=(v0,v2,v4,v7,Z,W,v5).
    numerator=(T(K(coefficients[0]+3*coefficients[1]))*w[0]
               +T(K(coefficients[2]))*w[1]+T(K(coefficients[4]))*w[2]
               +T(K(coefficients[7]))*w[3]-T(K(coefficients[1]))*w[4]
               +T(K(coefficients[3]))*w[5])
    v6_image=-numerator/K(g)
    images=[w[0],3*w[0]+v6_image-w[4],w[1],w[5],w[2],w[6],v6_image,w[3]]
    def restrict(poly):
        total=T.zero
        for mon,c in poly.items():
            term=T(K(c))
            for power,image in zip(mon,images):
                if power:term*=image**power
            total+=term
        return total
    assert restrict(ell)==0
    assert restrict(v[0])==w[0]and restrict(v[2])==w[1]and restrict(v[4])==w[2]
    assert restrict(v[7])==w[3]and restrict(3*v[0]-v[1]+v[6])==w[4]
    assert restrict(v[3])==w[5]and restrict(v[5])==w[6]
    def unit(poly,factors):
        assert poly
        for f in factors:
            while poly.rem(f)==0:poly=poly.exquo(f)
        return poly.is_ground
    assert all(unit(c.denom,[C,J])for image in images for c in image.values())
    restricted=[restrict(quadrics[i])for i in [0,1,3,5]]
    assert all(not mon[6]for q in restricted for mon in q)
    assert all(unit(c.denom,[C,J])for q in restricted for c in q.values())
    cleared=[];multipliers=[]
    for q in restricted:
        den=R.one
        for c in q.values():den=den.lcm(c.denom)
        assert unit(den,[C,J])
        out={m:(c*den).numer.exquo((c*den).denom)for m,c in q.items()}
        factor=R.one
        for f in [a,F,B,C,J]:
            while all(c.rem(f)==0 for c in out.values()):
                out={m:c.exquo(f)for m,c in out.items()};factor*=f
        multiplier=K(den)/K(factor)
        assert T.from_dict({m:K(c)for m,c in out.items()})==q*multiplier
        cleared.append(out);multipliers.append(multiplier)
    return dict(R=R,K=K,old_ring=P,old_variables=v,new_ring=T,new_variables=w,
                old_quadrics=quadrics,ell=ell,ell_coefficients=coefficients,
                images=images,restrict=restrict,restricted_qfour=restricted,
                cleared_qfour=cleared,clearing_unit_multipliers=multipliers,
                F=F,B=B,C=C,J=J,g=g,input_paths=[qp,pp])


def verify():
    started=time.monotonic();c=load_chart();R=c['R']
    rows=[];sample={R.symbols[0]:1,R.symbols[1]:2,R.symbols[2]:3}
    T=c['new_ring'];w=c['new_variables']
    products=[q*w[j]for q in c['restricted_qfour']for j in range(6)]
    monomials=sorted(set().union(*products));assert len(products)==24 and len(monomials)==56
    matrix=s.Matrix([[s.cancel(q.get(mon,c['K'].zero).as_expr()).subs(sample)for mon in monomials]for q in products])
    assert matrix.rank()==24
    note=ROOT/'research/notes/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.md'
    files=[Path(__file__),note,*c['input_paths']]
    return dict(id='S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62',all_passed=True,
                field='QQ(a,e,d)',known_open='a*F*B*C*J!=0, J=5ad-6e-4',
                light_coordinate_order=['v0','v2','v4','v7','3v0-v1+v6','v3'],heavy_coordinate='H=v5, weight2',
                exact_chart_coefficient=str(c['g']),
                images_of_v0_through_v7=[str(p)for p in c['images']],
                ell_maps_identically_to_zero=True,inverse_light_and_heavy_coordinate_identities_verified=True,
                raw_coordinate_and_qfour_denominators_only_C_J_powers=True,
                unscaled_restricted_qfour=[str(p)for p in c['restricted_qfour']],
                cleared_qfour_sparse=[[[list(mon),[[list(ex),str(co)]for ex,co in coefficient.items()]]
                                      for mon,coefficient in q.items()]for q in c['cleared_qfour']],
                exact_clearing_unit_multipliers=list(map(str,c['clearing_unit_multipliers'])),
                retained_four_quadric_source_indices=[0,1,3,5],
                linear_product_matrix_shape=[24,56],generic_sample=[1,2,3],generic_sample_rank=24,
                all_core_rank24_or_rank41_claimed=False,
                ambient_nilpotence_or_quartic_coverage_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},
                seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS',result['id'],result['seconds'])
