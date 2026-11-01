#!/usr/bin/env python3
"""Explicit toric-cubic plus two sections algebra of the known-five plane.

The sample identities are exact over QQ. The attached note proves their
rational generic form using the already certified universal plane and degrees.
This is not the full S61 fiber algebra.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61
ROOT=Path(__file__).resolve().parents[1]


def enc(matrix):return [[str(v) for v in row] for row in matrix.tolist()]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    plane_path=ROOT/'reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json';five_path=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';toric_path=ROOT/'reports/TORIC_ALGEBRA_DEGREE_20261031.json';section_path=ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json'
    plane=json.loads(plane_path.read_text());five=json.loads(five_path.read_text());toric=json.loads(toric_path.read_text());section=json.loads(section_path.read_text());model=S61()
    one=plane['weight1_u_indices_zero_based'];two=plane['weight2_u_indices_zero_based'];u=list(map(s.Symbol,plane['beta_variables']));known=list(map(s.sympify,five['toric_S61_coordinates']));sub=dict(zip(u,[known[i] for i in one],strict=True))
    center=s.Matrix(plane['q_center']).subs(sub);N=s.Matrix(plane['q_linear_matrix_N']).subs(sub);T=s.Matrix(plane['generic_plane_recipe']['sample_kernel_T']);d=s.Matrix(section['quadratic_section_q']).subs(sub)-center
    x,y,sigma=s.symbols('x y sigma');variables=(x,y,sigma);q=center+N*T*s.Matrix([x,y])+d*sigma
    values=list(known)
    for i,f in zip(two,q,strict=True):values[i]=f
    restricted=[s.expand(sum(s.Rational(n,den)*s.prod(values[i] for i in mon) for mon,n,den in equation)) for equation in model.reduced['equations']]
    c18,c19=s.symbols('c18 c19');nf=[s.sympify(f).subs({c18:x,c19:y}) for f in toric['quadratic_normal_forms']]
    conics=[s.expand(mon-f) for mon,f in zip((x*x,x*y,y*y),nf,strict=True)];constants=[f.subs({x:0,y:0}) for f in conics];assert any(constants)
    generators=[sigma*x,sigma*y,*[s.expand(f+k*(sigma-2*sigma*sigma)) for f,k in zip(conics,constants,strict=True)]]
    monomials=[s.S.One,x,y,sigma,x*x,x*y,x*sigma,y*y,y*sigma,sigma*sigma]
    G=s.Matrix([[s.Poly(f,*variables).coeff_monomial(m) for m in monomials] for f in generators]);F=s.Matrix([[s.Poly(f,*variables).coeff_monomial(m) for m in monomials] for f in restricted])
    columns=list(G.rref()[1]);assert len(columns)==5
    combinations=F[:,columns]*G[:,columns].inv();assert combinations*G==F
    rows=list(F.T.rref()[1]);assert len(rows)==5
    reverse=G[:,columns]*F[rows,columns].inv();assert reverse*F[rows,:]==G
    # Exact multiplication table in basis1,x,y,sigma,sigma².
    basis=[s.S.One,x,y,sigma,sigma*sigma]
    def vector(f):
        pol=s.Poly(f,*variables);return s.Matrix([pol.coeff_monomial(b) for b in basis])
    products=[s.expand(f-k*sigma+2*k*sigma*sigma) for f,k in zip(nf,constants,strict=True)]
    xx,xy,yy=map(vector,products);zero=s.zeros(5,1)
    mx=s.Matrix.hstack(vector(x),xx,xy,zero,zero)
    my=s.Matrix.hstack(vector(y),xy,yy,zero,zero)
    ms=s.Matrix.hstack(vector(sigma),zero,zero,vector(sigma*sigma),vector((sigma+sigma*sigma)/2))
    assert all(a*b==b*a for a,b in ((mx,my),(mx,ms),(my,ms)))
    identity=s.eye(5)
    def evaluate(f):
        return sum((co*mx**mon[0]*my**mon[1]*ms**mon[2] for mon,co in s.Poly(f,*variables).terms()),s.zeros(5))
    assert all(evaluate(f)==s.zeros(5) for f in generators)
    sigma_relation=sigma*(1+sigma-2*sigma*sigma);assert evaluate(sigma_relation)==s.zeros(5)
    idempotents=[1+sigma-2*sigma*sigma,(4*sigma*sigma-4*sigma)/3,(2*sigma*sigma+sigma)/3]
    matrices=list(map(evaluate,idempotents));assert sum(matrices,s.zeros(5))==identity
    assert all(a*a==a for a in matrices)
    assert all(matrices[i]*matrices[j]==s.zeros(5) for i in range(3) for j in range(i+1,3))
    assert [a.rank() for a in matrices]==[3,1,1]
    # The cyclic unit proves the displayed table has dimension5, rather than
    # merely giving an upper bound through formal reduction rules.
    unit=vector(1);assert s.Matrix.hstack(unit,mx*unit,my*unit,ms*unit,ms**2*unit)==identity
    report={'id':'S61_FIVE_COVER_PRODUCT_20261101','all_passed':True,'sample_field':'QQ','plane_coordinates':list(map(str,variables)),
        'toric_conics':list(map(str,conics)),'toric_conic_constants':list(map(str,constants)),
        'five_generators':list(map(str,generators)),'all134_generator_combinations':enc(combinations),
        'independent_original_rows':rows,'reverse_combinations':enc(reverse),'quadratic_rowspace_rank':5,
        'quotient_basis':list(map(str,basis)),'multiplication_by_x':enc(mx),'multiplication_by_y':enc(my),'multiplication_by_sigma':enc(ms),
        'sigma_relation':str(sigma_relation),'idempotents_toric_q51_q52':list(map(str,idempotents)),'idempotent_ranks':[3,1,1],
        'generic_product_formula':'K[x,y,sigma]/(sigma*x,sigma*y,R_i(x,y)+R_i(0)*(sigma−2sigma²)) = A_toric × K × K',
        'generic_hypotheses':['Restrict to the established open with three reduced toric points spanning their2-plane.','The origin of the c-plane is not a toric point (some R_i(0) != 0).','The five known points impose five independent conditions on polynomials of degree<=2, and the134 restricted quadrics have rank5. All are open conditions witnessed at the sample.'],
        'scope':'The algebra of the intersection with the known-five affine3-plane only; not the full S61 fiber algebra.',
        'proof_note':'research/notes/S61_FIVE_COVER_PRODUCT_20261101.md'}
    sources=[Path(__file__),ROOT/'code/s61_model.py',plane_path,five_path,toric_path,section_path];report['source_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'quadratic_rowspace_rank':5,'product_dimensions':[3,1,1],'output':str(args.output)}))

if __name__=='__main__':main()
