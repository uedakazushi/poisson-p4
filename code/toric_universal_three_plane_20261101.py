#!/usr/bin/env python3
"""Universal affine-linear toric reconstruction and the known five-lift plane.

All16 first-jet variables remain symbolic. This gives a rational affine
three-plane containing the five known generic lifts, not a containment of the
full S61 fiber. The attached note proves generic length five of its intersection.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from itertools import combinations,combinations_with_replacement
from pathlib import Path
import sympy as s
from toric_algebra_fiber import TRIPLES,multiplication_from_jet,associativity_quadrics
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J,COORDS
from s61_model import S61,X,quotient_coefficients
ROOT=Path(__file__).resolve().parents[1]


def encode(matrix):
    return [[str(f) for f in row] for row in matrix.tolist()]


def derive(model):
    z,u,beta=primitive_beta(model);y=s.Matrix(z);c=s.symbols('c0:20')
    vector=s.Matrix([s.diff(beta[i,j],z[k]) for (i,j),k in COORDS])
    m=multiplication_from_jet(vector,dict(zip(TRIPLES,c,strict=True)))
    form,equations=associativity_quadrics(m)
    my=sum((y[i]*m[i] for i in range(4)),s.zeros(4))
    quadratic=(my*J*my.T-y*(y.T*form)*J-J*form*y*y.T).applyfunc(s.expand)
    pairs=list(combinations(range(4),2));ym=list(combinations_with_replacement(range(4),2))
    cm=list(combinations_with_replacement(range(20),2))
    coefficients=[s.Poly(quadratic[i,j],*y).coeff_monomial(y[a]*y[b])
                  for i,j in pairs for a,b in ym]
    polynomials=[s.Poly(f,*c) for f in equations]
    leading=s.Matrix([[f.coeff_monomial(c[a]*c[b]) for a,b in cm] for f in polynomials])
    assert not leading.free_symbols
    rows=list(leading.T.rref()[1]);assert len(rows)==45
    independent=leading[rows,:];columns=list(independent.rref()[1])
    quadratic_symbol=s.Matrix([[s.Poly(f,*c).coeff_monomial(c[a]*c[b]) for a,b in cm]
                              for f in coefficients])
    reduction=quadratic_symbol[:,columns]*independent[:,columns].inv()
    assert reduction*independent==quadratic_symbol
    assert quadratic_symbol.rank()==45
    reduced=(s.Matrix(coefficients)-reduction*s.Matrix([equations[i] for i in rows])).applyfunc(s.expand)
    assert all(s.Poly(f,*c).total_degree()<=1 for f in reduced)
    # All constant leading-symbol relations give homogeneous beta*c equations.
    syzygies=leading.T.nullspace();assert len(syzygies)==51
    linear=[]
    for relation in syzygies:
        f=s.expand(sum(value*equations[i] for i,value in enumerate(relation) if value))
        assert s.Poly(f,*c).total_degree()<=1
        assert f.subs(dict.fromkeys(c,0))==0
        linear.append(f)
    constraint=s.Matrix(linear).jacobian(c)
    assert not constraint.free_symbols.intersection(c)
    assert all(s.Poly(f,*u).total_degree()==1 for f in constraint if f)
    # Keep the affine two-jet; arbitrary homogeneous B_0i quadratics contribute
    # only the affine cubic tail and cannot change the weight-two coordinates.
    p=s.zeros(4);first=-J*beta*J
    for k,(i,j) in enumerate(pairs):
        p[i,j]=J[i,j]+first[i,j]+sum(reduced[k*10+l]*y[a]*y[b] for l,(a,b) in enumerate(ym))
        p[j,i]=-p[i,j]
    homogeneous=s.zeros(5)
    for i,j in pairs:
        homogeneous[i+1,j+1]=sum(value*X[0]**(2-sum(mon))*s.prod(X[k+1]**n for k,n in enumerate(mon))
                                for mon,value in s.Poly(p[i,j],*y).terms())
        homogeneous[j+1,i+1]=-homogeneous[i+1,j+1]
    normalized=model.normalize(homogeneous)
    weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    assert [s.expand(normalized['parameters'][i]-a) for i,a in zip(one,u,strict=True)]==[0]*16
    for i in range(1,5):
        for a,b in combinations_with_replacement(range(1,5),2):
            tail=s.zeros(5);tail[0,i]=X[a]*X[b];tail[i,0]=-tail[0,i]
            raw=quotient_coefficients(tail)
            assert all(raw[model.certificate['raw_free_indices'][j]]==0 for j in one+two)
    q=s.Matrix([s.expand(normalized['parameters'][i]) for i in two])
    center=q.subs(dict.fromkeys(c,0));matrix=q.jacobian(c)
    assert (q-center-matrix*s.Matrix(c)).applyfunc(s.expand)==s.zeros(45,1)
    assert not matrix.free_symbols.intersection(c)
    assert all(s.Poly(f,*u).total_degree()==1 for f in matrix if f)
    assert all(s.Poly(f,*u).total_degree()==2 for f in center if f)
    return u,c,one,two,center,matrix,constraint,syzygies,rows,columns,reduction


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();start=time.time()
    model=S61();u,c,one,two,center,matrix,constraint,syzygies,hrows,hcols,reduction=derive(model)
    section_paths=[ROOT/f'reports/{name}_UNIVERSAL_SECTION_20261101.json' for name in ('Q51','Q52')]
    q51,q52=[s.Matrix(list(map(s.sympify,json.loads(path.read_text())['quadratic_section_q']))) for path in section_paths]
    assert (3*center-2*q51-q52).applyfunc(s.expand)==s.zeros(45,1)
    source=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';five=json.loads(source.read_text())
    toric=s.Matrix(list(map(s.sympify,five['toric_S61_coordinates'])))
    substitution=dict(zip(u,[toric[i] for i in one],strict=True))
    cbase=center.subs(substitution);nbase=matrix.subs(substitution);lbase=constraint.subs(substitution)
    columns=list(lbase.rref()[1]);rows=list(lbase[:,columns].T.rref()[1]);assert len(columns)==len(rows)==18
    free=[j for j in range(20) if j not in columns];assert free==[18,19]
    block=lbase[rows,columns];kernel=s.zeros(20,2)
    solved=-block.inv()*lbase[rows,free]
    for row,index in enumerate(columns):kernel[index,:]=solved[row,:]
    for row,index in enumerate(free):kernel[index,row]=1
    assert lbase*kernel==s.zeros(51,2)
    toric_directions=nbase*kernel;elliptic_direction=(q52-center).subs(substitution)
    directions=toric_directions.row_join(elliptic_direction)
    assert toric_directions.rank()==2 and directions.rank()==3
    direction_rows=list(directions.T.rref()[1]);direction_minor=directions[direction_rows,:].det();assert direction_minor
    # Compare the entire toric cubic algebra with the affine-linear formula.
    algebra_path=ROOT/'reports/TORIC_ALGEBRA_DEGREE_20261031.json';algebra=json.loads(algebra_path.read_text())
    t=s.Symbol('t');csub={c[18]:s.Rational(13,3430)*t*t-t/7-s.Rational(14,3),c[19]:t}
    cvalues=s.Matrix([s.sympify(f).subs(csub,simultaneous=True) for f in algebra['eliminated_c0_through_c19']])
    assert (cbase+nbase*cvalues-s.Matrix([toric[i] for i in two])).applyfunc(s.expand)==s.zeros(45,1)
    assert (lbase*cvalues).applyfunc(s.expand)==s.zeros(51,1)
    # Every degree-three monomial in the plane variables is in the ideal of
    # the quadratic leading forms at the sample. This supplies relative
    # finiteness after localization in beta; no full-fiber finiteness follows.
    tau,rho,sigma=s.symbols('tau rho sigma');variables=(tau,rho,sigma)
    A=s.Rational(24010,13)
    quadrics=[tau*sigma,tau*(rho-1),rho*(rho-1),sigma*(3*rho+2*sigma-2),tau*tau-A*rho*(1+2*sigma)]
    tops=[sum(co*s.prod(v**n for v,n in zip(variables,mon,strict=True))
              for mon,co in s.Poly(f,*variables).terms() if sum(mon)==2) for f in quadrics]
    cubics=list(combinations_with_replacement(range(3),3));cubic_monomials=[s.prod(variables[i] for i in mon) for mon in cubics]
    macaulay=s.Matrix([[s.Poly(s.expand(v*f),*variables).coeff_monomial(mon) for mon in cubic_monomials]
                      for f in tops for v in variables])
    cubic_rows=list(macaulay.T.rref()[1]);assert len(cubic_rows)==10
    cubic_minor=macaulay[cubic_rows,:].det();assert cubic_minor
    report={'id':'TORIC_UNIVERSAL_THREE_PLANE_20261101','field':'QQ[all16 first-jet variables]','all_passed':True,
        'weight1_u_indices_zero_based':one,'weight2_u_indices_zero_based':two,'beta_variables':list(map(str,u)),'c_variables':list(map(str,c)),
        'leading_associativity_rank':45,'quadratic_frame_symbol_rank':45,'H_pivot_rows':hrows,'H_pivot_columns':hcols,
        'quadratic_frame_symbol_reduction_matrix':encode(reduction),
        'all51_homogeneous_linear_consequences':encode(constraint),
        'all51_constant_combinations_of96':[[[i,str(v)] for i,v in enumerate(relation) if v] for relation in syzygies],
        'q_center':list(map(str,center)),'q_linear_matrix_N':encode(matrix),
        'universal_center_identity':'3*q_center = 2*Q51 + Q52','universal_center_identity_checked':True,
        'cubic_tail_basis_checks':40,
        'generic_plane_recipe':{'L_rows':rows,'L_pivot_columns':columns,'L_free_columns':free,
            'formula':'On det(L[rows,pivots]) != 0, let T_free=I2 and T_pivots=-L_pivots^(-1)*L_free. Plane q=q_center+N*T*(r,s)^t+(Q52-q_center)*sigma.',
            'additional_open':'The displayed three direction vectors have rank3.',
            'sample_linear_minor':str(block.det()),'sample_kernel_T':encode(kernel),
            'sample_plane_direction_rows':direction_rows,'sample_plane_direction_minor':str(direction_minor)},
        'sample_toric_cubic_algebra_linear_reconstruction_checked':True,
        'sample_plane_finiteness':{'variables':list(map(str,variables)),'five_quadrics':list(map(str,quadrics)),
            'degree3_monomials':list(map(str,cubic_monomials)),'leading_cubic_macaulay_matrix':encode(macaulay),
            'pivot_rows':cubic_rows,'minor_determinant':str(cubic_minor),'rank':10,'sample_intersection_length':5},
        'generic_intersection_length_with_this_plane':5,
        'generic_intersection_justification':'Universal plane contains the known3+1+1 lifts. The nonzero leading-cubic minor gives a finite plane intersection near the saved beta; its reduced length5 sample yields an upper bound5 on a smaller open by Nakayama. Known generic5 distinct lifts attain the bound.',
        'limitations':['No full S61 fiber containment in the three-plane is proved.','No full generic affine fiber upper bound follows.','The center identity is not asserted to be the trace of the three toric lifts without a separate trace calculation.'],
        'proof_note':'research/notes/TORIC_UNIVERSAL_THREE_PLANE_20261101.md'}
    inputs=[Path(__file__),ROOT/'code/toric_algebra_fiber.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/s61_model.py',*section_paths,source,algebra_path]
    report['source_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'universal_center_identity':True,'generic_plane_dimension':3,'generic_plane_intersection_length':5,'seconds':round(time.time()-start,2),'output':str(args.output)}))

if __name__=='__main__':main()
