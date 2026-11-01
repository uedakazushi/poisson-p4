#!/usr/bin/env python3
"""Exact five known lifts and their three-dimensional affine span.

The full S61 intersection with this explicitly imposed three-plane has length
five and is reduced. This does NOT prove that the whole45-variable fiber lies
in the plane or exclude additional components/points outside it.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61,X
from toric_algebra_fiber import multiplication_from_jet,trace_form,TRIPLES
from weight1_structure import J
ROOT=Path(__file__).resolve().parents[1]


def toric_coordinates(model):
    algebra=json.loads((ROOT/'reports/TORIC_ALGEBRA_FIBER_20261031.json').read_text())
    degree=json.loads((ROOT/'reports/TORIC_ALGEBRA_DEGREE_20261031.json').read_text())
    t=s.Symbol('t');cubic=t**3-s.Rational(24010,13)*t;red=lambda f:s.rem(s.expand(f),cubic,t)
    c=s.symbols('c0:20')
    substitution={c[18]:s.Rational(13,3430)*t*t-t/7-s.Rational(14,3),c[19]:t}
    solution=[s.sympify(f).subs(substitution,simultaneous=True) for f in degree['eliminated_c0_through_c19']]
    matrices=multiplication_from_jet(s.Matrix(algebra['fixed_first_jet']),dict(zip(TRIPLES,solution,strict=True)))
    matrices=[m.applyfunc(red) for m in matrices];form=trace_form(matrices).applyfunc(red)
    y=s.Matrix(s.symbols('y0:4'))
    frame=s.eye(4)+sum((y[i]*matrices[i] for i in range(4)),s.zeros(4))-y*(y.T*form)
    p=(frame*J*frame.T).applyfunc(red)
    cubic_part=p.applyfunc(lambda f:sum(co*s.prod(y[i]**n for i,n in enumerate(mon))
                          for mon,co in s.Poly(f,*y).terms() if sum(mon)==3))
    vector=s.Matrix([sum(s.diff(cubic_part[j,i],y[j]) for j in range(4))/5 for i in range(4)])
    assert (cubic_part-y*vector.T+vector*y.T).applyfunc(red)==s.zeros(4)
    homogeneous=s.zeros(5)
    for i in range(4):
        homogeneous[0,i+1]=-vector[i].subs(dict(zip(y,X[1:])),simultaneous=True)
        homogeneous[i+1,0]=-homogeneous[0,i+1]
        for j in range(i+1,4):
            value=sum(co*X[0]**(2-sum(mon))*s.prod(X[k+1]**n for k,n in enumerate(mon))
                      for mon,co in s.Poly(p[i,j]-cubic_part[i,j],*y).terms())
            homogeneous[i+1,j+1]=value;homogeneous[j+1,i+1]=-value
    normalized=model.normalize(homogeneous)
    return t,cubic,[red(v) for v in normalized['parameters']]


def verify_mod_polynomial(equations,values,variable,modulus):
    base,x=ring(str(variable),QQ);mod=base.from_expr(modulus);val=[base.from_expr(v) for v in values]
    for equation in equations:
        total=base.zero
        for mon,n,d in equation:
            term=base(QQ(n,d))
            for i in mon:term=(term*val[i]).rem(mod)
            total+=term
        assert not total.rem(mod)


def modular_jacobian_witness(model,point,prime):
    columns=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2];index={old:new for new,old in enumerate(columns)}
    def mod(value):
        value=s.Rational(value);assert int(value.q)%prime
        return int(value.p)*pow(int(value.q),-1,prime)%prime
    values=[mod(v) for v in point];matrix=[]
    for equation in model.reduced['equations']:
        row=[0]*45
        for mon,n,d in equation:
            for position,old in enumerate(mon):
                if old not in index:continue
                value=n*pow(d,-1,prime)%prime
                for j,k in enumerate(mon):
                    if j!=position:value=value*values[k]%prime
                row[index[old]]=(row[index[old]]+value)%prime
        matrix.append(row)
    basis={};selected=[];determinant=1
    for i,source in enumerate(matrix):
        row=list(source)
        for pivot,old in sorted(basis.items()):
            factor=row[pivot]
            if factor:row=[(x-factor*y)%prime for x,y in zip(row,old,strict=True)]
        pivot=next((j for j,v in enumerate(row) if v),None)
        if pivot is None:continue
        coefficient=row[pivot];determinant=determinant*coefficient%prime
        inverse=pow(coefficient,-1,prime);basis[pivot]=[v*inverse%prime for v in row];selected.append(i)
        if len(basis)==45:break
    assert len(basis)==45
    # Independent determinant on the actual45x45 minor fixes any pivot order sign.
    work=[list(matrix[i]) for i in selected];determinant=1
    for k in range(45):
        rr=next(i for i in range(k,45) if work[i][k])
        if rr!=k:work[k],work[rr]=work[rr],work[k];determinant=-determinant
        coefficient=work[k][k];determinant=determinant*coefficient%prime;inverse=pow(coefficient,-1,prime)
        for i in range(k+1,45):
            factor=work[i][k]*inverse%prime
            if factor:work[i]=[(x-factor*y)%prime for x,y in zip(work[i],work[k],strict=True)]
    assert determinant%prime
    return {'prime':prime,'rows_zero_based':selected,'columns':'all45 weight2 coordinates',
            'minor_determinant_mod_prime':determinant%prime}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    t,cubic,toric=toric_coordinates(model);verify_mod_polynomial(model.reduced['equations'],toric,t,cubic)
    saved=next(p for p in json.loads((ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json').read_text())['points'] if p['component_id']=='T01')
    assert [v.subs(t,0) for v in toric]==list(map(s.Rational,saved['S61_coordinates_u1_through_u61']))
    assert all(not toric[i].has(t) for i in one)
    u=s.symbols('u1:62');substitution={u[i]:toric[i] for i in one};elliptic=[]
    for name in ('Q51','Q52'):
        data=json.loads((ROOT/f'reports/{name}_UNIVERSAL_SECTION_20261101.json').read_text());point=list(toric)
        for i,f in zip(two,data['quadratic_section_q'],strict=True):point[i]=s.sympify(f).subs(substitution)
        assert all(model.evaluate(eq,point)==0 for eq in model.reduced['equations']);elliptic.append(point)
    q51,q52=elliptic;z=s.Symbol('z');primitive=toric[two[0]]
    toric_polynomial=s.Poly(s.resultant(cubic,z-primitive,t),z).monic().as_expr()
    a,b=q51[two[0]],q52[two[0]];assert a!=b and toric_polynomial.subs(z,a)*toric_polynomial.subs(z,b)!=0
    quintic=s.expand(toric_polynomial*(z-a)*(z-b));assert s.gcd(quintic,s.diff(quintic,z))==1
    cyclic=s.Matrix.hstack(*(s.Matrix([s.Poly(s.rem(primitive**k,cubic,t),t).coeff_monomial(t**j)
                                     for j in range(3)]) for k in range(3)));assert cyclic.det()!=0
    inverse=cyclic.inv();interpolation=[]
    for i in two:
        c=inverse*s.Matrix([s.Poly(toric[i],t).coeff_monomial(t**j) for j in range(3)])
        f=sum(c[k]*z**k for k in range(3));aa=(q51[i]-f.subs(z,a))/toric_polynomial.subs(z,a);bb=(q52[i]-f.subs(z,b))/toric_polynomial.subs(z,b)
        interpolation.append(s.expand(f+toric_polynomial*(aa+(bb-aa)*(z-a)/(b-a))))
    assert interpolation[0]==z
    interpolated=list(toric)
    for i,f in zip(two,interpolation,strict=True):interpolated[i]=f
    verify_mod_polynomial(model.reduced['equations'],interpolated,z,quintic)
    coefficient=s.Matrix([[s.Poly(f,z).coeff_monomial(z**k) for f in [s.S.One,*interpolation]] for k in range(5)])
    pivots=list(coefficient.rref()[1]);assert pivots==[0,1,2,3]
    basis=coefficient[:,pivots];rows=list(basis.T.rref()[1]);left_inverse=basis[rows,:].inv();v=s.symbols('v0:3');plane=[]
    for j in range(45):
        c=left_inverse*coefficient[rows,j+1];assert basis*c==coefficient[:,j+1]
        plane.append(c[0]+sum(c[k+1]*v[k] for k in range(3)))
    assert plane[:3]==list(v)
    values=list(toric)
    for i,f in zip(two,plane,strict=True):values[i]=f
    restricted=[]
    for equation in model.reduced['equations']:
        restricted.append(s.expand(sum(s.Rational(n,d)*s.prod(values[j] for j in mon) for mon,n,d in equation)))
    monomials=[s.S.One,*v,*[v[i]*v[j] for i in range(3) for j in range(i,3)]]
    full=s.Matrix([[s.Poly(f,*v).coeff_monomial(mon) for mon in monomials] for f in restricted])
    selected=list(full.T.rref()[1]);assert len(selected)==5
    equations=[restricted[i] for i in selected];basis5=full[selected,:];columns=list(basis5.rref()[1]);combination=full[:,columns]*basis5[:,columns].inv()
    assert combination*basis5==full
    shape=[v[1]-interpolation[1].subs(z,v[0]),v[2]-interpolation[2].subs(z,v[0]),quintic.subs(z,v[0])]
    gb=s.groebner(equations,v[2],v[1],v[0],order='lex');known=s.groebner(shape,v[2],v[1],v[0],order='lex')
    assert all(gb.reduce(f)[1]==0 for f in shape) and all(known.reduce(f)[1]==0 for f in equations)
    # Special trace relation explains why five points have only affine span3.
    A=s.Rational(24010,13)
    trace=[3*s.Poly(f,t).coeff_monomial(1)+2*A*s.Poly(f,t).coeff_monomial(t*t) for f in toric]
    assert s.Matrix(trace)==2*s.Matrix(q51)+s.Matrix(q52)
    assert s.Poly(t*t-A,t).is_irreducible
    geometric=s.symbols('tau rho sigma');tau,rho,sigma=geometric
    origin=s.Matrix([f.subs(t,0) for f in toric]);linear=s.Matrix([s.Poly(f,t).coeff_monomial(t) for f in toric]);quadratic=s.Matrix([s.Poly(f,t).coeff_monomial(t*t) for f in toric])
    embedding=origin+tau*linear+A*rho*quadratic+sigma*(s.Matrix(q52)-origin)
    substitution=dict(zip(v,[embedding[i] for i in two[:3]],strict=True))
    assert all(s.expand(plane[j].subs(substitution,simultaneous=True)-embedding[i])==0 for j,i in enumerate(two))
    simple=[tau*sigma,tau*(rho-1),rho*(rho-1),sigma*(3*rho+2*sigma-2),tau*tau-A*rho*(1+2*sigma)]
    geo_monomials=[s.S.One,*geometric,*[geometric[i]*geometric[j] for i in range(3) for j in range(i,3)]]
    simple_matrix=s.Matrix([[s.Poly(f,*geometric).coeff_monomial(mon) for mon in geo_monomials] for f in simple])
    transformed=s.Matrix([[s.Poly(s.expand(f.subs(substitution,simultaneous=True)),*geometric).coeff_monomial(mon) for mon in geo_monomials] for f in equations])
    selected_columns=list(simple_matrix.rref()[1]);simple_combination=transformed[:,selected_columns]*simple_matrix[:,selected_columns].inv()
    assert simple_combination*simple_matrix==transformed and simple_combination.det()!=0
    prime=32003;Amod=int(A.p)*pow(int(A.q),-1,prime)%prime;root=int(s.sqrt_mod(Amod,prime))
    conjugate_point=[s.Poly(f,t).eval(root) for f in toric]
    conjugate_witness=modular_jacobian_witness(model,conjugate_point,prime)
    q51_witness=modular_jacobian_witness(model,q51,prime)
    sources=['TORIC_ALGEBRA_FIBER_20261031.json','TORIC_ALGEBRA_DEGREE_20261031.json','Q51_UNIVERSAL_SECTION_20261101.json','Q52_UNIVERSAL_SECTION_20261101.json','S61_KNOWN_FIBER_AUDIT_20261031.json']
    report={'id':'S61_KNOWN_FIVE_POINT_SECTION_20261101','all_passed':True,'field':'QQ',
        'q_to_u_indices_zero_based':two,'toric_cubic':str(cubic),'toric_S61_coordinates':list(map(str,toric)),
        'q51_S61_coordinates':list(map(str,q51)),'q52_S61_coordinates':list(map(str,q52)),
        'primitive_q_coordinate_old_index_one_based':two[0]+1,'primitive_toric_polynomial':str(toric_polynomial),
        'q51_primitive_value':str(a),'q52_primitive_value':str(b),'known_fiber_quintic':str(quintic),
        'quintic_discriminant':str(s.discriminant(quintic,z)),
        'q_interpolation_polynomials':list(map(str,interpolation)),
        'affine_span_dimension':3,'affine_plane_free_q_indices':[0,1,2],
        'affine_plane_substitution':list(map(str,plane)),
        'restricted_equation_indices':selected,'restricted_five_quadrics':list(map(str,equations)),
        'all134_restricted_span_coefficients':[[str(x) for x in row] for row in combination.tolist()],
        'lex_groebner_basis':list(map(str,[p.as_expr() for p in gb.polys])),
        'intersection_length':5,'intersection_reduced':True,
        'geometric_plane_embedding_S61_coordinates':list(map(str,embedding)),
        'geometric_variables':list(map(str,geometric)),'five_simple_quadrics':list(map(str,simple)),
        'all134_geometric_quadric_span_coefficients':[[str(x) for x in row] for row in (combination*simple_combination).tolist()],
        'trace_identity':'sum(three toric S61 points)=2*q51+q52 at this fixed jet',
        'q51_full_fiber_jacobian_witness':q51_witness,
        'toric_quadratic_conjugate_witness':conjugate_witness|{'t_mod_prime':root,'irreducible_QQ_factor':str(t*t-A),
            'reason':'A rational minor polynomial nonzero after this modular specialization cannot vanish modulo the irreducible quadratic factor overQQ; hence both characteristic-zero conjugate points have rank45.'},
        'scope':['This is the full S61 intersection with the imposed affine3-plane.','Membership of its42 affine-linear equations in the fullfiber radical is unproved.','No additional points or positive-dimensional components outside the plane are excluded.','The displayed trace identity has only been checked at this fixed first jet.'],
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [Path(__file__),ROOT/'code/s61_model.py',*[ROOT/'reports'/name for name in sources]]}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'affine_span':3,'intersection_length':5,'not_full_fiber_exhaustion':True,'output':str(args.output)}))

if __name__=='__main__':main()
