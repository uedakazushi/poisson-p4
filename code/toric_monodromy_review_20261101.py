#!/usr/bin/env python3
"""Independent arithmetic audit of the toric S3 certificate.

No model/monodromy modules are imported. Reconstructs the line multiplication
from the old audited affine20 matrices and a direct beta tensor correction;
checks all96 commutators against the reported rank3 algebra, then its trace
discriminant, cyclic cubic and coprimality certificates.
"""
import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    input_path=ROOT/'reports/TORIC_MONODROMY_20261101.json'
    prior_path=ROOT/'reports/TORIC_ALGEBRA_FIBER_20261031.json'
    data=json.loads(input_path.read_text());prior=json.loads(prior_path.read_text())
    t,z=s.symbols('t z');cs=s.symbols('c0:20');x0,y0=s.symbols('x y')
    K=QQ.frac_field(t);R,x,y=ring('x,y',K)
    parse=lambda v:s.sympify(v,locals={'t':t,'z':z,'x':x0,'y':y0})
    kr=lambda v:K.from_sympy(parse(v))
    pr=lambda v:R.from_expr(parse(v))
    lin=[pr(v) for v in data['linear_solutions_c0_through_c19']]
    assert lin[-2:]==[x,y]
    old=[s.Matrix(block).applyfunc(s.sympify) for block in prior['multiplication_matrices']]
    J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    # beta(t)-beta*=t*y1*dy1 wedge dy3: only beta_0[0,2] changes.
    db=[s.zeros(4) for _ in range(4)];db[0][0,2]=1;db[0][2,0]=-1
    delta=[s.Matrix.hstack(*(J*s.Matrix([(db[i][j,k]+db[j][i,k])/3 for k in range(4)])
                             for j in range(4))) for i in range(4)]
    assert [J*dm+dm.T*J for dm in delta]==db
    assert all(s.trace(dm)==0 for dm in delta)
    assert all(delta[i][:,j]==delta[j][:,i] for i in range(4) for j in range(4))
    m=[]
    for i in range(4):
        block=[]
        for a in range(4):
            row=[]
            for b in range(4):
                poly=s.Poly(old[i][a,b],*cs)
                entry=R(kr(poly.coeff_monomial(1)+t*delta[i][a,b]))
                entry+=sum((R(kr(poly.coeff_monomial(c)))*value for c,value in zip(cs,lin)),R.zero)
                row.append(entry)
            block.append(row)
        m.append(block)
    q=[[sum((m[i][a][b]*m[j][b][a] for a in range(4) for b in range(4)),R.zero)/3
        for j in range(4)] for i in range(4)]
    equations=[]
    for i,j in combinations(range(4),2):
        for a in range(4):
            for b in range(4):
                equations.append(sum((m[i][a][r]*m[j][r][b]-m[j][a][r]*m[i][r][b]
                                      for r in range(4)),R.zero)
                                 +(q[j][b] if a==i else R.zero)-(q[i][b] if a==j else R.zero))
    nf=[pr(v) for v in data['nf']];quadratics=[x*x,x*y,y*y]
    assert all(max(map(sum,v.keys()),default=0)<=1 for v in nf)
    relations=[a-b for a,b in zip(quadratics,nf)]
    for f in equations:
        assert f-sum((R(f.coeff(mon))*rel for mon,rel in zip(quadratics,relations)),R.zero)==0
    leading=DomainMatrix([[f.coeff(mon) for mon in quadratics] for f in equations],(96,3),K)
    assert leading.rank()==3
    def vec(f):return [f.coeff(R.one),f.coeff(x),f.coeff(y)]
    ux,uxy,uy=map(vec,nf)
    MX=DomainMatrix([[K.zero,ux[0],uxy[0]],[K.one,ux[1],uxy[1]],[K.zero,ux[2],uxy[2]]],(3,3),K)
    MY=DomainMatrix([[K.zero,uxy[0],uy[0]],[K.zero,uxy[1],uy[1]],[K.one,uxy[2],uy[2]]],(3,3),K)
    assert MX==DomainMatrix([[kr(v) for v in row] for row in data['Mx']],(3,3),K)
    assert MY==DomainMatrix([[kr(v) for v in row] for row in data['My']],(3,3),K)
    identity=DomainMatrix.eye(3,K).to_dense()
    for lhs,row in [(MX.matmul(MX),ux),(MX.matmul(MY),uxy),(MY.matmul(MY),uy)]:
        assert lhs==identity.scalarmul(row[0])+MX.scalarmul(row[1])+MY.scalarmul(row[2])
    assert MX.matmul(MY)==MY.matmul(MX)
    unit=DomainMatrix([[K.one],[K.zero],[K.zero]],(3,1),K)
    cyclic=unit.hstack(MY.matmul(unit)).hstack(MY.matmul(MY).matmul(unit))
    assert cyclic.det()==kr(data['cyclic_determinant']) and cyclic.det()
    # Derive the cubic via the determinant, independently of charpoly().
    entries=MY.to_list()
    trace=sum(entries[i][i] for i in range(3))
    second=sum(entries[i][i]*entries[j][j]-entries[i][j]*entries[j][i] for i,j in combinations(range(3),2))
    a,b,c=-trace,second,-MY.det()
    cubic=z**3+K.to_sympy(a)*z**2+K.to_sympy(b)*z+K.to_sympy(c)
    assert s.cancel(cubic-parse(data['cubic']))==0
    disc=a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
    assert disc==kr(data['discriminant'])
    # Trace pairing gives the same discriminant square class in basis1,x,y.
    basis=[identity,MX,MY]
    tr=lambda mat:sum(mat.to_list()[i][i] for i in range(3))
    gram=DomainMatrix([[tr(u.matmul(v)) for v in basis] for u in basis],(3,3),K)
    assert disc==gram.det()*cyclic.det()**2
    D=s.Poly(parse(data['discriminant_denominator_factor']),t,domain=s.QQ)
    Q=s.Poly(parse(data['discriminant_numerator_factor']),t,domain=s.QQ)
    assert (D.degree(),Q.degree())==(10,12)
    assert disc==kr(7474295088540*Q.as_expr()**2/D.as_expr()**3)
    u,v=[s.Poly(parse(w),t,domain=s.QQ) for w in data['squarefree_denominator_bezout']]
    assert u*D+v*D.diff()==s.Poly(1,t,domain=s.QQ)
    u,v=[s.Poly(parse(w),t,domain=s.QQ) for w in data['coprime_num_den_bezout']]
    assert u*D+v*Q==s.Poly(1,t,domain=s.QQ)
    assert s.cancel(cubic.subs(t,0))==z**3-s.Rational(24010,13)*z
    assert parse(data['discriminant']).subs(t,0)==s.Rational(55365148804000,2197)
    for name in ('linear_determinant','leading_determinant','cyclic_determinant'):
        assert parse(data[name]).subs(t,0)!=0
    report={'id':'TORIC_MONODROMY_REVIEW_20261101','all_passed':True,
        'field':'QQ(t)','model_module_imports':False,'all96_commutators_independently_reconstructed':True,
        'residual_ideal_equals_three_quadratic_ideal':True,'rank3_quotient_relations_verified':True,
        'cyclic_cubic_determinant_verified':True,'trace_discriminant_change_of_basis_verified':True,
        'denominator_degree':10,'numerator_degree':12,'squarefree_and_coprime_bezout_checked':True,
        'complex_discriminant_valuation_at_each_denominator_root':-3,
        'accepted_geometric_monodromy':'S3',
        'geometric_dependencies':['Prior degree3 theorem and geometrically irreducible GL4/S5 source.',
                                  'Normality of the finite-etale base makes a rational square root of a unit regular, so it restricts to the chosen line.'],
        'scope':'T01 toric cover only; no full S61 degree, finiteness or exhaustion.',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),input_path,prior_path)}}
    with args.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({'all_passed':True,'monodromy':'S3','independent_residuals_checked':96}))

if __name__=='__main__':main()
