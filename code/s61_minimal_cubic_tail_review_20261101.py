#!/usr/bin/env python3
"""Independent all-equation and truncated-series check of the M1 cubic tail."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json'


def run():
    r=json.loads(INPUT.read_text());m=S61();un=r['universal'];pt=r['point'];ld=pt['local_data']
    matrix=lambda a:s.Matrix(a).applyfunc(s.sympify)
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2]
    index={i:j for j,i in enumerate(two)}
    source=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    def equations(beta):
        result=[]
        for eq in source:
            terms={}
            for mon,n,d in eq:
                qmon=tuple(index[i] for i in mon if i in index)
                value=s.Rational(n,d)*s.prod(beta[one.index(i)] for i in mon if i in one)
                terms.setdefault(qmon,[]).append(value)
            result.append({mon:s.expand(s.Add(*values)) for mon,values in terms.items()})
        return result
    b=s.symbols('b0:10');K=matrix(un['beta_kernel']);assert K.shape==(16,10) and K.rank()==10
    q=s.Matrix(list(map(s.Rational,un['q_point'])));cols=un['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
    beta_variables=s.symbols('beta0:16');beta_rows=[]
    for eq in source[:134]:
        terms=[]
        for mon,n,d in eq:
            qmon=[index[i] for i in mon if i in index]
            if len(qmon)!=2:continue
            coefficient=s.Rational(n,d)*s.prod(q[i] for i in qmon)
            if coefficient:terms.append(coefficient*s.prod(beta_variables[one.index(i)] for i in mon if i in one))
        value=s.expand(s.Add(*terms));row=[value.coeff(v) for v in beta_variables]
        assert s.expand(value-sum(a*v for a,v in zip(row,beta_variables)))==0
        beta_rows.append(row)
    beta_matrix=s.Matrix(beta_rows)
    assert beta_matrix.rank()==6 and beta_matrix*K==s.zeros(134,10)
    F=equations(K*s.Matrix(b));D=s.zeros(135,45)
    for row,f in enumerate(F):
        value=s.S.Zero
        for mon,c in f.items():
            if len(mon)!=2:continue
            i,j=mon;value+=c*q[i]*q[j];D[row,i]+=c*q[j];D[row,j]+=c*q[i]
        assert s.expand(value)==0
    D=D[:,cols].applyfunc(s.expand)
    assert (D-matrix(un['leading_jacobian'])).applyfunc(s.expand)==s.zeros(135,44)
    T=matrix(un['tangent']);L=matrix(un['quadratic_lift'])
    assert (D.to_DM()*T.to_DM()).is_zero_matrix
    assert T.to_DM().rank()==2
    sample={v:i+1 for i,v in enumerate(b)}
    assert D.subs(sample).rank()==42
    H=s.zeros(135,3)
    fullT=s.zeros(45,2)
    for j,i in enumerate(cols):fullT[i,:]=T[j,:]
    for row,f in enumerate(F):
        for mon,c in f.items():
            if len(mon)!=2:continue
            i,j=mon
            H[row,0]+=c*fullT[i,0]*fullT[j,0]
            H[row,1]+=c*(fullT[i,0]*fullT[j,1]+fullT[i,1]*fullT[j,0])
            H[row,2]+=c*fullT[i,1]*fullT[j,1]
    H=H.applyfunc(s.cancel)
    assert (H-matrix(un['quadratic_matrix'])).applyfunc(s.cancel)==s.zeros(135,3)
    assert (D.to_DM()*L.to_DM()-H.to_DM()).is_zero_matrix

    # The concrete full homogenization and its transverse series are rebuilt
    # directly, including lower q-degree terms multiplied by the missing t powers.
    point=list(map(s.Rational,ld['point61']));beta=[point[i] for i in one]
    F=equations(beta);qp=s.Matrix([point[i] for i in two]);cols=ld['chart_columns']
    fixed=next(i for i in range(45) if i not in cols)
    Jac=s.zeros(135,45)
    for row,f in enumerate(F):
        for mon,c in f.items():
            if len(mon)==2:
                for j,i in enumerate(cols):Jac[row,j]+=c*sum(qp[mon[1-a]] for a in range(2) if mon[a]==i)
            elif len(mon)==1:Jac[row,44]+=c*qp[mon[0]]
    assert Jac==matrix(ld['jacobian']) and Jac.rank()==43 and Jac[:,:44].rank()==42
    tangent=matrix(ld['tangent_basis']);assert Jac*tangent==s.zeros(135,2)
    assert tangent.rank()==2
    normal=pt['normal_kernel_coordinate'];curve_coordinate=pt['parameter_kernel_coordinate']
    assert normal!=curve_coordinate
    series=[s.zeros(45,1),tangent[:,normal],s.Matrix(list(map(s.Rational,pt['transverse_expansion'][0]['lift']))),
            s.Matrix(list(map(s.Rational,pt['transverse_expansion'][1]['lift'])))]
    rows=pt['jacobian_pivot_rows'];piv=pt['jacobian_pivot_columns']
    determinant=Jac.extract(rows,piv).det(method='domain-ge');assert determinant!=0
    free=[i for i in range(45) if i not in piv]
    assert all(series[n][i]==0 for n in (2,3) for i in free)
    def multiply(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(4)]
    coordinates=[[qp[i],0,0,0] for i in range(45)]
    for j,i in enumerate(cols):
        for n in (1,2,3):coordinates[i][n]=series[n][j]
    tseries=[0,*[series[n][-1] for n in (1,2,3)]]
    residual=s.zeros(135,4)
    for row,f in enumerate(F):
        for mon,c in f.items():
            term=[c,0,0,0]
            for i in mon:term=multiply(term,coordinates[i])
            for _ in range(2-len(mon)):term=multiply(term,tseries)
            for n in range(4):residual[row,n]+=term[n]
    assert residual[:,:3]==s.zeros(135,3)
    assert residual[:,3]==s.Matrix(list(map(s.Rational,pt['transverse_expansion'][1]['residual'])))
    assert residual[66,3]==s.Rational(84035,389017)
    assert residual[rows,3]==s.zeros(43,1)

    z=s.Symbol('z');curve=s.Matrix(list(map(s.sympify,pt['tail_q_coordinates'])))
    assert curve.subs(z,1)==qp
    assert all(s.expand(sum(c*s.prod(curve[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
    derivative=curve.diff(z).subs(z,1)-qp*curve.diff(z).subs(z,1)[fixed]/qp[fixed]
    curve_tangent=derivative[cols,:].col_join(s.zeros(1,1))
    curve_coordinates=s.Matrix(list(map(s.Rational,pt['tail_tangent_in_kernel'])))
    assert tangent*curve_coordinates==curve_tangent and curve_coordinates[curve_coordinate]!=0

    # Finite normalization, not just infinitesimal covariance.
    a,bb,d=s.symbols('a b d');x,y,zz,w=X[1:];V=s.Matrix([x,y,zz,w])
    def tensor(aa,bbb,dd):
        Q=s.zeros(4);Q[0,1]=x*x/2;Q[1,2]=-x*zz-bbb*x*x
        Q[1,3]=aa*x*x+2*dd*x*zz;Q[2,3]=x*x/2
        return Q-Q.T
    change=matrix(r['normalization']['matrix'])
    assert (change*J[1:,1:]*change.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
    normalized=change.inv()*tensor(a,bb,d).subs(dict(zip(V,change*V)),simultaneous=True)*change.inv().T
    assert (normalized-tensor(1,0,0)).applyfunc(s.cancel)==s.zeros(4)
    assert (a-2*bb*d).subs({a:3+13*z,bb:3,d:z}).expand()==3+7*z
    paths=[Path(__file__),INPUT,ROOT/'code/s61_minimal_cubic_tail_20261101.py',
           ROOT/'research/notes/S61_MINIMAL_CUBIC_TAIL_20261101.md',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
           m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_CUBIC_TAIL_REVIEW_20261101',all_passed=True,
        canonical_beta_equation_rank=6,universal_leading_rank=42,universal_quadratic_lift_verified=True,
        full_homogenized_rank=43,transverse_cubic_coefficient=str(residual[66,3]),
        full_Jacobian_minor=str(determinant),whole_curve_leading_zero=True,
        finite_symplectic_normalization_verified=True,
        accepted_conclusion='At the specified point, completion C[[u,e]]/(e^3) and t=0; extension only to the justified open.',
        scope='Independent all135 coefficient extraction and direct order3 substitution; finite-free/unique-branch proof reviewed separately in the note.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    args=p.parse_args();result=run()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'output':str(args.output)}))
