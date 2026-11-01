#!/usr/bin/env python3
"""Independent affine-tensor and invariant audit of the four flag endpoint strata."""
import argparse, hashlib, json, time
from pathlib import Path
from itertools import combinations, combinations_with_replacement
import sympy as s
from s61_model import S61, J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, sp_basis, infinitesimal_action

ROOT=Path(__file__).resolve().parents[1]


def run():
    start=time.time();table_path=ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json'
    invariant_path=ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101.json'
    table=json.loads(table_path.read_text());invariants=json.loads(invariant_path.read_text())
    for report in [table,invariants]:
        assert report['all_passed']
        for path,h in report['input_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    model=S61();z,u,beta=primitive_beta(model);x,y,zz,w=z;v=s.Matrix(z)
    tensor=s.Matrix([beta[i,j].diff(z[k])for(i,j),k in COORDS]);frame=tensor.jacobian(u)
    frame_rows=list(frame.T.to_DM().rref()[1]);inverse=frame[frame_rows,:].inv()
    action=s.Matrix.hstack(*[inverse*infinitesimal_action(tensor,A)[frame_rows,:]for A in sp_basis()])
    Q=s.zeros(4);Q[1,3]=-x*x;Q-=Q.T
    pairs=list(combinations(range(4),2));monomials=[s.prod(z[j]for j in m)for m in combinations_with_replacement(range(4),2)]
    def coefficients(P):
        return s.Matrix([s.expand(P[i,j]).coeff(m)for i,j in pairs for m in monomials])
    # Use all60 ordinary affine coefficients, independently of the S61 coordinate adapter.
    q=coefficients(Q);fixed=next(i for i,c in enumerate(q)if c)
    orbit=[]
    for A in sp_basis():
        dQ=sum(((A*v)[i]*Q.diff(z[i])for i in range(4)),s.zeros(4))-A*Q-Q*A.T
        d=coefficients(dQ);orbit.append(d-q*d[fixed]/q[fixed])
    O=s.Matrix.hstack(*orbit);stabilizer=s.Matrix.hstack(*O.nullspace())
    assert O.rank()==4 and stabilizer.cols==6
    dimensions={}
    for label,record in table['strata'].items():
        H=s.Matrix(record['beta_kernel']).applyfunc(s.sympify)
        source=s.Matrix(record['beta']).applyfunc(s.sympify)
        ann=s.Matrix.hstack(*H.T.nullspace()).T
        A=action.subs(dict(zip(u,source)))
        normal=ann*A*stabilizer;nr=normal.to_DM().rank()
        expected=2 if label=='pure_cubic'else 1
        assert nr==expected and 4+H.cols+nr==16
        sample={s.Symbol(k):s.sympify(value)for k,value in record['sample'].items()}
        D=H.row_join(A);joint=D.col_join(s.zeros(60,H.cols).row_join(O))
        assert D.subs(sample).rank()==joint.subs(sample).rank()==16
        dimensions[label]=dict(source_dimension=H.cols,normal_rank=nr,beta_rank=16,joint_rank=16)
    k,K,L,U,V=s.symbols('k K L U V');n=s.symbols('n0:6')
    source=s.Matrix([0,0,0,0,k,n[0],U,n[1],n[2],V,n[3],n[4],n[5],2*k+1,K-5*V,5*U-L])
    source_tensor=beta.subs(dict(zip(u,source)))
    checked=[]
    for record in table['presentations']:
        T=s.Matrix(record['matrix']).applyfunc(s.sympify);kv=s.Rational(record['k'])
        assert (T*J[1:,1:]*T.T-J[1:,1:]).applyfunc(s.expand)==s.zeros(4)
        G=s.sympify(record.get('G_infinity',record.get('residue_G')),
                    locals={f'x{i+1}':z[i]for i in range(4)})
        P=s.zeros(4);P[1,2]=-s.diff(G,w);P[1,3]=s.diff(G,zz);P-=P.T
        transformed=T.inv()*P.subs(dict(zip(z,T*v)),simultaneous=True)*T.inv().T
        assert (transformed-s.sympify(record['projective_factor'])*Q).applyfunc(s.factor)==s.zeros(4)
        C=T.T*source_tensor.subs(k,kv).subs(dict(zip(z,T*v)),simultaneous=True)*T
        transformed_beta=(inverse*s.Matrix([C[i,j].diff(z[h])for(i,j),h in COORDS])[frame_rows,:]).applyfunc(s.expand)
        assert (transformed_beta-s.Matrix(record['normalized_beta']).applyfunc(s.sympify)).applyfunc(s.expand)==s.zeros(16,1)
        ann=s.Matrix(table['strata'][record['target']]['annihilator']).applyfunc(s.sympify)
        assert (ann*transformed_beta).applyfunc(s.expand)==s.zeros(ann.rows,1)
        checked.append([record['kind'],record['source'],record['axis'],record['k']])
    assert len(checked)==12
    b=s.symbols('b0:15');bb=s.Matrix([0,*b]);Db=(action*stabilizer).subs(dict(zip(u,bb)))
    assert Db[0,:]==s.zeros(1,6);Db=Db[1:,:]
    ratios=[-s.S.One,-s.Integer(4),s.Rational(7,5)]
    for r in ratios:
        f=b[3]-r*b[12];der=s.Matrix([[s.diff(f,a)for a in b]])*Db
        assert der.subs({b[0]:0,b[1]:0,b[3]:r*b[12]}).applyfunc(s.expand)==s.zeros(1,6)
    F=s.sympify(invariants['pure_cubic']['quartic']);delta=5*b[3]-4*b[12]
    pure={b[0]:0,b[3]:s.Rational(4,5)*b[12]}
    assert Db[0,:].subs(pure).applyfunc(s.expand)==s.zeros(1,6)
    assert (s.Matrix([[s.diff(delta,a)for a in b]])*Db).subs(pure).applyfunc(s.expand)==s.zeros(1,6)
    derivatives=(s.Matrix([[s.diff(F,a)for a in b]])*Db).subs(pure)
    weights=[s.cancel(d/F)for d in derivatives];assert all(c.is_Rational for c in weights)
    t=s.Symbol('t');P=b[1]*t*t-s.Rational(6,5)*b[12]*t+b[14]-3*b[5]
    R=b[6]-(b[14]+b[5])*t+s.Rational(9,5)*b[12]*t*t-b[1]*t**3
    assert s.expand(125*s.resultant(P,R,t)-b[1]*F)==0
    disc=s.discriminant(F,b[6]);last=-25*b[1]*b[14]+75*b[1]*b[5]+9*b[12]**2
    assert s.expand(disc-16*(50*b[1]*b[5]-9*b[12]**2)**2*last)==0
    assert s.gcd_list(s.Poly(F,b[6]).all_coeffs())==1
    assert s.gcd_list(s.Poly(last,b[14]).all_coeffs())==1 and s.degree(last,b[14])==1
    paths=[Path(__file__),table_path,invariant_path,ROOT/'research/notes/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.md',
           ROOT/'research/notes/S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101.md',ROOT/'code/s61_model.py',
           ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_MINIMAL_FLAG_ENDPOINT_REVIEW_20261101',all_passed=True,
                independent_projective_orbit_dimension=4,independent_stabilizer_dimension=6,
                strata_dimensions=dimensions,twelve_finite_tensor_checks=checked,
                pure_invariant_weights_in_independent_basis=list(map(str,weights)),
                resultant_discriminant_content_verified=True,
                accepted=['Finite-source one-pole-coordinate DVR argument, including zero residue numerator.',
                          'Connected Borel invariance and geometrically irreducible16D incidence closure comparison.',
                          'Four distinct dominant16D strata; no t-saturation conclusion.'],
                scope=['Independent60-coefficient affine action, all12 saved finite tensor transformations, source dimensions and invariant algebra.',
                       'Archived normal-form classification and source residue formulas are inputs, not reclassified here.',
                       'Point Hessians and local completions are not independently recomputed.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
                seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=run()
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
