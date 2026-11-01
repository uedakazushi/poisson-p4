#!/usr/bin/env python3
"""Finite flag-boundary normalizations and exact linear-chart exceptions."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sympy as s
from s61_model import S61, X, J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, sp_basis, infinitesimal_action

ROOT = Path(__file__).resolve().parents[1]
M = lambda a: s.Matrix(a).applyfunc(s.sympify)
enc = lambda A: [[str(v) for v in row] for row in A.tolist()]


def run():
    start = time.time()
    source = ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'
    old = json.loads(source.read_text())
    assert old['all_passed']
    for p,digest in old['input_sha256'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest,p
    m = S61(); b = s.symbols('b0:15')
    H = M(old['flag']['beta_kernel']); O = M(old['flag']['q_orbit_differential'])
    stab = M(old['frame_change']['stabilizer']); db = M(old['frame_change']['beta_derivative'])
    z,u,beta = primitive_beta(m)
    v = s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS])
    basis = v.jacobian(u); rows = list(basis.T.rref()[1]); inv = basis[rows,:].inv()
    action = s.Matrix.hstack(*[inv*infinitesimal_action(v,A)[rows,:] for A in sp_basis()])
    action = action.subs(dict(zip(u,H*s.Matrix(b)))).applyfunc(s.expand)
    assert (action*stab-H*db).applyfunc(s.expand) == s.zeros(16,6)
    def finite_pull(T):
        B = beta.subs(dict(zip(u,H*s.Matrix(b))))
        C = T.T*B.subs(dict(zip(z,T*s.Matrix(z))),simultaneous=True)*T
        return (inv*s.Matrix([C[i,j].diff(z[k]) for (i,j),k in COORDS])[rows,:]).applyfunc(s.expand)
    t = s.Symbol('t'); flows = {}
    x,y,zz,w = X[1:]; spatial = s.Matrix([x,y,zz,w])
    Q = s.zeros(4);Q[1,3] = -x*x;Q -= Q.T
    for j in [0,2,3,4]:
        A = sum((stab[k,j]*sp_basis()[k] for k in range(10)),s.zeros(4))
        g = s.eye(4)+t*A
        assert A*A == s.zeros(4) and g*J[1:,1:]*g.T == J[1:,1:]
        assert (g.inv()*Q.subs(dict(zip(spatial,g*spatial)),simultaneous=True)*g.inv().T-Q).applyfunc(s.expand) == s.zeros(4)
        D = db[:,j].jacobian(b)
        assert D**4 == s.zeros(15)
        F = sum((t**k/s.factorial(k)*D**k for k in range(4)),s.zeros(15))
        assert (finite_pull(g)-H*F*s.Matrix(b)).applyfunc(s.expand) == s.zeros(16,1)
        flows[j] = F
    a,c = s.symbols('a c',nonzero=True)
    diagonal = s.diag(a,1/a,c,1/c)
    weights = [(s.cancel(db[i,1]/b[i]),s.cancel(db[i,5]/b[i])) for i in range(15)]
    torus = s.diag(*[a**p*c**q for p,q in weights])
    assert (finite_pull(diagonal)-H*torus*s.Matrix(b)).applyfunc(s.cancel) == s.zeros(16,1)
    def move(values,j,parameter):
        return (flows[j].subs(t,parameter)*values).applyfunc(s.cancel)
    def normalized(zero12=False):
        values = s.Matrix(b).subs(b[1],0)
        if zero12: values = values.subs(b[12],0)
        stages = []
        for j,choose in ([(4,lambda v:-v[3]/v[0]),(2,lambda v:-2*v[5]/v[0]),(3,lambda v:2*v[13]/(5*v[2])),(0,lambda v:-v[6]/v[0])] if zero12 else [(4,lambda v:-v[3]/v[0]),(2,lambda v:(v[14]-v[5])/(3*v[0])),(3,lambda v:v[5]/v[12]),(0,lambda v:-v[6]/v[0])]):
            parameter = s.cancel(choose(values))
            values = move(values,j,parameter)
            stages.append(dict(column=j,parameter=str(parameter)))
        return values,stages
    regular,stages = normalized()
    I = b[0]*b[2]-b[3]*b[12]+b[3]**2
    assert all(regular[j]==0 for j in [1,3,5,6,14])
    assert regular[0] == b[0] and regular[12] == b[12]
    assert s.cancel(regular[2]-I/b[0]) == 0
    r = s.Symbol('r')
    invariant = I-r*b[12]**2
    dinvariant = (s.Matrix([[s.diff(invariant,bb) for bb in b]])*db).subs(b[1],0).applyfunc(s.expand)
    assert dinvariant == s.Matrix([[0,-2*invariant,0,0,0,0]])
    zero,zero_stages = normalized(True)
    assert all(zero[j]==0 for j in [1,3,5,6,12,13])
    assert zero[0] == b[0] and s.cancel(zero[2]-(b[0]*b[2]+b[3]**2)/b[0]) == 0
    # Torus: a=b12,c^2=b0/b12 on the regular chart; a^2=I,c^2=b0/a on b12=0.
    assert s.cancel((torus*regular)[0]-b[0]/(a*c*c)) == 0
    assert s.cancel((torus*regular)[12]-b[12]/a) == 0
    assert s.cancel((torus*regular)[2]-c*c*I/(a*b[0])) == 0

    normalized_sub = {b[0]:1,b[1]:0,b[3]:0,b[5]:0,b[6]:0,b[12]:1,b[14]:0}
    cores = list(map(s.Rational,[-6,'-12/5','-119/20','-4/25']))
    image_cases = []
    for rr in cores:
        equation = I-rr*b[12]**2
        ann = s.zeros(3,16);ann[0,0] = ann[1,2] = 1
        for j,bb in enumerate(b): ann[2,j+1] = s.diff(equation,bb)
        sub = dict(normalized_sub);sub[b[2]] = rr
        N = (ann*action).subs(sub).applyfunc(s.factor)
        assert N.to_DM().rank() == 3
        sample = {bb:j+2 for j,bb in enumerate(b)};sample.update(sub)
        an = ann.subs(sample);tangent = s.Matrix.hstack(*an.nullspace())
        differential = tangent.row_join(action.subs(sample))
        joint = differential.col_join(s.zeros(44,13).row_join(O))
        assert differential.rank() == 16 and joint.rank() == 17
        assert N.subs(sample)[:,[4,5,6]].det() != 0
        image_cases.append(dict(core=str(rr),normal_matrix=enc(N),normal_minor=str(N.subs(sample)[:,[4,5,6]].det()),beta_rank=16,joint_rank=17))
    ann = s.zeros(3,16);ann[0,0] = ann[1,2] = ann[2,13] = 1
    N = (ann*action).subs({b[1]:0,b[12]:0}).applyfunc(s.factor)
    assert N.to_DM().rank() == 3
    sample = {bb:j+2 for j,bb in enumerate(b)};sample.update({b[1]:0,b[12]:0})
    tangent = s.Matrix.hstack(*ann.nullspace())
    differential = tangent.row_join(action.subs(sample))
    assert differential.rank() == 16 and differential.col_join(s.zeros(44,13).row_join(O)).rank() == 17
    zero_minor = N.subs(sample)[:,[4,5,6]].det()
    assert zero_minor == -4795
    ann4 = ann.col_join(s.zeros(1,16))
    for j,bb in enumerate(b):ann4[3,j+1] = s.diff(I,bb)
    exceptional_sub = {b[1]:0,b[12]:0,b[2]:-b[3]**2/b[0]}
    rank3 = (ann4*action).subs(exceptional_sub).applyfunc(s.factor)
    assert ann4.subs(exceptional_sub).to_DM().rank() == 4
    assert rank3.to_DM().rank() == 3 and all(rank3[:,j]==s.zeros(4,1) for j in range(10) if j not in [4,5,6])

    # Exact all-parameter ranks directly from the original134 plus G6 equations.
    one = [i for i,wg in enumerate(m.reduced['parameter_weights']) if wg==1]
    two = [i for i,wg in enumerate(m.reduced['parameter_weights']) if wg==2]
    ui,qi = {i:j for j,i in enumerate(one)},{i:j for j,i in enumerate(two)}
    weight6 = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    equations = m.reduced['equations']+[json.loads(weight6.read_text())['universal_new_equation']]
    q0 = M(old['flag']['point61'])[two,:]
    fixed = two.index(37);chart = [i for i in range(45) if i!=fixed]
    beta_normal = (H*s.Matrix(b)).subs(normalized_sub)
    rawJ = s.zeros(135,45);cache = {}
    for row,eq in enumerate(equations):
        for mon,aa,bb in eq:
            qs = tuple(qi[i] for i in mon if i in qi)
            if fixed not in qs:continue
            us = tuple(ui[i] for i in mon if i in ui)
            if us not in cache:cache[us] = s.prod(beta_normal[j] for j in us)
            co = s.Rational(aa,bb)*cache[us]
            if len(qs)==2:
                i,j = qs
                if i!=fixed:rawJ[row,chart.index(i)] += co*q0[j]
                if j!=fixed:rawJ[row,chart.index(j)] += co*q0[i]
            elif len(qs)==1:rawJ[row,44] += co*q0[fixed]
    core_ranks = []
    for rr in cores:
        A = rawJ.subs(b[2],rr)
        _,cols = A.to_DM().rref();_,rows = A[:,list(cols)].T.to_DM().rref()
        block = A[list(rows),list(cols)].to_DM()
        det = s.factor(block.domain.to_sympy(block.det()))
        leading = A[:,:44].to_DM().rank()
        assert det.is_Rational and det != 0
        assert (len(cols),leading) == ((35,34) if rr==-s.Rational(4,25) else (36,35))
        core_ranks.append(dict(core=str(rr),full_rank=len(cols),leading_rank=leading,minor=str(det),pivot_rows=list(rows),pivot_columns=list(cols)))
    paths = [Path(__file__),source,ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',weight6,ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101',all_passed=True,regular_chart=dict(condition='b0*b12 !=0, b1=0',unipotent_stages=stages,before_torus=list(map(str,regular)),torus='a=b12; c^2=b0/b12',finite_etale_cover_degree=2,core='(b0*b2-b3*b12+b3^2)/b12^2',remaining_variables=['b2','b4','b7','b8','b9','b10','b11','b13']),core_image_cases=image_cases,core_linear_ranks=core_ranks,zero12_chart=dict(condition='b0*(b0*b2+b3^2) !=0, b1=b12=0',unipotent_stages=zero_stages,before_torus=list(map(str,zero)),torus='a^2=b0*b2+b3^2; c^2=b0/a',finite_etale_cover_degree=4,remaining_variables=['b4','b7','b8','b9','b10','b11','b14'],beta_rank=16,joint_rank=17,normal_matrix=enc(N),normal_minor=str(zero_minor)),zero12_exception=dict(condition='b0 !=0, b1=b12=b0*b2+b3^2=0',source_dimension=12,normal_rank=3,beta_image_dimension_bound=15,normal_matrix=enc(rank3)),limitations=['This is a normalization and linear-rank theorem, not t-saturation removal.','The four core loci and b12=0 have dominant beta images.','Only b0!=0 inside the b1=0 boundary is treated; b0=b1=0 remains separate.'],proof_note='research/notes/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=run()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
