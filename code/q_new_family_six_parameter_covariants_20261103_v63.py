#!/usr/bin/env python3
"""Exact six-parameter actual NEW25 jet and compact N12/m circuit; no evaluations."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J, COORDS

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63'


def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(M):return [[str(x) for x in row] for row in M.tolist()]


class DAG:
    def __init__(self,names):
        self.names,self.nodes,self.memo=list(names),[],{}
        self.zero=self.q(0);self.one=self.q(1)
    def node(self,op,arg):
        key=(op,tuple(arg) if isinstance(arg,list) else arg)
        if key not in self.memo:self.memo[key]=len(self.nodes);self.nodes.append([op,arg])
        return self.memo[key]
    def q(self,c):return self.node('q',str(c))
    def add(self,*args):
        args=[x for x in args if x!=self.zero]
        if not args:return self.zero
        return args[0] if len(args)==1 else self.node('+',sorted(args))
    def mul(self,*args):
        if self.zero in args:return self.zero
        args=[x for x in args if x!=self.one]
        if not args:return self.one
        return args[0] if len(args)==1 else self.node('*',sorted(args))
    def inv(self,x):return self.node('inv',x)
    def expr(self,expr):
        if expr.is_Rational:return self.q(expr)
        if expr.is_Symbol:return self.node('v',self.names.index(str(expr)))
        if expr.is_Add:return self.add(*(self.expr(a) for a in expr.args))
        if expr.is_Mul:return self.mul(*(self.expr(a) for a in expr.args))
        if expr.is_Pow and expr.exp.is_Integer:
            z=self.expr(expr.base)
            if expr.exp<0:z=self.inv(z)
            return self.mul(*([z]*abs(int(expr.exp))))
        raise ValueError(expr)


def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    familypath=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz'
    geometrypath=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    family,geometry,spin,info=map(read,[familypath,geometrypath,spinpath,invpath]);assert geometry['all_passed']
    u,v,c,X,W,B=s.symbols('u v c X W B');variables=[u,v,c,X,W,B]
    xx=s.Matrix(s.symbols('x y z w'));q=s.Matrix(s.symbols('q0:4'))
    P=s.Matrix(family['full_affine_P']).subs(dict(zip(s.symbols('n0:3'),[u,1,v])))
    Z=-s.Rational(3,4)/u*(B+W/s.Integer(1536)-W**2+X**2/s.Integer(64)+s.Rational(2,9))
    point={xx[0]:X,xx[1]:0,xx[2]:Z,xx[3]:W}
    value=P.subs(point,simultaneous=True).applyfunc(s.cancel)
    linear=sum((P.diff(xx[k]).subs(point,simultaneous=True)*q[k] for k in range(4)),s.zeros(4)).applyfunc(s.cancel)
    A=s.Rational(1,36)-s.Rational(256,3)*W
    expected=s.Matrix([[0,A,3*X*A/4,s.Rational(8,9)],[-A,0,B,4*X/3],
                       [-3*X*A/4,-B,0,X**2],[-s.Rational(8,9),-4*X/3,-X**2,0]])
    assert (value-expected).applyfunc(s.cancel)==s.zeros(4)
    assert s.cancel(value[0,1]*value[2,3]-value[0,2]*value[1,3]+value[0,3]*value[1,2]-8*B/9)==0
    R=s.Matrix([[1,0,0,0],[0,0,0,s.Rational(9,8)/c],
                [-3*X/2,1,0,-9*A/8],[-9*X**2/(8*c*B),0,1/(c*B),-27*X*A/(32*c*B)]])
    assert (R*(c*value)*R.T-J).applyfunc(s.cancel)==s.zeros(4)
    inverse=R.inv()
    shifted=linear.subs(dict(zip(q,inverse*q)),simultaneous=True)
    normalized=(c*R*shifted*R.T).applyfunc(s.cancel)
    tau=normalized[0,1]+normalized[2,3]
    alpha=s.Matrix([s.diff(tau,z)/5 for z in q]);vv=J*alpha
    primitive=(normalized-2*alpha.dot(q)*J+q*vv.T-vv*q.T).applyfunc(s.cancel)
    beta=(-J*primitive*J).applyfunc(s.cancel)
    assert s.cancel(beta[0,1]+beta[2,3])==0
    for i,j,k in itertools.combinations(range(4),3):
        assert s.cancel(s.diff(beta[j,k],q[i])-s.diff(beta[i,k],q[j])+s.diff(beta[i,j],q[k]))==0
    beta24=s.Matrix([s.diff(beta[i,j],q[k]) for (i,j),k in COORDS])
    basis=s.Matrix(info['primitive24_by16_basis']);br=list(basis.T.rref()[1])
    beta16=(basis[br,:].inv()*beta24[br,:]).applyfunc(s.cancel)
    assert (basis*beta16-beta24).applyfunc(s.cancel)==s.zeros(24,1)
    weights=info['Cartan_weights'];cweights=[(1-sum(w))//2 for w in weights]
    assert all(1-sum(w)==2*e for w,e in zip(weights,cweights))
    assert all(s.cancel(b-c**e*b.subs(c,1))==0 for b,e in zip(beta16,cweights))
    dag=DAG(list(map(str,variables)));beta_nodes=[dag.expr(s.factor(z)) for z in beta16]
    monomials=spin['quadratic_monomials']
    mon_nodes=[dag.mul(beta_nodes[i],beta_nodes[j]) for i,j in monomials]
    def covariant(key):
        return [[dag.add(*(dag.mul(dag.q(s.Rational(co)),mon_nodes[k]) for k,co in enumerate(row) if co!='0'))
                 for row in spin[key][4*i:4*i+4]] for i in range(4)]
    C,D=map(covariant,['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136'])
    def mm(E,F):return [[dag.add(*(dag.mul(E[i][k],F[k][j]) for k in range(4))) for j in range(4)] for i in range(4)]
    def tr(E):return dag.add(*(E[i][i] for i in range(4)))
    def scale(q,E):return [[dag.mul(dag.q(s.Rational(q)),z) for z in row] for row in E]
    def add(E,F):return [[dag.add(E[i][j],F[i][j]) for j in range(4)] for i in range(4)]
    CC,CD,DD=mm(C,C),mm(C,D),mm(D,D)
    qC,qD=tr(CC),tr(DD)
    trD4=tr(mm(DD,DD));Snode=dag.mul(dag.q(s.Rational(1,2)),qD)
    Pnode=dag.add(dag.mul(dag.q(s.Rational(1,8)),qD,qD),dag.mul(dag.q(-s.Rational(1,4)),trD4))
    Delta8=dag.add(dag.mul(Snode,Snode),dag.mul(dag.q(-4),Pnode))
    N12=dag.add(tr(mm(mm(C,DD),mm(C,DD))),dag.mul(dag.q(-1),Snode,tr(mm(CC,DD))),dag.mul(Pnode,qC))
    mnode=dag.mul(N12,dag.q(s.Rational(1,2)),dag.inv(Delta8))
    moment=add(scale(-1,C),scale(2,D));F4=tr(mm(moment,moment))
    I8=[tr(mm(E,F)) for E,F in [(CC,CC),(CC,CD),(CC,DD),(CD,CD),(CD,DD)]]
    change=s.Matrix([[0,0,-s.Rational(1,2),0,-1],[-s.Rational(1,2),-1,1,0,2],
                     [s.Rational(1,2),0,-2,-s.Rational(1,2),-2],[-s.Rational(1,4),0,s.Rational(1,2),0,0],
                     [s.Rational(1,32),0,s.Rational(1,16),0,s.Rational(1,8)]])
    pnodes=[dag.add(*(dag.mul(dag.q(z),I8[j]) for j,z in enumerate(row) if z)) for row in change.tolist()]
    degree=[]
    for op,arg in dag.nodes:
        if op=='q':d0=None if arg=='0' else 0
        elif op=='v':d0=1 if arg==2 else 0
        elif op=='inv':d0=-degree[arg] if degree[arg] is not None else None
        elif op=='*':d0=None if any(degree[j] is None for j in arg) else sum(degree[j] for j in arg)
        else:
            degrees={degree[j] for j in arg if degree[j] is not None};assert len(degrees)<=1
            d0=next(iter(degrees),None)
        degree.append(d0)
    assert all(degree[i] in (None,e) for i,e in zip(beta_nodes,cweights))
    assert [degree[i] for i in [F4,Delta8,N12,mnode,*I8,*pnodes]]==[2,4,6,2,*([4]*10)]
    payload=dict(input_names=dag.names,nodes=dag.nodes,outputs=dict(beta16=beta_nodes,
                  C4=C,D4=D,F4=F4,trace8_basis=I8,determinant_pencil5=pnodes,
                  spectral_discriminant8=Delta8,N12=N12,m=mnode),
                  primitive_domain='u*v*(32*u+v)*c*B !=0; n1!=0 before the parameter gauge',
                  m_additional_domain='spectral_discriminant8 !=0',
                  coefficient_field='QQ',no_point_evaluation_performed=True)
    with gzip.open(out/'source-six-circuit.json.gz','wt') as stream:json.dump(payload,stream);stream.write('\n')
    cert=dict(source_variables=list(map(str,variables)),y_zero_value_z=str(Z),raw_value4=rows(value),
              raw_linear4=rows(linear),rational_Darboux_frame4=rows(R),normalized_linear4=rows(normalized),
              primitive_perspective_alpha=list(map(str,alpha)),primitive_beta4=rows(beta),
              primitive_beta24=list(map(str,beta24)),primitive_beta16=list(map(str,beta16)),
              exact_c_exponents_beta16=cweights,invariant_c_exponents=dict(F4=2,Inv8=4,Delta8=4,N12=6,m=2),
              circuit_nodes=len(dag.nodes),no_five_parameter_quotient_claimed=True,
              remaining_one_dimensional_source_redundancy_not_integrated=True)
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    files=[Path(__file__),familypath,geometrypath,spinpath,invpath,ROOT/'code/weight1_structure.py']
    result=dict(id=STEM,all_passed=True,six_parameter_actual_NEW25_jet_map=True,
                all_value_linear_primitive_coefficients_exact=True,explicit_rational_Darboux_frame=True,
                symbolic_m_N12_binding_to_actual_family=True,circuit_nodes=len(dag.nodes),
                overall_scalar_weights_verified_without_invariant_expansion=True,
                no_new_beta_point_no_high_degree_fit_no_GB=True,five_parameter_quotient_or_integrability_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','circuit_nodes','seconds']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())
