#!/usr/bin/env python3
"""Freeze full period/edge/Sp4 rational S61 maps as finite arithmetic DAGs.

No rational expression is expanded unnecessarily. A circuit is executable exact
arithmetic (Q constants, input variables, +, *, inverse); its61 outputs are the
actual free S61 coordinates after full projective/value/perspective normalization.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s

from mps_s61_catalogue_20261103_v62 import raw_family,DIAGRAMS,ROOT
from s61_model import S61,X,J,PAIRS,coefficients,pushforward,vertical_rref

MONS=list(itertools.combinations_with_replacement(range(5),2))
EXPS=[tuple(ids.count(i) for i in range(5)) for ids in MONS]


class DAG:
    def __init__(self,names):
        self.names=names;self.nodes=[];self.lookup={};self.zero=self.q(0);self.one=self.q(1)
        self.inputs={name:self.node(('v',i)) for i,name in enumerate(names)}
    def node(self,data):
        if data not in self.lookup:self.lookup[data]=len(self.nodes);self.nodes.append(data)
        return self.lookup[data]
    def q(self,value):return self.node(('q',str(Fraction(value))))
    def add(self,*values):
        pending=list(values);coefficients=defaultdict(Fraction);constant=Fraction(0)
        while pending:
            value=pending.pop();node=self.nodes[value]
            if node[0]=='+':pending.extend(node[1]);continue
            if node[0]=='q':constant+=Fraction(node[1]);continue
            coefficient=Fraction(1);base=value
            if node[0]=='*':
                factors=[]
                for j in node[1]:
                    if self.nodes[j][0]=='q':coefficient*=Fraction(self.nodes[j][1])
                    else:factors.append(j)
                base=factors[0] if len(factors)==1 else self.node(('*',tuple(factors)))
            coefficients[base]+=coefficient
        result=[self.mul(self.q(c),base) for base,c in sorted(coefficients.items()) if c]
        if constant:result.append(self.q(constant))
        if not result:return self.zero
        if len(result)==1:return result[0]
        return self.node(('+',tuple(sorted(result))))
    def mul(self,*values):
        pending=list(values);factors=[];constant=Fraction(1)
        while pending:
            value=pending.pop();node=self.nodes[value]
            if value==self.zero:return self.zero
            if node[0]=='q':constant*=Fraction(node[1])
            elif node[0]=='*':pending.extend(node[1])
            else:factors.append(value)
        if not constant:return self.zero
        if constant!=1:factors.append(self.q(constant))
        if not factors:return self.one
        if len(factors)==1:return factors[0]
        return self.node(('*',tuple(sorted(factors))))
    def neg(self,value):return self.mul(self.q(-1),value)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def inv(self,value):
        node=self.nodes[value]
        if node[0]=='q':return self.q(1/Fraction(node[1]))
        if node[0]=='inv':return node[1]
        return self.node(('inv',value))
    def power(self,value,n):
        if n<0:return self.inv(self.power(value,-n))
        result=self.one
        while n:
            if n&1:result=self.mul(result,value)
            n//=2
            if n:value=self.mul(value,value)
        return result
    def expression(self,value):
        value=s.sympify(value)
        if value.is_Rational:return self.q(Fraction(int(value.p),int(value.q)))
        if value.is_Symbol:return self.inputs[str(value)]
        if value.is_Add:return self.add(*(self.expression(v) for v in value.args))
        if value.is_Mul:return self.mul(*(self.expression(v) for v in value.args))
        if value.is_Pow and value.exp.is_Integer:return self.power(self.expression(value.base),int(value.exp))
        raise ValueError(value)
    def evaluate(self,point):
        values=[]
        for node in self.nodes:
            op=node[0]
            if op=='q':value=Fraction(node[1])
            elif op=='v':value=Fraction(point[self.names[node[1]]])
            elif op=='+':value=sum((values[i] for i in node[1]),Fraction(0))
            elif op=='*':
                value=Fraction(1)
                for i in node[1]:value*=values[i]
            elif op=='inv':value=1/values[node[1]]
            else:raise ValueError(op)
            values.append(value)
        return values


def identity(dag,n):return [[dag.one if i==j else dag.zero for j in range(n)] for i in range(n)]
def matmul(dag,A,B):
    return [[dag.add(*(dag.mul(A[i][k],B[k][j]) for k in range(len(B))))
             for j in range(len(B[0]))] for i in range(len(A))]
def transpose(A):return list(map(list,zip(*A)))
def block(dag,A):
    M=identity(dag,5)
    for i in range(4):
        for j in range(4):M[i+1][j+1]=A[i][j]
    return M


def sp4(dag,names):
    z=[dag.inputs[n] for n in names];a,b,c,d=z[:4]
    determinant=dag.sub(dag.mul(a,d),dag.mul(b,c));inverse=dag.inv(determinant)
    A=[[a,b],[c,d]]
    invT=[[dag.mul(d,inverse),dag.mul(dag.q(-1),c,inverse)],
          [dag.mul(dag.q(-1),b,inverse),dag.mul(a,inverse)]]
    S=[[z[4],z[5]],[z[5],z[6]]];T=[[z[7],z[8]],[z[8],z[9]]]
    AS=matmul(dag,A,S);TA=matmul(dag,T,A);TAS=matmul(dag,TA,S)
    Q=[A[i]+AS[i] for i in range(2)]+[
        TA[i]+[dag.add(TAS[i][j],invT[i][j]) for j in range(2)] for i in range(2)]
    # Convert (q1,q2,p1,p2) back to (q1,p1,q2,p2).
    order=[0,2,1,3];G=[[Q[order[i]][order[j]] for j in range(4)] for i in range(4)]
    j4=[[dag.q(J[i+1,j+1]) for j in range(4)] for i in range(4)]
    Gi=matmul(dag,matmul(dag,j4,transpose(G)),j4)
    Gi=[[dag.neg(v) for v in row] for row in Gi]
    return G,Gi,determinant


def transform(dag,values,M,Mi):
    substitutions=[]
    for k,l in MONS:
        row=[]
        for a,b in MONS:
            value=dag.mul(Mi[k][a],Mi[l][b])
            if a!=b:value=dag.add(value,dag.mul(Mi[k][b],Mi[l][a]))
            row.append(value)
        substitutions.append(row)
    result=[]
    for i,j in PAIRS:
        output=[[] for _ in MONS]
        for h,(r,t) in enumerate(PAIRS):
            factor=dag.sub(dag.mul(M[i][r],M[j][t]),dag.mul(M[i][t],M[j][r]))
            if factor==dag.zero:continue
            for k,coefficient in enumerate(values[h*15:(h+1)*15]):
                if coefficient==dag.zero:continue
                for column,value in enumerate(substitutions[k]):
                    if value!=dag.zero:output[column].append(dag.mul(factor,coefficient,value))
        result.extend(dag.add(*terms) for terms in output)
    return result


def construct(record,P,parameters):
    group_names=[f'g{i}' for i in range(10)];dag=DAG([*map(str,parameters),*group_names])
    raw=[]
    for i,j in PAIRS:
        poly=s.Poly(P[i,j],*X)
        raw.extend(dag.expression(poly.coeff_monomial(mon)) for mon in EXPS)
    at_one=P.subs(dict.fromkeys(X,1))
    C=s.Matrix(4,4,lambda i,j:at_one[i+1,j+1]-at_one[0,j+1]+at_one[0,i+1])
    pivot=next((i,j) for i,j in itertools.combinations(range(4),2) if s.cancel(C[i,j])!=0)
    perm=[*pivot,*[i for i in range(4) if i not in pivot]]
    C=C[perm,perm];a,b,c,d,e,f=[dag.expression(v) for v in [C[i,j] for i,j in itertools.combinations(range(4),2)]]
    delta=dag.add(dag.mul(a,f),dag.neg(dag.mul(b,e)),dag.mul(c,d));ai=dag.inv(a);di=dag.inv(delta)
    F0=[[dag.one,dag.zero,dag.zero,dag.zero],[dag.zero,ai,dag.zero,dag.zero],
        [dag.mul(d,ai),dag.neg(dag.mul(b,ai)),dag.one,dag.zero],
        [dag.mul(e,di),dag.neg(dag.mul(c,di)),dag.zero,dag.mul(a,di)]]
    Fi0=[[dag.one,dag.zero,dag.zero,dag.zero],[dag.zero,a,dag.zero,dag.zero],
         [dag.neg(dag.mul(d,ai)),b,dag.one,dag.zero],
         [dag.neg(dag.mul(e,ai)),c,dag.zero,dag.mul(delta,ai)]]
    reorder=[[dag.one if j==perm[i] else dag.zero for j in range(4)] for i in range(4)]
    F=matmul(dag,F0,reorder);Fi=matmul(dag,transpose(reorder),Fi0)
    G,Gi,group_det=sp4(dag,group_names)
    move=identity(dag,5);move_inverse=identity(dag,5)
    for i in range(1,5):move[i][0]=dag.q(-1);move_inverse[i][0]=dag.one
    M=matmul(dag,block(dag,G),matmul(dag,block(dag,F),move))
    Mi=matmul(dag,move_inverse,matmul(dag,block(dag,Fi),block(dag,Gi)))
    based=transform(dag,raw,M,Mi)
    model=S61();rref,pivots=vertical_rref()
    def quotient(values):
        return [dag.add(value,*(dag.mul(dag.q(-rref[i,k]),values[p]) for i,p in enumerate(pivots)
                               if rref[i,k])) for k,value in enumerate(values)]
    initial=quotient(based);response=model.perspective_matrix().inv()
    selected=[initial[i] for i in model.certificate['perspective_orbit_pivots_raw']]
    alpha=[dag.add(*(dag.mul(dag.q(-response[i,j]),selected[j]) for j in range(4))) for i in range(4)]
    perspective=identity(dag,5);perspective_inverse=identity(dag,5)
    for i in range(1,5):perspective[0][i]=alpha[i-1];perspective_inverse[0][i]=dag.neg(alpha[i-1])
    normal=quotient(transform(dag,based,perspective,perspective_inverse))
    outputs=[normal[i] for i in model.certificate['raw_free_indices']]
    return dag,outputs,dict(value_pivot_indices=list(pivot),value_pivot_node=a,value_pfaffian_node=delta,
                           sp4_big_cell_determinant_node=group_det,linear_normalizing_matrix_nodes=M,
                           source_parameters=list(map(str,parameters)),sp4_parameters=group_names)


def verify_sample(dag,outputs,metadata,P,parameters):
    group=[2,1,1,1,1,2,3,2,1,4]
    for height in range(2,16):
        point={str(p):height**(i+1)+3*i+1 for i,p in enumerate(parameters)}
        point.update({f'g{i}':v for i,v in enumerate(group)})
        try:values=dag.evaluate(point)
        except ZeroDivisionError:continue
        matrix=s.Matrix([[s.Rational(values[j].numerator,values[j].denominator) for j in row]
                         for row in metadata['linear_normalizing_matrix_nodes']])
        source=P.subs({p:point[str(p)] for p in parameters})
        old=S61();normalized=old.normalize(pushforward(source,matrix))
        expected=[s.Rational(values[j].numerator,values[j].denominator) for j in outputs]
        assert expected==normalized['parameters']
        checks=old.verify_point(normalized)
        return point,checks
    raise RuntimeError('No bounded generic rational chart was located')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--families',nargs='+',default=[f'T{i:02d}' for i in range(1,39)])
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    diagrams={r['id']:r for r in json.loads(DIAGRAMS.read_text())['diagrams']};records=[]
    for name in args.families:
        start=time.monotonic();record=diagrams[name];p,t,P,pf,_=raw_family(record)
        completed=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62'/f'{name}.json'
        if completed.is_file():
            body=json.loads(completed.read_text());P=s.zeros(5)
            for pair,expression in body['homogeneous_bivector'].items():
                i,j=map(int,pair.split(','));P[i,j]=s.sympify(expression);P[j,i]=-P[i,j]
        parameters=(*p,*t);dag,outputs,metadata=construct(record,P,parameters)
        point,checks=verify_sample(dag,outputs,metadata,P,parameters)
        payload=dict(component_id=name,input_names=dag.names,nodes=dag.nodes,output_u1_through_u61=outputs,
                     normalization=metadata,source_period_pfaffian_open=str(pf),
                     source_coefficient_denominators_required_nonzero=sorted({str(s.factor(s.denom(c)))
                         for i,j in PAIRS for c in s.Poly(P[i,j],*X).coeffs() if s.denom(c)!=1}),
                     residual_sp4_included=True,exact_rational_sample=point,sample_checks=checks,
                     general_identity_proof_dependency='Exact generic source Jacobi + finite rational projective transformations + existing S61 reduction theorem; sample verifies implementation/conventions only.',
                     expected_s61_dimension=record['projective_component_dimension']-9,
                     full_source_dimension=len(parameters)+10,seconds=time.monotonic()-start)
        encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode();path=args.output/f'{name}.json.gz'
        with path.open('wb') as stream:
            with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as archive:archive.write(encoded)
        records.append(dict(component_id=name,nodes=len(dag.nodes),bytes=path.stat().st_size,
                            circuit_sha256=hashlib.sha256(encoded).hexdigest(),archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            exact_sample_checked=True,seconds=payload['seconds']))
        print(name,'FULL SP4 CIRCUIT PASS',len(dag.nodes),payload['seconds'],flush=True)
    inputs=[Path(__file__),DIAGRAMS,ROOT/'code/mps_s61_catalogue_20261103_v62.py',ROOT/'code/s61_model.py',
            ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/report.json']
    (args.output/'report.json').write_text(json.dumps(dict(id='MPS_S61_ARITHMETIC_CIRCUITS_20261103_v62',all_passed=True,
        records=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        full_residual_sp4_included=True,global_exhaustion_claimed=False),indent=2)+'\n')
