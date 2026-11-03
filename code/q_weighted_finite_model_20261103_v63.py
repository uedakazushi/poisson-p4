#!/usr/bin/env python3
"""Build an exact finite weighted 35-variable transverse model at a=2/9.

All arithmetic selecting the slice and implicit blocks is over QQ(a), not a
finite field. Positive weights make the implicit solution polynomial and
triangular; no finite-order jet truncation is used.
"""
import argparse
from collections import defaultdict,Counter
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

ROOT=Path(__file__).resolve().parents[1]
X=s.symbols('x0:5');A=s.Symbol('a');ASTAR=s.Rational(2,9)
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))
BASIS=list(itertools.product(PAIRS,MONS))
WEIGHTS=[i+j-sum(mon)for(i,j),mon in BASIS]
TRIPLES=list(itertools.combinations(range(5),3))
MON3=list(itertools.combinations_with_replacement(range(5),3))


def fraction(value):
    return Fraction(int(s.numer(value)),int(s.denom(value)))


def coefficients(P):
    return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in mon))
                     for(i,j),mon in BASIS])


def point(a):
    P=s.zeros(5)
    for(i,j),v in {(1,4):a*X[0]**2,(2,3):-s.Rational(2,9)*X[0]**2,
                  (2,4):s.Rational(4,3)*X[0]*X[1],
                  (3,4):s.Rational(2,3)*X[0]*X[2]+X[1]**2}.items():
        P[i,j]=v;P[j,i]=-v
    return P


def frames():
    divergence=s.zeros(25,150)
    for n,((i,j),mon)in enumerate(BASIS):
        for k,target,sign in ((i,j,1),(j,i,-1)):
            if k in mon:
                rem=list(mon);rem.remove(k)
                divergence[5*target+rem[0],n]+=sign*mon.count(k)
    reduced,pivots=divergence.rref();assert len(pivots)==24
    free=[j for j in range(150)if j not in pivots];C=s.zeros(150,126)
    for k,j in enumerate(free):
        C[j,k]=1
        for row,p in enumerate(pivots):C[p,k]=-reduced[row,j]
    assert divergence*C==s.zeros(25,126)
    assert all(WEIGHTS[i]==WEIGHTS[free[j]]for i in range(150)for j in range(126)if C[i,j])
    Q=point(A);orbit=[];labels=[]
    for i,j in itertools.product(range(5),repeat=2):
        g=s.zeros(5);g[i,j]=1;v=g*s.Matrix(X)
        R=g*Q+Q*g.T-sum((v[k]*Q.diff(X[k])for k in range(5)),s.zeros(5))
        r=coefficients(R);assert divergence*r==s.zeros(25,1)
        orbit.append(r[free,:]);labels.append((i,j))
    O=s.Matrix.hstack(*orbit);Ostar=O.subs(A,ASTAR)
    columns=list(DomainMatrix.from_Matrix(Ostar).rref()[1]);assert len(columns)==17
    O1=Ostar[:,columns];aindex=BASIS.index(((1,4),(0,0)))
    cuts=[free.index(j)for j in free if WEIGHTS[j]>5 or WEIGHTS[j]==5 and j!=aindex]
    assert len(cuts)==7;rank=O1[cuts,:].rank();assert rank==7
    for k in sorted(range(126),key=lambda k:(-WEIGHTS[free[k]],free[k])):
        if k in cuts or free[k]==aindex:continue
        if O1[cuts+[k],:].rank()>rank:cuts.append(k);rank+=1
        if rank==17:break
    assert rank==17
    remaining=[k for k in range(126)if k not in cuts]
    positive=[k for k in remaining if free[k]!=aindex]
    degrees=[5-WEIGHTS[free[k]]for k in positive]
    assert len(positive)==108 and min(degrees)>0 and max(degrees)==12
    groupdet=s.factor(O[cuts,columns].det());assert groupdet.subs(A,ASTAR)!=0
    return C,free,cuts,positive,degrees,groupdet,[labels[i]for i in columns]


def equations(C,positive,degrees):
    # Coefficients of P are linear forms in a (index0), positive variables
    # (indices1..108), and the constant (index-1).
    forms=[defaultdict(Fraction)for _ in range(150)]
    Q0=coefficients(point(0));Qa=coefficients(point(1)-point(0))
    for i in range(150):
        if Q0[i]:forms[i][-1]=fraction(Q0[i])
        if Qa[i]:forms[i][0]=fraction(Qa[i])
        for j,k in enumerate(positive,1):
            if C[i,k]:forms[i][j]=fraction(C[i,k])
    entries={pair:{}for pair in PAIRS}
    for((pair,mon),form)in zip(BASIS,forms):
        if form:entries[pair][mon]=dict(form)
    cache={}
    for i,j in itertools.product(range(5),repeat=2):
        cache[i,j]=entries[i,j]if i<j else({m:{k:-v for k,v in f.items()}for m,f in entries[j,i].items()}if j<i else{})
    monindex={mon:i for i,mon in enumerate(MON3)}
    rows=[defaultdict(Fraction)for _ in range(350)]
    for n,(i,j,k)in enumerate(TRIPLES):
        for u,v,w in ((i,j,k),(j,k,i),(k,i,j)):
            for l in range(5):
                for m,f in cache[u,l].items():
                    for t,g in cache[v,w].items():
                        if l not in t:continue
                        rem=list(t);rem.remove(l);multiplicity=t.count(l)
                        row=rows[35*n+monindex[tuple(sorted((*m,*rem)))]]
                        for p,c in f.items():
                            for q,d in g.items():row[tuple(sorted((p,q)))]+=multiplicity*c*d
    rows=[{mon:c for mon,c in row.items()if c}for row in rows]
    weights={-1:0,0:0}|{i:d for i,d in enumerate(degrees,1)}
    target=[10-(sum(triple)-sum(mon))for triple in TRIPLES for mon in MON3]
    assert max(target)==19 and min(target)==1
    assert all(any(i>0 for i in mon)and sum(weights[i]for i in mon)==target[n]
               for n,row in enumerate(rows)for mon in row)
    linear=s.zeros(350,108)
    for i,row in enumerate(rows):
        for mon,c in row.items():
            if mon[0]in (-1,0)and mon[1]>0:
                linear[i,mon[1]-1]+=s.Rational(c.numerator,c.denominator)*(1 if mon[0]==-1 else A)
    return forms,rows,target,linear


def build(output):
    started=time.monotonic();C,free,cuts,positive,degrees,groupdet,groupcolumns=frames()
    print('Canonical126, weighted slice109, positive108',flush=True)
    forms,rows,target,linear=equations(C,positive,degrees)
    specialized=DomainMatrix.from_Matrix(linear.subs(A,ASTAR)).convert_to(s.QQ)
    pivots=list(specialized.rref()[1]);assert len(pivots)==74
    selected=list(specialized.extract(list(range(350)),pivots).transpose().rref()[1]);assert len(selected)==74
    remaining=[j for j in range(108)if j not in pivots];assert len(remaining)==34
    dag=DAG(['a']+[f'z{i}'for i in range(34)])
    values={-1:dag.one,0:dag.inputs['a']}
    for n,j in enumerate(remaining):values[j+1]=dag.inputs[f'z{n}']
    blocks=[];open_polynomial=groupdet
    for degree in sorted(set(degrees[j]for j in pivots)):
        cols=[j for j in pivots if degrees[j]==degree]
        rr=[i for i in selected if target[i]==degree]
        assert len(cols)==len(rr)
        matrix=linear[rr,cols]
        domain=DomainMatrix.from_Matrix(matrix).convert_to(s.QQ.frac_field(A))
        inverse=domain.inv().to_Matrix();assert(matrix*inverse).applyfunc(s.cancel)==s.eye(len(cols))
        determinant=s.factor(matrix.det());assert determinant.subs(A,ASTAR)!=0
        open_polynomial=s.factor(open_polynomial*determinant)
        rhs=[]
        for i in rr:
            terms=[]
            for j in remaining:
                if linear[i,j]:terms.append(dag.mul(dag.expression(linear[i,j]),values[j+1]))
            for mon,c in rows[i].items():
                if min(mon)>0:
                    assert all(degrees[j-1]<degree for j in mon)
                    terms.append(dag.mul(dag.q(c),*(values[j]for j in mon)))
            rhs.append(dag.neg(dag.add(*terms)))
        for n,j in enumerate(cols):
            values[j+1]=dag.add(*(dag.mul(dag.expression(inverse[n,k]),rhs[k])for k in range(len(rr))))
        blocks.append(dict(weight=degree,raw_positive_columns=cols,equation_rows=rr,
                           determinant=str(determinant),inverse=[[str(v)for v in row]for row in inverse.tolist()]))
        print('Exact implicit block',degree,len(cols),'nodes',len(dag.nodes),flush=True)
    assert len(values)==110
    reconstruction=[dag.add(*(dag.mul(dag.q(c),values[j])for j,c in form.items()))for form in forms]
    residual=[dag.add(*(dag.mul(dag.q(c),*(values[j]for j in mon))for mon,c in row.items()))for row in rows]
    print('All350 residual equations, nodes',len(dag.nodes),flush=True)
    # Verify weighted homogeneity of every circuit operation. All inverses have
    # weight0 and depend only on a; positive variables occur polynomially.
    inputweights=[0]+[degrees[j]for j in remaining];nodeweights=[]
    for op,*data in dag.nodes:
        if op=='q':w=None if Fraction(data[0])==0 else 0
        elif op=='v':w=inputweights[data[0]]
        elif op=='+':
            ww={nodeweights[j]for j in data[0]if nodeweights[j]is not None};assert len(ww)<=1
            w=next(iter(ww),None)
        elif op=='*':w=None if any(nodeweights[j]is None for j in data[0])else sum(nodeweights[j]for j in data[0])
        elif op=='inv':assert nodeweights[data[0]]==0;w=0
        else:raise ValueError(op)
        nodeweights.append(w)
    assert all(nodeweights[node]in (None,weight)for node,weight in zip(residual,target))
    assert all(nodeweights[values[j+1]]in(None,degrees[j])for j in range(108))
    numerator=s.fraction(open_polynomial)[0]
    reduced_open=s.prod(factor for factor,exponent in s.factor_list(numerator,A)[1])
    payload=dict(id='Q_WEIGHTED_FINITE_MODEL_20261103_V63',input_names=dag.names,input_weights=inputweights,
                 nodes=dag.nodes,all350_residual_nodes=residual,raw150_reconstruction_nodes=reconstruction,
                 canonical_free_raw_indices=free,group_cut_raw_indices=[free[k]for k in cuts],
                 group_lie_columns=groupcolumns,group_minor=str(groupdet),positive_raw_indices=[free[k]for k in positive],
                 positive_weights=degrees,free_positive_columns=remaining,implicit_blocks=blocks,
                 localized_a_open=str(open_polynomial),reduced_a_open=str(reduced_open),residual_weight_bounds=target,
                 raw_quadratic_equations=[[[*mon,str(c)]for mon,c in row.items()]for row in rows],
                 canonical_linear_forms=[[[j,str(c)]for j,c in form.items()]for form in forms])
    encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode()
    archive=output/'model.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as z:z.write(encoded)
    report=dict(id=payload['id'],construction_completed=True,construction_checks_passed=True,
                canonical_dimension=126,orbit_smooth_factor_dimension_affine=17,
                slice_dimension=109,implicit_variables=74,model_variables=35,
                positive_model_variables=34,zero_weight_parameter='a',center_a='2/9',
                maximum_implicit_weight=12,maximum_residual_weight=19,
                model_variable_weights=inputweights,circuit_nodes=len(dag.nodes),all350_equations_retained=True,
                exact_QQ_a_block_inverses_verified=True,positive_polynomiality_by_weight_verified=True,
                no_jet_truncation_used=True,global_cross_section_or_stabilizer_quotient_claimed=False,
                phase='construction; final validation and scope are recorded separately',
                localized_a_open=str(open_polynomial),reduced_a_open=str(reduced_open),
                decoded_sha256=hashlib.sha256(encoded).hexdigest(),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py']},
                seconds=time.monotonic()-started)
    (output/'construction.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Constructed exact finite model',report['seconds'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);build(args.output)
