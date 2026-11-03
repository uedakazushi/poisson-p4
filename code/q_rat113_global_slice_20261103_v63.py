#!/usr/bin/env python3
"""Exact two-Casimir Schur scheme in global37: polynomial A20 Rat113 branch."""
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
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from mps_s61_arithmetic_circuits_20261103_v62 import DAG
from q_smith_global_finite_model_20261103_v63 import closure

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2));PAIR_INDEX={p:i for i,p in enumerate(PAIRS)}
MON2=list(itertools.combinations_with_replacement(range(5),2))
MON3=list(itertools.combinations_with_replacement(range(5),3))


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def encode(M):return [[str(v)for v in row]for row in M.tolist()]


def jet_values(nodes,point):
    values=[]
    for op,*data in nodes:
        if op=='q':value=Fraction(data[0]);derivative={}
        elif op=='v':value=point[data[0]];derivative={data[0]:Fraction(1)}
        elif op=='+':
            value=sum((values[i][0]for i in data[0]),Fraction());derivative=defaultdict(Fraction)
            for i in data[0]:
                for j,c in values[i][1].items():derivative[j]+=c
        elif op=='*':
            value=Fraction(1);derivative={}
            for i in data[0]:
                v,d=values[i];new=defaultdict(Fraction)
                for j,c in derivative.items():new[j]+=c*v
                for j,c in d.items():new[j]+=c*value
                value*=v;derivative=new
        else:raise ValueError(op)
        values.append((value,{j:c for j,c in derivative.items()if c}))
    return values


def schur_source(parent):
    dag=DAG(parent['input_names'])
    dag.nodes=[(op,tuple(data[0]))if op in ('+','*')else(op,data[0])for op,*data in parent['nodes']]
    dag.lookup={node:i for i,node in enumerate(dag.nodes)};dag.zero=dag.lookup['q','0'];dag.one=dag.lookup['q','1']
    dag.inputs={name:dag.lookup['v',i]for i,name in enumerate(dag.names)};raw=parent['raw150_reconstruction_nodes']
    def entry(i,j,m):
        if i==j:return dag.zero
        return raw[15*PAIR_INDEX[i,j]+m]if i<j else dag.neg(raw[15*PAIR_INDEX[j,i]+m])
    contraction=[[entry(i,j,m)for j in range(5)]for i in range(5)for m in range(15)]
    pivot_rows=[45,30,50];pivot_columns=[2,3,4]
    pivot=[[contraction[i][j]for j in pivot_columns]for i in pivot_rows]
    assert pivot[0]==[dag.q(Fraction(2,9)),dag.zero,dag.zero]
    assert pivot[1]==[dag.zero,dag.q(Fraction(-2,9)),dag.zero]
    assert pivot[2][1]==dag.zero and pivot[2][2]==dag.one
    c=pivot[2][0];inverse=[[dag.q(Fraction(9,2)),dag.zero,dag.zero],
                         [dag.zero,dag.q(Fraction(-9,2)),dag.zero],
                         [dag.neg(dag.mul(dag.q(Fraction(9,2)),c)),dag.zero,dag.one]]
    assert [[dag.add(*(dag.mul(pivot[i][k],inverse[k][j])for k in range(3)))for j in range(3)]for i in range(3)]==[[dag.one if i==j else dag.zero for j in range(3)]for i in range(3)]
    E=[[contraction[i][j]for j in [0,1]]for i in pivot_rows]
    kernel=[[dag.one,dag.zero],[dag.zero,dag.one],
            [dag.neg(dag.mul(dag.q(Fraction(9,2)),E[0][j]))for j in range(2)],
            [dag.mul(dag.q(Fraction(9,2)),E[1][j])for j in range(2)],
            [dag.sub(dag.mul(dag.q(Fraction(9,2)),c,E[0][j]),E[2][j])for j in range(2)]]
    schur=[[dag.add(*(dag.mul(row[k],kernel[k][j])for k in range(5)))for j in range(2)]for row in contraction]
    cc=s.Symbol('c');EE=s.Matrix(3,2,s.symbols('e0:6'))
    DD=s.Matrix([[s.Rational(2,9),0,0],[0,-s.Rational(2,9),0],[cc,0,1]])
    KK=s.Matrix.vstack(-s.Rational(9,2)*EE[0,:],s.Rational(9,2)*EE[1,:],s.Rational(9,2)*cc*EE[0,:]-EE[2,:])
    assert (DD*KK+EE).applyfunc(s.expand)==s.zeros(3,2)
    roots=[i for row in schur for i in row]
    assert roots[121]==dag.neg(dag.inputs['a'])
    root_weights=[5-i-j+sum(mon)for i in range(5)for mon in MON2 for j in range(2)]
    weights=[]
    for op,*data in dag.nodes:
        if op=='q':w=None if not Fraction(data[0])else 0
        elif op=='v':w=parent['input_weights'][data[0]]
        elif op=='+':
            options={weights[i]for i in data[0]if weights[i]is not None};assert len(options)<=1;w=next(iter(options),None)
        elif op=='*':w=None if any(weights[i]is None for i in data[0])else sum(weights[i]for i in data[0])
        else:raise ValueError(op)
        weights.append(w)
    assert all(weights[i]in(None,w)for i,w in zip(roots,root_weights))
    return dag,contraction,pivot_rows,pivot_columns,pivot,inverse,kernel,roots,root_weights


def parameterize(source,parent,roots,root_weights,jac,rows,columns,kernel):
    free=[i for i in range(37)if i not in columns];dag=DAG([parent['input_names'][i]for i in free]);assigned=[dag.zero]*37
    for i,name in zip(free,dag.names):assigned[i]=dag.inputs[name]
    def clone(outputs):
        values={}
        for i in sorted(closure(source.nodes,outputs)):
            op,*data=source.nodes[i]
            if op=='q':value=dag.q(data[0])
            elif op=='v':value=assigned[data[0]]
            elif op=='+':value=dag.add(*(values[j]for j in data[0]))
            elif op=='*':value=dag.mul(*(values[j]for j in data[0]))
            else:raise ValueError(op)
            values[i]=value
        return values
    blocks=[]
    for d in range(1,9):
        rr=[i for i in rows if root_weights[i]==d];cc=[i for i in columns if parent['input_weights'][i]==d]
        assert len(rr)==len(cc)
        if not rr:continue
        matrix=jac[rr,cc];inverse=matrix.inv();values=clone([roots[i]for i in rr]);rhs=[values[roots[i]]for i in rr]
        answers=[dag.neg(dag.add(*(dag.mul(dag.q(inverse[j,k]),rhs[k])for k in range(len(rr)))))for j in range(len(cc))]
        for column,answer in zip(cc,answers):assigned[column]=answer
        blocks.append(dict(weight=d,rows=rr,columns=cc,matrix=encode(matrix),inverse=encode(inverse),determinant=str(matrix.det())))
    assert [len(b['rows'])for b in blocks]==[1,2,3,4,4,2]
    source_outputs=[*parent['raw150_reconstruction_nodes'],*parent['all350_residual_nodes'],*roots,*[i for row in kernel for i in row]]
    values=clone(source_outputs);outputs=[*assigned,*[values[i]for i in source_outputs]]
    used=closure(dag.nodes,outputs);used.update([dag.zero,dag.one,*dag.inputs.values()]);renumber={old:new for new,old in enumerate(sorted(used))};nodes=[]
    for old in sorted(used):
        op,*data=dag.nodes[old];nodes.append((op,tuple(renumber[i]for i in data[0]))if op in ('+','*')else(op,data[0]))
    outputs=[renumber[i]for i in outputs]
    R,*variables=ring(','.join(dag.names),s.QQ);polynomials=[]
    for op,*data in nodes:
        if op=='q':value=R(s.Rational(data[0]))
        elif op=='v':value=variables[data[0]]
        elif op=='+':value=sum((polynomials[i]for i in data[0]),R.zero)
        elif op=='*':
            value=R.one
            for i in data[0]:value*=polynomials[i]
        polynomials.append(value)
    parent_nodes=outputs[:37];raw=outputs[37:187];equations=outputs[187:537];schur=outputs[537:687];kernel_nodes=outputs[687:]
    assert all(not polynomials[i]for i in [*equations,*schur])
    assert polynomials[parent_nodes[0]]==R.zero
    for j,index in enumerate(free):assert polynomials[parent_nodes[index]]==variables[j]
    return dict(id='Q_RAT113_GLOBAL_SLICE_20261103_V63',input_names=dag.names,
                input_weights=[parent['input_weights'][i]for i in free],parent37_input_indices=free,
                nodes=nodes,parent37_input_nodes=parent_nodes,raw150_reconstruction_nodes=raw,
                all350_residual_nodes=equations,all150_Schur_nodes=schur,kernel5x2_nodes=kernel_nodes,
                triangular_blocks=blocks),dict(parent37_polynomials=[str(polynomials[i])for i in parent_nodes],
                raw150_polynomials=[str(polynomials[i])for i in raw],
                kernel5x2_polynomials=[str(polynomials[i])for i in kernel_nodes],
                raw150_total_terms=sum(len(polynomials[i])for i in raw))


def fixed_casimirs(quadratics):
    cubics=[m for m in MON3 if any(i>=2 for i in m)];ci={m:i for i,m in enumerate(cubics)}
    T=s.zeros(150,31);L=s.zeros(31,150)
    for pair,active,sign in [((2,3),4,1),((2,4),3,-1),((3,4),2,1)]:
        for j,mon in enumerate(MON2):
            raw=15*PAIR_INDEX[pair]+j;extended=tuple(sorted([*mon,active]))
            T[raw,ci[extended]]=sign*extended.count(active)
            L[ci[extended],raw]+=s.Rational(sign,1+sum(i>=2 for i in mon))
    assert L*T==s.eye(31)and DomainMatrix.from_Matrix(T).rank()==31
    divergence=s.zeros(25,150)
    for n,((i,j),mon)in enumerate(itertools.product(PAIRS,MON2)):
        for k,target,sign in [(i,j,1),(j,i,-1)]:
            if k in mon:
                changed=list(mon);changed.remove(k);divergence[5*target+changed[0],n]+=sign*mon.count(k)
    active=[15*PAIR_INDEX[p]+j for p in [(2,3),(2,4),(3,4)]for j in range(15)]
    assert divergence*T==s.zeros(25,31)and DomainMatrix.from_Matrix(divergence[:,active]).rank()==14
    row_data=[[(j,T[i,j])for j in range(31)if T[i,j]]for i in range(150)]
    for polynomial in quadratics:
        result=defaultdict(lambda:s.S.Zero)
        for i,j,c in polynomial:
            for k,b in row_data[i]:
                for ell,d in row_data[j]:result[tuple(sorted((k,ell)))]+=c*b*d
        assert all(not v for v in result.values())
    return dict(cubic_monomials=cubics,H31_to_raw150=encode(T),raw150_to_H31_relative_integration=encode(L),
                active45_raw_indices=active,divergence25_by45=encode(divergence[:,active]),
                fixed_W_divergence_rank=14,fixed_W_kernel_dimension=31,
                exact_integration_left_inverse=True,all350_fixed_W_Jacobi_polynomials_identically_zero=True)


def archive(path,payload):
    with path.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(payload,separators=(',',':'))+'\n').encode())


def build(output):
    started=time.monotonic();folder=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63'
    parent_path=folder/'model37.json.gz';frame_path=folder/'frames.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()));frames=json.loads(gzip.decompress(frame_path.read_bytes()))
    source,M,pivot_rows,pivot_columns,pivot,inverse,kernel,roots,root_weights=schur_source(parent)
    jets=jet_values(source.nodes,[Fraction(0)]*37)
    assert all(not jets[i][0]for i in roots)
    jac=s.Matrix([[jets[i][1].get(j,Fraction())for j in range(37)]for i in roots])
    rows=list(DomainMatrix.from_Matrix(jac.T).rref()[1]);columns=list(DomainMatrix.from_Matrix(jac[rows,:]).rref()[1])
    assert len(rows)==len(columns)==17
    minor=DomainMatrix.from_Matrix(jac[rows,columns]).det();assert str(minor)=='-7334/23914845'
    parameter,polynomials=parameterize(source,parent,roots,root_weights,jac,rows,columns,kernel)
    fixed=fixed_casimirs(frames['physical_raw350_quadratics'])
    dag=DAG(parameter['input_names']);dag.nodes=parameter['nodes'];values=dag.evaluate(dict.fromkeys(dag.names,Fraction(1)))
    point=[values[i]for i in parameter['parent37_input_nodes']]
    dense_jets=jet_values(parent['nodes'],point)
    assert all(not dense_jets[i][0]for i in parent['all350_residual_nodes'])
    full_jac=s.Matrix([[dense_jets[i][1].get(j,Fraction())for j in range(37)]for i in parent['all350_residual_nodes']])
    assert DomainMatrix.from_Matrix(full_jac).rank()==17
    dependent=full_jac[:,columns];dense_rows=list(DomainMatrix.from_Matrix(dependent.T).rref()[1]);assert len(dense_rows)==17
    dense_minor=DomainMatrix.from_Matrix(dependent[dense_rows,:]).det();assert dense_minor
    output.mkdir(parents=True,exist_ok=False)
    archive(output/'source-Schur.json.gz',dict(input_names=source.names,input_weights=parent['input_weights'],nodes=source.nodes,
        unchanged_parent_prefix=len(parent['nodes']),raw150_parent_nodes=parent['raw150_reconstruction_nodes'],
        parent350_nodes=parent['all350_residual_nodes'],contraction75x5_nodes=M,
        pivot_rows=pivot_rows,pivot_columns=pivot_columns,pivot3_nodes=pivot,pivot3_inverse_nodes=inverse,
        kernel5x2_nodes=kernel,Schur150_nodes=roots,Schur150_weights=root_weights))
    archive(output/'model20.json.gz',parameter)
    archive(output/'fixed-Casimir-linear.json.gz',fixed)
    certificate=dict(Schur150x37_Jacobian_at_Q0=encode(jac),Schur_selected_rows=rows,dependent_parent_input_columns=columns,
                     Schur17_minor=str(minor),dense_twenty_inputs=['1']*20,dense_parent37_point=list(map(str,point)),
                     dense_full350x37_Jacobian=encode(full_jac),dense_selected_equation_rows=dense_rows,
                     dense_dependent17_minor=str(dense_minor),**polynomials)
    (output/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    note=ROOT/'research/notes/Q_RAT113_GLOBAL_SLICE_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,frame_path,ROOT/'code/q_smith_global_finite_model_20261103_v63.py',
            ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',folder.parent/'Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW.json',
            ROOT/'legacy/v14/sections/rank2_global_eight.tex']
    report=dict(id=parameter['id'],all_passed=True,constant_contraction3_minor='-4/81',
                exact_Schur_rank_at_Q0=17,exact_Schur17_minor=str(minor),literal_a_equation=True,
                finite_triangular_block_sizes=[len(b['rows'])for b in parameter['triangular_blocks']],
                actual_global_Schur_quotient='QQ[20 parameters], with a=0',
                all350_Jacobi_and150_Schur_exact_QQ20_zero=True,
                raw150_total_polynomial_terms=polynomials['raw150_total_terms'],
                fixed_W_cubic_mod_Sym3W_rank=31,fixed_W_relative_divergence_rank=14,
                fixed_W_linear_integration_inverse_and_all350_Jacobi_verified=True,
                recognition_over_arbitrary_QQ_algebras_including_nilpotents=True,
                known_component='Rat(1,1,3)',known_component_newly_counted=False,
                branch_smooth_at_Q0_dimension=20,
                full_parent350x37_Jacobian_rank_at_dense_point=17,
                exact_dense17_minor=str(dense_minor),full_parent_regular_at_dense_point_dimension=20,
                prime_is_actual_minimal_component_of_full_R37=True,
                full_R37_local_and_completed_nonflat_over_a_at_Q0=True,
                nonflatness_reason='a belongs to the vertical minimal prime N113; minimal primes of a Noetherian ring are associated, so a is a zero divisor. This uses the full20 branch and actual full-Jacobian generic rank, not the closed2-plane alone.',
                full_R37_regular_at_Q0_claimed=False,full_remaining_branch_classification_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS actual Rat113 Schur=A20, universal350+150 zero, full-Jacobian17;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
