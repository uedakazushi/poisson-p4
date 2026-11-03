#!/usr/bin/env python3
"""All-a37 and a-unit35 exact finite Q models in polynomial Smith frames."""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

ROOT=Path(__file__).resolve().parents[1]
A=s.Symbol('a');DOMAIN=s.QQ.poly_ring(A)
PAIRS=list(itertools.combinations(range(5),2))
MON2=list(itertools.combinations_with_replacement(range(5),2))
MON3=list(itertools.combinations_with_replacement(range(5),3))
TRIPLES=list(itertools.combinations(range(5),3))
BASIS=list(itertools.product(PAIRS,MON2))
EXPS2=[tuple(mon.count(i)for i in range(5))for mon in MON2]
EXPS3=[tuple(mon.count(i)for i in range(5))for mon in MON3]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
@lru_cache(None)
def expression(text):return s.sympify(text,locals={'a':A})
def matrix(data):return s.Matrix([[expression(v)for v in row]for row in data])
def encode(M):return [[str(v)for v in row]for row in M.tolist()]
def dm(M):return DomainMatrix.from_Matrix(M).convert_to(DOMAIN).to_dense()


def inverse_unit(M):
    data=dm(M);numerator,denominator=data.inv_den(method='rref')
    denominator=DOMAIN.to_sympy(denominator)
    assert denominator and not denominator.has(A)
    inverse=numerator.to_Matrix().applyfunc(lambda v:s.expand(v/denominator))
    assert dm(M)*dm(inverse)==dm(s.eye(M.rows))
    return inverse


def saved_smith(block):
    B,U,V=matrix(block['matrix']),matrix(block['left']),matrix(block['right'])
    D=s.zeros(B.rows,B.cols)
    for i,text in enumerate(block['diagonal']):D[i,i]=expression(text)
    assert dm(U)*dm(B)*dm(V)==dm(D)
    for M in [U,V]:
        determinant=DOMAIN.to_sympy(dm(M).det());assert determinant and not determinant.has(A)
    return U,V


def smith(B,rows,columns,weight):
    value=dm(B);D,U,V=smith_normal_decomp(value)
    assert U*value*V==D
    determinants=[DOMAIN.to_sympy(M.det())for M in [U,V]]
    assert all(v and not v.has(A)for v in determinants)
    D,U,V=D.to_Matrix(),U.to_Matrix(),V.to_Matrix()
    diagonal=[s.factor(D[i,i])for i in range(min(D.shape))]
    assert all(not D[i,j]for i in range(D.rows)for j in range(D.cols)if i!=j)
    return dict(weight=weight,rows=rows,columns=columns,matrix=encode(B),left=encode(U),right=encode(V),
                diagonal=list(map(str,diagonal)),rank=sum(bool(v)for v in diagonal),
                left_determinant=str(determinants[0]),right_determinant=str(determinants[1]))


def raw_quadratics():
    positions={pair:i for i,pair in enumerate(PAIRS)};mon_index={m:i for i,m in enumerate(MON3)}
    out=[defaultdict(int)for _ in range(350)]
    def entry(i,j):
        if i==j:return[]
        pair,sign=((i,j),1)if i<j else((j,i),-1)
        return[(15*positions[pair]+m,mon,sign)for m,mon in enumerate(MON2)]
    for row,(i,j,k)in enumerate(TRIPLES):
        for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
            for ell in range(5):
                for p,mon,sign in entry(b,ell):
                    for q,other,second_sign in entry(c,d):
                        if ell not in other:continue
                        rem=list(other);rem.remove(ell)
                        target=35*row+mon_index[tuple(sorted([*mon,*rem]))]
                        out[target][tuple(sorted((p,q)))]+=sign*second_sign*other.count(ell)
    return[{mon:c for mon,c in row.items()if c}for row in out]


def read_frames(input_path):
    data=json.loads(gzip.decompress(input_path.read_bytes()))
    C=matrix(data['canonical_basis']);q=s.Matrix(list(map(expression,data['raw_Qa150'])))
    J=matrix(data['full350_by126_Jacobi']);O=matrix(data['full150_by25_orbit'])
    free=data['canonical_free_raw_indices'];canonical_index={v:i for i,v in enumerate(free)}
    degree=data['canonical_relative_weights'];target=data['Jacobi_row_relative_weights']
    orbit={-int(block['label'].rsplit('_',1)[1]):block for block in data['orbit_blocks']}
    group_rows=[];group_columns=[];group_row_weights=[];group_lie_weights=[]
    orbit_frames=[];B_parts=[];back_parts=[];weights=[];blocks=[]
    for block in data['initial_nonpositive_group_blocks']:
        S,T=saved_smith(block);rank=block['rank'];ell=int(block['label'].rsplit('_',1)[1])
        assert len(block['rows'])==rank
        assert all(not expression(v).has(A)for v in block['diagonal'])
        for raw in block['rows']:
            row=s.zeros(1,150);row[raw]=1;group_rows.append(row);group_row_weights.append(-ell)
        for k in range(rank):
            column=s.zeros(25,1)
            for j,index in enumerate(block['columns']):column[index]=T[j,k]
            group_columns.append(column);group_lie_weights.append(ell)
    assert len(group_rows)==len(group_columns)==7
    offset=0
    for d in sorted(v for v in set(degree)if v>0):
        canonical_columns=[j for j,v in enumerate(degree)if v==d]
        raw_rows=[free[j]for j in canonical_columns];n=len(raw_rows)
        if d in orbit:
            block=orbit[d];assert block['rows']==raw_rows
            S,T=saved_smith(block);rank=block['rank']
            # The saved rectangular Smith product is checked without padding assumptions.
            product=(S*O[raw_rows,block['columns']]*T).applyfunc(s.expand)
            assert all(product[i,j]==0 for i in range(product.rows)for j in range(product.cols)if i!=j)
            assert all(s.cancel(product[i,i]-expression(block['diagonal'][i]))==0 for i in range(min(product.shape)))
            assert all(not expression(v).has(A)for v in block['diagonal'])
            Sinv=inverse_unit(S);complement=Sinv[:,rank:]
            for k in range(rank):
                row=s.zeros(1,150)
                for j,raw in enumerate(raw_rows):row[raw]=S[k,j]
                group_rows.append(row);group_row_weights.append(d);column=s.zeros(25,1)
                for j,index in enumerate(block['columns']):column[index]=T[j,k]
                group_columns.append(column);group_lie_weights.append(-d)
        else:S=s.eye(n);rank=0;complement=S
        raw_basis=C[:,canonical_columns]*complement
        rr=[i for i,v in enumerate(target)if v==d]
        linear=J[rr,canonical_columns]*complement
        record=smith(linear,rr,list(range(offset,offset+complement.cols)),d)
        V=matrix(record['right']);Vinv=inverse_unit(V)
        B_parts.append((raw_basis*V).applyfunc(s.expand));back=s.zeros(complement.cols,150)
        transformed=Vinv*S[rank:,:]
        for i in range(back.rows):
            for j,raw in enumerate(raw_rows):back[i,raw]=s.expand(transformed[i,j])
        back_parts.append(back);weights.extend([d]*complement.cols);blocks.append(record)
        orbit_frames.append(dict(weight=d,canonical_columns=canonical_columns,raw_rows=raw_rows,
                                 removed_group_columns=rank,left=encode(S),complement=encode(complement)))
        offset+=complement.cols
    B=s.Matrix.hstack(*B_parts);back=s.Matrix.vstack(*back_parts)
    H=s.Matrix.vstack(*group_rows);G=s.Matrix.hstack(*group_columns)
    assert B.shape==(150,108)and back.shape==(108,150)and H.shape==(17,150)and G.shape==(25,17)
    assert dm(back)*dm(B)==dm(s.eye(108))and dm(back)*dm(q)==dm(s.zeros(108,1))
    assert dm(H)*dm(B)==dm(s.zeros(17,108))
    group_response=H*O*G;group_det=DOMAIN.to_sympy(dm(group_response).det())
    assert group_det and not group_det.has(A)
    assert sum(group_row_weights)+sum(group_lie_weights)==0
    assert all(not v.has(A)for v in H*q)
    assert dm(C[free,:])==dm(s.eye(126))
    raw_degrees=[5-sum(pair)+sum(mon)for pair,mon in BASIS]
    assert all(not q[i]or raw_degrees[i]==0 for i in range(150))
    assert all(not B[i,j]or raw_degrees[i]==weights[j]for i in range(150)for j in range(108))
    assert all(not G[i,j]or data['Lie_positions'][i][0]-data['Lie_positions'][i][1]==group_lie_weights[j]
               for i in range(25)for j in range(17))
    # At a general slice point the actual a-column is q'+B'w. The complete
    # action Jacobian is [q'+B'w, B, Lie_P G], rather than just H O_P G.
    column_weights=[0,*weights,*[-ell for ell in group_lie_weights]]
    assert sum(degree)-sum(column_weights)==0
    zero_action=s.Matrix.hstack(q.diff(A),B,O*G)[free,:]
    assert all(not s.expand(zero_action[i,j])for i in range(126)for j in range(126)
               if degree[i]!=column_weights[j])
    action_blocks=[];reordered_unit=s.S.One
    for d in sorted(set(degree)):
        rows=[i for i,w in enumerate(degree)if w==d];columns=[i for i,w in enumerate(column_weights)if w==d]
        assert len(rows)==len(columns)
        determinant=DOMAIN.to_sympy(dm(zero_action[rows,columns]).det())
        assert determinant and not determinant.has(A)
        reordered_unit*=determinant
        action_blocks.append(dict(weight=d,rows=rows,columns=columns,determinant=str(determinant)))
    return data,C,q,J,target,B,back,weights,blocks,dict(orbit_frames=orbit_frames,
        group_cut_matrix=encode(H),group_Lie_coefficient_matrix=encode(G),
        group_cut_targets=list(map(str,H*q)),group_cut_relative_weights=group_row_weights,
        group_Lie_weights=group_lie_weights,group_response=encode(group_response),group_determinant=str(group_det),
        full126_action_row_weights=degree,full126_action_column_weights=column_weights,
        full126_action_Qa_matrix=encode(zero_action),full126_action_Qa_blocks=action_blocks,
        full126_action_reordered_determinant=str(reordered_unit),
        full126_action_general_a_column="q'(a)+B'(a)*w",
        full126_action_general_other_columns='[B(a),Lie_P*G(a)]')


def source_model(q,B,target,weights,blocks,quadratics):
    dag=DAG(['a']+[f'w{i}'for i in range(108)])
    raw=[dag.add(dag.expression(q[i]),*(dag.mul(dag.expression(B[i,j]),dag.inputs[f'w{j}'])for j in range(108)if B[i,j]))for i in range(150)]
    equations=[dag.add(*(dag.mul(dag.q(c),raw[i],raw[j])for(i,j),c in sorted(row.items())))for row in quadratics]
    pivots=[]
    for block in blocks:
        U=matrix(block['left']);rr=block['rows']
        for i,diagonal_text in enumerate(block['diagonal']):
            diagonal=expression(diagonal_text)
            if not diagonal:continue
            kind='unit'if not diagonal.has(A)else'a'if s.degree(diagonal,A)==1 and diagonal.subs(A,0)==0 else'u'
            if kind=='u':assert s.degree(diagonal,A)==1 and diagonal.subs(A,s.Rational(2,9))==0
            node=dag.add(*(dag.mul(dag.expression(U[i,j]),equations[r])for j,r in enumerate(rr)if U[i,j]))
            pivots.append(dict(weight=block['weight'],column=block['columns'][i],diagonal=str(diagonal),kind=kind,source_node=node))
    assert [sum(p['kind']==kind for p in pivots)for kind in ['unit','a','u']]==[72,2,1]
    return dag,raw,equations,pivots


def closure(nodes,outputs):
    used=set(outputs);pending=list(used)
    while pending:
        op,*data=nodes[pending.pop()]
        children=data[0]if op in ('+','*')else[data[0]]if op=='inv'else[]
        for j in children:
            if j not in used:used.add(j);pending.append(j)
    return used


def eliminate(source,raw,equations,pivots,weights,a_open):
    chosen=[p for p in pivots if p['kind']=='unit'or(a_open and p['kind']=='a')]
    eliminated={p['column']for p in chosen};free=[i for i in range(108)if i not in eliminated]
    dag=DAG(['a']+[f'z{i}'for i in range(len(free))]);assigned=[dag.zero]*108
    for i,column in enumerate(free):assigned[column]=dag.inputs[f'z{i}']
    def clone(outputs):
        needed=closure(source.nodes,outputs);values={}
        for index in sorted(needed):
            op,*data=source.nodes[index]
            if op=='q':value=dag.q(data[0])
            elif op=='v':value=dag.inputs['a']if data[0]==0 else assigned[data[0]-1]
            elif op=='+':value=dag.add(*(values[j]for j in data[0]))
            elif op=='*':value=dag.mul(*(values[j]for j in data[0]))
            else:raise ValueError(op)
            values[index]=value
        return values
    for d in range(1,13):
        current=[p for p in chosen if p['weight']==d]
        if not current:continue
        values=clone([p['source_node']for p in current]);answers=[]
        for p in current:
            diagonal=expression(p['diagonal'])
            if p['kind']=='unit':inverse=dag.q(1/diagonal)
            else:inverse=dag.mul(dag.q(1/s.cancel(diagonal/A)),dag.inv(dag.inputs['a']))
            answers.append(dag.neg(dag.mul(inverse,values[p['source_node']])))
        for p,answer in zip(current,answers):assigned[p['column']]=answer
    values=clone([*raw,*equations,*[p['source_node']for p in pivots]])
    outputs=[*assigned,*[values[i]for i in raw],*[values[i]for i in equations],*[values[p['source_node']]for p in pivots]]
    keep=closure(dag.nodes,outputs);keep.update([dag.zero,dag.one,*dag.inputs.values()])
    numbering={old:new for new,old in enumerate(sorted(keep))};nodes=[]
    for old in sorted(keep):
        op,*data=dag.nodes[old]
        nodes.append((op,tuple(numbering[j]for j in data[0]))if op in ('+','*')else(op,numbering[data[0]])if op=='inv'else(op,data[0]))
    outputs=[numbering[i]for i in outputs]
    model=dict(id=f'Q_SMITH_GLOBAL_FINITE_{35 if a_open else 37}_20261103_V63',
               input_names=dag.names,input_weights=[0]+[weights[i]for i in free],
               localized_base='QQ[a,a^-1]'if a_open else'QQ[a]',nodes=nodes,
               free_weighted_columns=free,positive108_reconstruction_nodes=outputs[:108],
               raw150_reconstruction_nodes=outputs[108:258],all350_residual_nodes=outputs[258:608],
               Smith_transformed75_residual_nodes=outputs[608:],
               eliminated_Smith_pivots=[{k:v for k,v in p.items()if k!='source_node'}for p in chosen])
    assert len(model['input_names'])==(35 if a_open else 37)
    inverse_count=sum(node[0]=='inv'for node in nodes)
    assert inverse_count==(1 if a_open else 0)
    assert all(nodes[node[1]]==('v',0)for node in nodes if node[0]=='inv')
    node_weights=[];deps=[]
    for op,*data in nodes:
        if op=='q':weight=None if not Fraction(data[0])else 0;dep=0
        elif op=='v':weight=model['input_weights'][data[0]];dep=1<<data[0]
        elif op in ('+','*'):
            dep=0
            for j in data[0]:dep|=deps[j]
            if op=='+':
                options={node_weights[j]for j in data[0]if node_weights[j]is not None};assert len(options)<=1
                weight=next(iter(options),None)
            else:weight=None if any(node_weights[j]is None for j in data[0])else sum(node_weights[j]for j in data[0])
        elif op=='inv':assert deps[data[0]]==1 and node_weights[data[0]]==0;dep=1;weight=0
        node_weights.append(weight);deps.append(dep)
    for node,w in zip(model['positive108_reconstruction_nodes'],weights):assert node_weights[node]in(None,w)
    model['inverse_count']=inverse_count
    return model


def physical(raw):
    P={pair:{exp:raw[15*k+m]for m,exp in enumerate(EXPS2)if raw[15*k+m]}for k,pair in enumerate(PAIRS)}
    def entry(i,j):return P.get((i,j),{})if i<j else{m:-c for m,c in P.get((j,i),{}).items()}if i>j else{}
    derivatives={}
    for i in range(5):
        for j in range(5):
            for ell in range(5):
                value={}
                for mon,c in entry(i,j).items():
                    if mon[ell]:changed=list(mon);changed[ell]-=1;value[tuple(changed)]=c*mon[ell]
                derivatives[i,j,ell]=value
    for j in range(5):
        divergence=defaultdict(Fraction)
        for i in range(5):
            for mon,c in derivatives[i,j,i].items():divergence[mon]+=c
        assert all(not c for c in divergence.values())
    result=[]
    for i,j,k in TRIPLES:
        value=defaultdict(Fraction)
        for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
            for ell in range(5):
                for mon,co in entry(b,ell).items():
                    for other,second in derivatives[c,d,ell].items():
                        value[tuple(x+y for x,y in zip(mon,other))]+=co*second
        result.extend(value.get(mon,Fraction())for mon in EXPS3)
    return result


def validate(model,source,pivots,target,points):
    dag=DAG(model['input_names']);dag.nodes=model['nodes'];tests=[]
    eliminated={p['column']for p in model['eliminated_Smith_pivots']}
    for a_value in points:
        point={name:Fraction((-1)**i*(1+(7*i+3)%11),2+i%3)for i,name in enumerate(dag.names)};point['a']=a_value
        values=dag.evaluate(point);raw=[values[i]for i in model['raw150_reconstruction_nodes']]
        actual=physical(raw);assert actual==[values[i]for i in model['all350_residual_nodes']]
        assert all(values[node]==0 for node,p in zip(model['Smith_transformed75_residual_nodes'],pivots)if p['column']in eliminated)
        source_point={'a':a_value,**{f'w{i}':values[node]for i,node in enumerate(model['positive108_reconstruction_nodes'])}}
        source_values=source.evaluate(source_point)
        assert all(source_values[node]==v for node,v in zip(source.raw_nodes,raw))
        for i,column in enumerate(model['free_weighted_columns']):assert source_point[f'w{column}']==point[f'z{i}']
        tests.append(dict(a=str(a_value),all150_and350_match=True,divergence_zero=True,
                          selected_transformed_equations_zero=True,nonzero_full_residuals=sum(bool(v)for v in actual)))
    weights=[]
    for op,*data in model['nodes']:
        if op=='q':w=None if not Fraction(data[0])else 0
        elif op=='v':w=model['input_weights'][data[0]]
        elif op=='+':
            options={weights[i]for i in data[0]if weights[i]is not None};assert len(options)<=1;w=next(iter(options),None)
        elif op=='*':w=None if any(weights[i]is None for i in data[0])else sum(weights[i]for i in data[0])
        elif op=='inv':w=0
        weights.append(w)
    assert all(weights[node]in(None,d)for node,d in zip(model['all350_residual_nodes'],target))
    return tests


def archive(path,payload):
    with path.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(payload,separators=(',',':'))+'\n').encode())


def build(output):
    started=time.monotonic();input_path=ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz'
    data,C,q,J,target,B,back,weights,blocks,group=read_frames(input_path)
    quadratics=raw_quadratics();raw_linear=s.zeros(350,150)
    for r,row in enumerate(quadratics):
        for(i,j),c in row.items():raw_linear[r,i]+=c*q[j];raw_linear[r,j]+=c*q[i]
    assert dm(raw_linear)*dm(C)==dm(J)
    source,raw,equations,pivots=source_model(q,B,target,weights,blocks,quadratics);source.raw_nodes=raw
    print('globalframes+physicalsource ready',time.monotonic()-started,'nodes',len(source.nodes),flush=True)
    global37=eliminate(source,raw,equations,pivots,weights,False)
    open35=eliminate(source,raw,equations,pivots,weights,True)
    print('finite37/35 ready',time.monotonic()-started,'nodes',len(global37['nodes']),len(open35['nodes']),flush=True)
    tests37=validate(global37,source,pivots,target,[Fraction(0),Fraction(2,9),Fraction(-4,9)])
    tests35=validate(open35,source,pivots,target,[Fraction(2,9),Fraction(-4,9),Fraction(1,3),Fraction(44,81),Fraction(1)])
    output.mkdir(parents=True,exist_ok=False)
    frame_payload=dict(id='Q_SMITH_GLOBAL_FINITE_MODEL_20261103_V63',canonical_basis=encode(C),raw_Qa150=list(map(str,q)),
                       physical_raw350_quadratics=[[[i,j,c]for(i,j),c in sorted(row.items())]for row in quadratics],
                       raw150_positive108_basis=encode(B),positive108_from_raw150=encode(back),
                       positive108_weights=weights,Jacobi_row_weights=target,Jacobi_complement_Smith_blocks=blocks,
                       source_Smith_pivots=pivots,**group)
    archive(output/'frames.json.gz',frame_payload)
    archive(output/'source109.json.gz',dict(input_names=source.names,input_weights=[0]+weights,nodes=source.nodes,
                                          raw150_reconstruction_nodes=raw,all350_residual_nodes=equations,
                                          transformed_Smith_pivots=pivots))
    archive(output/'model37.json.gz',global37);archive(output/'model35.json.gz',open35)
    note=ROOT/'research/notes/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63.md'
    inputs=[Path(__file__),note,input_path,ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63_REVIEW.json',
            ROOT/'code/q_full_family_linear_smith_20261103_v63.py',ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py']
    report=dict(id=frame_payload['id'],all_passed=True,all_a_polynomial_model_variables=37,
                a_unit_model_variables=35,global37_new_inverses=0,a_unit35_only_inverse='a',
                unit_Smith_pivots_eliminated=72,additional_a_pivots_eliminated=2,true_u_pivot_retained='9a-2',
                source_nodes=len(source.nodes),model37_nodes=len(global37['nodes']),model35_nodes=len(open35['nodes']),
                all_raw350_physical_quadratics_saved=True,physical_linearization_matches_full126_input=True,
                all_weighted_Smith_identity_and_unit_determinants_checked=True,
                zero_section_group17_response_nonzero_QQ_constant=group['group_determinant'],
                full126_action_Jacobian_reordered_nonzero_QQ_constant=group['full126_action_reordered_determinant'],
                full126_action_Jacobian_includes_a_frame_derivatives=True,
                polynomial_global_positive_coordinate_inverse_checked=True,
                global_group_etale_along_identity_times_entire_slice=True,
                exact_ideal_elimination_by_positive_weight_triangularity=True,
                all_original350_outputs_retained=True,all_raw150_outputs_retained=True,
                dense_fraction_tests37=tests37,dense_fraction_tests35=tests35,
                frozen_parent_models_modified=False,global_group_cross_section_claimed=False,
                full_branch_or_primary_classification_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in output.glob('*.gz')},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact global37/all-a and35/a-unit finite models;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
