#!/usr/bin/env python3
"""Exact closed (a,weight7,weight8) slice: ideal (a*s,a*t), rank2 plane."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output):
    started=time.monotonic();folder=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63'
    parent_path=folder/'model37.json.gz';frame_path=folder/'frames.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()));frames=json.loads(gzip.decompress(frame_path.read_bytes()))
    records=[(i,p)for i,p in enumerate(frames['source_Smith_pivots'])if p['kind']=='a']
    assert [p['weight']for i,p in records]==[7,8]
    inputs=[parent['free_weighted_columns'].index(p['column'])+1 for i,p in records]
    assert inputs==[35,36]and[parent['input_names'][i]for i in inputs]==['z34','z35']
    R,a,u,v=ring('a,s,t',s.QQ);point=[a]+[R.zero]*36
    for i,z in zip(inputs,[u,v]):point[i]=z
    values=[]
    for op,*data in parent['nodes']:
        if op=='q':value=R(s.Rational(data[0]))
        elif op=='v':value=point[data[0]]
        elif op=='+':value=sum((values[i]for i in data[0]),R.zero)
        elif op=='*':
            value=R.one
            for i in data[0]:value*=values[i]
        else:raise ValueError(op)
        values.append(value)
    residuals=[values[i]for i in parent['all350_residual_nodes']]
    raw=[values[i]for i in parent['raw150_reconstruction_nodes']]
    generators=[a*u,a*v];forward=[]
    for polynomial in residuals:
        coefficient=[R.zero,R.zero]
        for mon,c in polynomial.items():
            assert mon[0]>=1 and mon[1]+mon[2]>=1
            if mon[1]:new=(mon[0]-1,mon[1]-1,mon[2]);column=0
            else:new=(mon[0]-1,mon[1],mon[2]-1);column=1
            coefficient[column]+=R({new:c})
        assert polynomial==sum((q*g for q,g in zip(coefficient,generators)),R.zero)
        forward.append(coefficient)
    reverse=[]
    for ordinal,(source_index,pivot)in enumerate(records):
        block=next(b for b in frames['Jacobi_complement_Smith_blocks']if b['weight']==pivot['weight'])
        row=block['columns'].index(pivot['column']);coefficients=[R.zero]*350
        for equation,text in zip(block['rows'],block['left'][row]):coefficients[equation]=R.from_expr(s.sympify(text)/6)
        assert values[parent['Smith_transformed75_residual_nodes'][source_index]]==6*generators[ordinal]
        assert sum((q*f for q,f in zip(coefficients,residuals)),R.zero)==generators[ordinal]
        reverse.append(coefficients)
    assert sum(bool(f)for f in residuals)==15
    assert all(mon[1]+mon[2]<=1 for polynomial in raw for mon in polynomial)
    X=s.symbols('x0:5');ss,tt=s.symbols('s t')
    H=X[0]*X[2]**2/3+X[1]**2*X[2]-s.Rational(4,3)*X[0]*X[1]*X[3]-s.Rational(2,9)*X[0]**2*X[4]-ss*X[3]*X[4]**2+tt*X[4]**3/3
    pairs=[(i,j)for i in range(5)for j in range(i+1,5)]
    monomials=[(i,j)for i in range(5)for j in range(i,5)]
    tensor={pair:s.S.Zero for pair in pairs}
    for pair_index,pair in enumerate(pairs):
        tensor[pair]=s.expand(sum(raw[15*pair_index+k].as_expr().subs({s.Symbol('a'):0})*X[i]*X[j]
                                   for k,(i,j)in enumerate(monomials)))
    assert all(not tensor[pair]for pair in pairs if pair not in [(2,3),(2,4),(3,4)])
    assert s.expand(tensor[2,3]-s.diff(H,X[4]))==0
    assert s.expand(tensor[2,4]+s.diff(H,X[3]))==0
    assert s.expand(tensor[3,4]-s.diff(H,X[2]))==0
    output.mkdir(parents=True,exist_ok=False)
    evidence=dict(parent_kept_input_indices=inputs,closed_variables=['a','s','t'],weights=[0,7,8],
                  all350_residuals=list(map(str,residuals)),all150_raw_coefficients=list(map(str,raw)),
                  ideal_generators=list(map(str,generators)),
                  forward350_by2=[[str(c)for c in row]for row in forward],
                  reverse2_by350=[[str(c)for c in row]for row in reverse],
                  a_zero_cubic_H=str(s.expand(H)),
                  a_zero_tensor_entries={str(pair):str(poly)for pair,poly in tensor.items()if poly})
    archive=output/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(evidence,separators=(',',':'))+'\n').encode())
    note=ROOT/'research/notes/Q_GLOBAL_HIGH_PIVOT_CLOSED_SLICE_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,frame_path,
            ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW.json',
            ROOT/'legacy/v14/sections/rank2_global_eight.tex',
            ROOT/'legacy/v14/sections/stable_rank2.tex']
    report=dict(id='Q_GLOBAL_HIGH_PIVOT_CLOSED_SLICE_20261103_V63',all_passed=True,
                exact_closed_slice_ideal='(a*s,a*t)',all350_forward_and_reverse_certificates_verified=True,
                exact_primary_decomposition='(a) intersection (s,t)',closed_total_scheme_reduced=True,
                no_embedded_prime_in_closed_slice=True,nonzero_residual_count=15,
                nonzero_a_fiber='one reduced point',a_zero_fiber='entire reduced affine plane A2',
                a_zero_plane_actual_generic_Poisson_rank=2,
                a_zero_plane_exact_Jacobian_cubic_presentation=True,
                a_zero_plane_in_existing_Rat113_component_by_explicit_Jacobian_map=True,
                closed_slice_nonflat_over_a=True,
                full_parent_nonflatness_at_a_zero_inferred=False,
                full_Qa_zero_local_branches_classified=False,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={'certificate.json.gz':sha(archive)},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact high-pivot closed slice (a*s,a*t), actual rank2 a0-plane;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
