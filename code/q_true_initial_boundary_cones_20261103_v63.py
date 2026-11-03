#!/usr/bin/env python3
"""Two explicit 15-dimensional cones in (57Q,17G), not the full initial ideal."""
import argparse,gzip,hashlib,json,time
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
QP=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
GP=ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/generators.json.gz'

def add(a,b,f=F(1)):
    c=a.copy()
    for m,x in b.items():
        z=c.get(m,0)+f*x
        if z:c[m]=z
        elif m in c:del c[m]
    return c

def mul(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            p=tuple(sorted(m+n));z=c.get(p,0)+x*y
            if z:c[p]=z
            elif p in c:del c[p]
    return c

def substitute(p,i,r):
    out={}
    for m,x in p.items():
        e=m.count(i);term={tuple(j for j in m if j!=i):x}
        for _ in range(e):term=mul(term,r)
        out=add(out,term)
    return out

def encode(p):return [[list(m),str(c)]for m,c in sorted(p.items())]

def main(output):
    start=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    q=json.load(gzip.open(QP,'rt'));g=json.load(gzip.open(GP,'rt'))
    names=q['variables'];idx={n:i for i,n in enumerate(names)}
    allQ=[{(i,j):F(c)for i,j,c in row}for row in q['all350_quadratic_rows']]
    cs=[{tuple(m):F(c)for m,c in row['terms']}for row in g['cubic_generators']]
    original=[allQ[i]for i in g['quadratic_generators']]+cs
    packets=[]
    for branch in [0,1]:
        zero=['z11','z12','z23','z33','z0','z2','z8','z10','z21','z32','z30','z18','z7']
        assignments=[(idx[n],{})for n in zero]
        assignments += [(idx['z22'],{():F(1)}),(idx['z20'],{(idx['z31'],):F(-2)}),
            (idx['z6'],{}if branch==0 else{(idx['z27'],):F(5,7)})]
        rows=[p.copy()for p in original];assigned=set();operations=[]
        for i,rhs in assignments:
            rows=[substitute(p,i,rhs)for p in rows];assigned.add(i)
        # Each step uses an actual equation with constant nonzero coefficient.
        while True:
            options=[]
            for j,p in enumerate(rows):
                for m,c in p.items():
                    if len(m)!=1:continue
                    i=m[0]
                    if i in assigned or any(i in n for n in p if n!=m):continue
                    options.append((len(p),i,j,c))
            if not options:break
            _,i,j,c=min(options)
            rhs={m:-x/c for m,x in rows[j].items()if m!=(i,)}
            operations.append(dict(variable=names[i],equation=j,coefficient=str(c),rhs=encode(rhs)))
            assignments.append((i,rhs));assigned.add(i)
            rows=[substitute(p,i,rhs)for p in rows]
        assert not any(rows)
        free=[i for i in range(35)if i not in assigned]
        assert len(free)==14
        graph={i:{(i,):F(1)}for i in free}
        for i,rhs in reversed(assignments):
            r=rhs.copy()
            for k,expression in graph.items():
                if any(k in m for m in r):r=substitute(r,k,expression)
            assert all(set(m)<=set(free)for m in r)
            graph[i]=r
        assert graph[idx['z22']]=={():F(1)}
        for p in allQ+cs:
            r=p.copy()
            for i in range(35):
                if i not in free:r=substitute(r,i,graph[i])
            assert not r
        packets.append(dict(branch=branch,equation='z6=0'if branch==0 else'7*z6=5*z27',
            free_variables=[names[i]for i in free],affine_chart_dimension=14,
            cone_dimension=15,linear_eliminations=operations,
            all35_graph=[encode(graph[i])for i in range(35)],
            all350_original_quadrics_and17_true_cubics_vanish_exactly=True))
    assert packets[0]['free_variables']==packets[1]['free_variables']
    proof=dict(id='Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63',variables=names,branches=packets,
        cone_map='x_i=h*graph_i(p), h invertible; inverse h=z22, p_i=x_i/z22 for listed14free',
        coefficient_field='QQ',original_J_generators=74)
    with(output/'certificate.json.gz').open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0)as z:z.write((json.dumps(proof,sort_keys=True,separators=(',',':'))+'\n').encode())
    paths=[Path(__file__),QP,GP,ROOT/'research/notes/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63.md']
    report=dict(id=proof['id'],all_passed=True,exact_parameterized_cones=2,
        each_chart='A14, z22=1',each_cone_open='Gm x A14',
        truncated_initial_lower_ideal_J_dimension_lower_bound=15,
        every12_homogeneous_linear_section_dimension_lower_bound=3,
        fixed_prior12_section_Artinian_target_impossible=True,
        full_initial_ideal_or_actual_local_dimension_lower_bound15_claimed=False,
        full_initial_or_actual_model_dimension_upper12_ruled_out=False,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS two exact15D cones in J; 12-linear Artinian target impossible',report['seconds'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    main(p.parse_args().output)
