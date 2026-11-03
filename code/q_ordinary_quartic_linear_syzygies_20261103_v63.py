#!/usr/bin/env python3
"""Exact101 linear syzygies of57Q, carried to true quartics and15-planes."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import json
from math import gcd,lcm
from pathlib import Path
import resource
import shlex
import subprocess
import time

from q_true_initial_eta_open_20261103_v63 import compose,pack,unpack,mul
from q_ordinary_cubic_initial_20261103_v63 import subtract_scaled,rank_QQ

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--timeout',type=int,default=180);args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2500000000,2500000000))
    cpath=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/cubic-initials.json.gz';c=json.loads(gzip.decompress(cpath.read_bytes()))
    qpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz';q=json.loads(gzip.decompress(qpath.read_bytes()))
    mpath=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz';model=json.loads(gzip.decompress(mpath.read_bytes()))
    ppath=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz';plane=json.loads(gzip.decompress(ppath.read_bytes()))
    rpath=ROOT/'reports/Q_ORDINARY_CUBIC_RANK_QQ_20261103_v63/report.json'
    if not rpath.exists():
        candidates=list((ROOT/'reports').glob('Q_ORDINARY_CUBIC*/**/report.json'))
        rpath=next(p for p in candidates if 'quadratic_linear_multiple_rank_QQ'in json.loads(p.read_text()))
    ranks=json.loads(rpath.read_text());rr={r['old_weight']:r['QQ_quadratic_multiple_rank']for r in ranks['block_ranks']}
    quadrics=[{(i,j):F(v)for i,j,v in row}for row in q['all350_quadratic_rows']]
    cubics=[unpack(row)for row in c['raw350_cubic_terms']]
    weights=model['input_weights'];blocks=defaultdict(list)
    for source in c['quadratic_linear_multiples_recipe']['quadratic_original_rows']:
        for variable in range(35):
            p={tuple(sorted(m+(variable,))):v for m,v in quadrics[source].items()}
            blocks[model['residual_weight_bounds'][source]+weights[variable]].append({'source':source,'variable':variable,'polynomial':p})
    inp=args.output/'integer-blocks.in';metadata={}
    with inp.open('w')as out:
        out.write(str(len(blocks))+'\n')
        for weight,rows in sorted(blocks.items()):
            mons=sorted({m for row in rows for m in row['polynomial']});ix={m:i for i,m in enumerate(mons)}
            out.write(f'{weight} {len(rows)} {len(mons)} {rr[weight]}\n')
            for row in rows:
                p=row['polynomial'];den=lcm(*(v.denominator for v in p.values()))
                nums={ix[m]:v.numerator*(den//v.denominator)for m,v in p.items()};content=gcd(*nums.values())
                row['integer_row_scale']=F(den,content)
                out.write(str(len(nums))+' '+' '.join(f'{i} {v//content}'for i,v in sorted(nums.items()))+'\n')
            metadata[weight]=rows
    native=ROOT/'code/q_ordinary_quartic_linear_syzygies_20261103_v63.cpp';exe=args.output/'kernel'
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],check=True,capture_output=True,text=True).stdout)
    subprocess.run(['c++','-std=c++17','-O3',str(native),*flags,'-o',str(exe)],check=True)
    result=args.output/'kernels.txt'
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(2500000000,2500000000));resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    with (args.output/'native.log').open('w')as log:
        completed=subprocess.run([str(exe.resolve()),str(inp.resolve()),str(result.resolve())],stdout=log,stderr=subprocess.STDOUT,
                                 timeout=args.timeout,preexec_fn=limits)
    assert completed.returncode==0
    tokens=iter(result.read_text().split());syzygies=[]
    for weight,rows in sorted(blocks.items()):
        assert int(next(tokens))==weight and int(next(tokens))==len(rows);nullity=int(next(tokens))
        for _ in range(nullity):
            count=int(next(tokens));coeff={int(next(tokens)):F(int(next(tokens)))for _ in range(count)}
            exact={};recipe=[]
            for i,h in coeff.items():
                row=rows[i];h*=row['integer_row_scale']
                subtract_scaled(exact,row['polynomial'],-h)
                recipe.append([row['source'],row['variable'],str(h)])
            assert not exact
            syzygies.append({'weight':weight,'linear_actual350_cofactor_recipe':recipe})
    assert next(tokens)=='TOTAL'and int(next(tokens))==101 and len(syzygies)==101
    names=q['variables'];h=names.index('z22');branchresults=[]
    for b in plane['branches']:
        graph={i:{(h,)if not m else tuple(m):F(v)for m,v in data}for i,data in enumerate(b['all35_graph'])}
        restrictedC={source:compose(cubics[source],graph)for source in c['quadratic_linear_multiples_recipe']['quadratic_original_rows']}
        restrictions=[]
        for syz in syzygies:
            target={}
            for source,var,coef in syz['linear_actual350_cofactor_recipe']:
                part=mul(graph[var],restrictedC[source]);subtract_scaled(target,part,-F(coef))
            restrictions.append(target)
        basis,piv=rank_QQ(restrictions)
        branchresults.append({'branch':b['branch'],'all101_true_quartic_restrictions':list(map(pack,restrictions)),
                              'nonzero':sum(bool(p)for p in restrictions),'QQ_rank':len(basis),
                              'QQ_basis':[pack(p)for _,p in sorted(basis.items())],'selected_syzygies':piv})
    payload={'variables':names,'all101_exact_linear_syzygies':syzygies,'plane_results':branchresults}
    cp=args.output/'certificate.json.gz'
    with gzip.open(cp,'wt')as out:json.dump(payload,out,separators=(',',':'))
    inputs=[Path(__file__),native,cpath,qpath,mpath,ppath,rpath]
    report={'id':'Q_ORDINARY_QUARTIC_LINEAR_SYZYGIES_20261103_V63','all_passed':True,'coefficient_field':'QQ',
            'exact_linear_quadratic_syzygies':101,'every_global_degree3_cancellation_replayed':True,
            'actual350_cofactor_degree':1,'quartic_formula':'sum h_(source,var)*x_var*CubicTail(F_source)',
            'branches':[{k:r[k]for k in ['branch','nonzero','QQ_rank','selected_syzygies']}for r in branchresults],
            'complete_quartic_initial_space_or_dimension_bound_claimed':False,
            'seconds':time.monotonic()-started,'peak_rss_kb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'certificate_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),
            'input_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs}}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
