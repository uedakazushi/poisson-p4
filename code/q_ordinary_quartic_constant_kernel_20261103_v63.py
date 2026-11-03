#!/usr/bin/env python3
"""True quartics from exact constant dependencies among frozen cubic tails."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import resource
import time

from q_ordinary_cubic_initial_20261103_v63 import subtract_scaled,rank_QQ
from q_true_initial_eta_open_20261103_v63 import compose,pack,unpack

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--quartic',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_CPU,(180,181))
    qpath=args.quartic.resolve()/'quartic-tails.json.gz';j=json.loads(gzip.decompress(qpath.read_bytes()))
    modelpath=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz';model=json.loads(gzip.decompress(modelpath.read_bytes()))
    plane_path=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz';plane=json.loads(gzip.decompress(plane_path.read_bytes()))
    oldqpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz';oldq=json.loads(gzip.decompress(oldqpath.read_bytes()))
    oldcpath=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/cubic-initials.json.gz';oldc=json.loads(gzip.decompress(oldcpath.read_bytes()))
    rows={r['source']:{'C':unpack(r['cubic_terms']),'D':unpack(r['quartic_terms']),
                       'CF':{(k,):F(c)for k,c in r['constant_full350_cofactors']}}
          for r in j['all293_constant_quadratic_cancellation_rows']}
    blocks=defaultdict(list)
    for k,r in rows.items():
        if r['C']:blocks[model['residual_weight_bounds'][k]].append(k)
    kernels=[];block_report=[]
    for weight,sources in sorted(blocks.items()):
        basis={};relations=[]
        for source in sorted(sources,key=lambda k:(len(rows[k]['C']),k)):
            p=rows[source]['C'].copy();cf={(source,):F(1)}
            while p:
                pivot=min(p);c=p[pivot]
                if pivot not in basis:
                    basis[pivot]=({m:a/c for m,a in p.items()},{m:a/c for m,a in cf.items()});break
                subtract_scaled(p,basis[pivot][0],c);subtract_scaled(cf,basis[pivot][1],c)
            else:relations.append(cf)
        block_report.append({'weight':weight,'rows':len(sources),'rank':len(basis),'kernel':len(relations)})
        for relation in relations:
            C={};D={};CF={}
            for (source,),c in relation.items():
                subtract_scaled(C,rows[source]['C'],-c);subtract_scaled(D,rows[source]['D'],-c)
                subtract_scaled(CF,rows[source]['CF'],-c)
            assert not C
            kernels.append({'weight':weight,'cubic_tail_relation':relation,'quartic':D,'CF':CF})
        print('C_KERNEL_BLOCK',block_report[-1],time.monotonic()-started,flush=True)
    assert len(kernels)==55 and sum(b['rank']for b in block_report)==80
    # Replay each true-4 initial directly from the original actual350 tails,
    # checking degree2 and degree3 cancel EXACTLY before calling it quartic.
    original=[{(i,k):F(c)for i,k,c in row}for row in oldq['all350_quadratic_rows']]
    cubic=[unpack(row)for row in oldc['raw350_cubic_terms']]
    quartic=[unpack(row)for row in j['all350_quartic_terms']]
    for n,row in enumerate(kernels):
        for polys,target in [(original,{}),(cubic,{}),(quartic,row['quartic'])]:
            total={}
            for (k,),c in row['CF'].items():subtract_scaled(total,polys[k],-c)
            assert total==target
    names=j['variables'];h=names.index('z22');branch_results=[]
    for b in plane['branches']:
        graph={i:{(h,)if not m else tuple(m):F(c)for m,c in data}for i,data in enumerate(b['all35_graph'])}
        restrictions=[compose(row['quartic'],graph)for row in kernels]
        basis,piv=rank_QQ(restrictions)
        branch_results.append({'branch':b['branch'],'all55_restrictions':list(map(pack,restrictions)),
                               'nonzero':sum(bool(p)for p in restrictions),'QQ_rank':len(basis),
                               'QQ_basis':[pack(p)for _,p in sorted(basis.items())],
                               'selected_kernel_rows':piv})
    payload={'variables':names,'constant_cubic_kernel_quartics':[
        {'weight':r['weight'],'constant_cubic_tail_relation':[[m[0],str(c)]for m,c in sorted(r['cubic_tail_relation'].items())],
         'true_quartic':pack(r['quartic']),'constant_actual350_cofactors':[[m[0],str(c)]for m,c in sorted(r['CF'].items())]}
        for r in kernels],'plane_results':branch_results}
    cp=args.output/'certificate.json.gz'
    with gzip.open(cp,'wt')as handle:json.dump(payload,handle,separators=(',',':'))
    inputs=[Path(__file__),qpath,modelpath,plane_path,oldqpath,oldcpath]
    report={'id':'Q_ORDINARY_QUARTIC_CONSTANT_KERNEL_20261103_V63','all_passed':True,'coefficient_field':'QQ',
            '135_nonzero_cubic_tail_rank':80,'exact_constant_cubic_dependencies':55,'block_ranks':block_report,
            'nonzero_new_true_quartics':sum(bool(r['quartic'])for r in kernels),
            'all55_actual350_cofactor_lower2_lower3_zero_and_degree4_replay':True,
            'branch_results':[{k:r[k]for k in ['branch','nonzero','QQ_rank','selected_kernel_rows']}for r in branch_results],
            'complete_degree4_initial_or_dimension_bound_claimed':False,
            'seconds':time.monotonic()-started,'peak_rss_kb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'certificate_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),
            'input_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs}}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
