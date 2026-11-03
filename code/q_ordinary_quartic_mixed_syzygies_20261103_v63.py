#!/usr/bin/env python3
"""All219 cubic-space dependencies including the63 mixed Q/C relations."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip,hashlib,json
from math import gcd,lcm
from pathlib import Path
import resource,shlex,subprocess,time
from q_true_initial_eta_open_20261103_v63 import compose,pack,unpack,mul
from q_ordinary_cubic_initial_20261103_v63 import subtract_scaled,rank_QQ
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--quartic',required=True,type=Path);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--timeout',type=int,default=180);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2500000000,2500000000))
    cpath=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/cubic-initials.json.gz';qpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
    mpath=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz';ppath=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz'
    rpath=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/exact-rank/report.json';tpath=args.quartic.resolve()/'quartic-tails.json.gz'
    load=lambda p:json.loads(gzip.decompress(p.read_bytes()))
    c,q,model,planes,tails=map(load,[cpath,qpath,mpath,ppath,tpath]);ranks=json.loads(rpath.read_text())
    rr={r['old_weight']:r['QQ_cubic_initial_rank']for r in ranks['block_ranks']};weights=model['input_weights'];blocks=defaultdict(list)
    Q=[{(i,j):F(v)for i,j,v in row}for row in q['all350_quadratic_rows']];C=[unpack(row)for row in c['raw350_cubic_terms']]
    G={r['source']:{'C':unpack(r['cubic_terms']),'D':unpack(r['quartic_terms']),'CF':r['constant_full350_cofactors']}
       for r in tails['all293_constant_quadratic_cancellation_rows']}
    for source in c['quadratic_linear_multiples_recipe']['quadratic_original_rows']:
        for variable in range(35):
            p={tuple(sorted(m+(variable,))):v for m,v in Q[source].items()}
            blocks[model['residual_weight_bounds'][source]+weights[variable]].append({'type':'W','source':source,'variable':variable,'P':p})
    for source,row in G.items():
        if row['C']:blocks[model['residual_weight_bounds'][source]].append({'type':'G','source':source,'P':row['C']})
    inp=args.output/'integer-blocks.in'
    with inp.open('w')as out:
        out.write(str(len(blocks))+'\n')
        for weight,rows in sorted(blocks.items()):
            mons=sorted({m for row in rows for m in row['P']});ix={m:i for i,m in enumerate(mons)};out.write(f'{weight} {len(rows)} {len(mons)} {rr[weight]}\n')
            for row in rows:
                den=lcm(*(v.denominator for v in row['P'].values()));nums={ix[m]:v.numerator*(den//v.denominator)for m,v in row['P'].items()};content=gcd(*nums.values())
                row['scale']=F(den,content);out.write(str(len(nums))+' '+' '.join(f'{i} {v//content}'for i,v in sorted(nums.items()))+'\n')
    native=ROOT/'code/q_ordinary_quartic_linear_syzygies_20261103_v63.cpp';exe=args.output/'kernel'
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],check=True,capture_output=True,text=True).stdout)
    subprocess.run(['c++','-std=c++17','-O3',str(native),*flags,'-o',str(exe)],check=True);result=args.output/'kernels.txt'
    def limits():resource.setrlimit(resource.RLIMIT_AS,(2500000000,2500000000));resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
    with (args.output/'native.log').open('w')as log:
        done=subprocess.run([str(exe.resolve()),str(inp.resolve()),str(result.resolve())],stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout,preexec_fn=limits)
    assert done.returncode==0;tokens=iter(result.read_text().split());syzygies=[]
    for weight,rows in sorted(blocks.items()):
        assert int(next(tokens))==weight and int(next(tokens))==len(rows);nullity=int(next(tokens))
        for _ in range(nullity):
            count=int(next(tokens));coeff={int(next(tokens)):F(int(next(tokens)))for _ in range(count)};zero={};recipe=[];actual=defaultdict(dict)
            for i,h in coeff.items():
                row=rows[i];h*=row['scale'];subtract_scaled(zero,row['P'],-h)
                if row['type']=='W':
                    recipe.append(['W',row['source'],row['variable'],str(h)]);subtract_scaled(actual[row['source']],{(row['variable'],):h},F(-1))
                else:
                    recipe.append(['G',row['source'],str(h)])
                    for k,v in G[row['source']]['CF']:subtract_scaled(actual[k],{():h*F(v)},F(-1))
            assert not zero
            syzygies.append({'weight':weight,'row_recipe':recipe,'actual350_cofactors':[[k,pack(p)]for k,p in sorted(actual.items())if p]})
    assert next(tokens)=='TOTAL'and int(next(tokens))==219 and len(syzygies)==219
    # Verify the actualF cofactor expressions, rather than relying only on
    # the cubic-space row kernel.  Each actual coefficient is constant plus
    # linear; check all degree2 and degree3 cancellations in the original35
    # variables BEFORE restricting any quartic to a plane.
    for syz in syzygies:
        lower2={};lower3={}
        for source,data in syz['actual350_cofactors']:
            for mon,coef in unpack(data).items():
                if not mon:
                    subtract_scaled(lower2,Q[source],-coef);subtract_scaled(lower3,C[source],-coef)
                else:
                    assert len(mon)==1
                    product={tuple(sorted(m+mon)):v for m,v in Q[source].items()}
                    subtract_scaled(lower3,product,-coef)
        assert not lower2 and not lower3
    names=q['variables'];h=names.index('z22');branchresults=[]
    for b in planes['branches']:
        graph={i:{(h,)if not m else tuple(m):F(v)for m,v in data}for i,data in enumerate(b['all35_graph'])}
        restrictedC={s:compose(C[s],graph)for s in c['quadratic_linear_multiples_recipe']['quadratic_original_rows']}
        restrictedG={s:compose(r['D'],graph)for s,r in G.items()if r['C']};restrictions=[]
        for syz in syzygies:
            target={}
            for term in syz['row_recipe']:
                if term[0]=='W':_,source,var,coef=term;part=mul(graph[var],restrictedC[source])
                else:_,source,coef=term;part=restrictedG[source]
                subtract_scaled(target,part,-F(coef))
            restrictions.append(target)
        basis,piv=rank_QQ(restrictions);branchresults.append({'branch':b['branch'],'all219_restrictions':list(map(pack,restrictions)),
            'nonzero':sum(bool(p)for p in restrictions),'QQ_rank':len(basis),'QQ_basis':[pack(p)for _,p in sorted(basis.items())],'selected_syzygies':piv})
    cp=args.output/'certificate.json.gz'
    with gzip.open(cp,'wt')as out:json.dump({'variables':names,'all219_exact_syzygies':syzygies,'plane_results':branchresults},out,separators=(',',':'))
    inputs=[Path(__file__),native,cpath,qpath,mpath,ppath,rpath,tpath]
    report={'id':'Q_ORDINARY_QUARTIC_MIXED_SYZYGIES_20261103_V63','all_passed':True,'coefficient_field':'QQ','exact_full_nonzero_cubic_space_dependencies':219,
        'all_degree3_cancellations_replayed_globally':True,'actual350_cofactor_polynomials_saved':True,
        'all219_actual350_lower2_and_lower3_polynomial_cancellations_replayed':True,
        'branches':[{k:r[k]for k in ['branch','nonzero','QQ_rank','selected_syzygies']}for r in branchresults],
        'full_degree4_initial_space_claimed':False,'local_dimension_bound_claimed':False,'seconds':time.monotonic()-started,
        'peak_rss_kb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'certificate_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),
        'input_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs}}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
