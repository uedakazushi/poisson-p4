#!/usr/bin/env python3
"""One fixed quartic scalar membership target, split by the preserving torus.

No power search: only D(p0)^2 is tested. All Jacobi generators are retained;
multigrading merely omits Macaulay columns unable to occur in this target.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,time
from classification_projected_veronese_zero_review import monomials,plus,scale
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2))


def weight(mon,weights):return tuple(sum(weights[i][k]for i in mon)for k in range(2))


def compute(case,source):
    started=time.monotonic();name=case['name'];smooth=name=='smooth_scroll_12'
    coordinate_weights=[(1,0),(1,1),(0,0),(0,1),(0,2)]if smooth else[(0,0),(0,1),(0,2),(0,3),(1,0)]
    quadric_weights=[(1,1),(1,2),(0,2)]if smooth else[(0,2),(0,3),(0,4)]
    weights=[tuple(coordinate_weights[i][k]+coordinate_weights[j][k]-q[k]for k in range(2))for q in quadric_weights for i,j in PAIRS]
    quadcols=[tuple(v)for v in source['all465_quadratic_parameter_monomials']]
    selected=case['selected_independent_original350_rows'];generators={i:{quadcols[int(j)]:Q(c)for j,c in case['all350_Jacobi_rows'][i].items()}for i in selected}
    generator_weights={}
    for i,row in generators.items():
        chars={weight(mon,weights)for mon in row};assert len(chars)==1;generator_weights[i]=next(iter(chars))
    point=[0,0,0,1,0]if smooth else[1,0,0,1,0];D0={}
    for exponent,row in zip(monomials(5,5),case['all126_actual_D_rows']):
        coefficient=1
        for a,e in zip(point,exponent):coefficient*=a**e
        if coefficient:D0=plus(D0,{quadcols[int(j)]:coefficient*Q(c)for j,c in row.items()})
    target={}
    for a,c in D0.items():
        for b,d in D0.items():target=plus(target,{tuple(sorted(a+b)):c*d})
    targets={}
    for mon,c in target.items():targets.setdefault(weight(mon,weights),{})[mon]=c
    allmon=list(combinations_with_replacement(range(30),4));results=[]
    for character,subtarget in sorted(targets.items()):
        columns=[mon for mon in allmon if weight(mon,weights)==character];ci={mon:j for j,mon in enumerate(columns)}
        t={ci[mon]:c for mon,c in subtarget.items()};basis={};proof={};chosen=[]
        def reduce(v,w=None):
            w={}if w is None else w
            for p in sorted(basis):
                if p in v:
                    c=v[p];v=plus(v,scale(basis[p],-c));w=plus(w,scale(proof[p],-c))
            return v,w
        count=0;certificate=None;labels=[]
        for i,g in generators.items():
            required=tuple(character[k]-generator_weights[i][k]for k in range(2))
            for mon2 in quadcols:
                if weight(mon2,weights)!=required:continue
                row={ci[tuple(sorted(mon+mon2))]:c for mon,c in g.items()};assert len(row)==len(g)
                labels.append((i,mon2));v,w=reduce(row,{count:Q(1)});count+=1
                if v:
                    p=min(v);a=v[p];basis[p]=scale(v,1/a);proof[p]=scale(w,1/a);chosen.append(count-1)
        residual,minusfunctional=reduce(t)
        contained=not residual
        if contained:
            functional=scale(minusfunctional,-1);expanded={}
            for k,c in functional.items():
                i,mon2=labels[k]
                expanded=plus(expanded,{tuple(sorted(mon+mon2)):c*d for mon,d in generators[i].items()})
            assert expanded==subtarget
            certificate=[{'original350_Jacobi_row':labels[k][0],'quadratic_multiplier':labels[k][1],'coefficient':str(c)}for k,c in sorted(functional.items())]
        results.append({'target_character':character,'quartic_columns':columns,'Macaulay_row_count':count,'rank':len(basis),
            'D0_square_weight_component_in_Jacobi_ideal':contained,'exact_original350_membership_functional':certificate,
            'target_quartic_coefficients':{str(j):str(c)for j,c in t.items()},'exact_reduced_residual':{str(j):str(c)for j,c in residual.items()}})
        print(name,'weight',character,'Macaulay',count,len(columns),'rank',len(basis),'target',contained,flush=True)
    return {'name':name,'coordinate_weights':coordinate_weights,'quadric_weights':quadric_weights,'all30_parameter_weights':weights,
        'actual_D_at_p0':[{ 'parameter_monomial':mon,'coefficient':str(c)}for mon,c in sorted(D0.items())],
        'fixed_exponent':2,'blocks':results,'D0_squared_exactly_in_full_Jacobi_ideal':all(b['D0_square_weight_component_in_Jacobi_ideal']for b in results),
        'seconds':time.monotonic()-started,
        'scope':'A positive answer proves D(p0)=0 on all Poisson sources and, by the dense preserving point orbit, excludes all rank4 tensors. A negative answer is only one failed quartic scalar membership criterion. No higher powers are searched.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('surface',choices=['smooth_scroll_12','twisted_cubic_cone_03']);args=ap.parse_args()
    destination=ROOT/f'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_SQUARE_{args.surface.upper()}_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));signal.alarm(100)
    sourcepath=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_JACOBI_20261107.json';source=json.loads(sourcepath.read_text());case=next(c for c in source['cases']if c['name']==args.surface)
    result={'task_id':'G','input_report':str(sourcepath.relative_to(ROOT)),'input_sha256':hashlib.sha256(sourcepath.read_bytes()).hexdigest(),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'case':compute(case,source)}
    result['maximum_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;destination.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
