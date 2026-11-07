#!/usr/bin/env python3
"""Whole rank-zero sources on both integral degree-three surfaces in P4.

Exact rational row membership in all350 projective Jacobi equations, with all
30 source parameters and all465 parameter quadrics. No point/prime sampling.
"""
from pathlib import Path
from itertools import combinations, combinations_with_replacement
from fractions import Fraction as Q
import hashlib,json,resource,signal,time
from classification_projected_veronese_zero_review import Ring,monomials,plus,scale,sparse_rref

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));TRIPLES=list(combinations(range(5),3));QUADS=list(combinations(range(5),4))
COLS=list(combinations_with_replacement(range(30),2));COLINDEX={v:i for i,v in enumerate(COLS)}
CASES={
    'smooth_scroll_12':['x0*x3-x1*x2','x0*x4-x1*x3','x2*x4-x3**2'],
    'twisted_cubic_cone_03':['x0*x2-x1**2','x0*x3-x1*x2','x1*x3-x2**2'],
}


def add(target,key,value):
    target[key]=target.get(key,0)+value
    if not target[key]:del target[key]


def product(first,second,derivative=None):
    out={}
    for (p,e),a in first.items():
        for (q,f),b in second.items():
            if derivative is not None:
                if not f[derivative]:continue
                g=list(f);b*=g[derivative];g[derivative]-=1;f=tuple(g)
            add(out,(tuple(sorted((p,q))),tuple(i+j for i,j in zip(e,f))),a*b)
    return out


def compute(name,texts):
    X=Ring([f'x{i}'for i in range(5)]);quadrics=[X.parse(t)for t in texts]
    P={}
    for k,ij in enumerate(PAIRS):
        P[ij]={(10*g+k,e):c for g,q in enumerate(quadrics)for e,c in q.items()}
        P[ij[1],ij[0]]={v:-c for v,c in P[ij].items()}
    strict={tri:{}for tri in TRIPLES}
    for tri in TRIPLES:
        for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]:
            for l in range(5):
                for key,value in product(P.get((a,l),{}),P[(b,c)],l).items():add(strict[tri],key,value)
    J=[];rowlabels=[]
    for quad in QUADS:
        poly={}
        for pos,missing in enumerate(quad):
            tri=quad[:pos]+quad[pos+1:]
            for (pq,e),value in strict[tri].items():
                f=list(e);f[missing]+=1;add(poly,(pq,tuple(f)),(-1)**pos*value)
        for e in monomials(5,4):
            J.append({COLINDEX[pq]:value for (pq,f),value in poly.items()if e==f});rowlabels.append({'quad':quad,'exponents':e})
    D={}
    for missing in range(5):
        a,b,c,d=[i for i in range(5)if i!=missing]
        for ij,kl,sign in [((a,b),(c,d),1),((a,c),(b,d),-1),((a,d),(b,c),1)]:
            for (pq,e),value in product(P[ij],P[kl]).items():
                f=list(e);f[missing]+=1;add(D,(pq,tuple(f)),(-1)**missing*sign*value)
    Drows=[{COLINDEX[pq]:value for (pq,f),value in D.items()if e==f}for e in monomials(5,5)]
    basis,proof,chosen=sparse_rref(J,True);certificates=[];all_contained=True
    for e,target in zip(monomials(5,5),Drows):
        if not target:continue
        residual={j:Q(c)for j,c in target.items()};w={}
        for p in sorted(basis):
            if p in residual:
                coefficient=residual[p];residual=plus(residual,scale(basis[p],-coefficient));w=plus(w,scale(proof[p],coefficient))
        all_contained &=not residual
        reconstructed=plus(*(scale(J[i],c)for i,c in w.items()))
        assert plus(reconstructed,residual)==target
        certificates.append({'exponents':e,'original350_Jacobi_functional':w,'residual':residual})
    print(name,'Jac rank',len(basis),'actual D contained',all_contained,'nonzero D coefficients',len(certificates),flush=True)
    return {'name':name,'quadratic_ideal_generators':texts,'full_source_dimension':30,'all350_Jacobi_rows':J,
        'all126_actual_D_rows':Drows,'full350_Jacobi_row_labels':rowlabels,'Jacobi_row_rank':len(basis),
        'selected_independent_original350_rows':chosen,'exact_RREF_rows':basis,'exact_RREF_original350_functionals':proof,
        'all_actual_D_coefficients_in_full350_Jacobi_rowspace':all_contained,'actual_D_coefficient_certificates':certificates,
        'scope':'The entire30-dimensional rankzero section space on this whole projective surface orbit. Positive membership implies no rank4 Poisson tensor vanishes on this surface. A negative residual disproves only the linear membership test; it is not a Poisson point or a nonexistence proof.'}


def main():
    destination=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_JACOBI_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    cases=[compute(name,texts)for name,texts in CASES.items()]
    result={'task_id':'G','all30_source_parameters_retained':True,'all350_Jacobi_equations_retained':True,
        'all465_quadratic_parameter_monomials':COLS,'parameter_order':[[g,i,j]for g in range(3)for i,j in PAIRS],
        'cases':cases,'both_orbits_actual_D_zero_on_all_Poisson_sources':all(c['all_actual_D_coefficients_in_full350_Jacobi_rowspace']for c in cases),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_sha256':{'code/classification_projected_veronese_zero_review.py':hashlib.sha256((ROOT/'code/classification_projected_veronese_zero_review.py').read_bytes()).hexdigest()},
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(json.dumps({'both_orbits_no_rank4':result['both_orbits_actual_D_zero_on_all_Poisson_sources'],'seconds':result['seconds']}))


if __name__=='__main__':main()
