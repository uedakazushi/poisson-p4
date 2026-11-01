#!/usr/bin/env python3
"""Parent review: independent local data, saved quartic ideal, and proof audit."""
import argparse,hashlib,json,sys,time
from pathlib import Path
from itertools import combinations_with_replacement
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate

def review_case(a):
    start=time.time();p=32003
    model=S61();two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    point=s.zeros(61,1)
    for k,v in a['sample'].items():point[int(k[1:])-1]=s.Rational(v)
    point[29]=1
    r0=calculate(model,True,point)
    for i,v in zip(two,a['same_beta_regular_point_q']):point[i]=s.Rational(v)
    r1=calculate(model,True,point)
    stats=lambda r:[r['homogenized_jacobian_rank'],r['quadratic_span_rank'],len(r['tangent_basis'][0])-r['common_hessian_rank'],r['homogenizer_zero_on_tangent']]
    assert stats(r0)==[34,14,2,True] and stats(r1)==[35,5,3,True]
    # Independent polynomial-dictionary construction of the saved, NONHOMOGENEOUS ideal.
    n=8;mon=[tuple(v) for d in [2,3,4] for v in combinations_with_replacement(range(n),d)]
    idx={v:i for i,v in enumerate(mon)}
    base=np.concatenate([np.array(a[k],dtype=np.int64) for k in ['residual_quadratic_mod_prime','residual_cubic_mod_prime','residual_quartic_mod_prime']],axis=1)
    multipliers=[()]+[(i,) for i in range(n)]+list(combinations_with_replacement(range(n),2))
    rows=[]
    for f in base:
        support=[(mon[i],int(f[i])) for i in np.flatnonzero(f)]
        for v in multipliers:
            out=np.zeros(len(mon),dtype=np.int64)
            for w,c in support:
                if len(v)+len(w)<=4:out[idx[tuple(sorted(v+w))]]+=c
            if np.any(out):rows.append(out%p)
    M=np.array(rows,dtype=np.int64)
    # Reverse both orders and use forward elimination, independent of saved pivots.
    def rank(M):
        M=M[::-1,::-1].copy()%p;r=0
        for c in range(M.shape[1]):
            nonzero=np.flatnonzero(M[r:,c])
            if not len(nonzero):continue
            j=r+int(nonzero[-1]);M[[r,j]]=M[[j,r]]
            M[r]=M[r]*pow(int(M[r,c]),-1,p)%p
            affected=np.flatnonzero(M[r+1:,c])+r+1
            M[affected]=(M[affected]-M[affected,c,None]*M[r,None,:])%p
            r+=1
            if r==len(M):break
        return r
    t=np.concatenate([np.array(a[k],dtype=np.int64) for k in ['homogenizer_quadratic_mod_prime','homogenizer_cubic_mod_prime','homogenizer_quartic_mod_prime']])
    rk=rank(M);rkt=rank(np.vstack([M,t]))
    assert rk==430 and rkt==431 and 1+n+len(mon)-rk==65
    return dict(extra_zero=a['extra_zero'],specified_point_rank_quad_radical=stats(r0),same_beta_regular_rank_quad_radical=stats(r1),independent_saved_relation_rank=rk,independent_augmented_t_rank=rkt,seconds=time.time()-start)

def review():
    path=ROOT/'reports/S61_BINARY_ENDPOINTS_20261101.json';a=json.loads(path.read_text());assert a['all_passed']
    cases=[review_case(c) for c in a['cases']]
    paths=[Path(__file__),path,ROOT/'code/s61_binary_endpoint_engine_20261101.py',ROOT/'code/s61_binary_endpoints_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'research/notes/S61_BINARY_ENDPOINTS_20261101.md']
    return dict(id='S61_BINARY_ENDPOINTS_REVIEW_20261101',all_passed=True,cases=cases,proof_audit=['Both specified transverse algebras have length65 by quintic finiteness, a430 relation minor, and SAME-beta generic length lower65.','Relative binary support over the eleven-dimensional beta slice is smooth of dimension14; its special transverse algebra is the same length65 algebra.','The relative completed ring is free rank65 by the same generic lower65, hence t^2=0.','The Sp4-times-slice beta morphism is smooth at the nonzero16-rank minor; full equivariance and faithful flatness descend nilpotence to the full ambient homogenization.','Connected projective stabilizer preserves each11D slice, so each image is16D; the remaining bad subset has dimension<=15.'],accepted_conclusion='Both u14=0 andu49=0 endpoint strata are removed over sufficiently general beta; specified completed rings are free rank65, quotient(t) free64.',limitations=['Quartic coefficients were read from the main certificate; relation algebra and QQ local data were reconstructed independently.','The u4=0 binary/R0 junction is not covered; no full binary or affine exhaustion.'],input_sha256={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=review()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,cases=r['cases'])))
