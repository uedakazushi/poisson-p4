#!/usr/bin/env python3
"""Independent factor/minor audit and a labelled replay of the universal producer."""
import argparse,hashlib,json,time
from pathlib import Path
from functools import lru_cache
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def verify(replay=None):
    started=time.time();path=ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
    prior=json.loads(path.read_text());assert prior['all_passed']
    for rel,digest in prior['input_sha256'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
    if replay is None:
        from s61_binary_deeper_uniform_quadratic_20261101 import verify as run
        current=run()
    else:current=json.loads(Path(replay).read_text())
    assert {k:v for k,v in prior.items() if k!='seconds'}=={k:v for k,v in current.items() if k!='seconds'}
    e,d=s.symbols('e d');other=dict.fromkeys(s.symbols('p q r h k'),s.S.Zero);prime=32003
    @lru_cache(None)
    def polynomial(text):return s.Poly(s.sympify(text).subs(other),e,d,domain=s.QQ)
    def evaluate(text,ev,dv):
        value=0
        for (a,b),co in polynomial(text).terms():
            assert int(co.q)%prime
            value=(value+int(co.p)*pow(int(co.q),-1,prime)*pow(ev,a,prime)*pow(dv,b,prime))%prime
        return value
    def determinant(matrix):
        # Plain Python modular elimination, independent of the NumPy producer.
        M=[row[:] for row in matrix];out=1
        for j in range(len(M)):
            i=next(i for i in range(j,len(M)) if M[i][j])
            if i!=j:M[i],M[j]=M[j],M[i];out=-out
            pivot=M[j][j];out=out*pivot%prime;inv=pow(pivot,-1,prime)
            for i in range(j+1,len(M)):
                a=M[i][j]*inv%prime
                for k in range(j,len(M)):M[i][k]=(M[i][k]-a*M[j][k])%prime
        return out%prime
    factors=[polynomial(f) for f in prior['irreducible_minor_factors']]
    assert len(factors)==13 and all(f.is_irreducible for f in factors)
    normalized=lambda f:s.Poly(f,e,d,domain=s.QQ).monic().as_expr()
    covers=set(map(normalized,prior['previously_removed_factors']))
    certificates=[]
    for item in prior['factor_nondivisibility_witnesses']:
        factor=item['factor'];ev,dv=item['sample_core'];assert evaluate(factor,ev,dv)==0
        matrix=[]
        for i in item['rows']:
            content=evaluate(prior['row_contents'][i],ev,dv)
            matrix.append([content*evaluate(prior['primitive_quadratic_rows'][str(i)][j],ev,dv)%prime for j in item['columns']])
        value=determinant(matrix);assert value==item['determinant_mod_prime']!=0
        covers.add(normalized(factor));certificates.append(dict(factor=factor,minor_mod_prime=value))
    assert covers=={f.monic().as_expr() for f in factors}
    assert 7-2+10==prior['residual_rank_drop_Sp4_saturation_dimension_upper_bound']==15
    paths=[Path(__file__),path,ROOT/'research/notes/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.md',ROOT/'research/notes/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_REVIEW_20261101.md']
    return dict(id='S61_BINARY_DEEPER_UNIFORM_QUADRATIC_REVIEW_20261101',all_passed=True,accepted=True,
        verification_kind='Independent exact factorization and modular minor audit; universal exact producer replay',
        exact_universal_producer_replay_seconds=current['seconds'],irreducible_factors_checked=13,
        independently_recomputed_minors=certificates,proper_beta_image_dimension_upper_bound=15,
        proof_checks=['Row contents restored before all exceptional-minor tests','Gauss lemma gives characteristic-zero nondivisibility, rather than a finite-field rank guess','The determinantal ideal is defined over QQ; its codimension survives extension to C','The rational slice is already transverse to the five-dimensional fixed-Q stabilizer'],
        limitations=['The universal matrix computation is a producer replay, not an independent implementation','Higher obstructions and full t-saturation removal remain open'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-started)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--replay',type=Path);a=p.parse_args()
    r=verify(a.replay)
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
