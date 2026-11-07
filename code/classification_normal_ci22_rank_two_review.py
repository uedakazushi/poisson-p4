#!/usr/bin/env python3
"""Independent normal-CI22 tangent/normalizer and original350 factor proof.

Fraction/Laurent arithmetic only; no producer, SymPy or legacy imports.
"""
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
import ast,hashlib,json,resource,signal,time
import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_minimal_quartic_review import LaurentRing

ROOT=Path(__file__).resolve().parents[1]
PREF='reports/CLASSIFICATION_NORMAL_CI22_RANK_TWO_'
PARAMS=('u0','u1','u2','u3','u4','a','b','c','d')
COORDS=tuple('x'+str(i)for i in range(5))
PAIRS=list(combinations(range(5),2));QUADS=list(combinations(range(5),4))


def main():
    started=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    out=ROOT/(PREF+'REVIEW_20261107.json')
    if out.exists():raise FileExistsError(out)
    job=json.loads((ROOT/'research/classification/jobs/normal_ci22_rank_two_review.json').read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256']
    sources=json.loads((ROOT/(PREF+'SOURCE_20261107.json')).read_text())['cases']
    proofs={r['model']:r for r in json.loads((ROOT/(PREF+'FACTOR_PROOF_20261107.json')).read_text())['cases']}
    R=LaurentRing(PARAMS,COORDS);p=R.parse;x=[R.variable(n)for n in COORDS];z=(0,)*9
    add=lambda *Ps:{ij:R.add(*(P.get(ij,{})for P in Ps))for ij in PAIRS}
    scale=lambda P,c:{ij:R.scale(P.get(ij,{}),c)for ij in PAIRS}
    def value(P,i,j):return P.get((i,j),{})if i<j else R.scale(P.get((j,i),{}),-1)if i>j else{}
    def jacobi(P):
        J={ijk:R.add(*(R.mul(value(P,a,l),R.derivative(value(P,b,c),l))for a,b,c in (ijk,ijk[1:]+ijk[:1],ijk[2:]+ijk[:2])for l in range(5)))for ijk in combinations(range(5),3)}
        return {I:R.add(*(R.scale(R.mul(x[I[k]],J[I[:k]+I[k+1:]]),(-1)**k)for k in range(4)))for I in QUADS}
    def pf(P):
        out=[]
        for missing in range(5):
            a,b,c,d=[j for j in range(5)if j!=missing]
            out.append(R.add(R.mul(value(P,a,b),value(P,c,d)),R.scale(R.mul(value(P,a,c),value(P,b,d)),-1),R.mul(value(P,a,d),value(P,b,c))))
        return out
    def rows(poly,mons):return [{e[:9]+(0,)*5:c for e,c in poly.items()if e[9:]==m}for m in mons]
    results=[]
    for row in sources:
        name=row['model'];proof=proofs[name];f,g=p(row['q1']),p(row['q2'])
        JU=[]
        for missing in range(5):
            I=[j for j in range(5)if j!=missing];T={ij:{}for ij in PAIRS}
            for a,i in enumerate(I):
                rem=I[:a]+I[a+1:]
                for b,j in enumerate(rem):
                    ij=tuple(rem[:b]+rem[b+1:]);T[ij]=R.add(T[ij],R.scale(R.mul(R.derivative(f,i),R.derivative(g,j)),(-1)**(a+b)))
            JU.append(T)
        gens=[[F(v)for v in r]for r in row['full_two_constant_skew_generators']]
        explicit={'sem221':[(0,1,1),(2,3,1)],'mixed21_11':[(0,2,1),(3,4,1)],'nil41':[(0,4,1),None]}[name]
        expected=[[F(0)]*10 for _ in range(2)]
        for k,e in enumerate(explicit):
            if e is not None:i,j,v=e;expected[k][PAIRS.index((i,j))]=F(v)
            else:expected[k][PAIRS.index((0,3))]=F(3,4);expected[k][PAIRS.index((1,2))]=F(1,4)
        assert gens==expected
        # Rebuild the full constant-skew conormal equation, including both
        # unknown coefficients ofq1,q2; full kernel, not a chosen submatrix.
        polys=[R.add(R.mul(R.derivative(f,i),R.derivative(g,j)),R.scale(R.mul(R.derivative(f,j),R.derivative(g,i)),-1))for i,j in PAIRS]+[R.scale(f,-1),R.scale(g,-1)]
        M=[[P.get(z+e,F(0))for P in polys]for e in q.monomials(2)]
        assert M==[[F(v)for v in r]for r in row['full_conormal_matrix']]
        K=[[F(v)for v in r]for r in row['conormal_augmented_nullspace']]
        assert q.rank(M)[0]==10 and q.rank(K)[0]==2
        assert all(not sum(a*b for a,b in zip(r,k))for r in M for k in K)
        assert q.rank([k[:10]for k in K]+gens)[0]==2
        basis=JU+[{ij:R.scale(f,v)for ij,v in zip(PAIRS,gens[0])},{ij:R.scale(g,v)for ij,v in zip(PAIRS,gens[0])},{ij:R.scale(f,v)for ij,v in zip(PAIRS,gens[1])},{ij:R.scale(g,v)for ij,v in zip(PAIRS,gens[1])}]
        for B,saved in zip(basis,row['all150_source_tensors']):
            for key,data in saved.items():
                ij=ast.literal_eval(key);assert B.get(ij,{})=={z+tuple(e):F(c)for e,c in data}
        P={ij:R.add(*(R.mul(p(t),B.get(ij,{}))for t,B in zip(PARAMS,basis)))for ij in PAIRS}
        pure={ij:R.add(*(R.mul(p(t),B.get(ij,{}))for t,B in zip(PARAMS[:5],JU)))for ij in PAIRS}
        assert not any(jacobi(pure).values())and not any(pf(pure))
        for B in JU:
            for h in (f,g):
                assert all(not R.add(*(R.mul(value(B,i,j),R.derivative(h,j))for j in range(5)))for i in range(5))
        monparams=[p(PARAMS[i]+'*'+PARAMS[j])for i,j in row['full_parameter_pair_order']]
        decode_row=lambda values:R.add(*(R.scale(monparams[int(i)],F(v))for i,v in values.items()))
        J=jacobi(P);Jrows=[v for I in QUADS for v in rows(J[I],q.monomials(4))]
        assert Jrows==[decode_row(t)for t in row['full350_Jacobi_rows']]
        Fraw=R.add(*(R.scale(R.mul(x[i],t),(-1)**i)for i,t in enumerate(pf(P))))
        assert rows(Fraw,q.monomials(5))==[decode_row(t)for t in row['full126_actualPf_rows']]
        for cert in proof['required_equation_functionals']:
            total=R.add(*(R.scale(Jrows[int(i)],F(c))for i,c in cert['full_original350_functional'].items()))
            assert total==p(cert['required_equation'])
        factor=None
        if name!='nil41':
            keep=4 if name=='sem221'else 1;factor=4 if name=='sem221'else 1
            def project(poly):
                out={}
                for ex,c in poly.items():
                    if any(ex[i]for i in range(5)if i!=keep)or ex[7]:continue
                    e=list(ex);c*=(-1)**e[6];e[5]+=e[6];e[6]=0;e=tuple(e);out[e]=out.get(e,F(0))+c
                return{e:c for e,c in out.items()if c}
            reduced={ij:project(t)for ij,t in P.items()}
            assert not any(jacobi(reduced).values())
            Dred=project(Fraw);assert Dred and all(e[9+factor]>=1 for e in Dred)
            assert all(not value(reduced,factor,j)for j in range(5))
        results.append(dict(model=name,full25_tangent_source_and2_normalizer=True,full150_all350_all126_reconstructed=True,all_required_equations_original350_functionals=True,actual_fixed_linear_factor=None if factor is None else COORDS[factor],nil_reduced_coefficient_determinant_forcedzero=name=='nil41',rankone_coefficient_case='Common quadratic factor before any Jacobi restriction',whole_normal_branch_empty=True))
    result=dict(task_id='G',all_passed=True,no_producer_legacy_or_CAS_imports=True,cases=results,global_Euler_and_normal_pencil_classification='Separately independently read in review note; allregular remainingtypes andall-singular smoothellipticcone retained.',no_non_normal_CI22_claim=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    out.write_text(json.dumps(result,indent=2)+'\n');print('PASS',round(result['seconds'],3),'seconds',flush=True)


if __name__=='__main__':main()
