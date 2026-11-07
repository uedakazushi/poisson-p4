#!/usr/bin/env python3
"""Full rank-zero source on a smooth intersection of two quadrics.

All20 constant-skew coefficients are retained. Coordinate sign characters
split the210 parameter quadrics; every actual Pfaffian lies in one of five
12-column weight-four blocks. First mode is an exact rational diagnostic.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import hashlib,json,resource,signal,sys,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));QUADS=list(combinations(range(5),4));PQUADS=list(combinations_with_replacement(range(20),2))


def monomials(n,d):
    out=[]
    for terms in combinations_with_replacement(range(n),d):
        e=[0]*n
        for j in terms:e[j]+=1
        out.append(tuple(e))
    return sorted(out,reverse=True)


def build(lam):
    labels=[(quad,e) for quad in QUADS for e in monomials(5,4)];idx={v:i for i,v in enumerate(labels)};cols={v:i for i,v in enumerate(PQUADS)};J=[{} for row in labels]
    def pair(i,j):return None if i==j else (PAIRS.index(tuple(sorted((i,j)))),1 if i<j else -1)
    for quad in QUADS:
        for pos,missing in enumerate(quad):
            tri=quad[:pos]+quad[pos+1:]
            for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]:
                for l in range(5):
                    il,jk=pair(a,l),pair(b,c)
                    if il is None:continue
                    for first in range(2):
                        for second in range(2):
                            col=cols[tuple(sorted((10*first+il[0],10*second+jk[0])))];base=2*(-1)**pos*il[1]*jk[1]
                            for t in range(5):
                                e=[0]*5;e[missing]+=1;e[l]+=1;e[t]+=2;row=J[idx[quad,tuple(e)]]
                                value=base*(lam[t] if first else 1)*(lam[l] if second else 1)
                                row[col]=row.get(col,0)+value
    return labels,[{j:s.expand(v) for j,v in row.items() if s.expand(v)!=0} for row in J]


def targets(quad):
    a,b,c,d=quad;out=[[0]*len(PQUADS) for k in range(3)]
    for ij,kl,sign in [((a,b),(c,d),1),((a,c),(b,d),-1),((a,d),(b,c),1)]:
        p,q=PAIRS.index(ij),PAIRS.index(kl)
        for u,v,slot in [(p,q,0),(p,10+q,1),(10+p,q,1),(10+p,10+q,2)]:out[slot][PQUADS.index(tuple(sorted((u,v))))]+=sign
    return out



def symbolic():
    start=time.monotonic();signal.alarm(90);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    diagnostic_path=ROOT/'reports/CLASSIFICATION_SMOOTH_CI22_ZERO_PARITY_DIAGNOSTIC_20261107.json';diagnostic=json.loads(diagnostic_path.read_text());assert diagnostic['all_five_blocks_passed']
    lam=s.symbols('l0:5');labels,J=build(lam);case=diagnostic['cases'][0];quad=tuple(case['quad']);assert quad==(0,1,2,3)
    columns=case['parameter_monomial_columns'];selected=case['rank_minor']['rows'];pivot_cols=case['rank_minor']['columns']
    small=s.Matrix([[J[i].get(j,0) for j in columns] for i in selected]);U=s.Matrix([[J[i].get(j,0) for j in pivot_cols] for i in selected]);T=s.Matrix(targets(quad))[:,columns]
    det=s.factor(U.det(method='domain-ge'));assert det
    local_cols=[columns.index(j) for j in pivot_cols];L=(T[:,local_cols]*U.inv(method='DM')).applyfunc(s.cancel)
    assert all(s.cancel(v)==0 for v in L*small-T)
    denoms=sorted(set(s.factor(s.denom(v)) for v in L if v),key=str);factors=[];allowed=True
    for denominator in denoms:
        coefficient,parts=s.factor_list(denominator)
        packet={'denominator':str(denominator),'rational_unit':str(coefficient),'difference_factors':[],'unhandled_factors':[]}
        for factor,power in parts:
            match=None
            for i,j in combinations(range(5),2):
                ratio=s.cancel(factor/(lam[i]-lam[j]))
                if ratio.is_Rational:match={'i':i,'j':j,'power':int(power),'rational_unit':str(ratio)};break
            if match:packet['difference_factors'].append(match)
            else:packet['unhandled_factors'].append({'factor':str(factor),'power':int(power)});allowed=False
        factors.append(packet)
    functions=[{str(selected[j]):str(s.factor(L[k,j])) for j in range(len(selected)) if L[k,j]!=0} for k in range(3)]
    out={'task_id':'G','all_symbolic_target_identities_passed':True,'lambda_symbols':list(map(str,lam)),
        'input_diagnostic_sha256':hashlib.sha256(diagnostic_path.read_bytes()).hexdigest(),
        'representative_quad':list(quad),'parameter_monomial_columns':columns,'selected_original350_rows':selected,'selected_pivot_columns':pivot_cols,
        'selected_rank_minor':str(det),'symbolic_six_by_twelve_matrix':[[str(s.expand(v)) for v in row] for row in small.tolist()],
        'three_Pfaffian_target_rows':[[int(v) for v in row] for row in T.tolist()],
        'three_exact_full350_Jacobi_functionals':functions,'three_exact_full350_Jacobi_functional_fractions':[{str(selected[j]):{'numerator':str(s.expand(s.numer(L[k,j]))),'denominator':str(s.expand(s.denom(L[k,j])))} for j in range(len(selected)) if L[k,j]!=0} for k in range(3)],'all_functional_denominators':factors,
        'only_Vandermonde_differences_inverted':allowed,'all_smooth_CI22_rankzero_source_has_no_rank4':allowed,
        'all_five_quad_blocks_by_full_coordinate_permutation_symmetry':True,
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-start,
        'scope':'Exact symbolic full-source PfA/PfAB/PfB identities in one four-index block. Coordinate S5 symmetry carries this same polynomial construction to every other block. If all functional denominators are Vandermonde units, every smooth quadric-pencil eigenvalue tuple is covered and the actual anticanonical quintic is zero on every Poisson source point. Nonzero rank2 tensors remain; singular pencils and rank2 surface restrictions are separate.'}
    dest=ROOT/'reports/CLASSIFICATION_SMOOTH_CI22_ZERO_PARITY_SYMBOLIC_20261107.json'
    if dest.exists():raise FileExistsError(dest)
    dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'all_smooth_CI22_rankzero_no_rank4':allowed,'minor':str(det),'denominators':list(map(str,denoms)),'seconds':out['seconds']}))

def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    lam=list(map(s.Integer,[0,1,2,3,4]));labels,J=build(lam);records=[]
    dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ)
    for quad in QUADS:
        columns=[]
        for k,(i,j) in enumerate(PQUADS):
            parity=set(PAIRS[i%10])^set(PAIRS[j%10])
            if parity==set(quad):columns.append(k)
        assert len(columns)==12
        rows=[i for i,row in enumerate(J) if any(j in row for j in columns)]
        assert all(set(J[i])<=set(columns) for i in rows)
        M=s.Matrix([[J[i].get(j,0) for j in columns] for i in rows]);T=s.Matrix(targets(quad))[:,columns]
        _,piv=dm(M).rref();rank=len(piv);contained=len(dm(M.col_join(T)).rref()[1])==rank;functions=None;minor=None
        if contained:
            selected=list(dm(M.T).rref()[1]);small=M[selected,:];pc=list(dm(small).rref()[1]);inverse=small[:,pc].inv();L=T[:,pc]*inverse
            assert L*small==T;functions=[]
            for k in range(3):functions.append({str(rows[selected[j]]):str(L[k,j]) for j in range(len(selected)) if L[k,j]})
            minor={'rows':[rows[i] for i in selected],'columns':[columns[j] for j in pc],'determinant':str(small[:,pc].det())}
        records.append({'quad':list(quad),'parameter_monomial_columns':columns,'full_Jacobi_rows_in_block':rows,'block_shape':list(M.shape),'rank':rank,
            'all_three_Pfaffian_targets_in_full_Jacobi_rowspace':contained,'full350_Jacobi_functionals_for_three_targets':functions,'rank_minor':minor})
        print('FIXED LAMBDA BLOCK',quad,'rank',rank,'Pf targets',contained,flush=True)
    result={'task_id':'G','lambda_values':list(map(str,lam)),'all20_source_parameters_retained':True,'all350_Jacobi_rows_reconstructed':True,
        'parameter_quadratic_monomial_indices':[list(v) for v in PQUADS],
        'full350_by210_sparse_Jacobi_matrix':{str(i):{str(j):str(v) for j,v in row.items()} for i,row in enumerate(J) if row},
        'rows':[{'quad':list(quad),'exponents':list(e)} for quad,e in labels],'cases':records,
        'all_five_blocks_passed':all(v['all_three_Pfaffian_targets_in_full_Jacobi_rowspace'] for v in records),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-start,
        'scope':'Exact diagnostic at five distinct rational pencil eigenvalues only. No universal smooth-CI22 exclusion until symbolic parameter identities and every Vandermonde complement are handled. This is the full20 source/all350 equations, not a point in that source.'}
    dest=ROOT/'reports/CLASSIFICATION_SMOOTH_CI22_ZERO_PARITY_DIAGNOSTIC_20261107.json'
    if dest.exists():raise FileExistsError(dest)
    dest.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    if '--symbolic' in sys.argv:symbolic()
    else:main()
