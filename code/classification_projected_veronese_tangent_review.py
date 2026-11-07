#!/usr/bin/env python3
"""Independent full tangent source and lifted dimension bound for projected PV.

Only the previously independent standard-library Fraction polynomial/matrix
engine is reused. No producer, CAS, Singular, SymPy or NumPy is imported.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement,permutations,product
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,math,resource,sys,time
import classification_projected_veronese_zero_review as r


def main():
    start=time.monotonic();sys.setrecursionlimit(100000)
    resource.setrlimit(resource.RLIMIT_CPU,(235,236));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    root=Path(__file__).resolve().parents[1]
    paths=[root/'reports'/('CLASSIFICATION_PROJECTED_VERONESE_TANGENT'+suffix+'_20261107.json') for suffix in ['', '_DIMENSION']]
    saved,dimdata=[json.loads(p.read_text()) for p in paths]
    assert hashlib.sha256(saved['producer_bytes'].encode()).hexdigest()==saved['producer_sha256']
    assert hashlib.sha256(dimdata['producer_bytes'].encode()).hexdigest()==dimdata['producer_sha256']
    assert dimdata['input_sha256']==hashlib.sha256(paths[0].read_bytes()).hexdigest()
    X=r.Ring(['x'+str(i) for i in range(5)]);S=r.Ring(['a','b','c']);Z=r.Ring(['z'+str(i) for i in range(15)]);ZX=r.Ring(Z.names+X.names)
    q=list(map(S.parse,['a*a-b*b','b*b-c*c','a*b','a*c','b*c']))
    assert q==list(map(S.parse,saved['surface_parametrization']))
    cols=[(ij,e) for ij in r.PAIRS for e in r.monomials(5,2)]
    assert cols==[(tuple(p['pair']),tuple(p['exponents'])) for p in saved['ambient_coefficient_columns']]
    labels=[(tri,e) for tri in r.TRIPLES for e in r.monomials(3,6)];row_index={v:i for i,v in enumerate(labels)};M=[{} for row in labels]
    for j,(ij,e) in enumerate(cols):
        value=S.const(1)
        for a,power in zip(q,e):value=r.times(value,S.power(a,power))
        for tri in r.TRIPLES:
            if set(ij).issubset(tri):
                k=next(k for k in tri if k not in ij)
                for exponent,c in r.times(q[k],value).items():M[row_index[tri,exponent]][j]=(-1)**tri.index(k)*c
    tangent={}
    for tri in r.TRIPLES:
        value={}
        for perm in permutations(range(3)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
            value=r.plus(value,r.scale(r.times(*(r.derivative(q[tri[i]],perm[i]) for i in range(3))),sign))
        tangent[tri]=value
    assert tangent=={tri:S.parse(saved['tangent_three_vector'][str(tri)]) for tri in r.TRIPLES}
    mons3=r.monomials(3,3);assert mons3==[tuple(e) for e in saved['restricted_cubic_monomial_exponents']]
    for j,e in enumerate(mons3):
        mon={e:F(1)}
        for tri in r.TRIPLES:
            for exponent,c in r.times(mon,tangent[tri]).items():M[row_index[tri,exponent]][150+j]=-c
    assert M==[{int(j):F(c) for j,c in saved['full280_by160_sparse_matrix'].get(str(i),{}).items()} for i in range(280)]
    red,_,_=r.sparse_rref(M);assert len(red)==saved['exact_matrix_rank']==121
    free=[i for i in range(160) if i not in red];K=[]
    for j in free:
        v=[F(0)]*160;v[j]=F(1)
        for p,row in red.items():v[p]=-row.get(j,F(0))
        K.append(v)
    assert len(K)==39 and len(r.sparse_rref([v[:150] for v in K])[0])==39
    assert all(not sum(v[j]*c for j,c in row.items()) for v in K for row in M)
    original=[[F(c) for c in v] for v in saved['full_kernel_basis_columns']]
    assert len(original)==39 and len(r.sparse_rref(original)[0])==39 and all(not sum(v[j]*c for j,c in row.items()) for v in original for row in M)
    cert=saved['rank_certificate'];det=r.det([[M[i].get(j,F(0)) for j in cert['columns']] for i in cert['rows']])
    assert det==F(cert['determinant']) and det and len(cert['rows'])==121
    def decode(v,ring=X):
        P={ij:{} for ij in r.PAIRS};offset=len(ring.names)-5
        for c,(ij,e) in zip(v,cols):
            if c:P[ij][(0,)*offset+e]=c
        for i,j in r.PAIRS:P[j,i]=r.scale(P[i,j],-1)
        return P
    def vector(P):return [P.get(ij,{}).get(e,F(0)) for ij,e in cols]
    Euler=[]
    for a,b in product(range(5),repeat=2):
        P={}
        for i,j in r.PAIRS:P[i,j]=r.plus(r.times(X.var('x'+str(i)),X.var('x'+str(b))) if j==a else {},r.scale(r.times(X.var('x'+str(j)),X.var('x'+str(b))),-1) if i==a else {})
        Euler.append(vector(P))
    assert len(r.sparse_rref(Euler)[0])==24
    assert all(not sum(v[j]*c for j,c in row.items() if j<150) for v in Euler for row in M)
    canonical=[]
    for v in K:
        P=decode(v[:150]);div=[r.plus(*(r.derivative(P.get((i,j),{}),i) for i in range(5))) for j in range(5)]
        C={ij:r.plus(P[ij],r.scale(r.plus(r.times(X.var('x'+str(i)),div[j]),r.scale(r.times(X.var('x'+str(j)),div[i]),-1)),-F(1,5))) for ij in r.PAIRS for i,j in [ij]}
        canonical.append(vector(C))
    basis,_,_=r.sparse_rref(canonical);B=[[row.get(j,F(0)) for j in range(150)] for p,row in sorted(basis.items())]
    assert len(B)==15 and len(r.sparse_rref(Euler+B)[0])==39
    for v,packet in zip(B,saved['canonical_basis']):
        P=decode(v);assert all(P[ij]==X.parse(packet.get(str(ij[0])+','+str(ij[1]),'0')) for ij in r.PAIRS)
    P={ij:r.plus(*(r.times(ZX.var('z'+str(k)),decode(v,ZX)[ij]) for k,v in enumerate(B))) for ij in r.PAIRS}
    for i,j in r.PAIRS:P[j,i]=r.scale(P[i,j],-1)
    assert all(P[ij]==ZX.parse(saved['whole_canonical_source'].get(str(ij[0])+','+str(ij[1]),'0')) for ij in r.PAIRS)
    print('FULL280x160 / KERNEL39 / EULER24 / QUOTIENT15 PASS',flush=True)
    strict={tri:r.plus(*(r.times(P.get((a,l),{}),r.derivative(P.get((b,c),{}),15+l)) for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])] for l in range(5))) for tri in r.TRIPLES}
    jac={quad:r.plus(*(r.scale(r.times(ZX.var('x'+str(quad[k])),strict[quad[:k]+quad[k+1:]]),(-1)**k) for k in range(4))) for quad in r.QUADS}
    J=[]
    for quad in r.QUADS:
        for e in r.monomials(5,4):J.append({p[:15]:c for p,c in jac[quad].items() if p[15:]==e})
    assert J==list(map(Z.parse,saved['all350_projective_Jacobi']))
    mons2=r.monomials(15,2);JC=[[poly.get(e,F(0)) for e in mons2] for poly in J]
    assert JC==[[F(c) for c in row] for row in saved['full350_by_parameter_quadric_matrix']]
    assert len(r.sparse_rref(JC)[0])==54
    signed=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5) if a!=missing]
        signed.append(r.scale(r.plus(r.times(P[i,j],P[k,l]),r.scale(r.times(P[i,k],P[j,l]),-1),r.times(P[i,l],P[j,k])),(-1)**missing))
    D=r.plus(*(r.times(ZX.var('x'+str(i)),v) for i,v in enumerate(signed)))
    assert D==ZX.parse(saved['actual_anticanonical_quintic'])
    def primitive(poly):
        den=1;g=0
        for v in poly.values():den=math.lcm(den,v.denominator)
        for v in poly.values():g=math.gcd(g,int(v*den))
        return r.scale(poly,F(den,g))
    selected=dimdata['full54_independent_Jacobi_row_indices'];assert len(selected)==54 and len(r.sparse_rref([JC[i] for i in selected])[0])==54
    gens=[primitive(J[i]) for i in selected];assert gens==list(map(Z.parse,dimdata['original_ideal_generators']))
    print('ALL350 JACOBI / FULL54 GENERATORS / ACTUAL QUINTIC PASS',flush=True)
    G=list(map(Z.parse,dimdata['lifted_polynomials']));sums=[{} for g in G]
    for count,lift in enumerate(dimdata['exact_ideal_membership_lifts']):
        i,j=lift['generator'],lift['polynomial'];assert 0<=i<54 and 0<=j<len(G)
        sums[j]=r.plus(sums[j],r.times(gens[i],Z.parse(lift['coefficient'])))
        if count and count%2500==0:print('EXACT LIFTS',count,'/',len(dimdata['exact_ideal_membership_lifts']),flush=True)
    assert G==sums
    leading=[max(poly,key=lambda e:(sum(e),tuple(-a for a in reversed(e)))) for poly in G]
    assert leading==[tuple(e) for e in dimdata['leading_monomial_exponents']]
    edges=tuple(sorted(set(sum(1<<i for i,a in enumerate(e) if a) for e in leading)))
    @lru_cache(None)
    def cover(edges):
        if not edges:return 0
        edge=min(edges,key=int.bit_count);assert edge
        return 1+min(cover(tuple(e for e in edges if not e&(1<<i))) for i in range(15) if edge&(1<<i))
    height=cover(edges);bound=15-height
    assert height==dimdata['leading_support_minimum_vertex_cover']==13 and bound==dimdata['fixed_Poisson_source_dimension_upper']==2
    assert dimdata['whole_projective_PGL_image_dimension_upper']==24+bound-1-3==22
    inputs=paths+[Path(r.__file__)]
    result={'task_id':'G','all_passed':True,'no_producer_or_CAS_or_Singular_import':True,
        'full280_by160_tangency_matrix_reconstructed':True,'matrix_rank':121,'maximal_minor_determinant':str(det),
        'complete_kernel_dimension':39,'full_Euler_dimension':24,'complete_canonical_tangent_source_dimension':15,
        'all350_projective_Jacobi_reconstructed':True,'all54_generators_span_full_Jacobi_rowspace':True,'actual_anticanonical_quintic_reconstructed':True,
        'every_exact_membership_lift_reconstructed':True,'exact_lift_count':len(dimdata['exact_ideal_membership_lifts']),
        'lifted_polynomial_count':len(G),'every_grevlex_leading_monomial_reconstructed':True,'leading_minimum_vertex_cover':height,
        'fixed_whole_Poisson_source_dimension_upper':bound,'whole_projective_PGL_image_dimension_upper':22,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reviewer_source_bytes':Path(__file__).read_text(),
        'seconds':time.monotonic()-start,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':'Complete tangent source, full projective equations and every exact lifted leading certificate independently reconstructed using Fraction. Initial-ideal containment gives dimension<=2 without rerunning GB or asserting its completeness. The preserving-group3 and whole clean rank2 double-surface reduction are separate geometric proofs; no general irreducible-quintic exhaustion or E membership asserted.'}
    dest=root/'reports/CLASSIFICATION_PROJECTED_VERONESE_TANGENT_REVIEW_20261107.json'
    if dest.exists():raise FileExistsError(dest)
    dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_passed':True,'seconds':result['seconds'],'maximum_RSS_KiB':result['maximum_RSS_KiB']}))

if __name__=='__main__':main()
