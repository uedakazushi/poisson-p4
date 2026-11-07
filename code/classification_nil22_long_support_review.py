#!/usr/bin/env python3
"""Independent Fraction review of whole nil22/top2 source dimension <= 9.

No producer, CAS, Singular, NumPy or SymPy import. The complete old ninety
coefficient source is a prior independently checked input; all original Jacobi,
core generators, exact membership lifts, leading covers, uniform fibers and
cofactor orbit are reconstructed here.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,permutations
from pathlib import Path
import hashlib,json,math,resource,sys,time
from classification_quadratic_cofactor_review import Ring,PAIRS
import classification_quadratic_cofactor_review as helper
import classification_two_lines_plane_cubic_component_review as linear


def main():
    start=time.monotonic();sys.setrecursionlimit(100000)
    resource.setrlimit(resource.RLIMIT_CPU,(110,111));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    root=Path(__file__).resolve().parents[1];R=Ring(20);zero=F(0)
    paths={key:root/'reports'/('CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_'+suffix+'_20261107.json') for key,suffix in
        [('core','CORE'),('critical','CRITICAL_CORES'),('low','LOW_RANK_CORE'),('support','SUPPORT')]}
    source_path=root/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    source=json.loads(source_path.read_text())['full_nilpotent_sources'][6]
    payload={key:json.loads(path.read_text()) for key,path in paths.items()}
    support=payload['support']['nil22_long_support']
    core=[2,8,11,12,14,17,18,19];avars=[1,7,10,13,16];bvars=[0,3,4,5,6,9,15]
    variable=lambda i:R.variable('n'+str(i))
    def parse(t):return R.parse(t)
    def scale(a,c):return R.scale(a,F(c))
    def subtract(a,b):return R.add(a,scale(b,-1))
    def determinant(A):
        out={};size=len(A)
        for perm in permutations(range(size)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(size) for j in range(i+1,size))
            out=R.add(out,scale(R.mul(*(A[i][perm[i]] for i in range(size))),sign))
        return out
    def substitute(poly,mapping):
        out={}
        for e,c in poly.items():
            term=R.constant(c)
            for i,degree in enumerate(e):
                for k in range(degree):term=R.mul(term,mapping.get(i,R.variable(R.names[i])))
            out=R.add(out,term)
        return out
    def primitive(poly):
        den=1
        for v in poly.values():den=math.lcm(den,v.denominator)
        nums=[int(v*den) for v in poly.values()];g=0
        for v in nums:g=math.gcd(g,v)
        out=scale(poly,F(den,g));leading=max(out)
        return out if out[leading]>0 else scale(out,-1)
    def sparse_id(poly):return tuple(sorted(poly.items()))
    def linear_system(poly,indexes):
        row=[{} for i in indexes];constant={}
        for e,v in poly.items():
            degree=sum(e[i] for i in indexes);assert degree<=1
            if not degree:constant[e]=v;continue
            pos=next(i for i in indexes if e[i]);f=list(e);f[pos]-=1
            row[indexes.index(pos)][tuple(f)]=v
        return row,scale(constant,-1)
    # Full original twenty-variable source; no core equation is taken on faith.
    E=list(map(parse,source['brackets']));P=[[{} for j in range(4)] for i in range(4)]
    for (i,j),v in zip(PAIRS,E):P[i][j]=v;P[j][i]=scale(v,-1)
    Jac=[R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in combinations(range(4),3)]
    assert Jac==list(map(parse,source['all_four_jacobi']))
    actual_coeff=[]
    for poly in Jac:
        grouped={}
        for e,v in poly.items():
            tail=e[20:];head=e[:20]+(0,)*4;grouped.setdefault(tail,{})[head]=v
        actual_coeff += [v for _,v in sorted(grouped.items(),reverse=True) if v]
    J=list(map(parse,source['distinct_nonzero_jacobi_coefficients']))
    assert all(any(v==p or v==scale(p,-1) for p in J) for v in actual_coeff)
    assert all(any(p==v or p==scale(v,-1) for v in actual_coeff) for p in J)
    core_generators=[]
    for poly in J:
        p=primitive(poly)
        if any(e[i] for e in p for i in range(20) if i not in core):continue
        if p not in core_generators:core_generators.append(p)
    assert len(core_generators)==7
    # Exact matrices from the actual coefficient equations, including all RHS.
    A=[]
    for i in [2,8]:
        row,rhs=linear_system(J[i],avars);assert not rhs;A.append(row)
    B=[];rhs=[]
    for i in [13,23,9,14,15]:
        row,right=linear_system(J[i],bvars);B.append(row);rhs.append(right)
    assert A==[list(map(parse,row)) for row in support['A_matrix']]
    assert B==[list(map(parse,row)) for row in support['B_matrix']]
    assert rhs==list(map(parse,support['B_right_hand_side']))
    assert R.add(A[1][1],scale(A[0][4],6),scale(A[0][3],-2),scale(A[1][2],-2))==R.constant(-4)
    C=[[B[i][j] for j in [0,4,5,6]] for i in [1,2]]
    assert C==[list(map(parse,row)) for row in support['C_matrix']]
    assert not any(B[i][j] for i in [1,2] for j in [3,2,1])
    forcing_indices=[17,18,11,14,19,12];forcing=[];forcing_rhs=[]
    for p in [v for row in C for v in row]:
        a,b=linear_system(p,forcing_indices);forcing.append(a);forcing_rhs.append(b)
    forced={int(k[1:]):parse(v) for k,v in support['C_zero_forced_substitution'].items()}
    assert set(forced)==set(forcing_indices)
    assert all(not substitute(p,forced) for row in C for p in row)
    def constant(poly):assert set(poly)<={R.zero};return poly.get(R.zero,zero)
    rows=support['C_zero_forcing_minor']['rows']
    forced_minor=linear.determinant([[constant(forcing[i][j]) for j in range(6)] for i in rows])
    assert forced_minor==F(support['C_zero_forcing_minor']['value']) and forced_minor
    assert substitute(J[25],forced)==R.constant(3)
    regular=determinant([[B[i][j] for j in [3,2,1]] for i in [0,3,4]])
    assert regular==parse(support['B_regular_three_pivot_minor'])==parse('-(n14+1)*(2*n14+1)*(3*n14+1)')
    for packet in support['B_critical_two_pivot_minors']:
        det=determinant([[B[i][j] for j in packet['columns']] for i in packet['rows']])
        assert substitute(det,{14:R.constant(F(packet['n14']))})==R.constant(F(packet['value']))
    print('FULL20 JACOBI / A2x5 / B5x7 / C CONTRADICTION PASS',flush=True)
    # Five independently checked membership/leading certificates. No GB is run.
    packets=[payload['core']['nil22_long_core'],*payload['critical']['nil22_long_critical_cores'],payload['low']['nil22_long_low_rank_core']]
    certificates=[]
    for packet in packets:
        begin=time.monotonic();gens=list(core_generators)
        fixed=packet.get('fixed_n14_value')
        if fixed is not None:gens.append(subtract(variable(14),R.constant(F(fixed))))
        if packet.get('imposed_A_rank_at_most_one'):
            for j,k in combinations(range(5),2):
                minor=subtract(R.mul(A[0][j],A[1][k]),R.mul(A[0][k],A[1][j]))
                if minor:gens.append(minor)
        assert gens==list(map(parse,packet['primitive_full_Jacobi_generators']))
        G=list(map(parse,packet['lifted_polynomials']));sums=[{} for p in G]
        for lift in packet['exact_ideal_membership_lifts']:
            i,j=lift['generator'],lift['polynomial'];assert 0<=i<len(gens) and 0<=j<len(G)
            sums[j]=R.add(sums[j],R.mul(gens[i],parse(lift['coefficient'])))
        assert G==sums
        leading=[]
        for poly in G:
            assert poly and all(not any(e[20:]) for e in poly)
            e=max(poly,key=lambda e:(sum(e[:20]),tuple(-a for a in reversed(e[:20]))));leading.append(e[:20])
        assert leading==[tuple(e) for e in packet['leading_monomial_exponents']]
        edges={frozenset(i for i,a in enumerate(e) if a) for e in leading};assert all(edge and edge<=set(core) for edge in edges)
        cover=None;tried=0
        for size in range(9):
            for subset in combinations(core,size):
                tried+=1
                if all(set(subset)&edge for edge in edges):cover=subset;break
            if cover is not None:break
        height=len(cover);bound=8-height
        assert height==packet['leading_monomial_minimum_vertex_cover'] and bound==packet['core_dimension_upper']
        assert 20-height==packet['fixed_source_dimension_upper']
        certificates.append({'fixed_n14':fixed,'rank_A_at_most_one':bool(packet.get('imposed_A_rank_at_most_one')),
            'all_original_generators_reconstructed':True,'all_membership_lifts_passed':True,'lift_count':len(packet['exact_ideal_membership_lifts']),
            'all_grevlex_leading_monomials_reconstructed':True,'leading_polynomial_count':len(G),
            'minimum_vertex_cover':list(cover),'minimum_cover_size':height,'all_smaller_covers_exhausted':True,'cover_subsets_checked':tried,
            'core_dimension_upper':bound,'seconds':time.monotonic()-begin})
        print('CORE LIFTS/LEADING PASS',fixed,'rankA1',bool(packet.get('imposed_A_rank_at_most_one')),'dimension',bound,flush=True)
    assert [(p['fixed_n14'],p['rank_A_at_most_one'],p['core_dimension_upper']) for p in certificates]==[(None,False,2),('-1',False,1),('-1/2',False,2),('-1/3',False,1),('-1/2',True,1)]
    # Original affine cofactor orbit: construct every one of twenty columns.
    NA=[[F(v) for v in row] for row in source['field_linear_matrix']];nc=list(map(F,source['field_constant']));columns=[]
    for index in range(20):
        M=[[zero for j in range(4)] for i in range(4)];t=[zero]*4
        if index<16:M[index//4][index%4]=F(1)
        else:t[index-16]=F(1)
        change=[[sum(M[i][k]*NA[k][j]-NA[i][k]*M[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
        change_t=[sum(M[i][k]*nc[k]-NA[i][k]*t[k] for k in range(4)) for i in range(4)]
        columns.append([v for row in change for v in row]+change_t)
    orbit=linear.transpose(columns)
    assert orbit==[[F(v) for v in row] for row in support['exact_affine_orbit_tangent_matrix']]
    assert linear.rank(orbit)[0]==12 and linear.rank([row[:16] for row in orbit[:16]])[0]==8
    for p in support['complete_dimension_strata']:
        assert p['first_fiber_dimension_upper']==5-p['rank_A_at_least']
        assert p['second_fiber_dimension_upper']==7-p['rank_B_at_least']
        assert p['total_source_dimension_upper']==p['core_dimension']+5-p['rank_A_at_least']+7-p['rank_B_at_least']==9
    assert support['fixed_source_dimension_upper']==9 and support['saturated_projective_dimension_upper']==9+12+4-1==24
    inputs=[source_path,*paths.values(),Path(helper.__file__),Path(linear.__file__)]
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'no_producer_or_CAS_import':True,
        'all_original20_parameter_source_Jacobi_reconstructed':True,'all_seven_core_generators_reconstructed':True,
        'full_A_2_by_5_and_B_5_by_7_and_all_rhs_reconstructed':True,'A_rank_one_uniform_constant_identity':True,
        'C_zero_forcing_nonzero_minor':str(forced_minor),'C_zero_actual_core_contradiction':'3',
        'all_regular_and_critical_B_pivots_reconstructed':True,'all_five_core_certificates':certificates,
        'affine_cofactor_orbit_dimension':12,'linear_cofactor_orbit_dimension':8,'full_affine_orbit_matrix_reconstructed':True,
        'all_boundary_dimension_strata_retained':True,'fixed_full20_source_dimension_upper':9,'saturated_projective_dimension_upper':24,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reviewer_source_bytes':Path(__file__).read_text(),
        'seconds':time.monotonic()-start,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':'All source/Jacobi/core/lift/leading-monomial-cover/uniform-fiber/orbit arithmetic independently reconstructed by standard-library Fraction. Only ideal membership of displayed lifted polynomials is required, not a fresh GB or completeness of that GB. The whole-source fiber dimension proof and projective saturation argument are separately reviewed in the note. Entire reducible/zero-Pfaffian boundary retained; no E component membership inferred.'}
    (root/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds'],'maximum_RSS_KiB':out['maximum_RSS_KiB']}))


if __name__=='__main__':main()
