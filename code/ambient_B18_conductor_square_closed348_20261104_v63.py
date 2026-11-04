#!/usr/bin/env python3
"""Exact closed conductor-pair normalization and sparse 348-dimensional ring.

Only fixed polynomial source identities, fixed S-pairs and normal forms.
No Groebner basis search, source solve, kernel matrix or numeric point is used.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b,c=F(1)):
    out=dict(a)
    for m,x in b.items():
        y=out.get(m,F(0))+c*x
        if y:out[m]=y
        else:out.pop(m,None)
    return out
def times(a,b):
    out={}
    for m,c in a.items():
        for n,d in b.items():
            k=tuple(i+j for i,j in zip(m,n));v=out.get(k,F(0))+c*d
            if v:out[k]=v
            else:out.pop(k,None)
    return out
def mon(z,r,c=1):return {(z,r):F(c)} if c else {}
def encode(p):return [[list(m),str(c)] for m,c in sorted(p.items(),reverse=True)]
def divide(p,basis):
    original=dict(p);quot=[{} for _ in basis];rem={};steps=0
    leading=[max(g) for g in basis]
    while p:
        steps+=1;assert steps<1000
        m=max(p);c=p[m]
        k=next((i for i,n in enumerate(leading) if all(a>=b for a,b in zip(m,n))),None)
        if k is None:rem[m]=c;del p[m];continue
        n=tuple(a-b for a,b in zip(m,leading[k]));factor={n:c/basis[k][leading[k]]}
        quot[k]=plus(quot[k],factor);p=plus(p,times(factor,basis[k]),-1)
    check=dict(rem)
    for q,g in zip(quot,basis):check=plus(check,times(q,g))
    assert check==original
    return quot,rem
def apply(columns,value):
    out={}
    for i,c in value.items():out=plus(out,columns[i],c)
    return out
def action_sum(a,b,c=F(1)):return [plus(x,y,c) for x,y in zip(a,b)]
def records(columns):return [[[i,str(c)] for i,c in sorted(v.items())] for v in columns]


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):
        hashes[str(p.relative_to(ROOT))]=sha(p)
        return json.loads(p.read_text())
    def gate(stem):
        p=ROOT/'reports'/stem/'report.json';r=read(p);assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            q=p.parent/name;assert sha(q)==h,name;hashes[str(q.relative_to(ROOT))]=h
        for rec in r.get('input_hashes',{}).values():
            assert sha(ROOT/rec['path'])==rec['sha256'];hashes[rec['path']]=rec['sha256']
        return r
    gate('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW')
    gate('AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW')
    gate('AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63_REVIEW')
    b=ROOT/'reports'
    source=read(b/'AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json')
    assert source['variables']==[f'zA_{i}' for i in range(7)]+['r_CA']
    assert F(source['A0_determinant'])==F(-1,1024)
    c=[]
    for beta in source['beta2_after_ell']:
        assert not beta or len(beta)==1 and beta[0][0]==[0]*7+[2]
        c.append(F(beta[0][1]) if beta else F(0))
    assert c==list(map(F,['9835/864',0,0,0,'29875/3456','-4565/864',0]))
    assert source['gA_after_ell']==[[[0]*7+[3],'125/24']]
    assert source['gB_after_ell']==[[[0]*7+[3],'-15125/144']]
    assert c[0]
    # A rational constant GL7, not a new formal coordinate guess.
    forward=[[F(0) for _ in range(7)] for _ in range(7)]
    inverse=[[F(0) for _ in range(7)] for _ in range(7)]
    forward[0][0]=1/c[0];inverse[0][0]=c[0]
    for i in range(1,7):forward[i][i]=1;forward[i][0]=-c[i]/c[0];inverse[i][i]=1;inverse[i][0]=c[i]
    assert [sum(forward[i][j]*c[j] for j in range(7)) for i in range(7)]==[F(1)]+[F(0)]*6
    for i,j in itertools.product(range(7),repeat=2):
        assert sum(forward[i][k]*inverse[k][j] for k in range(7))==F(i==j)
        assert sum(inverse[i][k]*forward[k][j] for k in range(7))==F(i==j)
    (out/'preflight.json').write_text(json.dumps(dict(
        source_A0_unit='-1/1024',source_c0=str(c[0]),actual_gA_cubic='125/24',actual_gB_graph_cubic='-15125/144',
        fixed_two_variable_basis_sizes=[5,4,3,2],fixed_S_pairs=20,
        predicted_block_lengths=[26,14,6,2],W6_degree_multiplicities=[1,6,21,56],
        predicted_total348=True,sparse_action_columns=8*348,dense348_matrices_forbidden=True,
        no_source_solve_higher_jet_GB_decomposition_search_or_kernel_matrix=True,
        limits=dict(seconds=60,address_bytes=2_000_000_000,output_bytes=32*1024*1024,expanded_terms=150000)),indent=2)+'\n')

    Z=mon(1,0);R2=mon(0,2);R3=mon(0,3);R6=mon(0,6)
    Az=[times(Z,Z),times(Z,R3),R6]
    Zplus=plus(Z,R2)
    Bz=[times(Zplus,Zplus),times(Zplus,R3),R6]
    src0=[times(a,b) for a in Az for b in Bz]
    src1=[times(times(a,b),m) for a in [Z,R3] for b in [Zplus,R3] for m in [Z,R2]]
    src2=[times(Z,Z),times(Z,R2),times(R2,R2)]
    src3=[Z,R2]
    G0=[plus(plus(mon(4,0),mon(3,2,2)),mon(2,4)),
        plus(mon(3,3),mon(1,7,-1)),plus(mon(2,5),mon(1,7)),mon(1,8),mon(0,10)]
    G1=[plus(mon(3,0),mon(1,4,-1)),plus(mon(2,2),mon(1,4)),mon(1,5),mon(0,7)]
    sources=[src0,src1,src2,src3];bases=[G0,G1,src2,src3]
    # Literal constant source cofactor vectors: no linear solve or basis search.
    cofactors=[
        [{0:1},{1:2,3:-1},{3:1,1:-1},{4:1,2:-1},{6:1,2:1,4:-2}],
        [{0:1,1:-1},{1:1},{3:1},{5:1,3:-1}],
        [{i:1} for i in range(3)],[{i:1} for i in range(2)]]
    certificates=[];standards=[];spair_count=0
    for q,(src,gb,cofs) in enumerate(zip(sources,bases,cofactors)):
        for g,cof in zip(gb,cofs):
            value={}
            for i,c0 in cof.items():value=plus(value,src[i],F(c0))
            assert value==g
        sdiv=[]
        for p in src:
            quot,rem=divide(dict(p),gb);assert not rem
            sdiv.append([encode(v) for v in quot])
        pairs=[]
        for i,j in itertools.combinations(range(len(gb)),2):
            a,b=max(gb[i]),max(gb[j]);lcm=tuple(max(u,v) for u,v in zip(a,b))
            pi=times({tuple(x-y for x,y in zip(lcm,a)):F(1)},gb[i])
            pj=times({tuple(x-y for x,y in zip(lcm,b)):F(1)},gb[j])
            p=plus(pi,pj,-1);quot,rem=divide(dict(p),gb);assert not rem
            pairs.append(dict(pair=[i,j],S=encode(p),quotients=[encode(v) for v in quot],remainder=[]));spair_count+=1
        leading=[max(g) for g in gb]
        assert all(g[max(g)]==1 for g in gb)
        nz=max(z for z,r in leading);nr=max(r for z,r in leading)
        std=[(z,r) for z in range(nz) for r in range(nr)
            if not any(z>=a and r>=b for a,b in leading)]
        assert len(std)==[26,14,6,2][q]
        standards.append(std)
        certificates.append(dict(W_degree=q,source_generators=[encode(g) for g in src],
            fixed_lex_basis=[encode(g) for g in gb],source_constant_cofactors=[[[i,str(c)] for i,c in sorted(v.items())] for v in cofs],
            all_original_source_division_quotients=sdiv,all_S_pair_divisions=pairs,
            leading_monomials=[list(m) for m in leading],standard_monomials=[list(m) for m in std],length=len(std)))
    assert spair_count==20

    def wmons(q):
        return [tuple(comb.count(i) for i in range(6)) for comb in itertools.combinations_with_replacement(range(6),q)]
    normal_basis=[];block_counts=[]
    for q in range(4):
        before=len(normal_basis)
        for w in wmons(q):
            for z,r in standards[q]:normal_basis.append((z,*w,r))
        block_counts.append(len(normal_basis)-before)
    assert block_counts==[26,84,126,112] and len(normal_basis)==348
    index={m:i for i,m in enumerate(normal_basis)};assert normal_basis[0]==(0,)*8
    def normal_monomial(m):
        q=sum(m[1:7])
        if q>=4:return {}
        quot,rem=divide(mon(m[0],m[7]),bases[q])
        return {index[(z,*m[1:7],r)]:c0 for (z,r),c0 in rem.items()}
    canonical=[]
    for k in range(8):
        columns=[]
        for m in normal_basis:
            n=list(m);n[k]+=1;columns.append(normal_monomial(tuple(n)))
        canonical.append(columns)
    for A in canonical:
        for B in canonical:
            assert all(apply(A,B[j])==apply(B,A[j]) for j in range(348))
    for j,m in enumerate(normal_basis):
        value={0:F(1)}
        for k,e in enumerate(m):
            for _ in range(e):value=apply(canonical[k],value)
        assert value=={j:F(1)}

    def generator_apply(side,g,value):
        if g==7:
            for _ in range(3):value=apply(canonical[7],value)
            return value
        if side=='A' or g:
            return apply(canonical[g],value)
        r2=apply(canonical[7],apply(canonical[7],value))
        return plus(apply(canonical[0],value),r2)
    pairs=list(itertools.combinations_with_replacement(range(8),2));zero_products=[]
    for ai,aj in pairs:
        for bi,bj in pairs:
            value={0:F(1)}
            for side,g in [('A',ai),('A',aj),('B',bi),('B',bj)]:value=generator_apply(side,g,value)
            assert not value
            zero_products.append([ai,aj,bi,bj])
    assert len(zero_products)==1296
    # Constant GL readback gives the actual original CA-F9 coordinate actions.
    old=[action_sum([{} for _ in range(348)],canonical[0],c[0])]
    for i in range(1,7):old.append(action_sum(canonical[i],canonical[0],c[i]))
    old.append(canonical[7])
    for i in range(7):
        rebuilt=[{} for _ in range(348)]
        for j in range(7):rebuilt=action_sum(rebuilt,old[j],forward[i][j])
        assert rebuilt==canonical[i]
    # Weighted homogeneity gives a tiny exact nilpotence check, without
    # enumerating43758 ordinary monomials of degree10.
    for gb in bases:
        for g in gb:assert len({2*z+r for z,r in g})==1
    weights=[2*sum(m[:7])+m[7] for m in normal_basis]
    assert max(weights)==9 and normal_monomial((0,)*7+(9,))
    source_nf=dict(actual_base='Rbar=R0/(ell), actual CA-F9 zA7,r_CA with s_CA=r_CA',
        source_file='reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
        source_A0_determinant='-1/1024',gA_leading='125/24*r^3',gB_on_graph_leading='-15125/144*r^3',
        c7=list(map(str,c)),
        allorder_pair_proof='IA zero7 are actual CA-F9 coordinates. Its own-plane gA(r,r) is r^3 times a unit. For IB the actual zero7 Jacobian in zA7 is a unit, so formal IFT gives a unique graph z=phi(r) with phi=-c*r^2+O(r^3). Use the actual eighth generator after its own zero7 Hadamard adjustment; its leading cubic along this graph is -15125/144*r^3, unchanged by z=O(r^2) and higher coordinate corrections. It generates(r^3) in the graph quotient. Consequently IB=(z-phi(r),r^3)=(z+c*r^2,r^3) exactly. All higher phi terms are absorbed by r^3. No all-order conclusion is inferred merely from equality of a truncated length20 quotient.',
        actual_pair_before_constant_GL=['A=(zA7,r^3)','B=(zA7+c*r^2,r^3)'],
        canonical_coordinates=['Z=zA0/c0']+[f'W{i}=zA{i}-(c{i}/c0)*zA0' for i in range(1,7)]+['r=r_CA'],
        constant_GL7_forward=[[str(v) for v in row] for row in forward],
        constant_GL7_inverse=[[str(v) for v in row] for row in inverse],
        constant_GL7_determinant=str(1/c[0]),
        actual_source_readbacks=['zA0=c0*Z']+[f'zA{i}=W{i}+c{i}*Z' for i in range(1,7)]+['r_CA=r','s_CA=r (ell=0 only)'],
        canonical_pair=['A=(Z,W1..W6,r^3)','B=(Z+r^2,W1..W6,r^3)'],
        coefficient_ideal_decomposition=['W-degree0:A0^2*B0^2','W-degree1:A0*B0*(A0+B0)','W-degree2:(A0+B0)^2','W-degree3:A0+B0','W-degree>=4:whole'],
        A0='(Z,r^3)',B0='(Z+r^2,r^3)',M='A0+B0=(Z,r^2)',
        exact_ring='R0/(ell,(IA*IB)^2) = QQ[[Z,W1..W6,r]] / A^2*B^2',
        no_normalization_before_ell_or_physical_rho_substitution=True)
    evidence=dict(source_normalization=source_nf,two_variable_certificates=certificates,
        canonical_variable_order=['Z','W1','W2','W3','W4','W5','W6','r'],
        standard_basis348=[list(m) for m in normal_basis],block_dimensions=block_counts,
        block_W_degree_multiplicities=[1,6,21,56],two_variable_block_lengths=[26,14,6,2],
        canonical_coordinate_actions8_column_sparse=[records(A) for A in canonical],
        actual_CA_F9_coordinate_order=[f'zA_{i}' for i in range(7)]+['r_CA'],
        actual_CA_F9_actions8_column_sparse=[records(A) for A in old],
        all1296_A2_B2_generator_products_annihilate_unit=zero_products,
        generator_action_convention='A/B source indices0..6 denote Z,W1..W6; index7 is r^3; B index0 is Z+r^2',
        all8_actions_commute=True,all348_basis_monomials_generate_from_unit=True,
        whole_module_annihilation_proof='Each1296 source product kills the unit. Every basis vector is a monomial in the commuting coordinate actions applied to the unit, so each product kills every348 basis vector.',
        exact_length348=True,max_weight_standard_monomial=9,weight_vector=[2]*7+[1],
        maximal_ideal_nilpotence_index_exact10=True,
        nilpotence_proof='All fixed block polynomials are homogeneous in weights(Z,W)=2,r=1. Normal forms preserve weight, and all standard weights are at most9. Thus all ordinary degree10 products vanish, while r^9 is a nonzero standard monomial.',
        actual_T_or_K_matrix_Tor_image_or_torsion_not_computed=True)
    terms=sum(len(v) for A in canonical+old for v in A)
    assert terms<150000
    (out/'actual-conductor-square348-ring.json').write_text(json.dumps(evidence,indent=2)+'\n')
    proof=dict(all_passed=True,exact_actual_ell0_pair_normalization_proved=True,
        fixed_source_GB_cofactors_verified=True,all20_fixed_S_pair_remainders_zero=True,
        all_original_source_generators_reduce_zero=True,standard_block_lengths=[26,14,6,2],
        W6_multiplicities=[1,6,21,56],length348=True,
        all8_sparse_actions_commute_and_all348_basis_words_correct=True,
        all1296_original_canonical_source_products_annihilate_whole_ring=True,
        actual_CA_F9_readback8_actions_checked=True,ring_m10_zero_m9_nonzero=True,
        no_source_solve_jet_GB_search_kernel_T_matrix_Tor_or_K_torsion_claim=True,
        serialized_action_nonzero_coefficients=terms)
    (out/'finite-polynomial-and-action-checks.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,actual_R0_mod_ell_I_squared_length=348,
        exact_source_pair_normalization=True,sparse_canonical_and_actual_actions=True,
        no_actual_T_K_torsion_comparison_matrix_or_depth_claim=True,
        input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        frozen=True,independent_review_pending=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
