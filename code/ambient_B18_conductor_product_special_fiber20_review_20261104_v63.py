#!/usr/bin/env python3
"""Independent fixed-coefficient and pivot-receipt audit of actual H/ell.

No producer code, source solve, new jet, symbolic CAS or Groebner basis.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b,c=F(1)):
    z=dict(a)
    for m,v in b.items():
        q=z.get(m,F(0))+c*v
        if q:z[m]=q
        else:z.pop(m,None)
    return z
def scale(a,c):return {m:c*v for m,v in a.items() if c*v}
def mul(a,b,cap):
    z={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(i+j for i,j in zip(m,n))
            if sum(k)>cap:continue
            q=z.get(k,F(0))+x*y
            if q:z[k]=q
            else:z.pop(k,None)
    return z
def var(n,i):m=[0]*n;m[i]=1;return {tuple(m):F(1)}
def exponent(n,indices):
    m=[0]*n
    for j in indices:m[j]+=1
    return tuple(m)
def decode(p):return {tuple(m):F(c) for m,c in p if F(c)}
def vec(p):return {int(i):F(c) for i,c in p if F(c)}
def hom(p,n):return {m:c for m,c in p.items() if sum(m)==n}
def compose(p,variables,cap):
    n=len(next(iter(variables[0])));one={(0,)*n:F(1)}
    answer={};powers={}
    for m,c in p.items():
        q=scale(one,c)
        for i,k in enumerate(m):
            if not k:continue
            if (i,k) not in powers:
                z=one
                for _ in range(k):z=mul(z,variables[i],cap)
                powers[i,k]=z
            q=mul(q,powers[i,k],cap)
        answer=add(answer,q)
    return answer
def matmul(A,B):return [[sum((x*y for x,y in zip(r,c)),F(0)) for c in zip(*B)] for r in A]
def identity(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def linear(A,ps):
    out=[]
    for r in A:
        p={}
        for c,q in zip(r,ps):p=add(p,q,c)
        out.append(p)
    return out
def rank20(rows):
    # Only fixed vector-space images inside the final dimension20 algebra.
    piv={}
    for p in rows:
        p=dict(p)
        while p:
            j=min(p)
            if j in piv:p=add(p,piv[j],-p[j])
            else:piv[j]=scale(p,1/p[j]);break
    return len(piv),list(piv.values())
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    b=ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63';p=b/'report.json';r=read(p)
    assert sha(p)=='fd6eac26d90e73599f61b91043ab55ad8c71ce9625388c56008795373e45231a' and r['all_passed']
    ep=ROOT/r['evidence']['path'];assert sha(ep)==r['evidence']['sha256'];E=read(ep)
    inputs={}
    for k,q in r['input_hashes'].items():
        assert sha(ROOT/q['path'])==q['sha256'];inputs[k]=read(ROOT/q['path'])
    assert inputs['length20_peer']['all_passed']
    full,joint,cb=inputs['full9_source2'],inputs['joint_plane3'],inputs['CB_plane3']
    assert joint['all105_equations_zero_through_degree3'] and cb['all57_equations_verified_through_degree3']
    mons9=[exponent(9,m) for m in full['source2_monomials']];vs9=[var(9,i) for i in range(9)]
    maps={}
    for family in ['CA','CB']:
        c=full['both57_source2'][family];original=inputs['original_source']['families'][family]['actual_zero_odd_I0_generators']
        assert original==E['source_conductor_generators_rows_and_ids'][family]
        assert [q[1] for q in original[:7]]==c['F9_nodes'][:7] and original[7]==[63,c['f63']]
        assert c['all57_equations_all55_coefficients_verified'] and c['unknown57'][-9:]==list(range(2,11))
        cache=dict(c['all57_source2jets'])
        maps[family]=[{m:F(q) for m,q in zip(mons9,cache[i]) if F(q)} for i in range(2,11)]
        assert maps[family]==list(map(decode,E['physical_coordinate_maps2'][family]))
    L={f:[[p.get(exponent(9,[i]),F(0)) for i in range(9)] for p in maps[f]] for f in maps}
    inverse=[[F(c) for c in row] for row in E['physical_CB_linear_inverse']]
    assert matmul(inverse,L['CB'])==matmul(L['CB'],inverse)==identity(9)
    YB1=linear(inverse,[hom(q,1) for q in maps['CA']]);YB2=list(map(decode,E['YB_in_YA_quadratic']))
    expected2=linear(inverse,[add(hom(a,2),compose(hom(b,2),YB1,2),-1) for a,b in zip(maps['CA'],maps['CB'])])
    assert YB2==expected2
    assert [compose(q,[add(a,b) for a,b in zip(YB1,YB2)],2) for q in maps['CB']]==maps['CA']
    Ylinear=[[p.get(exponent(9,[i]),F(0)) for i in range(9)] for p in YB1]
    assert Ylinear==[[F(c) for c in row] for row in E['YB_in_YA_linear']]==[[F(c) for c in row] for row in full['YB_to_YA_linear']]
    assert all(Ylinear[i][j]==0 for i in range(7) for j in [7,8])
    assert all(Ylinear[i][j]==0 for i in [7,8] for j in range(7))
    A0=[[F(c) for c in row] for row in E['A0']];Ai=[[F(c) for c in row] for row in E['A0_inverse']]
    assert A0==[row[:7] for row in Ylinear[:7]] and matmul(A0,Ai)==matmul(Ai,A0)==identity(7)
    A1=[[decode(q) for q in row] for row in E['A1_ordered_prefix']]
    b2=[]
    for i,q in enumerate(YB2[:7]):
        reconstructed={}
        for j in range(7):
            for m,c in A1[i][j].items():
                n=list(m);n[j]+=1
                # This is exactly the declared MAX-index ordered-prefix rule.
                assert j==max(k for k in range(7) if n[k])
            reconstructed=add(reconstructed,mul(A1[i][j],vs9[j],2))
        plane={m:c for m,c in q.items() if not any(m[:7])};b2.append(plane)
        assert add(reconstructed,plane)==q
    plane=[]
    for rec in joint['all8_CB_zero_restrictions'][:7]:
        plane.append({(0,)*7+tuple(m):F(c) for m,c in zip(joint['monomials'],rec['all_degree_le3_coefficients']) if F(c)})
    assert plane==list(map(decode,E['b_plane3'])) and [hom(q,2) for q in plane]==b2
    beta2=linear(Ai,b2)
    delta=[]
    for i in range(7):
        q=hom(plane[i],3)
        for j in range(7):q=add(q,mul(A1[i][j],beta2[j],3),-1)
        delta.append(q)
    beta3=linear(Ai,delta)
    assert beta2==list(map(decode,E['beta2_before_ell'])) and beta3==list(map(decode,E['beta3_before_ell']))
    v8=[var(8,i) for i in range(8)];restrict=v8[:7]+[v8[7],v8[7]]
    b28=[compose(q,restrict,3) for q in beta2];b38=[compose(q,restrict,3) for q in beta3]
    assert b28==list(map(decode,E['beta2_after_ell'])) and b38==list(map(decode,E['beta3_after_ell']))
    gA={tuple(m):F(c) for m,c in zip(joint['monomials'],joint['g_CA_all10_coefficients']) if F(c)}
    gB={tuple(m):F(c) for m,c in zip(cb['ordered_monomials'],cb['g_CB_all10_coefficients']) if F(c)}
    assert all(sum(m)==3 for q in [gA,gB] for m in q)
    assert gA==decode(E['gA_plane3']) and gB==decode(E['gB_own_plane3'])
    gA8=compose(gA,[v8[7],v8[7]],4)
    gB8=compose(gB,[compose(YB1[j],restrict,1) for j in [7,8]],4)
    assert gA8==decode(E['gA_after_ell'])=={(0,)*7+(3,):F(125,24)}
    assert gB8==decode(E['gB_after_ell'])=={(0,)*7+(3,):F(-15125,144)}
    IA=v8[:7]+[gA8];IB=[add(v8[i],add(b28[i],b38[i])) for i in range(7)]+[gB8]
    assert IA==list(map(decode,E['IA_truncated_generators'])) and IB==list(map(decode,E['IB_truncated_generators']))
    products=[mul(a,b,4) for a in IA for b in IB];assert products==list(map(decode,E['all64_actual_pair_products_mod_m5']))
    assert [sum(1 for q in products if q and min(map(sum,q))==k) for k in [2,4]]==[49,14] and products[63]=={}
    monomials=[tuple(m) for m in E['monomial_order']];assert len(monomials)==len(set(monomials))==495
    expectedmons=set()
    for k in range(5):
        for m in itertools.combinations_with_replacement(range(8),k):expectedmons.add(exponent(8,m))
    assert set(monomials)==expectedmons
    index={m:i for i,m in enumerate(monomials)}
    rows=list(map(vec,E['all_relation_rows_sparse']));assert len(rows)==579
    seen=set()
    for row,(k,m) in zip(rows,E['relation_source_pair_and_multiplier']):
        key=(k,tuple(m));assert key not in seen;seen.add(key)
        q=mul(products[k],{tuple(m):F(1)},4);assert row=={index[a]:c for a,c in q.items()}
    for k,m,rid,c in E['duplicate_source_pair_multiplier_row_and_scale']:
        key=(k,tuple(m));assert key not in seen;seen.add(key)
        q=mul(products[k],{tuple(m):F(1)},4);assert {index[a]:v for a,v in q.items()}==scale(rows[rid],F(c))
    required={(k,m) for k,q in enumerate(products) if q for m in monomials if sum(m)<=4-min(map(sum,q))}
    assert seen==required
    # Replay only the supplied fixed receipts. No rank search/basis discovery.
    piv={}
    for rec in E['pivot_receipts']:
        q=dict(rows[rec['source_row']])
        for j,c in rec['subtract']:
            assert j in piv and min(q)==j and q[j]==F(c)
            q=add(q,piv[j],-F(c))
        j=rec['pivot'];assert q and min(q)==j and j not in piv and q[j]==F(rec['divide_by'])
        piv[j]=scale(q,1/F(rec['divide_by']))
    assert len(piv)==475 and piv=={j:vec(q) for j,q in E['echelon_rows']}
    free=[j for j in range(495) if j not in piv];assert len(free)==20
    assert [monomials[j] for j in free]==list(map(tuple,E['quotient_basis_monomials']))
    reductions=list(map(vec,E['all495_monomial_reductions']))
    for j in reversed(range(495)):
        if j in free:assert reductions[j]=={free.index(j):F(1)}
        else:
            q={}
            for k,c in piv[j].items():
                if k!=j:assert k>j;q=add(q,reductions[k],-c)
            assert reductions[j]==q
    for row in rows:
        q={}
        for j,c in row.items():q=add(q,reductions[j],c)
        assert not q
    tables=[[vec(q) for q in a] for a in E['multiplication8_column_sparse']]
    def action(A,v):
        q={}
        for j,c in v.items():q=add(q,A[j],c)
        return q
    for i in range(8):
        for k,j in enumerate(free):
            m=list(monomials[j]);m[i]+=1
            assert tables[i][k]==(reductions[index[tuple(m)]] if sum(m)<=4 else {})
    for i,j,k in itertools.product(range(8),range(8),range(20)):
        assert action(tables[i],tables[j][k])==action(tables[j],tables[i][k])
    # Different annihilation check: apply the actual matrices, rather than
    # multiplying and reducing polynomial records as in the producer.
    for q in products:
        for k in range(20):
            result={}
            for m,c in q.items():
                v={k:c}
                for i,e in enumerate(m):
                    for _ in range(e):v=action(tables[i],v)
                result=add(result,v)
            assert not result
    powerbasis=[{j:F(1)} for j in range(20)];dims=[20]
    for _ in range(5):
        dim,powerbasis=rank20([action(A,v) for A in tables for v in powerbasis]);dims.append(dim)
    assert dims==E['m_power_dimensions']==[20,19,11,3,1,0]
    hf=[dims[j]-dims[j+1] for j in range(5)];assert hf==[1,8,8,2,1]==E['hilbert_function']
    proof={'physical9_and_F9_source_ancestry_bound':True,'both_physical_maps2_and_exact_composition_verified':True,
      'both_off_diagonal_linear_blocks_zero':True,'A0_inverse_and_ordered_A1_exact':True,
      'actual_beta2_beta3_before_after_ell_all_coefficients_verified':True,
      'actual_gA_gB_source_cubic_images':['125/24*r^3','-15125/144*r^3'],
      'actual_all64_products_orders':[49,14,1],'all579_rows_and_all_aliases_and_complete_multiplier_coverage':True,
      'fixed_pivot_receipts_replayed':475,'all495_reductions_verified':True,
      'all8_tables_commute_and_all64_products_annihilate20_basis':True,
      'm_power_dimensions':dims,'hilbert_function':hf,
      'same_length_upgrade':'accepted actual length20 plus independently computed truncated length20 makes actual->truncated iso',
      'actual_m5_zero_and_tables_exact':True,
      'frame':'CA F9 z7,r with s=r; rawphysicalrho readback is not constructed',
      'no_K_torsion_or_D0_depth_conclusion':True}
    ev=out/'independent-source-pivots-and-actions.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',p,ep,ROOT/'code/ambient_B18_conductor_product_special_fiber20_20261104_v63.py',ROOT/'research/notes/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63.md']+[ROOT/q['path'] for q in r['input_hashes'].values()]
    rep={'id':STEM,'all_passed':True,'source_bound_actual_special_fiber20_and_tables_accepted':True,
      'actual_HF_1_8_8_2_1_m5_zero_accepted':True,'fixed_receipts_and_every_source_coefficient_independently_checked':True,
      'no_author_producer_newsource_jet_GB_or_large_rank_search':True,
      'no_raw_physical_readback_K_torsion_or_D0_depth_claim':True,
      'input_sha256':{str(q.relative_to(ROOT)):sha(q) for q in paths},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['id','all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
