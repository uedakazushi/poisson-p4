#!/usr/bin/env python3
"""Independent fixed Buchberger certificates and sparse348 action audit.

No author import, basis search, dense348 matrix, source solve or new jet.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b,c=F(1)):
    z=dict(a)
    for m,x in b.items():
        q=z.get(m,F(0))+c*x
        if q:z[m]=q
        else:z.pop(m,None)
    return z
def scale(a,c):return {m:x*c for m,x in a.items() if x*c}
def mul(a,b):
    z={}
    for m,c in a.items():
        for n,d in b.items():
            k=tuple(i+j for i,j in zip(m,n));q=z.get(k,F(0))+c*d
            if q:z[k]=q
            else:z.pop(k,None)
    return z
def decode(p):return {tuple(m):F(c) for m,c in p if F(c)}
def sparse(p):return {int(i):F(c) for i,c in p if F(c)}
def action(A,v):
    q={}
    for j,c in v.items():q=add(q,A[j],c)
    return q
def remainder(p,gb):
    # Opposite divisor priority from the producer, after fixed GB proof.
    p=dict(p);r={};n=0
    while p:
        n+=1;assert n<1000
        m=max(p);c=p[m];i=next((j for j in reversed(range(len(gb))) if all(a>=b for a,b in zip(m,max(gb[j])))),None)
        if i is None:r[m]=c;del p[m]
        else:
            leading=max(gb[i]);v={tuple(a-b for a,b in zip(m,leading)):c/gb[i][leading]}
            p=add(p,mul(v,gb[i]),-1)
    return r
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    b=ROOT/'reports/AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63';p=b/'report.json';R=read(p)
    assert sha(p)=='93d5f8bf19fb97200209fa7500fe6eac8c2d467e23805891fbb31f72512f0297' and R['all_passed']
    for q,h in R['input_sha256'].items():assert sha(ROOT/q)==h,q
    for q,h in R['evidence_sha256'].items():assert sha(b/q)==h,q
    ep=b/'actual-conductor-square348-ring.json';E=read(ep)
    source=ROOT/E['source_normalization']['source_file'];H=read(source)
    peer=ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json';assert read(peer)['all_passed']
    N=E['source_normalization'];c=[F(q[0][1]) if q else F(0) for q in H['beta2_after_ell']]
    assert c==list(map(F,N['c7'])) and c[0]
    assert F(H['A0_determinant'])==F(-1,1024)
    assert H['gA_after_ell']==[[[0]*7+[3],'125/24']] and H['gB_after_ell']==[[[0]*7+[3],'-15125/144']]
    fw=[[F(q) for q in r] for r in N['constant_GL7_forward']];iv=[[F(q) for q in r] for r in N['constant_GL7_inverse']]
    for i,j in itertools.product(range(7),repeat=2):
        assert sum(fw[i][k]*iv[k][j] for k in range(7))==F(i==j)
        assert sum(iv[i][k]*fw[k][j] for k in range(7))==F(i==j)
    assert [sum(fw[i][j]*c[j] for j in range(7)) for i in range(7)]==[F(1)]+[F(0)]*6
    Z={(1,0):F(1)};r2={(0,2):F(1)};r3={(0,3):F(1)};r6={(0,6):F(1)};Zp=add(Z,r2)
    Az=[mul(Z,Z),mul(Z,r3),r6];Bz=[mul(Zp,Zp),mul(Zp,r3),r6]
    sources=[[mul(a,b) for a in Az for b in Bz],
      [mul(mul(a,b),q) for a in [Z,r3] for b in [Zp,r3] for q in [Z,r2]],
      [mul(Z,Z),mul(Z,r2),mul(r2,r2)],[Z,r2]]
    bases=[];standards=[];spairs=0;source_checks=0
    for k,rec in enumerate(E['two_variable_certificates']):
        src=list(map(decode,rec['source_generators']));gb=list(map(decode,rec['fixed_lex_basis']))
        assert src==sources[k]
        assert all(g[max(g)]==1 for g in gb)
        for g,cof in zip(gb,rec['source_constant_cofactors']):
            q={}
            for j,c0 in cof:q=add(q,src[j],F(c0))
            assert q==g
        for q,qs in zip(src,rec['all_original_source_division_quotients']):
            v={}
            for t,g in zip(qs,gb):v=add(v,mul(decode(t),g))
            assert v==q and not remainder(q,gb);source_checks+=1
        assert {tuple(q['pair']) for q in rec['all_S_pair_divisions']}==set(itertools.combinations(range(len(gb)),2))
        for q in rec['all_S_pair_divisions']:
            i,j=q['pair'];a,b=max(gb[i]),max(gb[j]);lcm=tuple(max(x,y) for x,y in zip(a,b))
            S=add(mul({tuple(x-y for x,y in zip(lcm,a)):F(1)},gb[i]),mul({tuple(x-y for x,y in zip(lcm,b)):F(1)},gb[j]),-1)
            assert S==decode(q['S']) and not q['remainder']
            value={}
            for coeff,g in zip(q['quotients'],gb):
                coeff=decode(coeff)
                if coeff:
                    assert tuple(x+y for x,y in zip(max(coeff),max(g)))<lcm
                value=add(value,mul(coeff,g))
            assert value==S and not remainder(S,gb);spairs+=1
        leading=[max(g) for g in gb]
        standard=[(z,r) for z in range(max(x for x,y in leading)) for r in range(max(y for x,y in leading)) if not any(z>=a and r>=b for a,b in leading)]
        assert list(map(list,standard))==rec['standard_monomials'] and len(standard)==[26,14,6,2][k]
        for g in gb:assert len({2*z+r for z,r in g})==1
        bases.append(gb);standards.append(standard)
    assert spairs==20 and source_checks==22
    basis=[];counts=[];mult=[]
    for k in range(4):
        ws=[tuple(w.count(i) for i in range(6)) for w in itertools.combinations_with_replacement(range(6),k)]
        mult.append(len(ws));before=len(basis)
        for w in ws:
            for z,r in standards[k]:basis.append((z,*w,r))
        counts.append(len(basis)-before)
    assert mult==E['block_W_degree_multiplicities']==[1,6,21,56]
    assert counts==E['block_dimensions']==[26,84,126,112]
    assert len(basis)==348 and basis==list(map(tuple,E['standard_basis348']))
    idx={m:i for i,m in enumerate(basis)};canonical=[[sparse(q) for q in A] for A in E['canonical_coordinate_actions8_column_sparse']]
    for k in range(8):
        for j,m in enumerate(basis):
            n=list(m);n[k]+=1;q=sum(n[1:7])
            if q>=4:value={}
            else:value={idx[(z,*n[1:7],r)]:c for (z,r),c in remainder({(n[0],n[7]):F(1)},bases[q]).items()}
            assert canonical[k][j]==value
    for i,j,k in itertools.product(range(8),range(8),range(348)):
        assert action(canonical[i],canonical[j][k])==action(canonical[j],canonical[i][k])
    for j,m in enumerate(basis):
        q={0:F(1)}
        for k,e in enumerate(m):
            for _ in range(e):q=action(canonical[k],q)
        assert q=={j:F(1)}
    pairs=list(itertools.combinations_with_replacement(range(8),2));expected={(a,b,c,d) for a,b in pairs for c,d in pairs}
    assert set(map(tuple,E['all1296_A2_B2_generator_products_annihilate_unit']))==expected
    for ai,aj,bi,bj in expected:
        q={0:F(1)}
        # Reverse application order is independent of the producer; the
        # already checked commutation makes every source order equivalent.
        for side,k in [('B',bj),('B',bi),('A',aj),('A',ai)]:
            if k==7:
                for _ in range(3):q=action(canonical[7],q)
            elif k==0 and side=='B':q=add(action(canonical[0],q),action(canonical[7],action(canonical[7],q)))
            else:q=action(canonical[k],q)
        assert not q
    actual=[[sparse(q) for q in A] for A in E['actual_CA_F9_actions8_column_sparse']]
    for i in range(7):
        for j in range(348):
            target={}
            for k in range(7):target=add(target,canonical[k][j],iv[i][k])
            assert actual[i][j]==target
    assert actual[7]==canonical[7]
    weights=[2*sum(m[:7])+m[7] for m in basis];assert max(weights)==9 and (0,)*7+(9,) in idx
    proof={'own_engine':'sparse Fraction polynomials with opposite fixed-divisor priority and sparse vector actions',
      'exact_allorder_actual_pair_normalization_accepted':True,'source_c0_GL7_and_both_cubic_units_verified':True,
      'source_cofactor_and_original_generator_identities':source_checks,
      'all20_Spair_standard_representations_strict_below_LCM_verified':True,
      'fixed_lex_Buchberger_criterion_proved_without_basis_search':True,
      'coefficient_ideal_lengths':[26,14,6,2],'W_multiplicities':mult,'total_length':348,
      'all8_canonical_actions_and_all348_cyclic_words_verified':True,
      'all8_actual_CA_F9_readback_actions_verified':True,
      'all1296_source_products_kill_unit_and_whole_by_commutation_cyclicity':True,
      'base_m10_zero_r9_nonzero':True,
      'no_dense348_matrix_source_solve_jet_GBsearch_T_K_or_torsion_claim':True}
    ev=out/'independent-fixed-bases-and-sparse-actions.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',p,ep,source,peer,
      ROOT/'code/ambient_B18_conductor_square_closed348_20261104_v63.py',ROOT/'research/notes/AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63.md']
    rep={'id':STEM,'all_passed':True,'actual_ell0_conductor_product_square_length348_accepted':True,
      'allorder_source_pair_unit_normalization_proved':True,'fixed20Spairs22sourceidentities_and_sparse_actions_verified':True,
      'actual_base_m10_zero_r9_nonzero_accepted':True,'no_T_K_socle_torsion_D0_or_generic_DVR_claim':True,
      'no_author_producer_newsource_jet_GBsearch_or_dense348_matrices':True,
      'input_sha256':{str(q.relative_to(ROOT)):sha(q) for q in paths},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['id','all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
