#!/usr/bin/env python3
"""Independent fixed-table audit of the literal actual three-target Cbar resolution."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_proper_core_five_lines_review_20261104_v63 import add,mul,scale,variable,power,one
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63';STEM=PRIMARY+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def dec(p):return {tuple(m):F(c) for m,c in p if F(c)}
def psum(polys):
    ans={}
    for p in polys:ans=add(ans,p)
    return ans
def act(A,v):
    ans={}
    for j,c in v.items():
        for i,b in A[j].items():ans[i]=ans.get(i,F(0))+c*b
    return {i:c for i,c in ans.items() if c}
def paction(p,Ts,v):
    ans={}
    for m,c in p.items():
        w=dict(v)
        for i,k in enumerate(m):
            for _ in range(k):w=act(Ts[i],w)
        for i,b in w.items():ans[i]=ans.get(i,F(0))+c*b
    return {i:c for i,c in ans.items() if c}

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PRIMARY;rp=folder/'report.json';A=read(rp)
    assert sha(rp)=='8885adc6fb2cfeccef8fc1924616866c16e9f694288b8a596bc0f218e95aab77' and A['all_passed']
    ep=ROOT/A['evidence']['path'];assert sha(ep)==A['evidence']['sha256'];E=read(ep)
    D={}
    for k,rec in A['input_hashes'].items():
        p=ROOT/rec['path'];assert sha(p)==rec['sha256'];D[k]=read(p)
        if k.endswith('_peer'):assert D[k]['all_passed']
    assert E['input_hashes']==A['input_hashes']
    C=D['Cbar'];T=[[{i:F(c) for i,c in col} for col in M] for M in C['C_coordinate_actions8_column_sparse']]
    assert C['C_basis_raw26_indices']==[0,1,2,3,4,5,6,14]
    assert E['actual_three_targets']==['tauA','tauB','tauA*tauB']
    assert E['actual_three_target_Cbar_basis_indices']==E['actual_three_target_raw26_indices']==[0,3,6]
    gs=[{i:F(1)} for i in [0,3,6]];chains=[]
    for g,n in zip(gs,[3,3,2]):
        v=g
        for _ in range(n):chains.append(v);v=act(T[7],v)
        assert not v
    assert chains==[{i:F(1)} for i in range(8)]
    assert chains==[{i:F(c) for i,c in rec} for rec in E['Cbar8_actual_r_chain_basis_vectors']]
    cs=[-T[i][3].get(5,F(0)) for i in range(7)]
    assert cs==list(map(F,E['CB_shift_coefficients7']))==[F(9835,864),0,0,0,F(29875,3456),F(-4565,864),0]
    for i,c in enumerate(cs):
        expected=[{} for _ in range(8)]
        if c:expected[3]={5:-c}
        assert T[i]==expected
        beta=D['H20']['beta2_after_ell'][i]
        assert dec(beta)==scale(power(variable(8,7),2,8),c)
    seqs=[]
    for s,n in enumerate([3,3,2]):
        seq=[add(variable(8,i),scale(power(variable(8,7),2,8),cs[i])) if s==1 else variable(8,i) for i in range(7)]
        seq.append(power(variable(8,7),n,8));seqs.append(seq)
        assert all(not paction(p,T,gs[s]) for p in seq)
    assert [[dec(p) for p in seq] for seq in E['actual_cyclic_sector_regular_sequences8']]==seqs
    # Literal source augmentation: no constant change of the accepted eight basis.
    lift=E['Cbar8_lifts_in_actual_three_targets'];j=0
    for s,n in enumerate([3,3,2]):
        for k in range(n):
            assert len(lift[j])==1 and lift[j][0][0]==s
            p=dec(lift[j][0][1]);assert p==power(variable(8,7),k,8)
            assert paction(p,T,gs[s])=={j:F(1)};j+=1
    bases=[[(s,w) for s in range(3) for w in itertools.combinations(range(8),k)] for k in range(9)]
    assert [[[s,list(w)] for s,w in b] for b in bases]==E['resolution_basis_sector_and_wedge']
    Difs=[];nonzero=0;monomials=0
    for k,saved in enumerate(E['differentials8_column_sparse'],1):
        prior={(s,w):i for i,(s,w) in enumerate(bases[k-1])};assert len(saved)==len(bases[k])
        M=[]
        for (s,w),rec in zip(bases[k],saved):
            actual={i:dec(p) for i,p in rec}
            expected={prior[s,w[:q]+w[q+1:]]:scale(seqs[s][v],(-1)**q) for q,v in enumerate(w)}
            assert actual==expected
            assert all(not p.get((0,)*8,F(0)) for p in actual.values())
            M.append(actual);nonzero+=len(actual);monomials+=sum(map(len,actual.values()))
        Difs.append(M)
    chaincolumns=0
    for k in range(1,8):
        for c in Difs[k]:
            product={}
            for j,p in c.items():
                for i,q in Difs[k-1][j].items():product[i]=add(product.get(i,{}),mul(p,q))
            assert all(not p for p in product.values());chaincolumns+=1
    assert nonzero==3072==E['differential_nonzero_polynomial_entries'] and monomials<150000
    diag=list(map(dec,E['S_resolution_D_diagonal_polynomials']))
    assert diag==[power(variable(8,7),k,8) for k in [3,3,2]]
    lifts=[[[dec(p) for p in row] for row in M] for M in E['T7_lifted_z_actions3_by3']]
    for c,L in zip(cs,lifts):
        assert len(L)==3 and all(len(row)==3 for row in L)
        for i,j in itertools.product(range(3),repeat=2):
            assert L[i][j]==(scale(power(variable(8,7),2,8),-c) if i==j==1 else {})
            assert mul(L[i][j],diag[j])==mul(diag[i],L[i][j])
    for U,V in itertools.product(lifts,repeat=2):
        for i,j in itertools.product(range(3),repeat=2):
            assert psum(mul(U[i][k],V[k][j]) for k in range(3))==psum(mul(V[i][k],U[k][j]) for k in range(3))
    ranks=[len(b) for b in bases];assert ranks==[3,24,84,168,210,168,84,24,3]==E['ranks']
    ev={'all8_actual_source_actions_and_original_r_chains':True,'all8_augmentation_basis_lifts':True,
      'all3072_literal_polynomial_entries':True,'monomial_terms':monomials,
      'all_chain_compositions_zero_columns':chaincolumns,'all_differential_entries_maximalideal':True,
      'regularity_proof':'z_i -> z_i+c_i*r^2, r -> r is a triangular automorphism; each transformed sequence z7,r^n is regular',
      'source_module_isomorphism':'the three actual cyclic chains give a surjective map from quotients of dimensions3+3+2=8 to the accepted8basis; augmentation is the literal identity',
      'minimal_resolution_ranks':ranks,'S_actions_T_equal_S_and_commute':True,
      'H_ij_zero_all_coherences':True,'no_generic_DVR_rank_or_K_T_depth_computation':True}
    out_ev=out/'independent-literal-source-resolution.json';out_ev.write_text(json.dumps(ev,indent=2)+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md';files=[rp,ep,Path(__file__),note,ROOT/'code/ambient_B18_proper_core_five_lines_review_20261104_v63.py']+[ROOT/r['path'] for r in A['input_hashes'].values()]
    report={'id':STEM,'all_passed':True,'actual_Cbar_literal_three_sector_minimal_resolution_accepted':True,
      'all3072_entries_augmentation_and_chain_squares_independently_checked':True,'Betti_numbers':ranks,
      'no_authorworker_source_or_rank_search_jet_GB_point':True,'K_T_generic_rank_length_depth_computed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{out_ev.name:sha(out_ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
