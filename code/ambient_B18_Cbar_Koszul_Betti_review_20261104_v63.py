#!/usr/bin/env python3
"""Independent source-bound sparse Koszul and prescribed rank-receipt replay."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_CBAR_KOSZUL_BETTI_20261104_v63_REVIEW'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decode(v):return {i:F(c) for i,c in v if F(c)}
def add(a,b,c=F(1)):
    out=dict(a)
    for i,v in b.items():
        q=out.get(i,F(0))+c*v
        if q:out[i]=q
        else:out.pop(i,None)
    return out
def apply(A,v):
    out={}
    for j,c in v.items():out=add(out,A[j],c)
    return out


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def gate(stem):
        p=ROOT/'reports'/stem/'report.json';r=read(p);assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            q=p.parent/name;assert sha(q)==h,name;hashes[str(q.relative_to(ROOT))]=h
        return r
    gate('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW')
    primary=gate('AMBIENT_B18_CBAR_KOSZUL_BETTI_20261104_v63')
    b=ROOT/'reports'
    source=read(b/'AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-physical12-actions-and-C-special.json')
    raw=read(b/'AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63/actual-raw26-actions-and-quotients.json')
    ev=read(b/'AMBIENT_B18_CBAR_KOSZUL_BETTI_20261104_v63/sparse-Koszul-and-filtration.json')
    names=[f'zA_{i}' for i in range(7)]+['r_CA']
    assert raw['actual_CA_F9_coordinate_order9']==names+['s_CA']
    assert raw['ell_zero_is_exact_s_equals_r'] and raw['ninth_s_CA_action_equals_coordinate7_r_CA']
    assert names==ev['actual_coordinate_order8'] and ev['common_r_is_CA_not_CB']
    assert source['C_special_dimension']==8
    assert source['C_basis_raw26_indices']==ev['basis_raw26_indices']==[0,1,2,3,4,5,6,14]
    assert ev['basis_labels']==[raw['raw26_basis'][i] for i in ev['basis_raw26_indices']]
    assert ev['action_matrices_column_sparse']==source['C_coordinate_actions8_column_sparse']
    X=[[decode(c) for c in A] for A in source['C_coordinate_actions8_column_sparse']]
    assert len(X)==8 and all(len(A)==8 for A in X)
    for A,B in itertools.product(X,repeat=2):
        for j in range(8):assert apply(A,B[j])==apply(B,A[j])
    wedge=[list(itertools.combinations(range(8),k)) for k in range(9)]
    assert [[[i for i in s] for s in w] for w in wedge]==ev['Koszul_wedge_order']
    dims=[8*math.comb(8,k) for k in range(9)]
    assert dims==ev['chain_dimensions']==[8,64,224,448,560,448,224,64,8]
    D=[None];nonzeros=[]
    for k in range(1,9):
        position={w:j for j,w in enumerate(wedge[k-1])};cols=[]
        for w in wedge[k]:
            for j in range(8):
                value={}
                for p,i in enumerate(w):
                    offset=8*position[w[:p]+w[p+1:]]
                    term={offset+r:c for r,c in X[i][j].items()}
                    value=add(value,term,F((-1)**p))
                cols.append(value)
        D.append(cols);nonzeros.append(sum(map(len,cols)))
    assert sum(nonzeros)==primary['expanded_differential_nonzeros']==1024
    assert sum(nonzeros)<150000
    (out/'preflight.json').write_text(json.dumps(dict(
        source='accepted actual CA-F9 coordinate actions on Cbar8',
        dimensions=dims,differential_nonzeros=nonzeros,total_nonzeros=1024,
        prescribed_receipts_only_no_new_pivot_selection=True,sparse_only=True,
        no_source_solve_jet_point_prime_GB=True,
        limits=dict(seconds=60,address_bytes=2_000_000_000,output_bytes=32*1024*1024,expanded_coefficients=150000)),indent=2)+'\n')

    # A rank is proved in BOTH directions: independent prescribed pivots plus
    # zero reductions of every original column. No pivot or rank is searched.
    def replay(columns,nrows,receipt):
        assert nrows==receipt['nrows'] and len(columns)==receipt['ncols']
        assert columns==[decode(v) for v in receipt['source_columns']]
        assert len(receipt['columns'])==len(columns)
        piv=[];leading=[]
        for j,(column,rec) in enumerate(zip(columns,receipt['columns'])):
            assert rec['column']==j
            value=dict(column)
            for k,c in rec['subtract']:
                assert 0<=k<len(piv)
                # Every prescribed reduction really clears that pivot.
                assert value.get(leading[k],F(0))==F(c)
                value=add(value,piv[k],-F(c))
            if rec.get('zero_remainder'):
                assert not value
            else:
                r=rec['new_pivot_row'];c=F(rec['normalizing_coefficient'])
                assert c and value.get(r)==c and 0<=r<nrows and r not in leading
                assert all(i>=r for i in value)
                normalized={i:a/c for i,a in value.items()}
                assert normalized==decode(rec['normalized'])
                piv.append(normalized);leading.append(r)
        assert len(piv)==receipt['rank']
        # Every saved pivot is a nonzero combination of source columns and has
        # distinct leading row, so independence is exact over QQ.
        return piv

    ranks=[0]
    for k in range(1,9):
        piv=replay(D[k],dims[k-1],ev['differential_rank_receipts'][str(k)])
        ranks.append(len(piv))
    assert ranks==ev['differential_ranks']==[0,5,35,105,175,175,105,35,5]
    chain_columns=0
    for k in range(2,9):
        for column in D[k]:assert not apply(D[k-1],column);chain_columns+=1
    assert chain_columns==ev['all_chain_square_zero_checked_columns']
    betti=[dims[k]-ranks[k]-(ranks[k+1] if k<8 else 0) for k in range(9)]
    assert betti==ev['Koszul_Betti_numbers']==primary['minimal_Betti_numbers']==[3*math.comb(8,k) for k in range(9)]
    # Reconstruct the filtration column inputs from the preceding verified
    # basis, then replay the fixed receipts, rather than run a rank algorithm.
    current=[{i:F(1)} for i in range(8)];filtration=[8]
    for receipt in ev['filtration_rank_receipts']:
        columns=[apply(A,v) for A in X for v in current]
        current=replay(columns,8,receipt);filtration.append(len(current))
    assert not current and filtration==ev['maximal_ideal_filtration_dimensions']==[8,5,2,0]
    socle_columns=[]
    for j in range(8):socle_columns.append({8*a+i:c for a in range(8) for i,c in X[a][j].items()})
    socle_rank=len(replay(socle_columns,64,ev['socle_rank_receipt']))
    socle=8-socle_rank
    assert socle==ev['socle_dimension']==primary['socle_dimension']==3
    loewy=[filtration[k]-filtration[k+1] for k in range(3)]
    assert loewy==ev['Loewy_Hilbert_function']==[3,3,2]
    proof=dict(all_passed=True,actual_coordinate_labels_and_all8_actions_source_bound=True,
        actual_frame='Rbar=R0/(ell), regular8 with CA-F9 zA0..6,r_CA; not raw physical rho or r_CB',
        chain_dimensions=dims,source_differential_nonzeros=sum(nonzeros),
        all_differential_source_columns_reconstructed=True,
        all_fixed_elimination_receipts_replayed_no_new_pivot_search=True,
        differential_ranks=ranks[1:],all_d_squared_zero_columns=chain_columns,
        minimal_Rbar_Cbar_Betti=betti,mu=3,socle=3,Loewy_HF=loewy,
        maximal_ideal_filtration=filtration,
        minimal_resolution_proof='Rbar is the regular completed8-dimensional actual CA-F9 coefficient ring after ell=s-r. The Koszul complex on these parameters resolves its residue field. Tensor with actual Cbar8 computes Tor; dimensions of its homology therefore equal ranks in the minimal free Rbar-resolution. No explicit minimal differentials are claimed.',
        Cbar_is_module_not_a_unital_algebra=True,
        no_T_K_torsion_generic_DVR_rank_or_D0_depth_inference=True,
        no_source_solve_newjet_CAS_GB_point_prime_or_producer_execution=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,minimal_Rbar_Betti_numbers=betti,
        mu=3,socle_dimension=3,Loewy_HF=loewy,
        no_K_T_torsion_DVR_rank_or_D0_depth_claim=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
