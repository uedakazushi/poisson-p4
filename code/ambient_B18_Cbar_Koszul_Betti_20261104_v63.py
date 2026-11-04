#!/usr/bin/env python3
"""Exact sparse Koszul homology of the accepted actual eight-dimensional Cbar.

Only finite rational arithmetic on the eight frozen coordinate actions is used.
Every column elimination is recorded; no source solve, CAS, field point or prime.
"""
import argparse, hashlib, itertools, json, math, resource, signal, time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_CBAR_KOSZUL_BETTI_20261104_v63'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v): return [[i, str(c)] for i, c in sorted(v.items()) if c]
def addto(v, w, c=F(1)):
    for i, a in w.items():
        b = v.get(i, F(0)) + c*a
        if b: v[i] = b
        else: v.pop(i, None)

def rank_receipt(columns, nrows):
    pivots = []; byrow = {}; records = []
    for j, col in enumerate(columns):
        v = dict(col); steps = []
        while v and min(v) in byrow:
            p = byrow[min(v)]; c = v[min(v)]
            steps.append([p, str(c)]); addto(v, pivots[p], -c)
        r = {'column': j, 'subtract': steps}
        if v:
            i = min(v); c = v[i]; v = {k:a/c for k,a in v.items()}
            assert 0 <= i < nrows
            byrow[i] = len(pivots); pivots.append(v)
            r.update({'new_pivot_row':i, 'normalizing_coefficient':str(c), 'normalized':enc(v)})
        else: r['zero_remainder'] = True
        records.append(r)
    # Exact replay independently reconstructs each residual from its source.
    saved = []
    for col, r in zip(columns, records):
        v = dict(col)
        for p,c in r['subtract']: addto(v, saved[p], -F(c))
        if r.get('zero_remainder'): assert not v
        else:
            c = F(r['normalizing_coefficient']); w = {i:a/c for i,a in v.items()}
            assert enc(w) == r['normalized'] and min(w) == r['new_pivot_row']
            assert all(r['new_pivot_row'] != min(p) for p in saved)
            saved.append(w)
    return {'rank':len(pivots), 'nrows':nrows, 'ncols':len(columns),
            'source_columns':list(map(enc, columns)), 'columns':records}, pivots

def main(out):
    start=time.monotonic(); signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-physical12-actions-and-C-special.json'
    raw=ROOT/'reports/AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63/actual-raw26-actions-and-quotients.json'
    gate=ROOT/'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json'
    assert sha(source)=='50e3408db2fc0b1dd81c2cf653a3b0711ad9f6e501ca2322e2991f49df9f84f4'
    assert sha(gate)=='79dd7a7e8d5833a110a138b6625919e4cf83bdc0f184b8533a964eba5303a873'
    data=json.loads(source.read_text()); rd=json.loads(raw.read_text())
    assert json.loads(gate.read_text())['all_passed']
    assert data['C_special_dimension']==8 and data['C_basis_raw26_indices']==[0,1,2,3,4,5,6,14]
    names=rd['actual_CA_F9_coordinate_order9']; assert names==[f'zA_{i}' for i in range(7)]+['r_CA','s_CA']
    assert rd['ninth_s_CA_action_equals_coordinate7_r_CA'] and rd['ell_zero_is_exact_s_equals_r']
    T=[[{i:F(c) for i,c in col} for col in mat] for mat in data['C_coordinate_actions8_column_sparse']]
    assert len(T)==8 and all(len(mat)==8 for mat in T)
    for a,b in itertools.product(range(8),repeat=2):
        for j in range(8):
            v={}; w={}
            for k,c in T[b][j].items(): addto(v,T[a][k],c)
            for k,c in T[a][j].items(): addto(w,T[b][k],c)
            assert v==w
    wedges=[list(itertools.combinations(range(8),k)) for k in range(9)]
    dims=[8*len(w) for w in wedges]; D=[None]; ranks=[0]; receipts={}
    preflight=[]; total=0
    for k in range(1,9):
        ix={s:i for i,s in enumerate(wedges[k-1])}; cols=[]
        for S in wedges[k]:
            for j in range(8):
                v={}
                for p,a in enumerate(S):
                    base=8*ix[S[:p]+S[p+1:]]
                    addto(v,{base+i:c for i,c in T[a][j].items()},F((-1)**p))
                cols.append(v)
        nz=sum(map(len,cols)); total+=nz
        preflight.append({'degree':k,'rows':dims[k-1],'columns':dims[k],'nonzeros':nz})
        assert nz==8*math.comb(7,k-1)
        D.append(cols)
    assert total==1024 and total<150000
    (out/'preflight.json').write_text(json.dumps({'matrices':preflight,'total_nonzeros':total,'dense_fallback':False},indent=2)+'\n')
    for k in range(1,9):
        rec,_=rank_receipt(D[k],dims[k-1]); ranks.append(rec['rank']);receipts[str(k)]=rec
    zero_columns=0
    for k in range(2,9):
        for col in D[k]:
            v={}
            for j,c in col.items(): addto(v,D[k-1][j],c)
            assert not v; zero_columns+=1
    betti=[dims[k]-ranks[k]-(ranks[k+1] if k<8 else 0) for k in range(9)]
    # Full maximal-ideal filtration and socle, still using only the same actions.
    current=[{i:F(1)} for i in range(8)]; filtration=[8]; filtration_receipts=[]
    while current:
        columns=[]
        for mat in T:
            for v in current:
                w={}
                for j,c in v.items(): addto(w,mat[j],c)
                columns.append(w)
        r,current=rank_receipt(columns,8);filtration.append(r['rank']);filtration_receipts.append(r)
        assert len(filtration)<10
    stack=[{8*a+i:c for a in range(8) for i,c in T[a][j].items()} for j in range(8)]
    sr,_=rank_receipt(stack,64);socle=8-sr['rank']
    assert ranks[1:]==[5,35,105,175,175,105,35,5]
    assert betti==[3,24,84,168,210,168,84,24,3]
    assert filtration==[8,5,2,0] and socle==3
    evidence={'actual_coordinate_order8':names[:8], 'common_r_is_CA_not_CB':True,
      'basis_raw26_indices':data['C_basis_raw26_indices'],
      'basis_labels':[data['raw26_basis'][i] for i in data['C_basis_raw26_indices']],
      'action_matrices_column_sparse':data['C_coordinate_actions8_column_sparse'],
      'Koszul_wedge_order':wedges,'chain_dimensions':dims,'differential_ranks':ranks,
      'all_chain_square_zero_checked_columns':zero_columns,'all_action_commutators_zero':True,
      'Koszul_Betti_numbers':betti,'differential_rank_receipts':receipts,
      'maximal_ideal_filtration_dimensions':filtration,'filtration_rank_receipts':filtration_receipts,
      'Loewy_Hilbert_function':[filtration[i]-filtration[i+1] for i in range(3)],
      'minimal_generator_count':betti[0], 'socle_dimension':socle,'socle_rank_receipt':sr,
      'Cbar_is_module_no_unital_algebra_claim':True}
    ev=out/'sparse-Koszul-and-filtration.json';ev.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_Cbar_dimension':8,
      'regular_coefficient_ring_dimension':8,'minimal_Betti_numbers':betti,
      'differential_ranks':ranks[1:],'mu':3,'socle_dimension':socle,'Loewy_HF':[3,3,2],
      'Loewy_filtration_dimensions':filtration,'all_d_squared_zero':True,
      'sparse_only':True,'expanded_differential_nonzeros':total,'no_new_source_solve_or_GB_or_point_or_prime':True,
      'generic_C_rank_K_torsion_D0_depth_claimed':False,'independent_review_pending':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [source,raw,gate,Path(__file__),note]},
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','minimal_Betti_numbers','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
