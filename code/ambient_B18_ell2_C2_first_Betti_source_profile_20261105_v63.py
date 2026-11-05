#!/usr/bin/env python3
"""One fixed-source d2(0) sparse rank/factorization from actual C2 actions."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,product_columns,srec
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):hashes[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    r=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/report.json');assert r['all_passed']
    peer=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json');assert peer['all_passed']
    proof=read('reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/report.json');assert proof['all_passed']
    d=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    A=[[parse(v) for v in M] for M in d['C2_coordinate_actions9_sparse_columns']];L=A[8]
    d1=[scale(v,-1) for M in [L]+A[:8] for v in M]
    blocks1=[(1,())]+[(0,(i,)) for i in range(8)];rowblock={w:k for k,(p,w) in enumerate(blocks1)}
    blocks2=[(2-q,w) for q in range(3) for w in itertools.combinations(range(8),q)]
    cols=[];labels=[]
    for p,w in blocks2:
        for j in range(14):
            v={}
            if p:
                sign=-1 if p%2 else 1;off=14*rowblock[w]
                v=add(v,{off+i:sign*c for i,c in L[j].items()})
            for t,i in enumerate(w):
                target=w[:t]+w[t+1:];off=14*rowblock[target];sign=(-1)**(p+t+1)
                v=add(v,{off+k:sign*c for k,c in A[i][j].items()})
            cols.append(v);labels.append([p,list(w),j])
    assert len(cols)==518 and sum(map(len,cols))==304
    assert all(not mv(d1,v) for v in cols),'literal source d1(0)d2(0)=0'
    pre={'rows':126,'columns':518,'input_nonzeros':304,'dense_work_upper':126*518,
      'source_constant_d1_rank':11,'no_d3_rank_or_new_source_coefficient':True,
      'resource_limits':{'seconds':60,'bytes':2_000_000_000,'aggregate_artifact_bytes':32*1024*1024,'total_slots':150000}}
    (out/'preflight.json').write_text(json.dumps(pre,separators=(',',':'))+'\n')
    piv={};receipts=[]
    for j,orig in enumerate(cols):
        v=dict(orig);steps=[]
        while v:
            p=min(v)
            if p not in piv:
                c=v[p];piv[p]=scale(v,1/c);receipts.append({'source_column':j,'pivot':p,'divide_by':str(c),'subtract':[[i,str(a)] for i,a in steps],'normalized':srec(piv[p])});break
            c=v[p];steps.append((p,c));v=add(v,piv[p],-c)
        assert sum(map(len,piv.values()))+304+sum(len(q['subtract']) for q in receipts)<150000
    porder=sorted(piv);pidx={p:i for i,p in enumerate(porder)};factor=[]
    for orig in cols:
        v=dict(orig);c={}
        while v:
            p=min(v);assert p in piv;a=v[p];c[pidx[p]]=a;v=add(v,piv[p],-a)
        factor.append(c)
    B=[piv[p] for p in porder];assert all(mv(B,c)==orig for c,orig in zip(factor,cols))
    selected=[q['source_column'] for q in receipts];pivot_sequence=[q['pivot'] for q in receipts]
    inversions=sum(pivot_sequence[i]>pivot_sequence[j] for i in range(len(pivot_sequence)) for j in range(i+1,len(pivot_sequence)))
    det=F((-1)**inversions)
    for q in receipts:det*=F(q['divide_by'])
    assert det
    fixed_minor=[[str(cols[j].get(i,F(0))) for j in selected] for i in porder]
    beta1=126-11-len(piv);assert beta1>=0
    evidence={'blocks_F1':[[p,list(w)] for p,w in blocks1],'blocks_F2':[[p,list(w)] for p,w in blocks2],
      'column_labels':labels,'all_d2_constant_columns':[srec(v) for v in cols],'rank':len(piv),
      'fixed_pivot_receipts':receipts,'normalized_source_basis_columns':[srec(v) for v in B],
      'all518_column_factorization_in_normalized_basis':[srec(v) for v in factor],
      'fixed_original_minor_rows':porder,'fixed_original_minor_columns':selected,'fixed_original_minor':fixed_minor,
      'fixed_original_minor_determinant':str(det),
      'determinant_proof':'normalized pivot-row minor is triangular diagonal1; source-column change determinant is product inverse pivot normalizers and pivot-order permutation',
      'all_d1_d2_constant_products_zero':True,'actual_beta0_C2_over_B':3,'actual_beta1_C2_over_B':beta1,
      'source_minimalization_consequence':f'additional {len(piv)} unit pairs between G1/G2 reduce first relation rank115 to {beta1}',
      'actual_V_Hadamard_coefficients_not_computed':True,'no_T2_or_H0_length_computed':True}
    ep=out/'actual-d2-origin-rank-and-source-factorization.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    slots=sum(map(len,cols))+sum(map(len,B))+sum(map(len,factor))+len(piv)**2+sum(len(q['subtract']) for q in receipts)
    assert slots<150000
    report={'id':STEM,'all_passed':True,'actual_d2_origin_rank':len(piv),'actual_beta1_C2_over_B':beta1,
      'fixed_minor_nonzero':True,'all518_column_factorizations_exact':True,'stored_QQ_and_receipt_slots':slots,
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','actual_d2_origin_rank','actual_beta1_C2_over_B','stored_QQ_and_receipt_slots','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
