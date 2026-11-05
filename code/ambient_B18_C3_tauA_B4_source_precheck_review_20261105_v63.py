#!/usr/bin/env python3
"""Independent formal carry audit. Does not pair any cached source leaves."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_TAUA_B4_SOURCE_PRECHECK_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_TAUA_B4_SOURCE_20261105_v63'
PIN='30a20e0df9a353209fa7305b3d4ee2050b0710ca4d53679298e03d5772db8484'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def add(a,b,k=Q(1)):
    d=dict(a)
    for m,c in b.items():
        d[m]=d.get(m,Q(0))+k*c
        if not d[m]:del d[m]
    return d
def scaled(a,k):return {m:c*k for m,c in a.items() if c*k}
def main():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    p=ROOT/'reports'/PRI/'report.json';assert sha(p)==PIN
    report=read(p);pins={str(p.relative_to(ROOT)):PIN}
    for path,h in report['input_sha256'].items():assert sha(ROOT/path)==h;pins[path]=h
    for path,h in report['evidence_sha256'].items():
        q=ROOT/'reports'/PRI/path;assert sha(q)==h;pins[str(q.relative_to(ROOT))]=h
    e=read(ROOT/'reports'/PRI/'four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json')
    beta=[[Q(c) for c in row] for row in e['beta2_plane_r2_rell_ell2']]
    a2,b1,b2=[Q(e['curve_ratios'][k]) for k in ['a2','b1','b2']]
    active=[i for i,row in enumerate(beta) if any(row)];assert active==[0,4,5]
    # Independent linear algebra in formal delta30/delta21/delta12/zr_i/zell_i.
    d30={'delta30':Q(1)};d21={'delta21':Q(1)};d12={'delta12':Q(1)}
    carry=[{} for _ in range(3)]
    for i,row in enumerate(beta):
        for k in range(3):carry[k]=add(carry[k],{f'zr_{i}':row[k]/4})
        for k in range(2):carry[k+1]=add(carry[k+1],{f'zell_{i}':row[k]/4})
    B=[add(scaled(d,-Q(1,4)),c) for d,c in zip([d30,d21,d12],carry)]
    A2=add(d12,d30,-a2)
    chi=add(add(B[2],B[0],-b2),A2,Q(1,4))
    claimed=scaled(d30,(b2-a2)/4)
    for i,row in enumerate(beta):
        claimed=add(claimed,{f'zr_{i}':row[2]/4,f'zell_{i}':row[1]/4})
    defect=add(chi,claimed,-1)
    expected={f'zr_{i}':-b2*row[0]/4 for i,row in enumerate(beta) if row[0]}
    assert defect==expected and len(defect)==3 and all(defect.values())
    corrected=scaled(d30,(b2-a2)/4)
    for i,row in enumerate(beta):
        corrected=add(corrected,{f'zr_{i}':(row[2]-b2*row[0])/4,f'zell_{i}':row[1]/4})
    assert corrected==chi
    proof={'formal_cubic_B_coefficients':'B=-delta/4+carry; carry30=sum(delta_zi*r*beta_i0)/4',
      'literal_chi':'B_rell2+A_rell2/4',
      'normal_form_requires':'B_rell2 means B12-b2*B30, so the carry30 contribution must be retained',
      'actual_nonzero_beta_r2_indices':active,'actual_b2':str(b2),
      'formal_claimed_formula_defect':{m:str(c) for m,c in sorted(defect.items())},
      'corrected_zr_coefficients':{str(i):str((row[2]-b2*row[0])/4) for i,row in enumerate(beta) if any(row)},
      'unknown_source_leaves_remain_formal':True,'no_cached_leaf_products_or_dots_or_ranks_evaluated':True,
      'frozen_primary_unmodified':True,'GO_not_fulfilled_by_this_receipt':True,
      'scope':'The V1 formula is not an identity for the stated literal cover functional. This is independent of any later literal-zero evidence for the selected actual leaves.',
      'bounds':{'numeric_slots':2000,'hard_seconds':60,'hard_address_bytes':2<<30,'hard_RAW_bytes':32<<20}}
    ep=out/'missing-r3-carry-formal-certificate.json';save(ep,proof)
    result={'id':ID,'all_passed':False,'audit_blocker_reproduced':True,'source_values_or_rank_evaluated':False,
      'primary_report_sha256':PIN,'input_sha256':pins,'code_sha256':sha(__file__),
      'evidence_sha256':{ep.name:sha(ep)},'elapsed_seconds':time.monotonic()-started}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':False,'blocker_reproduced':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
