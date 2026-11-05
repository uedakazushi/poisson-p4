#!/usr/bin/env python3
"""Companion: literal two-source minor in the full no-cut cubic cover."""
import hashlib,json,pathlib
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def mul(a,b):
    c={}
    for i,x in a.items():
        for j,y in b.items():c[i+j]=c.get(i+j,Q(0))+x*y
    return {i:x for i,x in c.items() if x}
def main():
    r=ROOT/'reports'/ID;p=r/'report.json';assert sha(p)=='fe2efb42fcff3b90486c8f1a61b1f6dbe41d40aaf5313e1b3b687fd189fe4871'
    d=json.loads((r/'independent-full-C-finiteness-Nakayama-and-conditional-graph03-proof.json').read_text())
    N=[[Q(x) for x in row] for row in d['full_no_cut_cubic_NF']]
    A=[row[3] for row in N];B=[row[7] for row in N]
    assert A==[Q(0),Q(0),Q(1),Q(0),Q(0),Q(0)]
    assert B==[Q(0),Q(0),Q(0),Q(0),Q(0),Q(1)]
    ca=[{1:x} if x else {} for x in A];cb=[{1:x} if x else {} for x in B]
    first=mul(ca[2],cb[5]);second=mul(ca[5],cb[2])
    minor=dict(first)
    for i,x in second.items():minor[i]=minor.get(i,Q(0))-x
    minor={i:x for i,x in minor.items() if x};assert minor=={2:Q(1)}
    result={'all_passed':True,'primary_review_report_sha256':sha(p),'full_no_cut_NF_ell3_columns_are_literal_e2_e5':True,
       'formal_graph_sources':'d03*N(eA*ell3),d03*N(eB*ell3)',
       'selected_rows':[2,5],'ordinary_QQ_d_minor':[[i,str(x)] for i,x in minor.items()],
       'd03_value_unread':True,'code_sha256':sha(__file__)}
    target=r/'universal-graph03-minor-companion.json';assert not target.exists();target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'companion_sha256':sha(target)}))
if __name__=='__main__':main()
