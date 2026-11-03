#!/usr/bin/env python3
"""Independent all-label polynomial certificates for a limited arc source."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha
from q_mps_arc_membership_review_20261103_v63 import dual

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    folder=ROOT/'reports/Q_MPS_CYCLE_SOURCE_OBSTRUCTION_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'};paths=[Path(__file__),public,ROOT/'code/q_smith_global_model_review_20261103_v63.py',ROOT/'code/q_mps_arc_membership_review_20261103_v63.py']
    for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
    for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h and p.read_bytes()==(args.replay/name).read_bytes();paths.append(p)
    saved=read(ROOT/'code/data/q_mps_cycle_source_cofactors_20261103_v63.json');records=read(ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json')['diagrams'];R,z,x,y=ring('z,x,y',s.QQ);eigen=[R.one,z,x,y];V=[[v**k for k in range(1,5)]for v in eigen];Delta=z*x*y*(z-1)*(x-1)*(y-1)*(x-z)*(y-z)*(y-x);pairs=list(itertools.combinations(range(4),2));labels=0
    decode=lambda data:R.from_dict({tuple(exp):s.QQ.from_sympy(s.Rational(c))for exp,c in data})
    for name,a0 in[('T39',s.Rational(2,9)),('T40',-s.Rational(8,9))]:
        record=next(v for v in records if v['id']==name);assert record['period_dimension']==1 and len(record['period_basis'])==1
        B=s.zeros(5)
        for(i,j),v in zip(itertools.combinations(range(1,5),2),record['period_basis'][0]):B[i,j]=s.Rational(v);B[j,i]=-B[i,j]
        for i in range(1,5):B[0,i]=sum(B[i,j]for j in range(1,5));B[i,0]=-B[0,i]
        forms={}
        for perm in itertools.permutations(range(5)):
            period=[B[perm[i],perm[j]]for i,j in itertools.combinations(range(1,5),2)];first=next(v for v in period if v);key=tuple(str(v/first)for v in period);forms.setdefault(key,0);forms[key]+=1
        cases={tuple(v['normalized_target_period']):v for v in saved if v['component']==name};assert len(forms)==len(cases)==12 and set(forms)==set(cases)and sum(forms.values())==120;labels+=120
        for key,case in cases.items():
            A=dual(list(map(s.Rational,key)));H=[[R.zero for j in range(4)]for i in range(4)]
            for i,j in pairs:H[i][j]=decode(case['power_form_numerators'][str((i,j))]);H[j][i]=-H[i][j]
            # This inverse identity uniquely determines the full power form on
            # Delta!=0, independent of the producer's complementary-minor code.
            for n,(i,j)in enumerate(pairs):assert sum((V[i][k]*H[k][l]*V[j][l]for k in range(4)for l in range(4)),R.zero)==Delta*R(A[n])
            U=decode(case['cofactor_u']);W=decode(case['cofactor_v']);target=R(9*a0)*(H[0][3]**2-H[1][2]**2)-8*H[0][3]*H[1][2]
            assert target==U*H[1][3]+W*H[2][3]
    assert labels==240
    paths += [ROOT/'reports/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_MPS_CYCLE_SOURCE_OBSTRUCTION_20261103_v63_REVIEW.md']
    report=dict(id='Q_MPS_CYCLE_SOURCE_OBSTRUCTION_20261103_v63_REVIEW',all_passed=True,
        fresh_producer_and_all_evidence_match=True,independent_all240_label_partition_into24_period_lines=True,
        independent_all24_V_H_Vtranspose_equals_Delta_A_polynomial_identities=True,
        independent_all24_explicit_QQ_obstruction_cofactor_identities=True,
        necessary_a_on_this_admissible_extended_toric_arc_source={'T39':'2/9','T40':'-8/9'},
        full_component_nonmembership_a_zero_exclusion_or_other_arc_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print('PASS all240 labels/24 polynomial cofactor certificates; limited arc source ONLY')

if __name__=='__main__':main()
