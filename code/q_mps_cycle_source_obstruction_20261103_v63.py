#!/usr/bin/env python3
"""All-label exact cofactors for the T39/T40 extended toric-arc source.

This restricts this particular source, not membership in an entire component.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.rings import ring

from q_mps_toric_arc_memberships_20261103_v63 import PAIRS4,PERMS,extend,dual

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    R,z,x,y=ring('z,x,y',sp.QQ,order='grlex');ds=[R.one,z,x,y]
    pairs=list(itertools.combinations(range(4),2));V=[[d**k for k in range(1,5)]for d in ds]
    Delta=z*x*y*(z-1)*(x-1)*(y-1)*(x-z)*(y-z)*(y-x)
    source=ROOT/'code/data/q_mps_cycle_source_cofactors_20261103_v63.json'
    saved=json.loads(source.read_text());assert len(saved)==24
    catalogue=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json'
    records=json.loads(catalogue.read_text())['diagrams']
    certificates=[]
    for name,a0 in [('T39',sp.Rational(2,9)),('T40',-sp.Rational(8,9))]:
        rec=next(q for q in records if q['id']==name)
        b=list(map(sp.Rational,rec['period_basis'][0]));den=sp.ilcm(*[q.q for q in b])
        B=extend([int(q*den)for q in b]);forms={}
        for perm in PERMS:
            period=[int(B[perm[i],perm[j]])for i,j in PAIRS4]
            first=next(q for q in period if q);key=tuple(str(sp.Rational(q,first))for q in period)
            forms.setdefault(key,[]).append(list(perm))
        assert len(forms)==12 and sum(map(len,forms.values()))==120
        data={tuple(item['normalized_target_period']):item for item in saved if item['component']==name}
        assert set(data)==set(forms)
        for key,permutations in forms.items():
            A=dual(list(map(sp.Rational,key)));H={}
            for i,j in pairs:
                cp=[k for k in range(4)if k not in[i,j]];num=R.zero
                for(k,l),aa in zip(pairs,A):
                    cr=[t for t in range(4)if t not in[k,l]]
                    num+=aa*(-1)**(i+j+k+l)*(V[cr[0]][cp[0]]*V[cr[1]][cp[1]]-V[cr[0]][cp[1]]*V[cr[1]][cp[0]])
                H[i,j]=num
            decode=lambda rows:R.from_dict({tuple(e):sp.QQ.from_sympy(sp.Rational(c))for e,c in rows})
            item=data[key];assert all(H[pair]==decode(item['power_form_numerators'][str(pair)])for pair in pairs)
            # This is the full complementary-minor inverse formula, checked
            # without rational-function expansion or a sampled evaluation.
            HM=[[R.zero for _ in range(4)]for _ in range(4)]
            for(i,j),value in H.items():HM[i][j]=value;HM[j][i]=-value
            for i,j in pairs:
                rhs=sum((V[i][k]*HM[k][l]*V[j][l]for k in range(4)for l in range(4)),R.zero)
                assert rhs==Delta*A[pairs.index((i,j))]
            f,g=H[1,3],H[2,3]
            target=R(a0*9)*(H[0,3]**2-H[1,2]**2)-8*H[0,3]*H[1,2]
            u,v=decode(item['cofactor_u']),decode(item['cofactor_v'])
            assert target==u*f+v*g
            certificates.append(dict(component=name,normalized_target_period=list(key),
                all_label_permutations=permutations,source_required_special_a=str(a0),
                complementary_minor_identity_exact=True,polynomial_cofactor_identity_exact=True,
                cofactor_term_counts=[len(u),len(v)],
                identity='9*a0*(H14^2-H23^2)-8*H14*H23=U*H24+V*H34',
                scope='Necessary for this extended arc source with distinct eigenvalues and rho*(rho^2-1)!=0; not full-component nonmembership.'))
    (args.output/'all-label-certificates.json').write_text(json.dumps(certificates,indent=2)+'\n')
    inputs=[Path(__file__),source,catalogue,ROOT/'code/q_mps_toric_arc_memberships_20261103_v63.py',
        ROOT/'research/notes/Q_MPS_CYCLE_SOURCE_OBSTRUCTION_20261103_v63.md']
    report=dict(id='Q_MPS_CYCLE_SOURCE_OBSTRUCTION_20261103_v63',all_passed=True,
        targets=['T39','T40'],permutations_each=120,projective_period_forms_each=12,
        exact_QQ_complementary_minor_and_cofactor_identities=24,
        required_source_a={'T39':'2/9','T40':'-8/9'},
        full_component_nonmembership_claimed=False,a_zero_nonmembership_claimed=False,
        source_obstruction_not_global_exhaustion=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS all240label source cofactor identities; nofull-componentnonmembership;',report['seconds'])


if __name__=='__main__':main()
