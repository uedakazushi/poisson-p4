#!/usr/bin/env python3
"""Exact arithmetic supplement to the normal sextic geometric exclusion.

No source search or CAS. Geometric descent, CM fibers, relative conic
embedding and rational-singularity arguments are hand proofs, not certified
by these finite polynomial/integer checks.
"""
import argparse,hashlib,json,time
from pathlib import Path
from classification_cubic_scroll_abelian_review import add,scale,parse,Q

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    variables=['a','b','g','s1','s2']
    canonical=parse('a*(2*g-8)+6+s2-2*a**2*(b-1)-a*s1',variables)
    nef_identity=parse('6-a*(8-2*g)-2*a**2*(b-1)-(a*s1-s2)',variables)
    assert canonical==nef_identity
    # Since g<=3,b>=1 and every m(a-m)>=0, 6>=2a, hence a<=3.
    candidates=[]
    for g in range(4):
        for b in range(1,4):
            for a in range(2,4):
                if not(b<=g and b>=3-g and g-1>=a*(b-1)):continue
                contribution=6-a*(8-2*g)-2*a*a*(b-1)
                if contribution>=0:candidates.append(dict(g=g,b=b,a=a,sum_m_a_minus_m=contribution))
    assert candidates==[dict(g=3,b=1,a=2,sum_m_a_minus_m=2),dict(g=3,b=1,a=3,sum_m_a_minus_m=0)]
    terms={a:{str(m):m*(a-m)for m in range(a+1)}for a in (2,3)}
    assert terms[2]=={'0':0,'1':1,'2':0}and terms[3]=={'0':0,'1':2,'2':2,'3':0}
    # Projective-bundle intersection identities use quotient P(E) convention.
    names=['e','d'];adjunction=parse('-(2*e+d)+2*(e+d)',names);assert adjunction==parse('d',names)
    e,d=Q(4),Q(-2);assert 2*e+d==6 and d==-2 and 5-e==1
    notes=['research/classification/RANK_ZERO_NORMAL_SEXTIC_REDUCTION.md','research/classification/RANK_ZERO_NORMAL_SEXTIC_SEPARATE_REVIEW.md','research/classification/NORMAL_SEXTIC_RULED_NUMERICS.md','research/classification/NORMAL_SEXTIC_RANK_ZERO_EXCLUSION.md','research/classification/RANK_ZERO_SMOOTH_NONPLANE_SEPARATE_REVIEW.md']
    paths=[Path(__file__),ROOT/'code/classification_cubic_scroll_abelian_review.py',*[ROOT/p for p in notes]]
    result=dict(task_id='G',all_passed=True,certificate_scope='Exact polynomial/integer supplement only. Geometric hypotheses and all new closure/descent/conic arguments are separate hand proofs.',
        universal_ruled_intersection_identity_verified=True,admissible_integer_triples_after_nef_and_Hurwitz=candidates,
        all_possible_single_multiplicity_contributions=terms,a3_pattern='All m_i=0 or3; every vertical L-degree is divisible by3.',a2_pattern='Exactly two m_i=1; every other m_i=0 or2.',
        relative_conic_intersection_identity='K_S.H=-(2e+d)+2(e+d)=d',relative_conic_solution=dict(degree_E='4',degree_D='-2',h0_E=5,h1_E=1),
        geometric_proof_status='New producer proof; independent root review required before global promotion.',geometric_conclusion_in_note='All normal integral nondegenerate degree6 rankzero sources with irreducible reduced quintic D are empty.',
        not_machine_verified=['Ruling descends from rational resolution','Flat pure-CM genus0 fibers','a3 positive-component argument','Every degree2 fiber and the relative embedding are plane conics','Serre-dual quotient and the negative-degree defining equation section'],
        proof_snapshots={p:(ROOT/p).read_text()for p in notes[-3:-1]},source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6))
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))


if __name__=='__main__':main()
