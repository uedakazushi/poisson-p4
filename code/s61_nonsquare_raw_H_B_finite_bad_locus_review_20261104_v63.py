#!/usr/bin/env python3
"""Independent assembly of three fully audited finite covers; no coefficient job."""
import argparse,hashlib,json,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_FINITE_BAD_LOCUS_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hashes={};b=ROOT/'reports';pre='S61_NONSQUARE_RAW_H_B_';end='_20261104_v63'
    def audit(p):
        d=read(p);hashes[str(p.relative_to(ROOT))]=sha(p)
        for n,v in d.get('input_sha256',{}).items():assert sha(ROOT/n)==v,n;hashes[n]=v
        for n,v in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==v,n;hashes[str((p.parent/n).relative_to(ROOT))]=v
        return d
    primary=audit(b/(pre+'FINITE_BAD_LOCUS'+end)/'report.json');assert primary['actual_B_rankbad_finite_or_empty']
    gates=[pre+'ADDITIONAL_NORM_PENCIL'+end+'_REVIEW',pre+'CURVE_COMPONENTS'+end+'_REVIEW',pre+'Q_COMPONENT_FINITE'+end+'_REVIEW']
    gg=[audit(b/g/'report.json') for g in gates];assert all(g['all_passed'] for g in gg)
    assert gg[0]['source_monic_norm_pencil_and_residual_cover_accepted'] and gg[0]['finite_residual_geometric_point_bound']==1600
    residual=read(b/(pre+'PROJECTION_CURVE_RESIDUAL_INPUT'+end)/'certificate.json');factor=read(b/(pre+'COMMON_CURVE_FACTOR'+end)/'certificate.json')
    norm=read(b/(pre+'NORM_PENCIL'+end)/'native-certificate.json')
    asdict=lambda p:{tuple(m):F(c) for m,c in p if F(c)}
    assert asdict(residual['common_G'])==asdict(factor['source_G'])==asdict(norm['final_common_curve_gcd'])
    h1,k1=residual['oldh_quotient_h1'],residual['newK_quotient_k1']
    assert [max(sum(m) for m,c in x) for x in (h1,k1)]==[16,25]
    old=4*16*25;assert old==1600
    L=read(b/(pre+'COMPONENT_COMMON_ROOT_L'+end)/'certificate.json');Q=read(b/(pre+'COMPONENT_EXACT_QUOTIENT_Q'+end)/'certificate.json')
    A=read(b/(pre+'Q_ROOT_ALL_MINORS'+end)/'native-certificate.json')
    degree=lambda p:max(int(i) for i,c in p if F(c))
    le=L['finite_e_exception_polynomials'];qe=Q['finite_e_exception_polynomials']
    D=[[3,'486'],[2,'243'],[1,'-297'],[0,'-54']]
    first=A['all22_evaluations'][0]['remainder_quadratic_norm_e'];assert first and all(m[0]==m[2]==0 for m,c in first)
    ne=[[m[1],c] for m,c in first];assert degree(ne)==38
    assert sum(map(degree,le))==65 and sum(map(degree,qe))==320
    bounds=[old,4*65,2*4*(320+3+38)];assert bounds==[1600,260,2888] and sum(bounds)==4748
    assert primary['three_geometric_point_upper_bounds']==bounds and primary['geometric_point_upper_bound_total']==sum(bounds)
    c=read(b/(pre+'FINITE_BAD_LOCUS'+end)/'certificate.json');covers=c['all_B_bad_contained_in_three_covers'];assert len(covers)==3
    assert [r['geometric_point_bound'] for r in covers]==bounds
    assert covers[1]['e_polynomial_factors']==le and covers[2]['e_polynomial_factors']==qe+[D,ne]
    assert L['generic_common_roots_all_U_boundary'] and not Q['generic_common_roots_all_U_boundary'] and 0 in A['nonzero_remainder_indices']
    proof=dict(all_passed=True,independent_upstream_norm_pencil_and_coprime_residual_gate_closed=True,
        independent_L_and_Q_source_arithmetic_gates_closed=True,
        all_three_covers_refer_to_the_same_G49_and_actual_monic_H=True,
        entire_actual_B_rankbad_on_FULLU_finite_or_empty=True,
        geometric_point_upper_bounds=bounds,total_geometric_points_atmost=4748,
        old_finite_residual_and_all_L_Q_denominator_zeros_retained=True,
        U_root_zero_outside_field_exceptions_cannot_contribute_FULLU_bad_points=True,
        no_scheme_length_empty_global_rank3_raw_H7_p0_cubic57_quartic239_or_exhaustion_claim=True,
        no_new_GCD_coefficient_job_point_prime_or_GB=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md'):hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,actual_B_rankbad_finite_or_empty=True,three_geometric_point_upper_bounds=bounds,
        geometric_point_upper_bound_total=4748,scheme_length_or_empty_or_H7_full57_full239_claim=False,
        input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,report_sha256=sha(out/'report.json'),seconds=report['seconds'])))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
