#!/usr/bin/env python3
"""Parent proof audit and producer rerun of relative C3 stratum removal."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_binary_mixed_stratum_removal_20261101 import verify
ROOT=Path(__file__).resolve().parents[1]
def review():
 start=time.time();path=ROOT/'reports/S61_BINARY_MIXED_STRATUM_REMOVAL_20261101.json';saved=json.loads(path.read_text());assert saved['all_passed'];fresh=verify()
 for key in ['vertex_jacobian','vertex_minor','actual_binary_family_tangent','binary_base_minor','smooth_beta_base_change_minor','theorem_conclusions']:assert fresh[key]==saved[key]
 assert s.Rational(saved['vertex_minor'])==s.Rational(6244,5) and s.Rational(saved['binary_base_minor'])==s.Rational(99,3122)
 paths=[Path(__file__),path,ROOT/'code/s61_binary_mixed_stratum_removal_20261101.py',ROOT/'research/notes/S61_BINARY_MIXED_STRATUM_REMOVAL_20261101.md',ROOT/'reports/S61_BINARY_MIXED_FORMAL_REVIEW_20261101.json',ROOT/'reports/S61_BINARY_MIXED_RELATIVE_DOUBLE_LINE_20261101.json']
 return dict(id='S61_BINARY_MIXED_STRATUM_REMOVAL_REVIEW_20261101',all_passed=True,verification_kind='independent proof audit plus exact small producer rerun',numeric_reconstruction_independent=False,proof_audit=['Equal-rank torsion-free specialization needs no flatness assumption: for one parameter the torsion quotient is regular in that parameter, and localization at the height-one prime fixes the reduced rank; tensor exactness and Nakayama kill the torsion. Induct using presentation minors for several parameters.','The relative binary vertex equations are smooth of relative dimension1, and two quadratic-form coefficients give the same smooth relative3fold with the old three base coordinates.','Q0 and Q1 lie in one projective quadratic-form fiber, justifying application of the universal generic lower65 to the SAME relative component.','The actual relative double line is A-flat, so the kernel N specializes exactly to the reviewed torsion-free N0.','Special finite length<=67 gives finite generation overR12; specialization rank<=65 and actual binary lower>=65 force rank65.','The torsion-free specialization lemma applies to N, not to the nonfree whole S. Generic t^2=0 and t=0 onD give t^2 inN and then full relative zero.','The seven-direction beta minor makes the Sp4 base change smooth and faithfully flat; full135 equivariance descends t^2 to the ambient local ring.','The irreducible stratum has dimension16 and dominatesW16; its proper closed complement has dimension<=15, so the whole sufficiently general beta fiber is removed.'],conclusions=saved['theorem_conclusions'],limitations=['Small producer computations were rerun, not independently implemented.','No statement that t^2=0 at every special point of the stratum or that S/N is free65/67.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=review()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','seconds']}))
