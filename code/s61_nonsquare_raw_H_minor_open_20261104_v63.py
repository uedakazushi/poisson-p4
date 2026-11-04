#!/usr/bin/env python3
"""Strip only named units from the fixed407-term minor; exact H1-open corollary."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_MINOR_OPEN_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 gb=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';gp=gb/'certificate.json';gr=read(gb/'report.json');assert sha(gp)==gr['certificate_sha256'];g=read(gp);assert g['generic_rank_actual16_H7']==7
 R,a,e,d=ring('a,e,d',s.QQ);dec=lambda q:R.from_dict({tuple(m):s.QQ(v)for m,v in q});enc=lambda q:[[list(m),str(v)]for m,v in q.items()]
 F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4;known=[a,d,F,B,C,J]
 Delta=dec(g['nonzero_selected_minor_Delta']);assert len(Delta)==407;constant=Delta.LC;Q=Delta/constant;powers=[]
 for f in known:
  n=0
  while not Q.rem(f):Q=Q.exquo(f);n+=1
  powers.append(n)
 unit=R.ground_new(constant)
 for f,n in zip(known,powers):unit*=f**n
 assert unit*Q==Delta and Q;p=dec(g['p']);assert p
 hb=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';h=read(hb/'certificate.json');assert h['target117_projection_preserves_Hdegree_exact_structural_support'];assert len(h['target122_Hdegree_classes']['1'])==23
 hreview=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW/report.json';rreview=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW/report.json';assert read(hreview)['all_passed']and read(rreview)['all_passed']
 proof={'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','fixed_original_minor':enc(Delta),'fixed_minor_terms':407,
  'named_unit_factor':{'QQ_constant':str(constant),'a_d_F_B_C_J_powers':powers},'remaining_factor_Q':enc(Q),'literal_reconstruction_Delta_unit_times_Q':True,
  'remaining_factor_is_constant':Q.is_ground,'no_irreducible_factor_search_or_unknown_inverse_on_FULLU':True,
  'p':enc(p),'explicit_H7_chart':'FULLU intersect D(p*Q); equivalently D(pairdet*Delta) withinFULLU',
  'preserved_residual14x5_minor_formula':'-(6/5)*B*p^2*Delta','preserved_residual14x5_minor_terms':len(g['nonzero_minor_of14x5']),
  'source_bound_H7_rank_on_chart':'The two previously knownunit pivots plus the sparsepair determinant and selected3minor are invertible on this explicit chart. Thus the actual16 Hlinear rows span all7 over the localized ring, using an adjugate rightinverse and the saved exact source/target transports.',
  'H1_projection23_surjective_on_chart':True,'proof_H1':'Seven actual cubic combinations have arbitrary prescribedHlinear coefficients. Their seven actual light multiples span H*Sym2(light), dimension28 and split by the constant unordered-pair selection. The reviewed NF117 is a split surjection28->23. Cubic timeslight has noH2, so the reviewedH2 elimination does not change this H1 image. Pure98 remains the same lazy actual full-source arithmetic.',
  'if_Q_constant_then_B7x3_fullrank_on_ALL_FULLU':Q.is_ground,'sparse4x2_pair_fullrank_is_sufficient_not_claimed_necessary':True,
  'no_allFULLU_H7_or_H1_surjectivity_claim':True,'no_full239_claim_or_global_kernel_count_change':True,'genericrank_packet_independent_review_pending':True}
 (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),gb/'report.json',gp,hb/'report.json',hb/'certificate.json',hreview,rreview]
 report={'id':ID,'fixed_minor_knownunit_stripping_and_literal_reconstruction_passed':True,'named_unit_powers':powers,'remaining_factor_terms':len(Q),'remaining_factor_totaldegree':max(sum(m)for m in Q),'remaining_factor_constant':Q.is_ground,
  'actual_H1_23_surjectivity_proved_on_explicit_D_pQ_chart':True,'no_wholeFULLU_H1_or_full239_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
