#!/usr/bin/env python3
"""Known-array factorization into eight formal psi contractions; no values."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 d=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE_PROFILE_20261105_v63/actual-four-active-source-complements-and-unknown-derivative-profile.json')
 g=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json')
 assert g['remaining_top_linear_quotient_coordinates']==[25,26]
 Pi=list(map(H.parse,g['Pi_remaining_top_linear_columns']));assert len(Pi)==27
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json')
 nf=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json')
 memo={i:H.parse(v) for i,v in nf['all264_linear_NF_columns']};mi={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 def monnf(target,u,v):
  if u==v==8:return {}
  m=[0]*9;m[u]+=1;m[v]+=1;return memo[lin['all264_NF_root_nodes'][target][mi[tuple(m)]]]
 Q=[monnf(2,7,7),monnf(2,7,8)];weights=[];verified=0
 for rec in d['remaining_four_source_directions']:
  b=H.parse(rec['unique_selected99_linear_correction']);a={};c={}
  for j in range(99):
   r=b.get(j*9+7,F(0));e=b.get(j*9+8,F(0))
   if r:a[j*2]=r;c[j*2+1]=r
   if e:c[j*2]=e
  rebuild={}
  for q,w in zip(Q,[a,c]):
   for slot,x in w.items():
    j,psi=divmod(slot,2)
    for u,p in enumerate(Pi):
     y=p.get(psi,F(0))
     if y:
      for output,z in q.items():rebuild=H.guard(add(rebuild,{output*99*3*9+j*27+u:-x*y*z}))
  expected=H.parse(rec['unknown_post27_PiY_Ystar1_sensitivity6_by99by3by9']);assert rebuild==expected
  verified+=6*99*3*9
  weights.append({'source_extra_index':rec['source_extra_index'],'formal_A_weights_on99psi2':srec(a),'formal_B_weights_on99psi2':srec(c),
   'all16038_formal_derivative_slots_verified':True})
 data={'quotient_free_coordinates':[25,26],'free_monomials':['top2*rA','top2*ell'],
  'NF_top2_r_squared':srec(Q[0]),'NF_top2_r_ell':srec(Q[1]),'all_four_pairs_of_formal_psi_contractions':weights,
  'formal_derivative_slots_checked':verified,'selected99_source741_indices':d['selected99_unit_source741_indices'],
  'psi_definition':'psi=PiTop2(PiY0*ACTUALpost27Ystar1) inSAMEcoherentnormalizedcomparisonframeasL48',
  'factor_formula':'unknownquadraticNF_i= -NF(top2*r²)*A_i -NF(top2*r*ell)*B_i; A=sum b1r*psi0, B=sum(b1ell*psi0+b1r*psi1)',
  'post27_source_identity':'G2linear*Ystar1_j=M0*(Xj2-Xp1*a1j)+M1*(Xj1-Xp0*a1j)-G2quadratic*Yj0',
  'identity_scope':'a1j isfirst27 bottomunitcoordinate ofknownXj1; Xj0=0; unknowna2 dropsM0Xp0=0; Yp0carryabsorbedusingG2linearYp0=M0Xp1+M1Xp0. No newcoefficient values are asserted.',
  'next_source_dependencies':['eightweightedrawX2rell contractions; honestsecondsensmustbecheckedforNEWweights',
   'M1onknownXeff1 weightedsourcevectors; no commonWv assumption','contractedpureknownC2G2quadratic response'],
  'frame_caution':'M1/psi mustremaininthenormalizedcoherentcomparisonframe; canonicalother2muunknownandtypednormalizerA0recipecannotbedropped',
  'unknown_values_evaluated':False,'new_rank_search':False,'actual_T2_HF2_or_length_not_computed':True}
 ep=out/'actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json';H.dump(ep,data)
 slots=sum(len(q) for q in Q)+sum(len(r['formal_A_weights_on99psi2'])+len(r['formal_B_weights_on99psi2']) for r in weights)+25000;assert slots<60000
 report={'id':STEM,'all_passed':True,'formal_scalar_contraction_count':8,'all_formal_derivative_slots_verified':verified,
  'unknown_values_or_new_rank_computed':False,'owned_QQ_slot_bound':slots,'input_sha256':hs,
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','formal_scalar_contraction_count','all_formal_derivative_slots_verified','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
