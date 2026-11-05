#!/usr/bin/env python3
"""Publish Q1/J primary result and SOURCE-only independent d-Horner coefficient/budget recipe."""
from pathlib import Path
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_J_DHORNER_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 actual=load('S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ENTRY_20261105_v63');n=actual['result']
 assert actual['all_passed']and n['J_valuation']==0 and n['first_nonzero_terms']==606
 custody=load('S61_NONSQUARE_PREPOWER_Q1_J_HASSE_WORKER_CUSTODY_20261105_v63');assert custody['all_passed']
 old=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63');assert old['all_passed']
 literal=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt')
 assert sha(literal)==old['evidence_sha256'][literal.name]
 value=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ENTRY_20261105_v63/first-nonzero-J-Hasse.txt')
 assert sha(value)==actual['evidence_sha256'][value.name]
 zero=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ENTRY_20261105_v63/completed-zero-orders.txt');assert zero.read_bytes()==b''
 with value.open()as f:
  assert f.readline().strip()=='Q1JHV001';k,dd,am=map(int,f.readline().split());assert k==0 and am==18
  assert int(f.readline())==606
 stage=[sum(a+e<=68-2*l for a in range(35-l)for e in range(52-l))for l in range(19)]
 assert max(stage)==1667
 own=2*(4811+2*1667)+16000;assert own==32290
 proof=dict(SOURCE_ONLY_NO_DHORNER_VALUES=True,
  independently_reparse='original ROOT normalized4811 integer literal at D197; never producerimports/nativeengine or expected606 values until finalcompare',
  grouping='I_l(a,e), l0..18 by literal packedsource d exponent, 4811 totalcoefficients; not a numerical point',
  algorithm='H=I18; for l17..0: Hnew=(4+6e)*H+(5a)^(18-l)*I_l',
  proof='H_l=sum_{j>=l}I_j*(4+6e)^(j-l)*(5a)^(18-j); at l0 this is exactly numeric numeratorP0 of Q1(d0) with common D*5^18*a^18',
  no_binomial_Rpowers_or_forward_reverse_Hasse_engine=True,
  source_frame='no a^3/C4 redivision; all4811 body unchanged, source d exponents remain0..18',
  denominator='D197*5^18 times symbolic a^18, both must match independently',
  target='zeroorders are EMPTY since firstnonzerok0; compareALL606 integer coefficients and full denominator/externala18; require nonzero',
  H_l_support='a<=34-l,e<=51-l,a+e<=68-2l',all_stage_support_bounds=stage,
  max_H1667=1667,input4811=4811,
  lifetime='inputgroup4811 plus H and ONE temporaryH; freeH thenreplace; no expected606 map; compare only after releasing input',
  whole_owned32290=own,proposed_owned34376=34376,
  metadata16000='all19groupindex/pow5/scalar/exponentkey/cursor/sourcepins and comparison/publication scratch; coefficient and packedkey countedseparately',
  height='||H_l||1<=||I||1*max_{j>=l}10^(j-l)5^(18-j)<=||I||1*10^18<2^361; stage-awarepartials samebound; conservative beforecancel364<518',
  actual_value_scope_NOT_evaluated_by_this_source=True,stored_hard518=518,conservative_preop364=364,
  all_product_upper65000=65000,product_detail='2*18*1667+4811+18=64841; choose65k. Scalarshifts5^m represented by monomialkey, no sourcepolynomialpowers.',
  RAW_upper4000000=4000000,soft45_hard60_RAM2GB=True,
  no_Jquotient_GCD_newfactor_tail_or_badlocus=True)
 dump('independent-d-Horner-SOURCE.json',proof)
 dump('primary-result-binding.json',dict(J_valuation0_PRIMARY=True,firstlayer606=True,primary_sha256=sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ENTRY_20261105_v63/report.json'),
  independent_dHorner_peer_STILL_REQUIRED=True,no_otherfactor_or_quotient=True))
 pin(ROOT/'research/notes/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_RESULT_20261105_v63.md');pin(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,actual_primary_publication=True,
  dHorner_independent_values_NOT_computed=True,source_owned32290=32290,source_products65000=65000,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-t)
 dump('report.json',r);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
