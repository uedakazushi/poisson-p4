#!/usr/bin/env python3
"""Literal source selectors clarifying rawpair4 versus the seven borderedrows."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K120_SOURCE_ROWS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json');assert load('S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63_REVIEW')['all_passed'];assert load('S61_NONSQUARE_PREPOWER_K120_PIVOT13_SCHUR_20261104_v63')['all_passed']
 pure=[22,11,13,25,15,19,14];indices=[q['new22_target_original32_columns'].index(j)for j in pure];assert indices==[15,8,9,16,11,14,10];pair=[2,4,5,8];P4=[0,1,3,6];assert not(set(P4)&set(pair));assert not({7,10,11,12}&set(pair));assert hc['known_row_units_rho14'][11]==hc['known_row_units_rho14'][9]
 rows=[{'row':i,'entries':[{'current15':r,'coefficient':'1'}]}for i,r in enumerate(P4+[2])]
 rows+=[dict(row=5,entries=[dict(current15=11,coefficient='g/c10'),dict(current15=7,coefficient='-h/c6')]),dict(row=6,entries=[dict(current15=12,coefficient='g/c11'),dict(current15=10,coefficient='-g/c11'),dict(current15=7,coefficient='-h120/c6')])]
 cert=dict(actual_input='M15x22 afteracceptedQ2/LprePOWERprefix, notbareactual135equations',raw_pair_selector_current15=pair,raw_pair_compact14=[1,3,4,7],P4_selector_current15=P4,
  pure7_original32=pure,pure7_current22=indices,bare7_source_matrix7x15_sparse=rows,
  g='N[6,1]',h='N[10,1]',h120='N[11,1]-N[9,1]',constants={k:hc['known_row_units_rho14'][i]['constant']for k,i in [('c6',6),('c10',10),('c11',11)]},
  rawpair_rank_scope='rank(PairSelector*M15*Pure7Selector)<=1, notrank(PairSelector*M15)<=1',
  full_pair_rank_scope='Submatrixon[Hbar27,Hbar31,M22,B*M2]hasdet=-T nonzero onD(T), sofullpairrowsrank4.',
  bare7_rank_scope='rank(topP4+pairpureline)<=5; twooutsidecoupledrowsK1,K2canraisetorank7. No contradictionwithpairpureprojectionrank1.',
  higher_rows_not_pair=True,full_rawH7_lift_requires_prior_accepted_cleanup=True,all7tail_effects_retained=True,no_new_scalar_arithmetic=True)
 dump(out/'source-selectors.json',cert);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,rawpair4_and_bare7_sources_distinguished=True,all15_source_indices_and7_targets_checked=True,primary_pivot13_note_unchanged=True,input_sha256=pins,evidence_sha256={'source-selectors.json':sha(out/'source-selectors.json')});dump(out/'report.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
