#!/usr/bin/env python3
"""Prove actualoldhead*r³zero beforeallowingcorrectedM1quarticrelations."""
import hashlib,json
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_QUARTIC_ZERO3_SOURCE_20261105_v63'
PATHS={
 'oldGB':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
 'oldGBpeer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
 'oldLoewy':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'first':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/actual-three-first-columns-and-fixed-quotient-proof.json',
 'firstpeer':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW/report.json',
 'quartic':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_PREFLIGHT_20261105_v63/tiny-five-relations-and-constructive-quartic-dual-source.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def main():
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(D[k]['all_passed']for k in ['oldGBpeer','oldLoewy','firstpeer'])
 gs=[{(t,tuple(e)):Q(c)for t,e,c in p}for p in D['oldGB']['homogeneous_module_generators40']];order=lambda k:(-k[0],)+k[1];lead=lambda p:max(p,key=order);proof=[]
 for head in[1,2,3]:
  initial={(head,(0,)*7+(3,)):Q(1)};p=dict(initial);trace=[]
  while p:
   k=lead(p);j=next((j for j,g in enumerate(gs)if lead(g)[0]==k[0]and all(a<=b for a,b in zip(lead(g)[1],k[1]))),None);assert j is not None,('OLD_GR3_NONZERO',head,p)
   g=gs[j];l=lead(g);shift=tuple(a-b for a,b in zip(k[1],l[1]));c=p[k]/g[l];trace.append([j,shift,c])
   for(t,e),v in g.items():key=(t,tuple(a+b for a,b in zip(e,shift)));p[key]=p.get(key,Q(0))-c*v
   p={k:c for k,c in p.items()if c};assert len(trace)<1000
  back={}
  for j,shift,c in trace:
   for(t,e),v in gs[j].items():key=(t,tuple(a+b for a,b in zip(e,shift)));back[key]=back.get(key,Q(0))+c*v
  back={k:c for k,c in back.items()if c};assert back==initial;proof.append({'oldhead':head,'source_free56_label':[22,23,29,30,55][head],'homogeneous_degree3_zero_reduction_trace':trace,'multiplyback_exact':True,'actual_zero_reason':'initialmembershipimpliesclassinm4Told; acceptedm4Told=0makesliteralmonomialzero'})
 F=D['first'];assert F['all576_actual_first_Koszul_checks']and F['actual_old_free56_labels']==[23,29,30]and F['old_incorrect_cache_and_rank_diagnostic_withdrawn']
 M=D['quartic']['source_actual_natural_origin_columns5'];assert all(all(Q(c)==0 for c in M[j])for j in[1,2,3])
 data={'three_actual_old_head_r3_zero_source_proofs':proof,'three_corrected_actual_M1_columns40':F['actual_three_minimal_first_columns40'],
  'source_oldhead_indices':[1,2,3],'source_free56_labels':[23,29,30],'all_three_M0_columns_zero':True,
  'genuine_quartic_relation_proof':'Actualoldhead*r³=0. M0col=0 exactly, sopushforwardhasfirstpossibledegree4=M1col*r³. EveryM>=2termisdegree>=5. UsecorrectedfullactualM1, notwithdrawncache.',
  'future_augmented_relation_shape':[3,8],'future_first5':'q3Pure*M0','future_last3':'constructivep4*(correctedM1_columns40 timesr³)',
  'augmented_matrix_rank_or_quartic_dual_not_computed':True,'new_physical_source_gauge_or_jet':False}
 ep=out/'three-actual-zero-heads-and-corrected-first-source-relations.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_oldactualzero3_to_genuinequarticrelations','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'actual_quartic_rank_dual_or_cube_computed':False};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
