#!/usr/bin/env python3
"""Formal DVR length corollary; no comparison coefficient is evaluated."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 sources={
 'Told':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'ranks':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json',
 'regular':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/report.json',
 'regular_peer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/report.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'T2':'reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63/report.json',
 'T2_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json'}
 D={k:json.loads((ROOT/p).read_text())for k,p in sources.items()};assert all(x['all_passed']for x in D.values())
 assert D['Told']['actual_T_length']==11 and D['Told']['actual_mu_D_torsion_K']==3
 assert D['ranks']['rank_D_K']==4 and D['ranks']['rank_D_C']==4
 assert D['regular']['actual_ell_regular_whole_E4_and_Ered']
 assert D['C2']['actual_C2_dimension']==14
 torsC2=14-2*4;assert torsC2==6
 result={'base':'D=k[[ell]], actual same physical parameter ell',
 'actual_exact_sequence':'0 -> K/ell^2 K -> T2 -> C[ell^2] -> 0, from the original free R0 comparison complex',
 'K_DVR_structure':'D^4 plus three D/(ell^e_i), all e_i>=1; module decomposition only',
 'C_ell2_torsion_length':torsC2,
 'formula':'length(T2)=14+sum_i min(2,e_i)=17+# {i:e_i>=2}',
 'length_T2_lower':17,'length_T2_upper':20,
 'possible_length_and_count':[[17+j,j]for j in range(4)],
 'actual_T2_length_evaluated':False,
 'H0_identification':'H0_m(D0)=D-torsion(K): E4 is ell-regular, K is D-finite, so ell-torsion of D0 equals the finite-length torsion of K',
 'H0_length_formula':'length H0_m(D0)=sum_i e_i >=3',
 'equivalences':['length(T2)=17','all three e_i=1','ell*H0_m(D0)=0','length H0_m(D0)=3'],
 'exact_exponents_or_H0_length_decided':False,
 'not_implied':'ell annihilation does not imply m annihilation; no socle dimension, explicit nil generators, nil cube, or original nil index follows',
 'no_new_source_jet_gauge_rank_or_matrix_job':True}
 p=out/'formal-DVR-length-and-H0-criterion.json';p.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'formal_corollary_only':True,'actual_T2_length_interval':[17,20],
 'actual_T2_length_or_H0_evaluated':False,'input_sha256':{p:sha(ROOT/p)for p in sources.values()},
 'code_sha256':sha(__file__),'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
