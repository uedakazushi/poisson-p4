#!/usr/bin/env python3
"""Formal hereditary-DVR corollary; no source coefficients or ranks evaluated."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_HONEST_T2_DVR_TYPE_20261105_v63'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 refs=['reports/AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63/report.json','reports/AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json']
 assert all(json.loads((ROOT/p).read_text())['all_passed']for p in refs)
 types=[dict(K_torsion_exponents_ge2=j,A_free=6+j,k=5-j,total_length=17+j,ell_rank=6+j,D_minimal_generators=11)for j in range(4)]
 assert all(2*x['A_free']+x['k']==x['total_length']and x['A_free']+x['k']==11 for x in types)
 data={'base':'D=k[[ell]], A=D/(ell²)','noncanonical_D_module_split':'T2 = K/ell²K plus C[ell²]',
 'proof':'D is hereditary; the homology0,1 truncation splits because Ext_D²(C,K)=0. Higher homology cannot enter H1 after derived tensor because Tor-dimensionD(A)=1.',
 'C_ell_squared_torsion_type':{'A_free':2,'k':2},'possible_types':types,
 'actual_D_minimal_generators_T2':11,'actual_j_or_T2_length_evaluated':False,
 'not_claimed':['R0-module splitting','algebra or product splitting','natural T2 to Told equals ell quotient','actual torsion exponents','H0 length','nil cube']}
 ep=out/'formal-DVR-type-and-ell-rank-checks.json';ep.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'formal_corollary_only':True,'actual_D_minimal_generators_T2':11,
 'actual_T2_length_or_ell_rank_evaluated':False,'input_sha256':{p:sha(ROOT/p)for p in refs},
 'code_sha256':sha(__file__),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
