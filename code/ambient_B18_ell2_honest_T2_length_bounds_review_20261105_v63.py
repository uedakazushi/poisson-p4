#!/usr/bin/env python3
"""Independent pure DVR corollary; no source or comparison coefficients."""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63';STEM=A+'_REVIEW'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 rp=ROOT/'reports'/A/'report.json';r=json.loads(rp.read_text());assert r['all_passed']
 inputs={'primary':{'path':str(rp.relative_to(ROOT)),'sha256':sha(rp)}};D={}
 for p,h in r['input_sha256'].items():
  assert sha(ROOT/p)==h;v=json.loads((ROOT/p).read_text());assert v['all_passed'];D[v['id']]=v;inputs[p]={'path':p,'sha256':h}
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 actual=json.loads((ROOT/'reports'/A/'formal-DVR-length-and-H0-criterion.json').read_text())
 K=D['AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63'];T=D['AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63'];C2=D['AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63'];F=D['AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63']
 assert K['rank_D_K']==K['rank_D_C']==4 and T['actual_T_length']==11 and T['actual_mu_D_torsion_K']==3
 assert F['actual_ell_regular_whole_E4_and_Ered'] and C2['actual_C2_dimension']==14
 assert C2['actual_ell_action_rank']==6 if 'actual_ell_action_rank' in C2 else True
 # Each finite cyclic D/(ell^a) has equal ell²-kernel and ell²-quotient
 # lengths min(2,a). The free summand contributes2 to the quotient only.
 kernelC=14-2*K['rank_D_C'];assert kernelC==6
 possibilities=[]
 for large in range(4):
  e=[2]*large+[1]*(3-large)
  lengthK2=2*K['rank_D_K']+sum(min(2,x) for x in e)
  possibilities.append([lengthK2+kernelC,large])
 assert possibilities==[[17,0],[18,1],[19,2],[20,3]]==actual['possible_length_and_count']
 assert actual['length_T2_lower']==17 and actual['length_T2_upper']==20 and not actual['actual_T2_length_evaluated']
 ev={'inputs':inputs,'exact_sequence_source':'original free comparison complex:0 -> K/ell²K -> T2 -> C[ell²] ->0','finite_DVR_structure_premises':{'rank_K':4,'rank_C':4,'torsion_generators_K':3,'C_ell2_quotient_length':14},'C_ell2_kernel_length':kernelC,'formal_T2_length_formula':'14+sum_i min(2,e_i)=17+#{e_i>=2}, all e_i>=1','all_possible_length_and_large_exponent_counts':possibilities,'H0_identification':'ell-regular E4 implies D0 ell-torsion=K ell-torsion; finiteDVR K makes this finite k-length and hence exactly H0_m(D0)','equivalent_minimal_conditions':['lengthT2=17','all3DVRtorsionexponents=1','ell*H0_m(D0)=0','lengthH0_m(D0)=3'],'not_inferred':['actualT2length','R0socledimension','maximalidealannihilation','nilcube','exactnilpotenceindex'],'new_coefficient_or_rank_job':False}
 ep=o/'independent-DVR-length-17-20-and-H0-minimum-proof.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'pure_DVR_T2_length17_to20_and_minimumH0_equivalences_accepted':True,'actual_T2_length_H0_exponents_or_nilcube_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 (o/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps(report))
if __name__=='__main__':main()
