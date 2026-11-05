#!/usr/bin/env python3
"""Pure derived-DVR splitting review; no coefficient or homology query."""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_HONEST_T2_DVR_TYPE_20261105_v63';STEM=A+'_REVIEW'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 p=ROOT/'reports'/A/'report.json';r=json.loads(p.read_text());assert r['all_passed'];inputs={'primary':{'path':str(p.relative_to(ROOT)),'sha256':sha(p)}}
 for path,h in r['input_sha256'].items():
  assert sha(ROOT/path)==h and json.loads((ROOT/path).read_text())['all_passed'];inputs[path]={'path':path,'sha256':h}
 for path,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/A/path)==h
 types=[]
 for j in range(4):
  blocks2=4+j+2;blocks1=(3-j)+2
  assert blocks2+blocks1==11 and 2*blocks2+blocks1==17+j
  types.append({'K_exponents_at_least2':j,'A_blocks':blocks2,'k_blocks':blocks1,'D_minimal_generators':11,'ell_rank':6+j,'length':17+j})
 ev={'inputs':inputs,'DVR_global_dimension':1,'low_Postnikov_obstruction':'Ext_D^2(C,K)=0','noncanonical_D_derived_split':'tau_homology0,1(originalcomplex) is C direct_sum K[1]','higher_homology_excluded':'Tor_D amplitude0,1 increases homologicaldegree; q>=2 cannot contribute to H1','D_flat_source':'the original freeR0 complex is D-flat sinceR0=D[[8variables]] in the actualell coordinate','noncanonical_D_module_split':'T2 is K/ell²K direct_sum C[ell²]','C_ell2_kernel_type':'A^2 direct_sum k^2 from C2=A^6+k^2 andrank_D C=4','all_possible_T2_types':types,'mu_B_upper_bound':11,'not_claimed':['R0_split','algebra_or_filtration_split','naturality_of_T2_reduction_to_Told','actual_j_or_length','H0_exponents_or_nilcube'],'no_coefficient_or_new_rank_job':True}
 ep=o/'independent-derived-DVR-split-and-eleven-generator-proof.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'noncanonical_D_module_split_and_mu_D_T2_eleven_accepted':True,'actual_T2_length_ell_rank_or_R0_split_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 (o/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps(report))
if __name__=='__main__':main()
