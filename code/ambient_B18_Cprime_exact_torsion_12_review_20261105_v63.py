#!/usr/bin/env python3
"""Dependency-only independent module-theorem review; no value computation."""
import hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_EXACT_TORSION_12_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 st=time.monotonic();OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 primary='AMBIENT_B18_CPRIME_EXACT_TORSION_12_20261105_v63';r=load(primary,h='ecda032b8c247df78de85cad3793b1a47dd99a9067660623997f6f498c05104f')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 p=load(primary,'exact-exponents-full-Hilbert-and-genuine-target-corollary.json')
 right=load('AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63_REVIEW',h='3f4190e9cccc0ac8eafded0e18cbc6c0ce9f86f67c50d80b669df7770dd75301');assert right['positive_long_exponent2_proved']and right['actual_Cprime_ell3_length']==9 and right['actual_Cprime_DVR_torsion_exponents']==[1,2]
 quadratic=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW');assert quadratic['actual_Cprime_HF_through2']==[3,3,3]
 assert r['Cprime_DVR_torsion_exponents']==[1,2]and r['Cprime_torsion_length']==3 and r['Cprime_full_HF_prefix']==[3,3,3,2,2,2]and r['genuine_Cprime_ell_torsion_dimension']==2
 proof=dict(hypotheses='Same actual finite R0-module M=Cprime; D=QQ[[ell]], ell in maximalideal m, rank_D M=2, D-torsion exponents1,2; actualHF0..2=(3,3,3); m3(M/ell3M)=0.',no_ell_cancellation='m3M⊂ell3M fromquotient, ell3M⊂m3M becauseell∈m, soequal. Then ell3M⊂ell*m2M⊂m3M=ell3M givesm3M=ell*m2M withoutdividingell.',all_higher_powers='Multiplying m3M=ell*m2M bym^(j−2) givesm^(j+1)M=ell*m^jM forj>=2. Base m3M=ell3M andinduction yieldm^jM=ell^jM forj>=3.',tail_dimension='AllD-torsion diesundermultiplicationell². Thereforeell^jM isD-free rank2 forj>=2, andell^jM/ell^(j+1)M hasdimension2. Actualm-adicHF=(3,3,3,2,2,...).',Hilbert_series='(3−t³)/(1−t)',noncanonical_D_type='D²⊕D/(ell)⊕D/(ell²); notR0split andnotalgebra split.',torsion_scope='Length3 is D-torsionofCprime only. No lengthclaimfororiginalH0(D0),originalnilideal, orsocle.',genuine_target='The canonicalmap M[ell]→ker(ell:M/ell³→M/ell³)/imell² is anR0-linear isomorphism. ByDstructurefreeblocks contributeexactlyimell², both torsionblocks remaininjective because exponents<3. Quotientdimension2; no basis/actionscomputed.',extra_corollaries='For k>=3, length(M/m^kM)=2k+3. Sinceell²M⊂m²M, length(m²M/ell²M)=length(M/ell²M)−length(M/m²M)=7−6=1, andm kills it becausem³M⊂ell²M. Thusm²M/ell²M≅k.',new_value_job=False)
 ep=OUT/'independent-exact-torsion-full-HF-and-canonical-target-proof.json';ep.write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n')
 for q in[Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/f'code/ambient_B18_Cprime_exact_torsion_12_20261105_v63.py',ROOT/f'research/notes/{primary}.md']:pins[str(q.relative_to(ROOT))]=sha(q)
 assert pins[f'code/ambient_B18_Cprime_exact_torsion_12_20261105_v63.py']==r['code_sha256']and pins[f'research/notes/{primary}.md']==r['note_sha256']
 final=dict(id=ID,all_passed=True,frozen=True,pure_same_module_proof_accepted=True,Cprime_torsion_exponents=[1,2],Cprime_DVR_torsion_length=3,Cprime_full_HF='3,3,3,2,2,...',Cprime_Hilbert_series='(3-t^3)/(1-t)',Cprime_m2_over_ell2_length=1,genuine_target_kerell_mod_imell2_canonical=True,no_original_H0_nilindex_or_actionbasis_claim=True,no_value_source_gauge_jet_job=True,input_sha256=pins,evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-st)
 (OUT/'report.json').write_text(json.dumps(final,indent=2)+'\n');print(json.dumps({k:v for k,v in final.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
