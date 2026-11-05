#!/usr/bin/env python3
"""Pure source obstruction to full-C cubic dimension4; no graph readbacks."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63'
PATHS={
 'reduction':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/report.json',
 'reduction_evidence':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json',
 'reduction_peer':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/report.json',
 'complete':'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63/complete155-source-cubic-operator-with-four-formal-alphas.json',
 'complete_peer':'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json',
 'alpha':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63/report.json',
 'alpha_peer':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'cut':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63/one-CUT1-residual-and-formal-coordinate-rank-certificate.json',
 'cut_peer':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63_REVIEW/report.json',
 'cache_source':'reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_SOURCE_20261105_v63/report.json',
 'cache_custody':'reports/AMBIENT_B18_C_FULL_GRAPH03_VALUE_CUSTODY_20261105_v63/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['reduction']['unconditional_m3C_equals_ell_m2C_proved']
 assert d['alpha']['nonzero_count']==0 and d['alpha']['actual_C3_HF_degree3']==4
 assert d['cut']['actual_C3_length19_proved']and d['cut']['actual_C3_HF_if_positive']==[3,6,6,4,0]
 assert d['complete']['all_other151_columns_zero_in_known_cover']
 nf=list(map(lambda row:list(map(Q,row)),d['reduction_evidence']['full_no_cut_cubic_NF6_rows12_columns']))
 # Actual C3 gr3 is exactly the complete four-dimensional cover because
 # all four alpha values are0. The no-cut ell3 columns are discarded here
 # only to identify this ALREADY PROVED quotient's nonzero class.
 cutNF=[nf[i]for i in[0,1,3,4]]
 rell2A=[row[2]for row in cutNF];assert rell2A==[0,1,0,0]
 rell2B=[row[6]for row in cutNF];assert rell2B==[0,0,0,1]
 # Both r*ell² classes are nonzero in the actual gr3(C/ell³C), not just
 # nonzero in an upper cover: completeness+alpha0 is essential here.
 assert len(d['reduction_evidence']['five_graph_witness_source_bindings'])==5
 assert not(ROOT/'reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_VALUES_20261105_v63').exists()
 data={
  'scope':'FULL actual comparison C over the same regular9 R0; D=QQ[[ell]], m=m_R0. Not original H0(D0), not an algebra splitting.',
  'accepted_data':{'generic_D_rank':4,'full_HF_degrees0_1_2':[3,6,6],'unconditional_m3C_equals_ell_m2C':True,
   'C3_actual_HF':[3,6,6,4,0],'C3_actual_length':19,'C_DVR_torsion_current':[1,1,2,'e>=3']},
  'actual_C3_gr3_basis':['A*r²ell','A*r ell²','B*r²ell','B*r ell²'],
  'actual_nonzero_r_times_ell2_classes':{'A':list(map(str,rell2A)),'B':list(map(str,rell2B))},
  'general_lemma_hypothetical_HF3_4_steps':[
   'Put F=m³C. The reduction theorem gives mF=m4C=ellF, so mu_D(F)=HF3(full C). F has D-rank4.',
   'If HF3(full C)=4, the finite DVR module F is D-free4.',
   'C3 m4=0 gives ellF=m4C⊂ell³C, while ell³C⊂F. Actual length(m³C3)=4 equals length(F/ellF)=4. Hence ell³C=ellF=m4C.',
   'For every x∈C choose y∈F with ell³x=ell*y. Then t=ell²x−y∈C[ell].',
   'For each u∈m, u*t∈F since u*ell²x and u*y lie in F. It is ell-torsion, hence is0 because F is D-free.',
   'Therefore u*ell²x=u*y∈mF=m4C=ell³C. Thus m*ell²C⊂ell³C.',
   'Take u=r,x=the actual A target. Its r*ell² image in C3 has the proven nonzero coordinate(0,1,0,0), contradiction.'
  ],
  'unconditional_full_C_HF_degree3_values':[5,6],
  'unconditional_m3C_DVR_torsion_generator_count_values':[1,2],
  'five_graph03_source_zero_proof':'A nonzero graph d03 would give the independent two-column no-cut cubic minor d03², forcing full HF3<=4. The full-C reduction/rank lower is4 and the preceding lemma excludes4. Thus each of the SAME five actual graph delta_ell³ coefficients is0, without evaluating any root.',
  'five_graph03_actual_coefficients_source_forced_zero':[[r['physical_row'],'0']for r in d['reduction_evidence']['five_graph_witness_source_bindings']],
  'no_cache_value_worker_required_or_run':True,
  'unrun_cache_history':'The distinct graph03 SOURCE/custody/code remain frozen but values are unnecessary. Older graph r*ell² V1/V2 workers remain unrun separately. No fresh physical/gauge/IFT/jet.',
  'conditional_theorem_status':'The previous nonzero-graph=>1123 theorem remains logically valid, but its graph antecedent is now disproved. Do not infer exponent3 from that conditional.',
  'limitations':['The remaining exponent e is still only>=3. This does not prove e>=4 or exclude1123 by another mechanism.',
   'Full HF3 is5or6, not the quotient C3 HF3=4. Further fullC cubic classes are uncomputed.',
   'Tor_D(m³C) is nonzero; this is not an assertion that m³Tor_D(C) is nonzero.',
   'No original H0 length, nilcube, exact nilindex or total Poisson exhaustion claim.'],
  'source_profile':{'QQ_matrix_degree':3,'numeric_envelope':25000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'full-C-cubic-obstruction-and-five-graph03-source-zeros.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_FULL_C_HF3_5_OR6_AND_FIVE_GRAPH03_ZERO_NO_VALUES',
  'actual_full_C_HF_degree3_in_5_6_proved':True,'all_five_graph03_coefficients_zero_source_proved':True,'longer_C_torsion_exponent_determined':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
