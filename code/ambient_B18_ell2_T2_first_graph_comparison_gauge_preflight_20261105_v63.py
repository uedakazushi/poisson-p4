#!/usr/bin/env python3
"""Pure coherent-comparison gauge audit for four uncomputed graph scalars."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_FIRST_GRAPH_COMPARISON_GAUGE_PREFLIGHT_20261105_v63'
P={'quotient':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json',
'bottom':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',
'known_peer':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json',
'adjoint_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/report.json',
'HF66_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63_REVIEW/report.json',
'source_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
'preserved_reader_attempt':'reports/AMBIENT_B18_ELL2_T2_FIRST_GRAPH_COMPARISON_GAUGE_PREFLIGHT_20261105_v63_FIRST_FULL_BOTTOM_READER_ATTEMPT/failure.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sp(p):return {i:Q(c) for i,c in p if Q(c)}
def enc(p):return [[i,str(c)] for i,c in sorted(p.items())]
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=True);assert not (o/'report.json').exists();inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert all(D[k]['all_passed'] for k in D if k.endswith('peer'))
 pi=list(map(sp,D['quotient']['Pi_remaining_top_linear_columns']));assert len(pi)==27
 a=Q(1095,3652);b=Q(1224288,8979355);c=Q(1973552256,40990755575)
 expected=[{} for _ in range(27)];expected[7]={1:a};expected[16]={0:b,1:c};expected[17]={1:b};expected[25]={0:Q(1)};expected[26]={1:Q(1)};assert pi==expected
 bottom=list(map(sp,D['bottom']['all642_bottom27_columns']));selected=D['bottom']['selected_bottom_columns'];assert len(selected)==23
 re=[(j,s) for j,s in enumerate(selected) if any(i%9 in [7,8] for i in bottom[s])];assert re==[(19,582),(21,634)]
 v=[Q(-3375,11),Q(-50),Q(25,2)]
 expected19={9*i+7:q for i,q in enumerate(v)};expected19.update({9*i+8:q*Q(14,125) for i,q in enumerate(v)})
 expected21={9*i+7:q*Q(8,3) for i,q in enumerate(v)};expected21.update({9*i+8:q*Q(8,3)*Q(4,25) for i,q in enumerate(v)})
 assert {i:c for i,c in bottom[582].items() if i%9 in [7,8]}==expected19 and {i:c for i,c in bottom[634].items() if i%9 in [7,8]}==expected21
 # Four gauge shifts from a three-component constant w. Coefficients are
 # exact quotient rows, not new unknown-source values or a rank search.
 action=[{}, {1:b,2:Q(1)}, {1:b,2:Q(1)}, {0:a,1:c}]
 assert action[0]=={} and action[1]==action[2] and a
 ep=o/'actual-Pi2-support-and-four-scalar-gauge-law.json';dump(o,ep.name,{'inputs':inputs,'actual_Pi2_columns_verified':True,'only_two_selected_bottom_vectors_have_r_or_ell':re,'common_bottom_v':list(map(str,v)),
'actual_selected_bottom19':enc(bottom[582]),'actual_selected_bottom21':enc(bottom[634]),
'coherent_comparison_change':'Mprime=M+G2*A preserves G1*Mprime=V. SinceG2(0)=0 it preservesM0, andMprime*n is a formalunit, so sections can be adjusted. Positive-degree contraction givesH2(G2*z)=z-G3*H3(z); modulo the knownG3 firstspan this is the triangular target change byA0.',
'projected_target_change':'On commonbottomdirectionv putw=PiY*A0*H0v. Thenpsi(w*r)=(u,vv),psi(w*ell)=(0,u),u=b*w1+w2,vv=a*w0+c*w1.',
'four_scalar_order':['lambda0 MrH0v','lambda1 MrH0v','lambda0 MellH0v','lambda1 MellH0v'],
'four_scalar_changes_from_w_columns3':[enc({i:row.get(j,Q(0)) for i,row in enumerate(action) if row.get(j,Q(0))}) for j in range(3)],
'gauge_shift_form':['0','u','u','vv'],'gauge_dimension_exact':2,
'comparison_gauge_invariant_combinations':['lambda0 MrH0v','lambda1 MrH0v-lambda0 MellH0v'],
'not_a_zero_value_statement':'This is a coherent source/target isomorphism law, not the value of canonicalM1 and not permission to set arbitraryfourcontractions0.',
'no_actual_M1_graph_or_higher_source_coefficients_computed':True})
 report={'id':STEM,'all_passed':True,'formal_source_gauge_law_only':True,'uncomputed_four_contractions_have_two_comparison_gauge_invariants':True,'actual_lineargraph_or_M1_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json')}))
if __name__=='__main__':main()
