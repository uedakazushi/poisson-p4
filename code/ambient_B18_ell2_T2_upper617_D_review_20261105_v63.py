#!/usr/bin/env python3
"""Independent1120 known coefficient dots; no X/source response evaluation."""
from pathlib import Path
from fractions import Fraction as Q
import json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='fc5529eeec78a76c8f85b154b51ec05c2b8b7b0b010788f1eff1761f6b316e24'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 e=rd(b+'actual-D-coefficients-two-by10-by56.json')
 pr='reports/AMBIENT_B18_ELL2_T2_UPPER617_D_INPUT_GATE_20261105_v63_REVIEW/report.json';assert rd(pr)['all_passed'] and pins[pr]=='4bd3832e3f064b638cdc61414c6599fbb398bd82ebf3bb086bff9f863301aa39'
 C=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json',['Cfull_two_by2_by56'])['Cfull_two_by2_by56'];C=[[[Q(x) for x in row] for row in ch] for ch in C]
 rho=list(map(sp,fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63/actual-two-rho1144-covectors-and-full918-identities.json',['rho_two_by1144_sparse_rows'])['rho_two_by1144_sparse_rows']))
 M=list(map(sp,fld('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['Mnew0_columns'])['Mnew0_columns']))
 gp='reports/AMBIENT_B18_ELL2_T2_UPPER617_D_INPUT_GATE_20261105_v63/pruned-D-coefficient-functional-and-ten-pair-source-input.json';g=rd(gp);pairs=list(map(tuple,g['quadratic_pair_order44']));active=list(map(tuple,e['active_pairs10']));assert active==list(map(tuple,g['active_pairs10']));ix={v:i for i,v in enumerate(pairs)}
 R=[[{i:rho[t0][44*i+ix[p]] for i in range(26) if rho[t0].get(44*i+ix[p])} for p in active] for t0 in range(2)]
 prior=0
 for t0 in range(2):
  for k,pair in enumerate(active):
   for j,col in enumerate(M):
    terms=[(x,col.get(i,Q(0))) for i,x in R[t0][k].items()];c=C[t0][1][j] if pair==(7,7) else C[t0][0][j] if pair==(7,8) else Q(0);terms.append((c,Q(1)))
    terms=[(x,y) for x,y in terms if x and y]
    D=math.lcm(*(x.denominator*y.denominator for x,y in terms));U=sum(abs(x.numerator*y.numerator)*(D//(x.denominator*y.denominator)) for x,y in terms);prior=max(prior,D.bit_length(),U.bit_length());assert prior<=92
 put('phase.json',{'phase':'independent_shared_rho_pairrow_positive_before1120_dots','source_prior_bits':prior,'no_X_or_callback_response':True})
 actual=0;output=[]
 def ck(q):
  nonlocal actual
  actual=max(actual,abs(q.numerator).bit_length(),q.denominator.bit_length());assert actual<=92;return q
 for t0 in range(2):
  ch=[]
  for k,pair in enumerate(active):
   row=[]
   for j,col in enumerate(M):
    q=Q(0)
    for i,x in reversed(list(R[t0][k].items())):q=ck(q+ck(x*col.get(i,Q(0))))
    c=C[t0][1][j] if pair==(7,7) else C[t0][0][j] if pair==(7,8) else Q(0);d=ck(c-q);assert d==Q(e['D_two_by10_by56'][t0][k][j]);row.append(str(d))
   ch.append(row)
  output.append(ch)
 assert sum(Q(x)!=0 for ch in output for row in ch for x in row)==122
 proof={'all1120_D_coefficients_includingzeros':True,'original_Cfull_Psi_rho_Mnew0_source':True,'independent_order':'Extractshared26componentrho_pairrow,computeentiredotwithM0firstinreversesourceorder,thensubtractfromC. ProducersequentiallysubtractedeachtermfromC.','Drr_Cell_Drell_Cr_no_factor2':True,'positive_prior_bits':prior,'actual_bits':actual,'nonzero122':True,'no_Xc2_Xp2_honestgradient_firstcross_RHSZ2_b1_upper_or_higher_values':True,'outer_four_OLDg2_and_native_L2S0_frontier_still_retained':True}
 ep='independent1120-original-D-pairrow-dots-and-zero-entries.json';put(ep,proof)
 owned=count([r,e,C,rho,M,g,pairs,active,ix,R,output,proof,pins])+8000;assert owned<60000
 report={'id':STEM,'all_passed':True,'all1120_known_D_coefficients_checked':True,'no_source_callback_or_upper_response_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()
