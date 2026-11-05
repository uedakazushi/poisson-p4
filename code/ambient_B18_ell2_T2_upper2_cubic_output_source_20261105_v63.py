#!/usr/bin/env python3
"""Known-point upper2 extraction at ONE cubic output monomial; no cycles."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER2_CUBIC_OUTPUT_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 first=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_INPUT_GATE_V2_20261105_v63/pruned-known-firsts-and-positive-M1-four-invariant-test-input.json');G2=list(map(H.parse,first['G2linear_columns']));lam=list(map(H.parse,first['lambda_two']))
 pre='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';Pi=list(map(H.parse,field(pre,'PiY3_from102_columns')))
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json');nf=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json');memo={i:H.parse(v) for i,v in nf['all264_linear_NF_columns']};ix={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 rr=[0]*9;rr[7]=2;re=[0]*9;re[7]=re[8]=1;pairs=[rr,re];top=[[memo[lin['all264_NF_root_nodes'][t][ix[tuple(m)]]] for m in pairs] for t in range(3)];Q=[top[2][0],top[2][1]]
 assert read('reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/report.json')['all_passed']
 maximum=0
 def mark(b):
  nonlocal maximum
  maximum=max(maximum,b[0].bit_length(),b[1].bit_length());return b
 def exact(x):x=F(x);return mark((x.denominator,abs(x.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def bounddot(row,col):
  q=(1,0)
  for i,x in row.items():q=plus(q,mul(exact(x),exact(col.get(i,F(0)))))
  return q
 GR=[{i//9:x for i,x in v.items() if i%9==7} for v in G2];GE=[{i//9:x for i,x in v.items() if i%9==8} for v in G2];SB=[[bounddot(lam[0],GE[j]),bounddot(lam[1],GE[j])] for j in range(102)]
 for j in range(102):
  for row in lam:bounddot(row,GR[j])
  for k in range(6):
   a=(1,0)
   for t,c in Pi[j].items():a=plus(a,mul(exact(c),exact(top[t][0].get(k,F(0)))))
   b=plus(mul(exact(Q[0].get(k,F(0))),SB[j][0]),mul(exact(Q[1].get(k,F(0))),SB[j][1]));plus(a,b)
   a=(1,0)
   for t,c in Pi[j].items():a=plus(a,mul(exact(c),exact(top[t][1].get(k,F(0)))))
   plus(a,mul(exact(Q[1].get(k,F(0))),SB[j][0]))
 assert maximum<=384;H.dump(out/'phase.json',{'phase':'fixed_known_matrix_and_NF_before_signed_point_checks','positive_prior_bits':maximum,'no619_relation_or_higher_source_coefficient':True})
 bits=0
 def ck(x):
  nonlocal bits
  bits=max(bits,abs(x.numerator).bit_length(),x.denominator.bit_length());assert bits<=maximum;return x
 def dot(a,b):
  q=F(0)
  for i,x in a.items():q=ck(q+ck(x*b.get(i,F(0))))
  return q
 s0=[];s1=[]
 for j in range(102):
  assert dot(lam[0],GR[j])==0;a=dot(lam[0],GE[j]);assert dot(lam[1],GR[j])==a;s0.append(a);s1.append(dot(lam[1],GE[j]))
 for j in range(102):
  for mon in range(2):
   for k in range(6):
    actual=F(0)
    for t,c in Pi[j].items():actual=ck(actual+ck(c*top[t][mon].get(k,F(0))))
    expected=ck(Q[1].get(k,F(0))*s0[j]) if mon else ck(ck(Q[0].get(k,F(0))*s0[j])+ck(Q[1].get(k,F(0))*s1[j]));assert actual==expected
 # The omitted z-containing quadratic slots vanish in the known cover and
 # ell² is zero in the coefficient ring; no generic Taylor slot is omitted.
 for target in range(6):
  for z in range(7):
   for var in range(9):
    mon=[0]*9;mon[z]+=1;mon[var]+=1;assert not memo[lin['all264_NF_root_nodes'][target][ix[tuple(mon)]]]
 ep=out/'fixed-upper2-one-cubic-output-identity-and-source-frontier.json';H.dump(ep,{'s0_on_F2_source102':srec({j:x for j,x in enumerate(s0) if x}),'s1_on_F2_source102':srec({j:x for j,x in enumerate(s1) if x}),'lambda_two':[srec(v) for v in lam],'Q_two6':[srec(v) for v in Q],
  'all102_times2_times6_point_NF_identities':True,'all378_z_times_linear_NF_zero':True,'universal_polynomial_identity':'ForU2homogeneousquadratic: upperNF=Q0*s0(Urr)+Q1*(s0(Urell)+s1(Urr)); lambda0(G2lin U2)_{r²ell}=s0(Urr), lambda1(...)=s0(Urell)+s1(Urr). Ordinarycoefficients,nofactorials.',
  'application_condition':'ApplyonlytoanactualF2vectorafterALL27/99/G3clearingwhoseconstantandlinearpartsarezero. RawYbar1isNOTzeroandunknownPi1/Yp0/a2/M1carriesMUSTremainuntilsourceboundcancellationproved.',
  'uncomputed619_source_frontier':['literalcorrectedF2 leadingzero sourceoperator','M0*effectiveX3_r²ell','M1*effectiveX2 andpointpivotYp0/a2carry','G2quadratic*remainingY1 ifnotfullycleared','unit99/G3clearingcoefficientsandtheirsourceadjoint contractions','honestquotientthirdpureplane r²ell sourceancestorandcachetest'],
  'no_full619_operator_or_cycle_values_or_cache_claim':True,'no_new_physical_gauge_jet_or_rank_search':True})
 report={'id':STEM,'all_passed':True,'upper2_point_NF_and_cubic_output_identities_verified':True,'cubic_output_monomial_only':[7,7,8],'source_order_target_only_not_actual619_formula':True,'no_new_source_coefficients_or_cycle_values':True,'positive_prior_bits':maximum,'actual_known_point_bits':bits,'owned_numeric_bound':count([first,G2,Pi,lin,nf,memo,GR,GE,SB,s0,s1])+10000,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};assert report['owned_numeric_bound']<60000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_prior_bits','actual_known_point_bits','owned_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
