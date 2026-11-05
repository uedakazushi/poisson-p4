#!/usr/bin/env python3
"""Fixed point Psi correction and conditional complete b1 source proof."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);hs={};maximum=0
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 qpath='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json';Q=list(map(H.parse,field(qpath,'Q_two6_columns')));L=list(map(H.parse,field(qpath,'Q_dual_two6_rows')))
 lin='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json';mon=field(lin,'degree2_monomial_order44');std=field(lin,'degree2_linear_standard_module_monomials');roots=field(lin,'all264_NF_root_nodes');N=dict((i,H.parse(v)) for i,v in field('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json','all264_linear_NF_columns'))
 pre='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';Pi=list(map(H.parse,field(pre,'PiBottom3_from56_columns')));assert len(Pi)==56
 psi_rr=[{i:row.get(2*i,F(0)) for i in range(3) if row.get(2*i)} for row in L];psi_re=[{i:row.get(2*i+1,F(0)) for i in range(3) if row.get(2*i+1)} for row in L]
 def mark(q):
  nonlocal maximum
  maximum=max(maximum,q[0].bit_length(),q[1].bit_length());return q
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def bounddot(a,b):
  q=(1,0)
  for i,x in a.items():q=plus(q,mul(ex(x),ex(b.get(i,F(0)))))
  return q
 for row in L:
  for col in Q:bounddot(row,col)
 for row in psi_rr+psi_re:
  for col in Pi:bounddot(row,col)
 assert maximum<=384;H.dump(out/'phase.json',{'phase':'positive_point_bound_before_LQ_and_DeltaC','positive_prior_bits':maximum,'no_actual_CU_K_or_relation_values':True})
 bits=0
 def ck(x):
  nonlocal bits
  bits=max(bits,abs(x.numerator).bit_length(),x.denominator.bit_length());assert bits<=maximum;return x
 def dot(a,b):
  q=F(0)
  for i,x in a.items():q=ck(q+ck(x*b.get(i,F(0))))
  return q
 assert [[dot(row,col) for col in Q] for row in L]==[[F(1),F(0)],[F(0),F(1)]]
 rr=[0]*9;rr[7]=2;re=[0]*9;re[7]=re[8]=1;pairs=[rr,re]
 for i in range(3):
  for k,m in enumerate(pairs):
   target=3+i;coord=2*i+k;mi=mon.index(m);assert std[coord]==[target,m] and N[roots[target][mi]]=={coord:F(1)}
 correction=[]
 for ch in range(2):
  correction.append({'channel':ch,'r_row56':[[j,str(dot(psi_re[ch],col))] for j,col in enumerate(Pi) if dot(psi_re[ch],col)],'ell_row56':[[j,str(dot(psi_rr[ch],col))] for j,col in enumerate(Pi) if dot(psi_rr[ch],col)]})
 premises={
 'exact_fullsource_clearing':'reports/AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63_REVIEW/report.json',
 'pure_C_compression':'reports/AMBIENT_B18_ELL2_T2_UPPER617_C_COMPRESSION_SOURCE_20261105_v63/report.json',
 'coherent_free3_frame':'reports/AMBIENT_B18_ELL2_T2_UPPER_FREE3_GAUGE_SOURCE_20261105_v63_REVIEW/report.json',
 'complete393_and617_sourcecoverage':'reports/AMBIENT_B18_ELL2_T2_LOWER_RANK2_HF_TAIL_SOURCE_20261105_v63/report.json',
 'all393_actual_zero_active222':'reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63_REVIEW/report.json',
 'all393_actual_zero_nonactive171':'reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUES_20261105_v63_REVIEW/report.json'}
 # Exact existing report paths are resolved before any response calculation;
 # only all_passed metadata is read for these already adopted source premises.
 for p in premises.values():assert read(p)['all_passed']
 ep=out/'known-Psi-point-correction-and-full-b1-source-annihilation.json';data={'Q_two6':field(qpath,'Q_two6_columns'),'Q_left_inverse_two6_rows':field(qpath,'Q_dual_two6_rows'),'all4_LQ_identities':True,'all6_bottom_standard_NF_identity':True,'psi_rr_three_rows':[[[i,str(x)] for i,x in sorted(v.items())] for v in psi_rr],'psi_rell_three_rows':[[[i,str(x)] for i,x in sorted(v.items())] for v in psi_re],'DeltaC_known_point_rows':correction,
 'full_coordinate_identity':'actualupper_t=f_Ut+Psi_t(B2),Cfull=C_U+DeltaC,Kfull=GammaXstar2+CfullXstar1','do_not_use_KU_annihilation_or_B2zero':True,
 'full_K_b1_annihilates_actual_G3linear_image_source_proof':'applyexactclearingtoX=0,Y=G3v;subtractfixed25sourceprojection;remaining393truecompletequadraticclasses0 =>minusKfullB99C G3linv=0. Qleftinverseandfullbottomtermareessential.',
 'conditional_strand_choice_independence':'IfrankG2active190andknownzero3tensorproof,kerG2linear918=imG3linear242. ActualfullKcovectorannihilateskernel,soexistsrhoonallquadraticF1outputswithrhoG2linear=alphaK. Actualpairing=rhoM0Z2. Explicitfixedrightsolveandresource/valuesgateSTILLrequired.',
 'no_actual_CU_K_X2_X3_M1_or_relation_values':True,'source_premises':premises};H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'fixed_point_Psi_and_full_source_proof_only':True,'positive_prior_bits':maximum,'actual_fixed_point_bits':bits,'actual_K_or_upper_values_queried':False,'G2_strand_rank_not_assumed':True,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'owned_numeric_bound':count([Q,L,mon,std,roots,N,Pi,psi_rr,psi_re,correction,data])+14000,'seconds':time.monotonic()-start};assert report['owned_numeric_bound']<60000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','positive_prior_bits','actual_fixed_point_bits','owned_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
