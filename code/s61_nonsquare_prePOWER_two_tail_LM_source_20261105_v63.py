#!/usr/bin/env python3
"""Saved-support SOURCE only: ordinary and exactly-named FULLU nondivision; no new target products/GCD."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 gp=load('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63');assert gp['all_passed']
 gr=load('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW');assert gr['all_passed']
 unitpath=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63/unit-inventory.json')
 assert sha(unitpath)==gp['evidence_sha256'][unitpath.name]
 assert sha(unitpath)==gr['input_sha256'][str(unitpath.relative_to(ROOT))]
 unit=json.loads(unitpath.read_text());assert unit['ring']=='QQ[a,e,d,1/(a*d*F*B*C*J)]'
 assert unit['factor_order']==['a','d','F','B','C','J']
 factors=[{tuple(e):Q(c)for e,c in p}for p in unit['named_factors']]
 normal=[]
 for p in factors:
  lead=p[max(p)];normal.append(tuple(sorted((e,c/lead)for e,c in p.items())))
 assert len(set(normal))==6
 # Complete localization product, not any old B-rank/GB assertion.
 old=load('S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63','certificate.json')
 expected={tuple(e):Q(c)for e,c in old['FULLU_product']};product={(0,0,0):Q(1)};largest=1
 for fac in factors:
  z={}
  for e,c in product.items():
   for f,v in fac.items():k=tuple(e[i]+f[i]for i in range(3));z[k]=z.get(k,Q(0))+c*v
  product={e:c for e,c in z.items()if c};largest=max(largest,len(product))
 assert product==expected and len(product)==49
 localization=dict(exact_ring=unit['ring'],complete_inverted_irreducible_list=unit['factor_order'],
  accepted_GLOBAL_review_sha256=sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW/report.json'),
  accepted_inventory_sha256=sha(unitpath),
  original_FULLU_product49_literal_match=True,old_certificate_used_ONLY_for_FULLU_product=True,
  no_new_or_hidden_inverted_factors=True,NOT_inverted_in_base_ring=unit['NOT_inverted']+['T','D12'],
  later_T_D12_are_corollary_localizers_ONLY=unit['later_corollary_localizers_only'],
  no_conclusion_after_adding_D12_Delta_Saux_or_any_other_unknown_localizer=True,
  factor_irreducibility=dict(a='variable',d='variable',F='linear in a, leadingconstant500 a unit overQQ[e,d]',
   C='linear in a, leadingconstant100',J='linear in e, leadingconstant-6 overQQ[a,d]',
   B='linear in d andprimitive: onlypossible coefficientcommonfactor3e-2 doesnotdivideb1 sinceb1(e=2/3)=50a'),
  six_pairwise_nonassociate_by_literal_scalar_normal_forms=True,
  base_change_C='all primitive/linear/nonassociation arguments remain true overC; no new numeric sourcepoint',
  ordinary_UFD_unit_group='nonzeroQQconstants times integerpowers of exactly a,d,F,B,C,J',
  fixed_unit_product_only_max_intermediate_terms=largest)
 dump('complete-localization-and-factor-list.json',localization)
 r0=load('S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_20261105_v63');assert r0['all_passed']
 r2=load('S61_NONSQUARE_PREPOWER_Q1_B4_ROOT_MULTINOMIAL_REVIEW_20261105_v63');assert r2['all_passed']
 paths=[('N0','S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63','Q1unused','T0C4QV001',4500),
        ('N2','S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63','Q1unused','Q1B4QV001',2317)]
 saved=[]
 for name,st,unused,magic,n in paths:
  primary=load(st);assert primary['all_passed'];p=pin(ROOT/'reports'/st/'quotient-integer-common-D.txt')
  assert sha(p)==primary['evidence_sha256'][p.name]
  lm=None;prev=None;mins=[10**9]*3;maxs=[0]*3;td=0;count=0
  with p.open()as f:
   assert f.readline().strip()==magic;D=int(f.readline());assert D>0 and D.bit_length()==197;assert int(f.readline())==n
   for line in f:
    key,c=map(int,line.split());assert c and(prev is None or prev<key);prev=key
    a,x=divmod(key,1048576);e,d=divmod(x,1024);v=(a,e,d);lm=v if lm is None or v>lm else lm
    mins=[min(s,t)for s,t in zip(mins,v)];maxs=[max(s,t)for s,t in zip(maxs,v)];td=max(td,sum(v));count+=1
   assert count==n
  saved.append(dict(name=name,terms=n,original_D197=D,LM_a_e_d=lm,min_exponents=mins,max_exponents=maxs,totaldegree=td,literal_sha256=sha(p)))
 assert saved[0]['original_D197']==saved[1]['original_D197']
 assert saved[0]['LM_a_e_d']==(16,10,18)and saved[1]['LM_a_e_d']==(12,11,14)
 assert saved[1]['min_exponents']==[0,0,0]
 valdeps={}
 for name,st in [('C4exact','S61_NONSQUARE_PREPOWER_T1_C_HASSE_HORNER_REVIEW_20261105_v63'),
  ('F0','S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_20261105_v63'),
  ('J0','S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_20261105_v63'),
  ('B4exact','S61_NONSQUARE_PREPOWER_Q1_B_ROOT_JET_HORNER_REVIEW_20261105_v63')]:
  r=load(st);assert r['all_passed'];valdeps[name]=sha(ROOT/'reports'/st/'report.json')
 coprime=dict(complete_factor_list_covered=['a','d','F','B','C','J'],saved_valuation_dependencies=valdeps,
  a='min_a(N2)=0 from literal2317 support',d='min_d(N2)=0 from literal2317 support',
  F='vF(Q1)=0; F dividingN2 woulddivideQ1=B4*N2/D, contradiction',
  J='vJ(Q1)=0; J dividingN2 woulddivideQ1, contradiction',
  C='exactvC(T1)=4; T1=a3*C4*Q1 andC doesnotdividea implyvC(Q1)=0; C dividingN2 woulddivideQ1',
  B='exactvB(Q1)=4 plus acceptedfullidentityQ1=B4*N2/D givesvB(N2)=0',
  no_valuation_for_D12_Delta_Saux_or_other_chartfactor_claimed=True)
 dump('N2-coprime-to-ALL-inverted-factors.json',coprime)
 proof=dict(SAVED_SUPPORT_ONLY_NO_NEW_TARGET_COEFFICIENT_VALUES=True,polynomials=saved,
  fixed_order='lex a>e>d; packedkey is order-compatible over allprovedboxes',
  LM_quotient_difference=[4,-1,4],
  ordinary_nondivision='LM(N2) has e11 whileLM(N0) has e10, so LM(N2)doesnotdivideLM(N0); hence N2 doesnotdivideN0 overQQ orC',
  named_localization_nondivision='ifN0=N2*H/s withs producta,d,F,B,C,J, thenN2 divides s*N0. ALLinvertedirreducibles coprimeN2, henceUFDcancellation impliesN2|N0, contradiction.',
  exact_ring_ONLY=unit['ring'],Q0='N0/D=T0/(a4*C4)',Q2='N2/D=T1/(a3*C4*B4)',
  ideal_consequence='(Q2) is strictly contained in (Q0,Q2) in this SIX-factor base ring',
  scope='two-tail schemeideals only; NOT actualbadFitting ideal, no strictopen/radicalsupport/unitideal/regularsequence/emptycomplement assertion',
  other_localizations_NOT_covered='no claim forR[1/D12],oldOmega,Delta407,Saux,etc; strictness can disappear after furtherlocalization',
  no_GCD_factor_search_rank_or_new_value_job=True)
 dump('ordinary-and-named-localized-LM-certificate.json',proof)
 pin(ROOT/f'research/notes/{ID}.md');pin(Path(__file__))
 result=dict(id=ID,all_passed=True,frozen=True,source_only=True,ordinary_N2_nondivides_N0=True,
  named_SIX_factor_localized_nondivision_source_proof=True,complete_inverted_factor_list_pinned=True,
  no_other_chart_localization_or_radical_actualbad_claim=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',result);print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
