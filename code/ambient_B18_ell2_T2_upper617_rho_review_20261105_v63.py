#!/usr/bin/env python3
"""Provided rho: original56unit and1836source identities, no reverseAD rerun."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63_REVIEW'
def main(expected):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']==expected
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 e=rd(b+'actual-two-rho1144-covectors-and-full918-identities.json');rho=list(map(sp,e['rho_two_by1144_sparse_rows']));unitrecords=e['actual_transpose_unit28_two_records'];assert len(rho)==len(unitrecords)==2
 s='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_GATE_20261105_v63_REVIEW/report.json';assert rd(s)['all_passed'] and pins[s]=='434d61e4d29ddf52307bacb489e56d97d6b557dcffcad13ac9fd27be6c34022d'
 g='reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_INPUT_GATE_20261105_v63/pruned-transpose-R-unit28-source-and-positive-rho-input.json';d=rd(g)
 original='reports/AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63/actual-factored918-rightsolve-source-circuit-and-all-generator-identities.json';z=fld(original,['active_unit_block28','quadratic_pair_order44','actual_full_G2linear_columns']);U=z['active_unit_block28']['matrix_integer28_by28'];pairs=z['quadratic_pair_order44'];assert pairs==[list(v) for v in __import__('itertools').combinations_with_replacement(range(9),2) if v!=(8,8)];ix={tuple(v):i for i,v in enumerate(pairs)}
 F=list(map(sp,fld('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['full_G2_linear234_by102_columns'])['full_G2_linear234_by102_columns']));assert F==list(map(sp,z['actual_full_G2linear_columns']))
 a='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63/actual-two-full918-alpha-covectors-and836-G3-identities.json';alpha=list(map(sp,fld(a,['alpha_two_by918_sparse_rows'])['alpha_two_by918_sparse_rows']));assert alpha==list(map(sp,d['alpha_two_by918_sparse_rows']))
 bound=0
 def positive(terms):
  nonlocal bound
  if not terms:return
  D=math.lcm(*(x.denominator*y.denominator for x,y in terms));V=sum(abs(x.numerator*y.numerator)*(D//(x.denominator*y.denominator)) for x,y in terms);bound=max(bound,D.bit_length(),V.bit_length());assert bound<=384
 for t0,row in enumerate(rho):
  assert set(row)<=set(d['possible_rho_input_slots_by_channel'][t0]);env={i:(int(D),int(V)) for i,D,V in d['positive_rho_envelopes_by_channel'][t0]}
  for i in range(1144):
   q=row.get(i,Q(0));D,V=env.get(i,(1,0));assert D%q.denominator==0 and abs(q.numerator)*D<=V*q.denominator
  w=list(map(Q,unitrecords[t0]['rhs28']));y=list(map(Q,unitrecords[t0]['solution28']));assert len(w)==len(y)==28
  for i in range(28):positive([(Q(U[j][i]),y[j]) for j in range(28)]+[(Q(-1),w[i])])
  for j,c in enumerate(F):
   for var in range(9):
    terms=[]
    for k,x in c.items():
     target,u=divmod(k,9);pair=tuple(sorted((u,var)))
     if pair!=(8,8):terms.append((x,row.get(44*target+ix[pair],Q(0))))
    terms.append((Q(-1),alpha[t0].get(9*j+var,Q(0))));positive(terms)
 put('phase.json',{'phase':'independent_positive_bound_before56plus1836_literal_identity_dots','positive_prior_bits':bound,'no_reverseAD_or_newunit_solve':True})
 actual=0
 def dotterms(ts):
  nonlocal actual
  q=Q(0)
  for x,y in reversed(ts):
   v=x*y;q=q+v;actual=max(actual,abs(v.numerator).bit_length(),v.denominator.bit_length(),abs(q.numerator).bit_length(),q.denominator.bit_length());assert actual<=bound
  return q
 checks=[]
 for t0,row in enumerate(rho):
  rec=unitrecords[t0];w=list(map(Q,rec['rhs28']));y=list(map(Q,rec['solution28']))
  for i in range(28):assert dotterms([(Q(U[j][i]),y[j]) for j in range(28)]+[(Q(-1),w[i])])==0
  assert rec['literal_Utranspose_residuals28']==['0']*28
  zs=[]
  for j,c in enumerate(F):
   for var in range(9):
    terms=[]
    for k,x in c.items():
     target,u=divmod(k,9);pair=tuple(sorted((u,var)))
     if pair!=(8,8):terms.append((x,row.get(44*target+ix[pair],Q(0))))
    terms.append((Q(-1),alpha[t0].get(9*j+var,Q(0))));assert dotterms(terms)==0;zs.append('0')
  checks.append(zs)
 assert checks==e['all1836_original_rhoF_minus_alpha_literal_zeros'] and len(checks[0])==len(checks[1])==918
 proof={'all56_original_Utranspose_given_rhs_solution_equations':True,'all1836_original_full9_Fsource_equations_rhoF_equalalpha':True,'all2288_saved_rho_envelopes_and_structuralzeros':True,'rho_supports':list(map(len,rho)),'validity_scope':'Providedrhoisavalidfullsourceadjoint. ForanytrueimageRHSFv, rho(Fv)=alpha(v). Equationsdonotclaimrhoisuniqueoutsidetheimage;thereverseDAGandunitRHSconstructionarenotrerun.','source_slot_convention':'F=26targetsby102sourcesoflinear9variables;rho1144=44*target+ordinarypairwithonlyell²omitted;alpha918=9*source+var.','positive_prior_bits':bound,'actual_identity_bits':actual,'no_reverseAD_unit_solve_RHSZ2_b1_upper_or_physical_higher_values':True,'rho_two_by1144_sparse_rows':e['rho_two_by1144_sparse_rows']}
 ep='independent-provided-rho56-unit-and1836-original-source-certificates.json';put(ep,proof)
 owned=count([r,e,rho,unitrecords,d,z,U,pairs,ix,F,alpha,checks,proof,pins])+12000;assert owned<60000
 report={'id':STEM,'all_passed':True,'all56_unit_and1836_original_source_identities_checked':True,'valid_rho_adjoint_on_true_image':True,'no_reverseAD_rerun_or_RHS_response_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'positive_prior_bits':bound,'actual_bits':actual,'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--expected-primary-sha',required=True);a=p.parse_args();main(a.expected_primary_sha)
