#!/usr/bin/env python3
"""Four remaining active G3 source directions; no rank or unknown values."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE_PROFILE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def select(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 g=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json')
 c=read('reports/AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63/actual-G3-plane-condition-fixed-rank-and-kernel.json')
 active16=[i for i in g['selected_columns'] if i<238];assert len(active16)==16 and set(active16)<=set(c['selected_original_source_columns'])
 extra=[j for j in c['selected_original_source_columns'] if j not in active16];assert extra==[77,101,117,125]
 fullpath='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
 G3=list(map(H.parse,select(fullpath,'actual_full_G3linear102_by418_columns')));PiY=list(map(H.parse,select(fullpath,'PiY3_from102_columns')))
 mup='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
 basis=list(map(H.parse,select(mup,'actual_selected_Y0_columns102')))
 norm={p:H.parse(v) for p,v in select(mup,'normalized_basis_rows')}
 reps={p:H.parse(v) for p,v in select(mup,'normalized_basis_original_selected_source_representations')}
 sources=select(mup,'actual_selected_Y0_source_columns');assert len(basis)==len(sources)==len(norm)==99
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json')
 nf=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json')
 memo={i:H.parse(v) for i,v in nf['all264_linear_NF_columns']};mi={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 def monnf(target,u,v):
  if u==v==8:return {}
  m=[0]*9;m[u]+=1;m[v]+=1;return memo[lin['all264_NF_root_nodes'][target][mi[tuple(m)]]]
 allrecords=[];allunknown=set();slots=20000+sum(map(len,G3+PiY+basis))+sum(map(len,norm.values()))+sum(map(len,reps.values()))
 for source in extra:
  factors=H.parse(g['all418_original_basis_factors'][source]);assert all(i<16 for i in factors)
  comb=add({source:F(1)},{active16[k]:v for k,v in factors.items()},-1)
  first=H.guard(mv(G3,comb));b1={}
  for var in range(9):
   q={i//9:v for i,v in first.items() if i%9==var};original=dict(q);assert not mv(PiY,q);a={}
   while q:
    p=min(q);assert p in norm;x=q[p];q=H.guard(add(q,norm[p],-x));a=H.guard(add(a,reps[p],x))
   assert mv(basis,a)==original
   for j,x in a.items():b1[j*9+var]=x
  # Unknowns are actual post27 selected-Ystar first derivatives, thenPiY0.
  # Everyz-sensitive term is retained until its literal NF check is zero.
  sensitivity={};zchecks=0
  for k,x in b1.items():
   j,var=divmod(k,9)
   for target in range(3):
    for u in range(9):
     v=monnf(target,var,u)
     if var<7 or u<7:assert not v;zchecks+=1
     if v:
      for output,y in v.items():
       flat=output*99*3*9+j*3*9+target*9+u;sensitivity=H.guard(add(sensitivity,{flat:-x*y}))
  for flat in sensitivity:
   rem=flat%(99*3*9);j,rem=divmod(rem,27);target,u=divmod(rem,9);allunknown.add((j,target,u))
  allrecords.append({'source_extra_index':source,'actual_projected_zero_linear_source238_combination':srec(comb),
   'full102_by9_first_term':srec(first),'unique_selected99_linear_correction':srec(b1),
   'unknown_post27_PiY_Ystar1_sensitivity6_by99by3by9':srec(sensitivity),
   'all_full102_first_equals_selectedY0_b1_identities':True,'all_z_unknown_sensitivity_terms_NF_zero':True,'literal_z_NF_checks':zchecks})
  slots+=len(comb)+len(first)+len(b1)+len(sensitivity)
 assert slots<150000
 data={'active_projected_rank16_sources':active16,'conditionrank20_sources':c['selected_original_source_columns'],
  'remaining_four_source_directions':allrecords,'required_unknown_derivative_labels':[list(k) for k in sorted(allunknown)],
  'selected99_unit_source741_indices':sources,'full_active_projected_kernel_dimension':222,'known_plane_kernel_dimension':218,
  'no_new_rank_search':True,'no_unknown_values_queried':True,
  'exact_actual_quadratic_formula':{'top':'NF(PiY0*G3_2*c - PiY0*Ystar1*b1)','bottom':'NF(-Hstar1*b1)',
   'Ystar1':'ACTUAL selected99 source derivative AFTER27unitcancellation; it includesunknownYp0*pivotcarry. DoNOTreplaceit withrawYj1 ordeclareM1only.',
   'Hstar1':'actualpost27bottom3 selectedsource derivative, availableinprinciplefromfull9X1; not evaluatedhere',
   'unknowns':'onlyr/ell derivative slots surviveafterknown378zNF checks; formalindependence isretained',
   'source_complement_proof':'active16selectedsources subsetcondition20selectedsources; extra4 adjustedbyactiveprojectedfactors areprojectedzero andindependentmodconditionkernel218 byfixedconditionminor. No newrank.'},
  'owned_QQ_slot_bound':slots,'actual_T2_HF2_or_length_not_computed':True}
 ep=out/'actual-four-active-source-complements-and-unknown-derivative-profile.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'source_complement_count':4,'remaining_formal_derivative_labels':len(allunknown),
  'new_rank_search_or_coefficient_values':False,'actual_T2_HF2_or_length_computed':False,'owned_QQ_slot_bound':slots,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_complement_count','remaining_formal_derivative_labels','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
