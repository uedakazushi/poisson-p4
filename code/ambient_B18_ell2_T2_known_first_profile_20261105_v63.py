#!/usr/bin/env python3
"""Known G3-linear and bottom-three first relations; no X2/M1 query."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_first_Hilbert_source_20261105_v63 import project_linear,quotient
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63'
def parse(v):return H.parse(v)
def fixed_rank(cols):
 norms={};reps={};original=[];selected=[];factors=[];receipts=[]
 for j,v0 in enumerate(cols):
  v=dict(v0);a={};steps=[]
  while v and min(v) in norms:
   p=min(v);c=v[p];v=H.guard(add(v,norms[p],-c));a=H.guard(add(a,reps[p],c));steps.append([p,str(c)])
  if v:
   p=min(v);c=v[p];i=len(original);original.append(v0);selected.append(j);norms[p]=H.guard(scale(v,1/c));reps[p]=H.guard(scale(add({i:F(1)},a,-1),1/c));a={i:F(1)}
   receipts.append({'column':j,'basis_column':i,'pivot':p,'divide_by':str(c),'subtract':steps})
  assert mv(original,a)==v0;factors.append(a)
 return norms,reps,original,selected,factors,receipts
def main(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized_known_first_arrays_no_X2'})
 hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 gr=read(gate);assert gr['all_passed'] and not gr['rank_started'];ip=str((Path(gate).parent/'actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'));d=read(ip);assert gr['evidence_sha256'][Path(ip).name]==hs[ip]
 read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json')
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json')
 G3=list(map(parse,d['actual_minimal6_G3linear_top27_by418_columns']));gn,greps,gorig,gsel,gfactors,greceipts=fixed_rank(G3)
 gfree,PiTop=quotient(27,gn)
 H.dump(out/'phase1-known-G3linear-fixed-source-rank.json',{'actual_projected_G3linear_columns':d['actual_minimal6_G3linear_top27_by418_columns'],
  'rank':len(gn),'selected_columns':gsel,'fixed_rank_receipts':greceipts,
  'normalized_basis':[[p,srec(v)] for p,v in sorted(gn.items())],'all418_original_basis_factors':[srec(v) for v in gfactors],
  'remaining_top_linear_quotient_coordinates':gfree,'Pi_remaining_top_linear_columns':[srec(v) for v in PiTop]})
 old640=H.selected(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9');old128=H.selected(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9');hs[H.PATHS['old640']]=H.sha(ROOT/H.PATHS['old640']);hs[H.PATHS['old128']]=H.sha(ROOT/H.PATHS['old128'])
 old=[];old0=[]
 for col in old640+old128:
  q={};v0={}
  for row,p in col:
   for k,c in parse(p).items():
    if k:q[row*9+k-1]=c
    else:v0[row]=c
  old.append(q);old0.append(v0)
 del old640,old128
 ellpath='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json';ell=list(map(parse,read(ellpath)['actual_first640_ell_columns']))
 source=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json');L=list(map(parse,source['L0_Res4_Ntranspose_columns']));S0=list(map(parse,source['native_S0_columns']));S1=list(map(parse,source['native_S1ell_columns']));eff=[mv(L,v) for v in S1]
 Pbottom=list(map(parse,d['PiBottom3_from56_columns']));M=list(map(parse,source['Mnew0_columns']))
 hs[H.PATHS['origin']]=H.sha(ROOT/H.PATHS['origin']);piv=mini['original27_X_source_columns'];Xp=[old0[j] for j in piv];assert sum(map(len,old0))==97;del old0
 C=mini['free30_V_indices'];B30inv=list(map(parse,mini['fixed_bottom30_inverse_columns']));rawcomb=[parse(r['original768_combination']) for r in d['all642_constant_Y0_kernel_source_combinations']]
 bnorm={};breps={};borig=[];bselected=[];bfactors=[];breceipts=[];selected_full1=[];kernels=[];allbottom=[];peak=0
 fixed=sum(map(len,old))+sum(map(len,ell+eff))+sum(map(len,G3))+sum(map(len,M+Pbottom+PiTop+B30inv+Xp))+sum(map(len,rawcomb))+12000
 for j,c in enumerate(rawcomb):
  cn={i-640:a for i,a in c.items() if i>=640};assert not mv(S0,cn)
  x=mv(old,c);xe={}
  for k,a in c.items():xe=add(xe,ell[k] if k<640 else eff[k-640],a)
  for row,a in xe.items():x=add(x,{row*9+8:a})
  # Every642 source column hasY0=0, soM0X1=0 inallnine directions.
  for var in range(9):assert not mv(M,{row:x.get(row*9+var,F(0)) for row in range(56) if x.get(row*9+var,F(0))})
  b=project_linear(Pbottom,[x])[0];allbottom.append(b);v=dict(b);a={};steps=[]
  while v and min(v) in bnorm:
   p=min(v);z=v[p];v=H.guard(add(v,bnorm[p],-z));a=H.guard(add(a,breps[p],z));steps.append([p,str(z)])
  if v:
   p=min(v);z=v[p];k=len(borig);borig.append(b);bselected.append(j);selected_full1.append(x)
   bnorm[p]=H.guard(scale(v,1/z));breps[p]=H.guard(scale(add({k:F(1)},a,-1),1/z));a={k:F(1)}
   breceipts.append({'column':j,'basis_column':k,'pivot':p,'divide_by':str(z),'subtract':steps})
  else:
   # Exactzero-bottom source combination; onlyexistingfirstsareused.
   gamma=add({j:F(1)},{bselected[k]:z for k,z in a.items()},-1);z=add(x,mv(selected_full1,a),-1);orig={}
   for k,c0 in gamma.items():orig=add(orig,rawcomb[k],c0)
   A={}
   for var in range(9):
    coordinates=mv(B30inv,{i:z.get(row*9+var,F(0)) for i,row in enumerate(C) if z.get(row*9+var,F(0))})
    assert all(i<27 for i in coordinates)
    assert mv(Xp,coordinates)=={row:z.get(row*9+var,F(0)) for row in range(56) if z.get(row*9+var,F(0))}
    for i,c0 in coordinates.items():A[i*9+var]=c0
   H.guard(gamma);H.guard(orig);H.guard(A);H.guard(z)
   kernels.append({'source642_combination':srec(gamma),'original768_combination':srec(orig),'actual_a1_27_times9':srec(A),'full56_Xgamma1_equals_Xp0a1':True})
  assert mv(borig,a)==b;bfactors.append(a)
  live=fixed+sum(map(len,allbottom+borig))+sum(map(len,bnorm.values()))+sum(map(len,breps.values()))+sum(map(len,bfactors))+sum(map(len,selected_full1))+sum(len(r['source642_combination'])+len(r['original768_combination'])+len(r['actual_a1_27_times9']) for r in kernels)+len(x)*3+2048
  peak=max(peak,live);assert peak<=150000,('source_profile_live_guard',j,peak)
  if (j+1)%64==0:H.dump(out/'phase.json',{'phase':'known_bottom_first_profile','completed642_columns':j+1,'bottom_rank':len(borig),'peak_slots':peak})
 assert len(kernels)==642-len(borig)
 evidence={'all642_bottom27_columns':[srec(v) for v in allbottom],'bottom_rank':len(borig),'selected_bottom_columns':bselected,
  'fixed_bottom_rank_receipts':breceipts,'normalized_bottom_basis':[[p,srec(v)] for p,v in sorted(bnorm.items())],
  'all642_bottom_original_basis_factors':[srec(v) for v in bfactors],
  'all_zero_bottom_and_Y0_source_combinations':kernels,
  'all642_M0Xgamma1_zero_full9':True,'all_kernel_full56_Xgamma1_Xp0a1_identities_verified':True,
  'known_G3_top_linear_rank':len(gn),'known_bottom_linear_rank':len(borig),
  'remaining_top_linear_quotient_dimension':len(gfree),
  'first_Hilbert_rank_decomposition':'rank_total=rank(bottom)+rank(G3top)+rank(remainingprojectedX2functionals_onzero-bottomkernel)',
  'unqueried_X2_source_formula':'G2linear*Y1=M0*(Xgamma2-Xp1*a1); projecttop3andmodknownG3linear; no M1',
  'actual_HF1_or_T2_length_not_claimed':True,'new_X2_M1_physical_or_implicit_values':False,'peak_slots':peak}
 ep=out/'actual-known-bottom-first-rank-and-corrected-source-kernel.json';H.dump(ep,evidence)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'known_G3_top_linear_rank':len(gn),'known_bottom_linear_rank':len(borig),
  'remaining_top_linear_quotient_dimension':len(gfree),'zero_bottom_kernel_count':len(kernels),
  'actual_HF1_or_T2_length_computed':False,'new_X2_M1_physical_or_implicit_values':False,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{p.name:H.sha(p) for p in [out/'phase1-known-G3linear-fixed-source-rank.json',ep]},'peak_slots':peak,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','known_G3_top_linear_rank','known_bottom_linear_rank','remaining_top_linear_quotient_dimension','zero_bottom_kernel_count','peak_slots','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);a=p.parse_args();main(a.output.resolve(),a.gate)
