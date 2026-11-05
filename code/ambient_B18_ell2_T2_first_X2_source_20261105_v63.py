#!/usr/bin/env python3
"""Coherent rell source folds for the two T2 first-Hilbert callbacks.

Second coefficients remain formal. Native individual effective-X shortcuts
are forbidden; L1ell*Omega_r is retained explicitly.
"""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63'
def parse(v):return H.parse(v)
def dot(a,b):return sum((c*b.get(i,F(0)) for i,c in a.items()),F(0))
def keyrec(v):return [[list(k),str(c)] for k,c in sorted(v.items()) if c]
def put(d,k,c):
 a=d.get(k,F(0))+c
 if a:d[k]=a
 else:d.pop(k,None)
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'source_fold_started_no_second_value_queries'})
 hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 adj=read('reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63/actual-two-full918-source-adjoint-certificates-or-obstructions.json')
 read('reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/report.json')
 known=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json')
 read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json')
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json')
 gate=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json')
 M=list(map(parse,gate['Mnew0_columns']));Gamma=[]
 for result in adj['results']:
  assert result['exact_rr_rell_adjoint_exists'];lam=parse(result['lambda52']);assert all(i%2==1 for i in lam)
  Gamma.append({j:sum((lam.get(2*i+1,F(0))*c for i,c in row.items()),F(0)) for j,row in enumerate(M) if sum((lam.get(2*i+1,F(0))*c for i,c in row.items()),F(0))})
 assert all(sorted(g)==[1,2,4,17,18,21,22,28,55] for g in Gamma)
 chain=read(H.PATHS['chain']);ops=chain['operator_scalar_nodes'];B=[dict(v) for v in chain['B48_by640_sparse_columns']]
 hs[H.PATHS['C0']]=H.sha(ROOT/H.PATHS['C0']);P=H.selected(H.PATHS['C0'],'honest_D0_relation_columns')
 hs[H.PATHS['old640']]=H.sha(ROOT/H.PATHS['old640']);old=H.selected(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9');oldr=[]
 for col in old:oldr.append({row:parse(v).get(8,F(0)) for row,v in col if parse(v).get(8,F(0))})
 ell=list(map(parse,read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json')['actual_first640_ell_columns']))
 hs[H.PATHS['G0']]=H.sha(ROOT/H.PATHS['G0']);G=H.selected(H.PATHS['G0'],'G8');hs[H.PATHS['S']]=H.sha(ROOT/H.PATHS['S']);S=H.selected(H.PATHS['S'],'native_S12_by32')
 S0=list(map(parse,gate['native_S0_columns']));Se=list(map(parse,gate['native_S1ell_columns']));Sr=[{} for _ in range(128)]
 for row,col,k,sign in S['koszul28_in_first8_rows']:
  for inp in range(4):
   for mask,c in enumerate(map(F,G[k][32:36])):
    if c and not(mask&inp):Sr[4*col+inp]=add(Sr[4*col+inp],{4*row+(mask|inp):F(sign)*c})
 L=list(map(parse,gate['L0_Res4_Ntranspose_columns']));piv=mini['original27_X_source_columns'];S0p=[S0[j-640] if j>=640 else {} for j in piv];Srp=[Sr[j-640] if j>=640 else {} for j in piv];Sep=[Se[j-640] if j>=640 else {} for j in piv]
 GL0=[]
 for gamma in Gamma:GL0.append({j:dot(gamma,v) for j,v in enumerate(L) if dot(gamma,v)})
 folds=[];neededroots=set();total=0
 for number,record in enumerate(known['all_zero_bottom_and_Y0_source_combinations']):
  c=parse(record['original768_combination']);A=parse(record['actual_a1_27_times9']);ar={i:A.get(i*9+7,F(0)) for i in range(27) if A.get(i*9+7,F(0))};ae={i:A.get(i*9+8,F(0)) for i in range(27) if A.get(i*9+8,F(0))};cn={j-640:v for j,v in c.items() if j>=640};assert not mv(S0,cn)
  omr=add(mv(Sr,cn),mv(S0p,ar),-1);ome=add(mv(Se,cn),mv(S0p,ae),-1)
  intrinsic_correction=add(mv(Srp,ae),mv(Sep,ar))
  fs=[];offset=[]
  for channel,gamma in enumerate(Gamma):
   f={};k=F(0)
   for j,v in c.items():
    if j>=640:continue
    for row,a in gamma.items():
     if row<48:
      root=B[j].get(row)
      if root is None:continue
      n=ops[root]
      if n[0]=='c':continue
      if n[0]=='actual_C0_scalar_function':key=(1,n[1])
      else:key=(0,root);neededroots.add(root)
     else:key=(1,P[j][row-48+1])
     put(f,key,a*v)
   for i,j in enumerate(piv):
    if j<640:k-=dot(gamma,oldr[j])*ae.get(i,F(0))+dot(gamma,ell[j])*ar.get(i,F(0))
   k-=dot(GL0[channel],intrinsic_correction);H.guard({0:k} if k else {})
   fs.append(f);offset.append(k)
  folds.append({'index':number,'original768_combination':record['original768_combination'],'a1_r':srec(ar),'a1_ell':srec(ae),
   'two_old_typedB_or_originalP_rell_root_functionals':[keyrec(v) for v in fs],
   'native_intrinsic_S2_source_combination':srec(cn),'native_Omega_r':srec(omr),'native_Omega_ell':srec(ome),
   'two_known_oldpivot_and_L0S1_corrections':list(map(str,offset))})
  total+=sum(map(len,fs))+len(cn)+len(omr)+len(ome)+len(ar)+len(ae)+sum(bool(v) for v in offset)
  assert total+20000<150000
 closure=H.closure(ops,neededroots);atoms=sorted({ops[i][1] for i in closure if ops[i][0]=='actual_C0_scalar_function'})
 forecast=total+6*len(closure)+2*len(closure)*max(1,len(atoms))+2*619*max(1,len(atoms))+20000
 data={'source_pair_only':['rA','ell'],'two_Gamma_on56_actual_X2':[srec(v) for v in Gamma],
  'two_Gamma_L0_on48_intrinsic_native_S2':[srec(v) for v in GL0],
  'all619_source_folds':folds,'old_root_kind0_typedB_kind1_originalC0':True,
  'needed_typed_closure_nodes':[[i,ops[i]] for i in sorted(closure)],'needed_original_oldB_atom_labels':atoms,
  'native_formula':'GammaL0*S2rell*cN + GammaL1r*Omega_ell + GammaL1ell*Omega_r minusGammaL0*(Sp1r*a1ell+Sp1ell*a1r); noL2sinceS0cN0',
  'native_individual_effective_Xp_firsts_not_used':True,
  'all619_intrinsic_nativeS0_zero_verified':True,
  'all_second_source_values_remain_formal':True,
  'candidate_cached_native_intrinsic_second':'MScanonical30Artin33; requires exactsourceframeandMarkpeer; allsecondrellcoeffsdescendonlyafterordinarydegree2functionalcheck',
  'candidate_cached_L1ell':'reconstructfullfirst5telescope/R12/T/N fromsavedcanonicalfirstell andphysical105firstell; noimplicitnewsolve; notyetcomputed',
  'old_AD_scope':'preserveallhonestsecondfunctionsandcoherentnegativealiases; no newordinary/C0gradients invented',
  'folded_coefficient_count':total,'next_typed_AD_forecast_upper':forecast,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384},
  'new_physical_gauge_or_implicit_coefficients':False}
 ep=out/'actual-two-rell-source-folds-and-retained-native-first-operators.json';H.dump(ep,data)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'all619_source_functionals_retained':True,'actual_second_values_computed':False,
  'folded_coefficient_count':total,'old_typed_closure_nodes':len(closure),'old_atom_count':len(atoms),
  'next_typed_AD_forecast_upper':forecast,'native_Omega_r_nonzero_count':sum(bool(r['native_Omega_r']) for r in folds),
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','folded_coefficient_count','old_typed_closure_nodes','old_atom_count','next_typed_AD_forecast_upper','native_Omega_r_nonzero_count','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
