#!/usr/bin/env python3
"""Independent known-resolution/actual-X1 first-layer postcheck; no X2."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63'
P={
 'primary':f'reports/{A}/report.json','g_rank':f'reports/{A}/phase1-known-G3linear-fixed-source-rank.json',
 'bottom':f'reports/{A}/actual-known-bottom-first-rank-and-corrected-source-kernel.json',
 'gate_report':f'reports/{G}/report.json','gate':f'reports/{G}/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json',
 'mini':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',
 'mini_peer':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json',
 'mu_rank':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json',
 'mu_peer':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json',
 'mu_gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
 'frame':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
 'active':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
 'G2':'reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json',
 'old640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'old128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'ell640':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def quotient(n,basis):
 free=[i for i in range(n) if i not in basis];ix={i:j for j,i in enumerate(free)};out={}
 for i in reversed(range(n)):
  if i in ix:out[i]={ix[i]:Q(1)}
  else:out[i]=add(*(scale(out[j],-c) for j,c in basis[i].items() if j!=i))
 return free,[out[i] for i in range(n)]
def project(P,cols):
 result=[]
 for col in cols:
  p={}
  for k,c in col.items():
   r,v=divmod(k,9);p=add(p,{i*9+v:a*c for i,a in P[r].items()})
  result.append(p)
 return result
def rank_postcheck(cols,selected,receipts,factors,normalized):
 B=[cols[j] for j in selected];basis={};reprs={}
 for j,r in enumerate(receipts):
  assert r['basis_column']==j and r['column']==selected[j];p=dict(B[j]);w={j:Q(1)}
  for i,c in r['subtract']:
   assert p.get(i)==Q(c);p=add(p,scale(basis[i],-Q(c)));w=add(w,scale(reprs[i],-Q(c)))
  i=r['pivot'];c=Q(r['divide_by']);assert c and min(p)==i and p[i]==c and i not in basis
  basis[i]=scale(p,1/c);reprs[i]=scale(w,1/c);assert mv(B,reprs[i])==basis[i]
 assert [[i,enc(p)] for i,p in sorted(basis.items())]==normalized
 assert len(factors)==len(cols) and all(mv(B,sp(p))==c for p,c in zip(factors,cols))
 return basis
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'KNOWN_G3_AND_BOTTOM_X1_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert all(D[k]['all_passed'] for k in ['primary','gate_report','mini_peer','mu_peer'])
 for name,folder in [('primary',A),('gate_report',G)]:
  for p,h in D[name]['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in D[name]['evidence_sha256'].items():assert sha(ROOT/'reports'/folder/p)==h
 mini=D['mini'];frame=D['frame'];pref=D['active'];gate=D['gate'];active=[0,1,4,5,7];zero=[2,3,6]
 yn={i:sp(p) for i,p in D['mu_rank']['normalized_basis_rows']};freeY,PiY=quotient(102,yn)
 assert freeY==mini['remaining_free3_Y_rows'] and list(map(enc,PiY))==gate['PiY3_from102_columns']
 inv30=list(map(sp,mini['fixed_bottom30_inverse_columns']));Pb=[{} for _ in range(56)]
 for i,r in enumerate(mini['free30_V_indices']):Pb[r]={j-27:c for j,c in inv30[i].items() if j>=27}
 assert list(map(enc,Pb))==gate['PiBottom3_from56_columns']
 fixed=list(map(sp,frame['active_F2_fixed182_columns']));norm={};reps={}
 for r in frame['active_F2_fixed182_receipts']:
  w={r['source_column']:Q(1)}
  for i,c in r['subtract']:w=add(w,scale(reps[i],-Q(c)))
  i=r['pivot'];reps[i]=scale(w,1/Q(r['divide_by']));norm[i]=sp(r['normalized']);assert mv(fixed,reps[i])==norm[i]
 d3=list(map(sp,pref['active_d3_columns']));selected=frame['active126_source_columns'];comp=[i for i in range(364) if i not in selected];I3=[]
 for j in comp:
  p=dict(d3[j]);w={}
  while p:
   i=min(p);c=p[i];w=add(w,scale(reps[i],c));p=add(p,scale(norm[i],-c))
  assert all(i>=56 for i in w)
  z=add({j:Q(1)},{selected[i-56]:-c for i,c in w.items()});assert not mv(d3,z);I3.append(z)
 assert comp==gate['active_F3_complement238'] and list(map(enc,I3))==gate['actual_active_I3_origin_columns']
 rows=[(2-q,w) for q in range(3) for w in itertools.combinations(active,q)];ix={w:i for i,(p,w) in enumerate(rows)};raw=[]
 for p,w in [(3-q,w) for q in range(4) for w in itertools.combinations(active,q)]:
  for j in range(14):
   col={}
   if p:col[(14*ix[w]+j)*9+8]=Q(1)
   for s,i in enumerate(w):col=add(col,{(14*ix[w[:s]+w[s+1:]]+j)*9+i:Q((-1)**(p+s))})
   raw.append(col)
 Pi42=list(map(sp,frame['active_Pi42_columns']));ag3=project(Pi42,[mv(raw,p) for p in I3]);assert list(map(enc,ag3))==gate['actual_active_G3linear42_by238_columns']
 ag2=[{(i//6)*9+(active+[8])[i%6]:c for i,c in sp(p).items()} for p in D['G2']['actual_linear_G2_active102_by42_columns']]
 full=list(ag3)
 for iz,z in enumerate(zero):
  for j,p in enumerate(ag2):
   v={j*9+z:Q(1)}
   for k,c in p.items():r,var=divmod(k,9);v=add(v,{(42+17*iz+r)*9+var:c})
   full.append(v)
 # Derive activeG1 linear from the literal block coefficients of fullG2.
 g2full=list(map(sp,D['mu_gate']['full_G2_linear234_by102_columns']));g1=[]
 for j in range(17):g1.append({(r-17)*9+v:c for k,c in g2full[42+j].items() for r,v in [divmod(k,9)] if 17<=r<20})
 pairs=list(itertools.combinations(range(3),2))
 for ip,(a,b) in enumerate(pairs):
  for j,p in enumerate(g1):
   v={(42+17*b+j)*9+zero[a]:Q(-1),(42+17*a+j)*9+zero[b]:Q(1)}
   for k,c in p.items():r,var=divmod(k,9);v=add(v,{(93+3*ip+r)*9+var:c})
   full.append(v)
 for j in range(3):full.append({(93+3*pairs.index(tuple(k for k in range(3) if k!=i))+j)*9+zero[i]:Q((-1)**i) for i in range(3)})
 assert len(full)==418 and sum(map(len,full))==1010 and list(map(enc,full))==gate['actual_full_G3linear102_by418_columns']
 G3=project(PiY,full);assert sum(map(len,G3))==135 and list(map(enc,G3))==gate['actual_minimal6_G3linear_top27_by418_columns']==D['g_rank']['actual_projected_G3linear_columns']
 mon2=[m for m in itertools.combinations_with_replacement(range(9),2) if m!=(8,8)];mi={m:i for i,m in enumerate(mon2)}
 for col in full:
  image={}
  for k,c in col.items():
   src,v=divmod(k,9)
   for i,a in g2full[src].items():
    r,w=divmod(i,9);m=tuple(sorted((v,w)))
    if m!=(8,8):image=add(image,{r*44+mi[m]:c*a})
  assert not image
 gr=D['g_rank'];gb=rank_postcheck(G3,gr['selected_columns'],gr['fixed_rank_receipts'],gr['all418_original_basis_factors'],gr['normalized_basis']);assert len(gb)==25
 freeTop,PiTop=quotient(27,gb);assert freeTop==gr['remaining_top_linear_quotient_coordinates']==[25,26] and list(map(enc,PiTop))==gr['Pi_remaining_top_linear_columns']
 dump(o,'phase1-independent-G3-source-and-rank25-proof.json',{'active_I3_and_both_unit_derivative_cancellations_verified':True,'all1010_full_G3linear_coefficients_and418_quadratic_chain_identities_checked':True,'all135_projected_coefficients_rank25_receipts_and418factors_verified':True,'remaining_quotient_coordinate_labels':freeTop,'literal_unit_duals_not_assumed':True})
 old=[];orig=[]
 for col in D['old640']['source_columns640_rows56_sparse_dual9']+D['old128']['source_columns128_rows48_sparse_dual9']:
  p={};q={}
  for row,rawp in col:
   for k,c in sp(rawp).items():
    if k:p[row*9+k-1]=c
    else:q[row]=c
  old.append(p);orig.append(q)
 assert sum(map(len,orig))==97
 ell=list(map(sp,D['ell640']['actual_first640_ell_columns']));L=list(map(sp,D['mu_gate']['L0_Res4_Ntranspose_columns']));S0=list(map(sp,D['mu_gate']['native_S0_columns']));S1=list(map(sp,D['mu_gate']['native_S1ell_columns']));eff=[mv(L,p) for p in S1];M=list(map(sp,D['mu_gate']['Mnew0_columns']))
 rawcols=[sp(r['original768_combination']) for r in mini['all642_zero_constant_combinations']];assert gate['all642_constant_Y0_kernel_source_combinations']==mini['all642_zero_constant_combinations']
 X=[];bottom=[]
 for c in rawcols:
  assert not mv(S0,{j-640:a for j,a in c.items() if j>=640}) and not mv(orig,c)
  x=mv(old,c);xe={}
  for j,a in c.items():xe=add(xe,scale(ell[j] if j<640 else eff[j-640],a))
  x=add(x,{row*9+8:a for row,a in xe.items()})
  for v in range(9):assert not mv(M,{row:x[row*9+v] for row in range(56) if row*9+v in x})
  X.append(x);bottom.append(project(Pb,[x])[0])
 B=D['bottom'];assert list(map(enc,bottom))==B['all642_bottom27_columns']
 bn=rank_postcheck(bottom,B['selected_bottom_columns'],B['fixed_bottom_rank_receipts'],B['all642_bottom_original_basis_factors'],B['normalized_bottom_basis']);assert len(bn)==23
 nonpivot=[j for j in range(642) if j not in B['selected_bottom_columns']];assert len(nonpivot)==619==len(B['all_zero_bottom_and_Y0_source_combinations'])
 Xp=[orig[j] for j in mini['original27_X_source_columns']];count_raw=0;count_a=0;maxraw=0;maxa=0
 for j,record in zip(nonpivot,B['all_zero_bottom_and_Y0_source_combinations']):
  fac=sp(B['all642_bottom_original_basis_factors'][j]);combo=add({j:Q(1)},{B['selected_bottom_columns'][k]:-a for k,a in fac.items()})
  assert enc(combo)==record['source642_combination'] and not mv(bottom,combo)
  rawc=mv(rawcols,combo);z=mv(X,combo);a1=sp(record['actual_a1_27_times9'])
  assert enc(rawc)==record['original768_combination']
  for v in range(9):
   given={i:a1[i*9+v] for i in range(27) if i*9+v in a1};actual={r:z[r*9+v] for r in range(56) if r*9+v in z}
   assert mv(Xp,given)==actual
  count_raw+=len(rawc);count_a+=len(a1);maxraw=max(maxraw,len(rawc));maxa=max(maxa,len(a1))
 assert (count_raw,maxraw,count_a,maxa)==(992,35,1887,89)
 assert D['primary']['peak_slots']==B['peak_slots']==24500<98958<150000
 ep=o/'independent-known-first-source-ranks25-23-and619-corrections.json'
 dump(o,ep.name,{'inputs':inputs,'full_G3linear_and_projection_source_derived':True,'all418_G2linearG3linear_identities_modell2':True,'known_G3_top_rank':25,'all642_bottom_first_columns_source_derived_and_rank':23,'all619_kernel_originalsource_and_full56_X1_Xp0a1_checked':True,'kernel_original_support_total_max':[992,35],'a1_support_total_max':[1887,89],'remaining_top_linear_quotient_dimension':2,'actual_HF1_interval_from_complete_presentation':[4,6],'exact_HF1_not_computed':True,'M1_drop_after_actual_Xp_correction_proved':True,'projected_X2_only_sufficiency_still_requires':'exact adjoint lambda*G2linear=psi on all102*9 linear coefficient slots; formalexactnessalone does not establish associatedgraded exactness','resource_scope':'the earlier perkernel40 support bound was a forecast, not a proveduniversalbound; actualmax35 and actualhardguard24500 verified withoutcapextension','new_X2_M1_physical_or_implicit_values_computed':False})
 report={'id':STEM,'all_passed':True,'known_G3_top25_bottom23_and619_actual_source_corrections_accepted':True,'actual_HF1_interval':[4,6],'exact_HF1_T2_length_or_X2_values_computed':False,'projected_X2_only_method_pending_exact_rowspan_certificate':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep),'phase1-independent-G3-source-and-rank25-proof.json':sha(o/'phase1-independent-G3-source-and-rank25-proof.json')},'seconds':time.monotonic()-t}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
