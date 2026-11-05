#!/usr/bin/env python3
"""ONE508 source/finite-ring input binding, no source function response."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return H.selected(p,k)
 sourcep='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json';s=read(sourcep)
 a1p='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json';rec=field(a1p,'all_zero_bottom_and_Y0_source_combinations')[508]
 assert rec['original768_combination']==[[595,'1']] and not rec['actual_a1_27_times9']
 kr=next(v for v in s['kernel617_source_recipes'] if v['free_case619']==508);assert kr['subtract_actual_basis563609']==['0','0']
 assert s['original619_c0_source_combinations'][508]==[[595,'1']]
 hpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json';h=read(hpath)
 t=dict(h['typed_source_program']);u=dict(s['typed_program']);assert all(t[i]==n for i,n in u.items() if i in t);t.update(u);c=dict(h['literal_C0_program']);boundaries={i for i,k,v in h['literal_boundary_firsts_point_r_ell_and_membership']}
 Droots=dict(s['D_old_root_columns'])[595];Groots=dict(h['pruned_original_old_columns_roots'])[595]
 tr={i for row,i in Droots};cr=set()
 for row,(kind,i) in Groots:
  assert kind in ['typed','C0'];(tr if kind=='typed' else cr).add(i)
 keep=set();todo=list(tr)
 while todo:
  i=todo.pop()
  if i in keep:continue
  keep.add(i);n=t[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0]=='actual_C0_scalar_function';cr.add(n[1])
 kc=set();leaf=set();todo=list(cr)
 while todo:
  i=todo.pop()
  if i in kc or i in leaf:continue
  if i not in c:assert i in boundaries;leaf.add(i);continue
  kc.add(i);n=c[i]
  if n[0] in ['a','m']:todo.extend(n[1:3])
  elif n[0]=='n':todo.extend(n[1:2])
  else:assert n[0]=='c'
 assert len(keep)==86 and len(kc)==38 and len(leaf)==39
 # B21 contains all nine first coordinates, the actual ten nonzz pairs, and
 # the single plane cubic r²ell. It is a monomial quotient, not a collection
 # of unrelated scalar Taylor projections.
 z=(0,)*9;units=[tuple(int(i==j) for i in range(9)) for j in range(9)]
 pairs=[tuple(v) for v in s['pair_order10']]
 b=[z]+units+[tuple(int(i==a)+int(i==q) for i in range(9)) for a,q in pairs]+[tuple(2 if i==7 else 1 if i==8 else 0 for i in range(9))]
 assert len(set(b))==21;bi={v:i for i,v in enumerate(b)}
 rawp='reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json';raw33=read(rawp)
 fb='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';nf=read(fb)
 rp=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json');assert rp['all_passed'] and rp['actual_product_m5_zero']
 source_prearith=count([s,h,raw33,nf,pins])+60000;assert source_prearith<150000
 rb=[]
 for sector,r,e in raw33['actual33_basis']:
  m=[0]*9
  if sector:m[sector-1]=1
  m[7]=r;m[8]=e;rb.append(tuple(m))
 images=[bi.get(m) for m in rb];select=[images.index(i) for i in range(21)]
 assert select==[0,7,6,5,4,3,2,1,8,9,22,23,20,21,14,15,12,13,24,25,30]
 mon=[tuple(m) for m in nf['monomials660_degree_le4']];norm=[dict(v) for v in nf['all660_actual_raw33_normal_forms']];assert len(mon)==len(norm)==660
 def project(v):return {j:v[i] for j,i in enumerate(select) if i in v and v[i]!='0'}
 for m,v in zip(mon,norm,strict=True):assert project(v)==({bi[m]:'1'} if m in bi else {})
 lookup=dict(zip(mon,norm));products=[]
 for i,a in enumerate(rb):
  row=[]
  for j,q in enumerate(rb):
   m=tuple(x+y for x,y in zip(a,q,strict=True));rawimage={} if sum(m)>=5 or m[8]>=2 else project(lookup[m])
   outidx=bi.get(m) if images[i] is not None and images[j] is not None else None
   assert rawimage==({outidx:'1'} if outidx is not None else {})
   row.append(outidx)
  products.append(row)
 # Scope source proofs, without reading any function-second/third values.
 for p in ['reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json']:
  assert read(p)['all_passed']
 dform=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json');v=dict(dform['all617_combined_OLD2_independent_symbol_coefficients'])[508];assert {f for ch in v for p,f,q in ch}=={20144,40922,299650,505938}
 # all96 firsts are metadata/source labels here; values stay in later input.
 ev={'selected_case619':508,'literal_original768_source':[[595,'1']],'a1_all27times9_exact_empty':True,'native_c0_exact_empty':True,'same617_source_recipe':kr,'complete_actual_upper_recipe':'Gamma X595_3(r²ell) +D X595_2; allXp2*a1 andXp1*a1 exact0 fromfulla1=[], nativec0=0; fullPsi alreadyinD/Gamma','physical_variable_order9':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],'ordinary_basis_B21':b,'raw33_to_B21_selector21':select,'all660_NF_projection_identities_verified':True,'all1089_raw33_products_preserved':True,'basis33_images_index_B21_or_null':images,'all1089_product_target_indices_or_null':products,'full_typed_program86':[[i,t[i]] for i in sorted(keep)],'literal_C0_program38':[[i,c[i]] for i in sorted(kc)],'function_leaf39_ids':sorted(leaf),'D_old595_source_roots':Droots,'Gamma_old595_source_roots':Groots,'Gamma_two':h['Gamma_two'],'D_two_by10_by56':s['D_two_by10_by56'],'D_second_functional508_only_physical20144_40922_and_order3conductors':True,'source_IFT_firsts_and_actual_physical_higher_not_extracted_here':True,'future_ordinary_gA3_source':'505938=CAprefix0[7790]; raw33->planeP7 preserves r²ell, sameCA-F9zero7 source','future_ordinary_gB3_source':'299650=fCB28539-sumzB_i*h_i; zB0/plane1=0,h0=0; gB_r²ell=fCB30-sum(zB24*h_ell+zB25*h_r)','future_honest_higher_extension':'all only79 honest seconds/thirds may be arbitrarily0 ONLYinside this whole projectedfunctional by separatelypeered all617J1 andOLD2 certificates; actualfunctionjet0 isNOTclaimed','all_first_cross_products_retained_by_full124_nodes':True,'no_response_gradient_or_upper_value_evaluated':True}
 ev['source_fiber_product_proof']='P6=P7/(r³); B18=B17 times_P5 P6, dimensions17+6-5=18; B21=B18 times_B7 T10, dimensions18+10-7=21. Threeinactivefirsts are annihilatedbythemaximalideal. LiteralNF/productchecksbindthecommonraw33map.'
 ep=out/'ONE508-source-closure-B21-and-exact-carry-zeros.json';H.dump(ep,ev)
 owned=count([s,h,raw33,nf,keep,kc,leaf,mon,norm,lookup,products,ev,dform,pins])+10000;assert owned<150000
 report={'id':STEM,'all_passed':True,'source_known_ring_identities_only':True,'selected_case508_only':True,'typed_nodes86_C0nodes38_leaves39':True,'all1089_ring_and660_NF_checks_passed':True,'no_native_or_pivot_response_shortcut_all_a1_exact0':True,'source_owned_bound':owned,'whole_RAW_reads_upper':raw,'no_actual_function_response_or_upper_value':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 report['source_prearith_bound']=source_prearith
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_owned_bound','whole_RAW_reads_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
