#!/usr/bin/env python3
"""Independent ONE508 source/ring/eligibility audit, no function responses."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63';STEM=A+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 return int(isinstance(v,(int,float)) or isinstance(v,str) and v.lstrip('-').replace('/','').isdigit())
def fields(p,names):
 txt=p.read_text();dec=json.JSONDecoder();i=txt.index('{')+1;out={}
 while True:
  while txt[i].isspace() or txt[i]==',':i+=1
  if txt[i]=='}':break
  k,i=dec.raw_decode(txt,i)
  while txt[i].isspace() or txt[i]==':':i+=1
  v,i=dec.raw_decode(txt,i)
  if k in names:out[k]=v
  del v
 assert set(out)==set(names),(p,set(names)-set(out));return out
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;assert raw<32*1024*1024
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 r=read(f'reports/{A}/report.json');assert r['all_passed'] and pins[f'reports/{A}/report.json']=='7b2bd2588c1d7f1cb30a8c6c7822a80d1a770417cfd114d324531f4925a93c22'
 ep='ONE508-source-closure-B21-and-exact-carry-zeros.json';d=read(f'reports/{A}/{ep}');assert pins[f'reports/{A}/{ep}']==r['evidence_sha256'][ep]
 src='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json'
 s=read(src,['original619_c0_source_combinations','kernel617_source_recipes','typed_program','D_old_root_columns','D_two_by10_by56','pair_order10'])
 a=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'][508]
 assert a['original768_combination']==s['original619_c0_source_combinations'][508]==[[595,'1']] and a['actual_a1_27_times9']==[]
 kr=next(x for x in s['kernel617_source_recipes'] if x['free_case619']==508);assert kr==d['same617_source_recipe'] and kr['subtract_actual_basis563609']==['0','0']
 hp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json'
 h=read(hp,['typed_source_program','literal_C0_program','literal_boundary_firsts_point_r_ell_and_membership','pruned_original_old_columns_roots','Gamma_two'])
 t=dict(h['typed_source_program']);u=dict(s['typed_program']);assert all(t[i]==n for i,n in u.items() if i in t);t.update(u)
 c=dict(h['literal_C0_program']);boundary={i:k for i,k,v in h['literal_boundary_firsts_point_r_ell_and_membership']}
 dr=dict(s['D_old_root_columns'])[595];gr=dict(h['pruned_original_old_columns_roots'])[595]
 assert dr==d['D_old595_source_roots'] and gr==d['Gamma_old595_source_roots'] and h['Gamma_two']==d['Gamma_two'] and s['D_two_by10_by56']==d['D_two_by10_by56']
 todo=[i for row,i in dr]+[i for row,(kind,i) in gr if kind=='typed'];Ctodo={i for row,(kind,i) in gr if kind=='C0'};Tseen=set()
 while todo:
  i=todo.pop()
  if i in Tseen:continue
  Tseen.add(i);n=t[i]
  if n[0]=='actual_C0_scalar_function':Ctodo.add(n[1])
  elif n[0]=='sum':todo.extend(j for j,a,b in n[1])
  else:assert n[0]=='product';todo.extend(n[3])
 Cseen=set();leaves=set();todo=list(Ctodo)
 while todo:
  i=todo.pop()
  if i in Cseen or i in leaves:continue
  if i not in c:assert i in boundary;leaves.add(i);continue
  Cseen.add(i);n=c[i]
  if n[0] in ['a','m']:todo.extend(n[1:3])
  elif n[0]=='n':todo.append(n[1])
  else:assert n[0]=='c'
 assert [[i,t[i]] for i in sorted(Tseen)]==d['full_typed_program86'] and len(Tseen)==86
 assert [[i,c[i]] for i in sorted(Cseen)]==d['literal_C0_program38'] and len(Cseen)==38 and sorted(leaves)==d['function_leaf39_ids'] and len(leaves)==39
 # Known source quotient maps only. All660 nonhomogeneous NF columns are
 # projected, and all1089 products are checked directly using the original
 # raw33 multiplication table (not the producer's monomial-NF lookup).
 nf=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json',['monomials660_degree_le4','all660_actual_raw33_normal_forms'])
 tab=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json',['basis33','table_columns_for_each_left_basis','m5_zero'])
 owner0=count([r,d,s,a,h,t,u,c,boundary,dr,gr,Tseen,Cseen,leaves,nf,tab,pins])+20000
 assert owner0<150000 and raw+1_000_000<32*1024*1024
 dump('phase-before-ring-checks.json',{'phase':'literal_source_eligibility_fullgraph_owner_PASS','source_owner_upper':owner0,'copy_selector_only_no_signed_function_response':True})
 b=[tuple(v) for v in d['ordinary_basis_B21']];bi={v:i for i,v in enumerate(b)}
 assert len(bi)==21 and b[0]==(0,)*9 and b[1:10]==[tuple(int(i==j) for i in range(9)) for j in range(9)]
 assert b[10:20]==[tuple(int(i==a)+int(i==q) for i in range(9)) for a,q in s['pair_order10']]
 assert b[20]==(0,0,0,0,0,0,0,2,1)
 for m in b:
  for i,a0 in enumerate(m):
   if a0:assert tuple(v-int(k==i) for k,v in enumerate(m)) in bi
 rb=[]
 for z,r0,e in tab['basis33']:
  m=[0]*9
  if z:m[z-1]=1
  m[7]=r0;m[8]=e;rb.append(tuple(m))
 images=[bi.get(m) for m in rb];sel=[images.index(i) for i in range(21)];assert sel==d['raw33_to_B21_selector21']
 for m,v in zip(nf['monomials660_degree_le4'],nf['all660_actual_raw33_normal_forms'],strict=True):
  expected={} if tuple(m) not in bi else {bi[tuple(m)]:'1'}
  got={images[i]:q for i,q in v if images[i] is not None and q!='0'};assert got==expected
 for i in range(33):
  for j in range(33):
   got={images[k]:q for k,q in tab['table_columns_for_each_left_basis'][i][j] if images[k] is not None and q!='0'}
   m=tuple(a+b for a,b in zip(rb[i],rb[j]));ix=bi.get(m) if images[i] is not None and images[j] is not None else None
   assert got==({ix:'1'} if ix is not None else {}) and ix==d['all1089_product_target_indices_or_null'][i][j]
 # Full fiber-product shape, not a naive unrelated list of Taylor callbacks.
 inactive={2,3,6};assert all(tuple(int(i==j)+v for i,v in enumerate(m)) not in bi for j in inactive for m in b[1:])
 p6=[v for v in b if all(v[i]==0 for i in range(7))];assert len(p6)==6 and (0,0,0,0,0,0,0,3,0) not in bi
 # Scope proofs and literal formal support: no actual second/third response.
 for p0 in ['reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json']:assert read(p0)['all_passed']
 j1=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63/617-formal-direct-old-third-and-second-symbol-sensitivities.json',['all617_formal_directOLD_J0_Jr_Jell'])['all617_formal_directOLD_J0_Jr_Jell']
 rows=next(v['two_formal_J0_Jr_Jell_symbol_rows'] for v in j1 if v['free_case619']==508)
 for row in rows:
  for f,co in row:
   if boundary[f]!='source164':assert f in [299650,505938]
   if f in [299650,505938]:assert co[0]=='0' # actual thirds not needed in this case.
 old2=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json',['all617_combined_OLD2_independent_symbol_coefficients'])['all617_combined_OLD2_independent_symbol_coefficients']
 rr=next(ch for j,ch in old2 if j==508);assert {f for ch in rr for pair,f,q in ch}=={20144,40922,299650,505938}
 # All source primary input hashes checked without reading or evaluating their
 # actual physical function response arrays. Hash bytes are charged explicitly.
 for p0,ha in r['input_sha256'].items():
  if p0 not in pins:raw+=(ROOT/p0).stat().st_size;pins[p0]=sha(ROOT/p0)
  assert pins[p0]==ha,p0
 owner=count([r,d,s,a,h,t,u,c,boundary,dr,gr,Tseen,Cseen,leaves,nf,tab,b,bi,rb,images,j1,rows,old2,pins])+20000
 assert owner<150000 and raw+1_000_000<32*1024*1024
 ep2='independent-ONE508-eligibility124graph-B21-and-whole-extension-scope.json'
 ev={'single_originalOLD595_fulla1empty_nativeempty_basis_subtraction0':True,'full86typed38C0nodes39leaves_retained':True,
  'all660_nonhomogeneous_NF_and_direct_original1089_table_products':True,'B21_divisorclosed_and_fiberproduct_inactive_firsts_annihilated':True,
  'whole_fullPsi_formula_GammaX595third_plus_DX595second_not_pureU':True,
  'all_genuine_honest_higher_symbols_total_independence_only':True,'actual_conductor_seconds_order3zero_thirds_coefficient0_in_case508':True,
  'all_first_cross_products_retained_and_no_actual_higher_jet_zero_claim':True,
  'no_function_response_or_upper_value_evaluated':True,'whole_source_owned_upper':owner,'whole_RAW_upper':raw+1_000_000}
 dump(ep2,ev);rep={'id':STEM,'all_passed':True,'source_only':True,'selected_case508_only':True,'input_sha256':pins,'evidence_sha256':{ep2:sha(out/ep2)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'whole_source_owned_upper':owner,'whole_RAW_upper':raw+1_000_000,'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owner,'RAW':raw+1_000_000,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
