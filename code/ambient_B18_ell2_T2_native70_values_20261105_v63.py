#!/usr/bin/env python3
"""Prepared witness-first native70 lower replay from fixed finite caches only."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63'
def count(x):
 if isinstance(x,(int,Q)):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return 0
def main(out,peer_path,custody_path):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);pins={};computed=[];cursor={};maximum=0
 def read(p,keys=None):pins[p]=sha(ROOT/p);return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 gb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63';gp=gb+'/report.json';g=read(gp);assert g['all_passed'];ip=gb+'/complete70-cached-input-source-and-positive-whole-image-bounds.json';d=read(ip);assert pins[ip]==g['evidence_sha256'][Path(ip).name]
 peer=read(peer_path);assert peer['all_passed'] and pins[gp] in json.dumps(peer)
 cb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_LIVE_BOUND_20261105_v63';cr=read(cb+'/report.json');ce=read(cb+'/explicit-final-input-and-worker-overlap-count.json');assert cr['all_passed'] and cr['source_gate_sha256']==pins[gp] and cr['source_input_sha256']==pins[ip] and cr['actual_worker_owned_bound']==49734<51263
 custody=read(custody_path);assert custody['all_passed'] and custody['safe_help_literal_PASS'] and custody['worker_sha256']==sha(Path(__file__)) and custody['source_gate_sha256']==pins[gp]
 for p,h in custody['worker_and_helper_sha256'].items():assert sha(ROOT/p)==h
 for p,h in custody['literal_external_source_sha256'].items():assert sha(ROOT/p)==h
 for z in g['input_sha256'].values():assert sha(ROOT/z['path'])==z['sha256']
 assert count([g,peer,cr,ce,custody,pins])<900 and d['whole_source_committed_bound']==128 and d['whole_source_temporary_bound']==385 and count(d)==15878
 def dump(name,v):
  (out/name).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<d['future_whole_RAW_bound']==8305349
 dump('phase.json',{'phase':'authorized_cache_only_native70_before_any_values','source_gate_sha256':pins[gp],'source_peer_sha256':pins[peer_path],'live_bound_sha256':pins[cb+'/report.json'],'no_physical_gauge_or_implicit_response_solve':True})
 def ck(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length())
  if maximum>128:
   dump('guard-frontier.json',{'cursor':cursor,'completed_prefix':computed,'offending_uncommitted_value':str(x),'height':maximum,'no_retry_or_cap_raise':True});raise RuntimeError('source128guard')
  return x
 def vadd(a,b):return [ck(x+y) for x,y in zip(a,b,strict=True)]
 def vscale(a,c):return [ck(x*c) for x in a]
 table={(i,j):[(k,Q(c)) for k,c in v] for i,j,v in d['C72_multiplication_records']}
 def vmul(a,b):
  v=[Q(0)]*72
  for i,x in enumerate(a):
   if not x:continue
   for j,y in enumerate(b):
    if not y:continue
    for k,c in table.get((i,j),[]):v[k]=ck(v[k]+ck(ck(x*y)*c))
  return v
 one=[Q(0)]*72;one[0]=Q(1);kap=[list(map(Q,p)) for p in d['actual_kappa9_C72']];prog=d['regular_G8_source_program'];core=[kap[i] for i in prog['core_kappa_indices4']];powers=[]
 for i,a in enumerate(core):
  cursor.update(stage='finite_G8_core_powers',core=i);p=[one,a]
  for e in range(2,7):p.append(vmul(p[-1],a))
  assert not any(vmul(p[-1],a));powers.append(p)
 def poly(terms):
  z=[Q(0)]*72
  for exp,n,den in terms:
   if sum(exp)>=7:continue
   v=one
   for i,e in enumerate(exp):
    if e:v=vmul(v,powers[i][e])
   z=vadd(z,vscale(v,Q(n,den)))
  return z
 def inv(k,P):
  if k==0:return one
  z=[Q(0)]*72
  for e in range(4):z=vadd(z,vscale(powers[1][2*e],Q((-1)**e*math.comb(k+e-1,e),5**(k+e)) if P else Q(math.comb(k+e-1,e))))
  return z
 G=[]
 for name,idx in zip(prog['five_pivot_order'],prog['five_pivot_kappa_indices'],strict=True):
  cursor.update(stage='finite_regular_G8',root=name);r=prog['five_actual_rational_roots'][name];root=vmul(vmul(poly(r['numerator']),inv(r['D_power'],False)),inv(r['P_power'],True));G.append(vadd(kap[idx],vscale(root,Q(-1))))
 G.extend(poly(p) for p in prog['regular_CI_generators3']);assert all(v[:40]==list(map(Q,w)) for v,w in zip(G,d['actual_G8_point_and_first40'],strict=True))
 ellp='reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json';ell=read(ellp,['actualG8'])['actualG8'];assert all([v[4*i+m] for i in [0,8,9,11] for m in range(4)]==list(map(Q,w)) for v,w in zip(G,ell,strict=True))
 opath='reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json';op=read(opath,['native_S12_by32'])['native_S12_by32'];defs={j:[] for j in range(28)}
 for row,col,gi,sgn in op['koszul28_in_first8_rows']:defs[col].append((row,gi,Q(sgn)))
 nb=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_C72_POSITIVE_SOURCE_20261105_v63/report.json');assert nb['all_passed']
 assert nb['input_sha256'][opath]['sha256']==pins[opath]
 # No actual N2 is constructed: this exact N0/N1 extension is used only on
 # the certified weighted intrinsic point kernels.
 N={(i,j):list(map(Q,p))+[Q(0)]*32 for i,j,p in d['actual_needed_N_point_and_first40']}
 dump('phase1-cached-C72-regular-G8-and-source-overlaps.json',{'actual_regular_G8_C72':list(map(lambda p:list(map(str,p)),G)),'all_G8_point_and_nine_first_coefficients_match':True,'all_C16_G8_coefficients_including_rell_match':True,'no_new_source_or_gauge_solve':True,'N_second_extension_only_for_weighted_kernel':True})
 pairs=list(map(tuple,d['ordinary_pair_order8']));first={i:sp(v) for i,v in d['actual_C0_leaf_first10']};second={i:list(map(Q,p)) for i,p in d['actual_13_honest_second8']};ops=dict(d['typed_closure_nodes']);jets={}
 def jm(a,b):
  v=[ck(a[0]*b[0])]+[ck(ck(a[0]*b[k])+ck(b[0]*a[k])) for k in range(1,10)]
  for j,(u,w) in enumerate(pairs):
   x=ck(ck(a[0]*b[10+j])+ck(b[0]*a[10+j]));x=ck(x+ck(a[u+1]*b[w+1]))
   if u!=w:x=ck(x+ck(a[w+1]*b[u+1]))
   v.append(x)
  return v
 for k,n in sorted(ops.items()):
  if n[0]=='actual_C0_scalar_function':v=[first[n[1]].get(i,Q(0)) for i in range(10)]+second.get(n[1],[Q(0)]*8)
  elif n[0]=='c':v=[Q(n[1],n[2])]+[Q(0)]*17
  elif n[0]=='sum':
   v=[Q(0)]*18
   for j,a,b in n[1]:v=[ck(x+ck(Q(a,b)*y)) for x,y in zip(v,jets[j],strict=True)]
  else:
   assert n[0]=='product';v=[Q(n[1],n[2])]+[Q(0)]*17
   for child in n[3]:v=jm(v,jets[child])
  jets[k]=v
 Gamma=list(map(sp,d['Gamma4_source56_times8']));B={j:dict(v) for j,v in d['oldB_source_rows_pruned']};Xp={p:list(map(sp,v)) for p,v in d['actual_pivot_old8_first_pruned']};E={p:sp(v) for p,v in d['actual_old_pivot_ell_first_pruned']};L0=list(map(sp,d['native_L0_rows_pruned']));Le=list(map(sp,d['native_L1ell_rows_pruned']));S0={p:sp(v) for p,v in d['native_pivot_S0_columns']};Se={p:sp(v) for p,v in d['native_pivot_S1ell_columns']}
 for p,j in d['native_pivot_local128_ids']:
  v={}
  for A,b in [(L0,Se[p]),(Le,S0[p])]:
   for col,c in b.items():
    for row,x in A[col].items():
     z=ck(v.get(row,Q(0))+ck(c*x))
     if z:v[row]=z
     else:v.pop(row,None)
  E[p]=v
 witness=None
 for rec,bounds in zip(d['source70_records'],d['all70_four_component_positive_envelopes'],strict=True):
  case=rec['case619'];assert case==bounds['case619'];cursor={'stage':'actual_native70_witness_scan','case619':case,'completed':len(computed)};cn=sp(rec['native128_constant_projection']);assert max(cn)<112
  S=[[Q(0)]*72 for _ in range(8)]
  for idx,c in cn.items():
   col,mask=divmod(idx,4);tau=[Q(0)]*72;tau[mask]=c
   for row,gi,sign in defs[col]:S[row]=vadd(S[row],vscale(vmul(tau,G[gi]),sign))
  assert all(not any(p[:4]) for p in S) and rec['all48_intrinsic_native_S0_components_exact_zero']
  X={g:[Q(0)]*72 for g in [0,4,5,7]}
  for g in X:
   for i in range(8):X[g]=vadd(X[g],vmul(N[i,g],S[i]))
  old={i:Q(c) for i,c in rec['original768_constant_combination'] if i<640};a1=sp(rec['actual_a1_27_times9']);values=[];pieces=[]
  for channel,row in enumerate(Gamma):
   oldv=Q(0);nativev=Q(0);carry=Q(0)
   for slot,c in row.items():
    source,pair=divmod(slot,8);u,v=pairs[pair];g,mask=divmod(source,4);nativev=ck(nativev+ck(c*X[g][4*(10+pair)+mask]))
    for j,c0 in old.items():
     root=B[j].get(source)
     if root is not None:oldv=ck(oldv+ck(c*ck(c0*jets[root][10+pair])))
    for p in d['needed_pivot_indices']:
     au=a1.get(9*p+u,Q(0));av=a1.get(9*p+v,Q(0))
     if au:
      firstv=Xp[p][v].get(source,Q(0)) if v<8 else E[p].get(source,Q(0));carry=ck(carry+ck(c*ck(au*firstv)))
     if u!=v and av:
      firstv=Xp[p][u].get(source,Q(0)) if u<8 else E[p].get(source,Q(0));carry=ck(carry+ck(c*ck(av*firstv)))
   z=ck(ck(oldv+nativev)-carry);values.append(z);pieces.append(list(map(str,[oldv,nativev,carry])))
   for name,x in [('old',oldv),('native',nativev),('individual_pivot_carry',carry),('combined',z)]:
    D,U=map(int,bounds['four_component_positive_envelopes'][channel][name]);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  computed.append({'case619':case,'actual_lower4_values':list(map(str,values)),'old_native_and_subtracted_pivot_pieces':pieces})
  if any(values):witness={**computed[-1],'original768_constant_combination':rec['original768_constant_combination'],'actual_a1_27_times9':rec['actual_a1_27_times9']};break
  if len(computed)%10==0:dump('phase.json',{'phase':'actual_native70_prefix','computed_prefix_count':len(computed),'exact_completed_prefix':computed,'all_source_kernel_and_bounds_checks_so_far':True})
 owned=count([d,table,kap,powers,G,N,jets,Gamma,B,Xp,E,L0,Le,S0,Se,defs,first,second,S,X,ell,op,computed,witness,pins,g,peer,cr,ce,custody])+count(d)+3000
 assert owned<51263
 ep='actual-native70-lower-prefix-and-witness-if-any.json';dump(ep,{'computed_prefix':computed,'computed_prefix_count':len(computed),'nonzero_witness':witness,'all70_lower4_zero':witness is None and len(computed)==70,'scope':'ActualLOWERfourclassesonly. A nonzero witness gives a true quadratic initial relation. Allzero closes onlylower619;upper2/HF2/totalT2remain separate.','actual_N2_not_set_zero':True,'individual_pivot_ell_SUM_retained':True,'no_new_physical_gauge_implicit_or_jet_coefficient_solve':True})
 r={'id':STEM,'all_passed':True,'computed_prefix_count':len(computed),'nonzero_lower_witness_found':witness is not None,'all70_lower4_zero':witness is None and len(computed)==70,'source_committed_bits':128,'actual_bits':maximum,'whole_numeric_cap':51263,'actual_end_conservative_owned_count':owned,'pre_value_count_addendum_sha256':pins[cb+'/report.json'],'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{n:sha(out/n) for n in ['phase1-cached-C72-regular-G8-and-source-overlaps.json',ep]},'seconds':time.monotonic()-t};dump('report.json',r);print(json.dumps({'all_passed':True,'prefix':len(computed),'witness':witness,'actual_bits':maximum,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output.resolve(),a.source_peer,a.custody)
