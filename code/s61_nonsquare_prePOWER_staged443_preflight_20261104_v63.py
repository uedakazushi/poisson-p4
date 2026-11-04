#!/usr/bin/env python3
"""Lazy4+4+3 source inverse plus ONE G115/h9target19 representation-cost profile."""
import argparse,hashlib,importlib.util,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def rankone_proof():
 def atom(i):return{(i,):1}
 def add(a,b,s=1):
  z=a.copy()
  for m,c in b.items():
   z[m]=z.get(m,0)+s*c
   if not z[m]:del z[m]
  return z
 def mul(a,b):
  z={}
  for m,c in a.items():
   for n,d in b.items():
    k=tuple(sorted(m+n));z[k]=z.get(k,0)+c*d
    if not z[k]:del z[k]
  return z
 def det(a):
  z={}
  for p in itertools.permutations(range(3)):
   t={():1}
   for i in range(3):t=mul(t,a[i][p[i]])
   parity=sum(p[i]>p[j]for i in range(3)for j in range(i+1,3));z=add(z,t,-1 if parity%2 else 1)
  return z
 A=[[atom(4*i+j)for j in range(4)]for i in range(4)];u=[atom(16+i)for i in range(4)];v=[atom(20+j)for j in range(4)];k=atom(24)
 for i in range(4):
  for j in range(4):
   rr=[r for r in range(4)if r!=j];cc=[c for c in range(4)if c!=i];M=[[A[r][c]for c in cc]for r in rr];tilde=[[add(mul(k,A[r][c]),mul(u[r],v[c]),-1)for c in cc]for r in rr]
   reduced=mul(k,det(M))
   for s,c in enumerate(cc):
    replaced=[[u[r]if t==s else A[r][cc[t]]for t in range(3)]for r in rr];reduced=add(reduced,mul(v[c],det(replaced)),-1)
   assert not add(det(tilde),mul(mul(k,k),reduced),-1)
 return dict(all16_adjugate_entry_kappa_squared_factors_verified=True,formal_variables=25,no_parameter_polynomial_products=True)
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 ks='S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63';ko=load(ks,'operator.json');kp=load(ks,'proof.json');assert load(ks+'_REVIEW','report.json')['all_passed']
 lift='S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63';assert load(lift,'report.json')['all_passed']
 hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');rank=load('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63');entry='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63';eb=load(entry,'source-binding.json');assert load(entry+'_REVIEW','report.json')['all_passed']
 nodes=[]
 def n(op,sh,args=(),**kw):
  assert all(0<=i<len(nodes)for i in args);i=len(nodes);nodes.append(dict(id=i,op=op,shape=sh,args=list(args),**kw));return i
 def ref(st,name,ptr,sh):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return n('external_frozen_reference',sh,reference=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),json_pointer=ptr))
 def prod(a,b):assert nodes[a]['shape'][1]==nodes[b]['shape'][0];return n('matrix_product',[nodes[a]['shape'][0],nodes[b]['shape'][1]],[a,b])
 def sel(a,r,c):return n('submatrix',[len(r),len(c)],[a],rows=r,columns=c)
 def sub(a,b):return n('difference',nodes[a]['shape'],[a,b])
 def add(a,b):return n('sum',nodes[a]['shape'],[a,b])
 def neg(a):return n('negative',nodes[a]['shape'],[a])
 def scale(s,a):return n('central_scalar_multiple',nodes[a]['shape'],[s,a])
 def eye(k):return n('identity',[k,k])
 def zero(m,k):return n('zero',[m,k])
 def block(rr):
  h=[nodes[r[0]]['shape'][0]for r in rr];w=[nodes[x]['shape'][1]for x in rr[0]]
  for i,r in enumerate(rr):
   for j,x in enumerate(r):assert nodes[x]['shape']==[h[i],w[j]]
  return n('block_matrix',[sum(h),sum(w)],[x for r in rr for x in r],blocks=rr)
 M=ref(ks,'operator.json','/exports/input15x22',[15,22]);Dcol=ref(ks,'operator.json','/exports/normalized_column_units',[22,22]);Dcoli=n('diagonal_known_FULLU_unit_inverse',[22,22],[Dcol]);cols=[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31]
 rows=[2,4,5,8,0,1,3,6,7,13,10];rest=[9,11,12,14];targets=[27,31,22,2,11,13,25,15,28,30,29];pc=[cols.index(j)for j in targets];rc=[i for i in range(22)if i not in pc]
 Pr=n('row_permutation',[15,15],order=rows+rest);Pri=n('transpose',[15,15],[Pr]);Pc=n('column_permutation',[22,22],order=pc+rc);Pci=n('transpose',[22,22],[Pc]);MM=prod(prod(Pr,prod(M,Dcol)),Pc)
 A4=sel(MM,list(range(4)),list(range(4)));B44=sel(MM,list(range(4)),list(range(4,8)));C44=sel(MM,list(range(4,8)),list(range(4)));Qraw=sel(MM,list(range(4,8)),list(range(4,8)));u=sel(MM,list(range(4,8)),[2])
 T=ref(entry,'root239.json','/polynomial',[1,1]);iT=n('certified_chart_inverse',[1,1],[T],chart='D(T)');kap=[ref(entry,f'root{i}.json','/polynomial',[1,1])for i in eb['cofactor_roots']];kaprow=n('horizontal_concatenation',[1,4],kap);adjA=n('lazy_polynomial_adjugate_with_last_row_cache',[4,4],[A4,kaprow],last_row_equals='minus kappa',other12_entries_unexpanded=True);Ai=neg(scale(iT,adjA));Q4=sub(Qraw,prod(prod(C44,Ai),B44))
 v=[]
 for j in range(4):
  a=sel(MM,[0,1,2],[0,1,4+j]);v.append(n('negative_universal_three_determinant',[1,1],[a],expanded=False))
 vv=n('horizontal_concatenation',[1,4],v);red=n('rank_one_adjugate_reduced_polynomial',[4,4],[kap[3],Qraw,u,vv],formula='entry(i,j)=(-1)^(i+j)*(k*det(A_withoutrowj_coli)-v_withouti^t*Adj(A_withoutrowj_coli)*u_withoutj)',no_k_inverse_or_parameter_division=True,all16_entries_unexpanded=True)
 Sax=ref('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63','native.json','/exact_Saux_Laurent_polynomial',[1,1]);c=n('QQ_constant',[1,1],value=kp['cprod']);iCS=n('certified_chart_inverse',[1,1],[prod(c,Sax)],chart='D(Saux)',known_cprod_QQ_unit=True);Qi=neg(scale(iCS,red))
 AB=prod(Ai,B44);CA=prod(C44,Ai);A8=sel(MM,list(range(8)),list(range(8)));A8i=block([[add(Ai,prod(prod(AB,Qi),CA)),neg(prod(AB,Qi))],[neg(prod(Qi,CA)),Qi]])
 B3=sel(MM,list(range(8,11)),list(range(8,11)));C38=sel(MM,list(range(8,11)),list(range(8)));N3=ref('S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63','operator.json','/exports/H5_B3',[3,3]);D3=ref('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63','certificate.json','/nonzero_selected_minor_Delta',[1,1]);iD3=n('certified_chart_inverse',[1,1],[D3],chart='D(Delta407)');rhos=[ref('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json',f'/known_row_units_rho14/{i}',[1,1])for i in[6,12,9]];RD=n('diagonal',[3,3],rhos,diagonal=rhos);RDi=n('diagonal_known_FULLU_unit_inverse',[3,3],[RD]);B3i=prod(scale(iD3,n('universal_polynomial_adjugate',[3,3],[N3],expanded=False)),RDi)
 K=sel(MM,list(range(11)),list(range(11)));Ki=block([[A8i,zero(8,3)],[neg(prod(prod(B3i,C38),A8i)),B3i]])
 B=sel(MM,list(range(11)),list(range(11,22)));C=sel(MM,list(range(11,15)),list(range(11)));D=sel(MM,list(range(11,15)),list(range(11,22)));KB=prod(Ki,B);CK=prod(C,Ki);res=sub(D,prod(C,KB));LS=block([[Ki,zero(11,4)],[neg(CK),eye(4)]]);LSi=block([[K,zero(11,4)],[C,eye(4)]]);UT=block([[eye(11),neg(KB)],[zero(11,11),eye(11)]]);UTi=block([[eye(11),KB],[zero(11,11),eye(11)]])
 fs=prod(LS,Pr);fsi=prod(Pri,LSi);ft=prod(prod(Dcol,Pc),UT);fti=prod(prod(UTi,Pci),Dcoli)
 operator=dict(nodes=nodes,exports=dict(A4=A4,A4inverse=Ai,Q4=Q4,Q4inverse=Qi,reduced_polynomial_adjugate=red,pivot11=K,pivot11inverse=Ki,residual4x11=res,fullsource15=fs,fullsource15_inverse=fsi,fulltarget22=ft,fulltarget22_inverse=fti),source_pivot_order=rows,source_remaining_order=rest,target_pivot_original32=targets,target_remaining_original32=[cols[i]for i in rc],chart='FULLU D(Delta407*Saux*T)',all8_tail_effects_preserved_tail2_is_unit_pivot=True,full57_source_lift_reference=lift,new_parameter_coefficients=0)
 dump(out/'operator.json',operator);dump(out/'rank-one-universal-proof.json',rankone_proof())
 # Source-only cost of ONE remainingG115 row/target19, two exact representations.
 path=ROOT/'code/s61_nonsquare_prePOWER_h22_kernel9_support_preflight_20261104_v63.py';sp=importlib.util.spec_from_file_location('bounds',path);b=importlib.util.module_from_spec(sp);sp.loader.exec_module(b);pins[str(path.relative_to(ROOT))]=sha(path)
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');fixed=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');cached={x['old_DAG_root']:b.atom(x['Laurent_polynomial'])for x in fixed['node_values']};vals=[]
 for i,a in enumerate(q['scalar_DAG']['nodes']):
  if i in cached:x=cached[i]
  elif a[0]=='z':x=None
  elif a[0]=='one':x=b.ONE
  elif a[0]in['a','s','m']:x=b.mm(vals[a[1]],vals[a[2]])if a[0]=='m'else None if a[0]=='s'and a[1]==a[2]else b.aa(vals[a[1]],vals[a[2]])
  else:
   num,den=cp['matrix25x32'][a[1]][a[2]]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   x=b.atom(num)
   if len(den)==1:assert den[0][0][1:]==[0,0];x=b.shift(x,-den[0][0][0])
   else:x=None if x is None else(x[0],x[1],1)
  vals.append(x)
 del cp,pd,fixed
 bottom=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];pre=load('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63','operator-and-cost.json');cofs=[]
 for i in pre['source_coefficients3']:cofs.append(b.atom(load('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63',f'source-cofactor-node{i}.json')['polynomial']))
 labels6=[22,11,13,25,15,19];colidx=[q['new22_target_original32_columns'].index(j)for j in labels6];Pm=[[vals[q['remaining15x22_scalar_roots'][r][j]]for j in colidx]for r in[0,1,3,6]]
 for j,label in enumerate(labels6):Pm[0][j]=b.atom(h['source_h9'][h['original_target_labels'].index(label)])
 Bhat=[]
 for j in range(6):
  z=None
  for co,r in zip(cofs,[2,0,1]):z=b.aa(z,b.mm(co,b.atom(bottom[r][j])))
  Bhat.append(z)
 brow10=rank['cleared_upper_row_coordinates7xr'][3];assert len(brow10[0])==415 and not brow10[1]and not brow10[2]
 assert hc['known_row_units_rho14'][6]['a_d_F_B_C_J_powers']==hc['known_row_units_rho14'][10]['a_d_F_B_C_J_powers']==[1,0,1,0,2,0]
 delta=b.atom(rank['nonzero_selected_minor_Delta']);lam=b.atom(brow10[0]);g=b.atom(hc['raw_compact14x5'][6][1]);hh=b.atom(hc['raw_compact14x5'][10][1]);profiles=[]
 for kind,c1,c2 in [('Delta_lambda',delta,lam),('actual_coefficient_chart',g,hh)]:
  Kr=[]
  for j in colidx:
   a=vals[q['remaining15x22_scalar_roots'][11][j]];c0=vals[q['remaining15x22_scalar_roots'][7][j]];assert a[2]==c0[2]==0;Kr.append(b.aa(b.mm(c1,a),b.mm(c2,c0)))
  matrix=Pm+[Bhat,Kr];states={0:b.ONE}
  for row in matrix:
   nex={}
   for mask,x in states.items():
    for j,y in enumerate(row):
     if not(mask>>j&1):nex[mask|1<<j]=b.aa(nex.get(mask|1<<j),b.mm(x,y))
   states=nex
  outbound=states[63];profiles.append(dict(representation=kind,one_remaining_source_current15=11,actual135_label=115,pivot_source_current15=7,pivot_actual135=85,original32_target=19,det6_source_rows='P4,Bhat4,K',det6_source_columns=labels6,cleared_output_bound=b.enc(outbound),input_entry_support_sum=sum(b.size(x)for row in matrix for x in row),K_six_input_support=[b.size(x)for x in Kr],height_lifetime_whole_worker_NOT_evaluated=True,scalar_value_NOT_evaluated=True))
 cost=dict(one_source_row_only=True,one_target_only=True,profiles=profiles,lambda415_against_Brow6_ONLY=True,other_two_B_Cramer_coordinates_literal_zero=True,
  exact_large_row='K=Delta407*G115/c10-lambda415*G85/c6',exact_small_row='Ksmall=N6,1*G115/c10-N10,1*G85/c6',
  polynomial_relation='Delta407*Ksmall=N6,1*K from publishedBrowcoordinate identity anddomain cancellation',
  projected_residual_target19='c10*det6(P4,Bhat,K)/(Delta407*Saux); on extraD(N6,1), equivalentlyc10*det6(P4,Bhat,Ksmall)/(N6,1*Saux)',
  small_chart_complement_retained=True,three_middle_coefficient_charts_cover_DDelta='B3selectedrow6 isunimodular afterDelta inversion; adjugate suppliesunitcombination ofits3entries',
  no_new_scalar_GCD_or_entry_or_rank=True)
 dump(out/'one-row-source-cost.json',cost)
 proof=dict(detA4='-T',row_swap_count_pair_first=16,column_swap_count_A4first=16,detQ4='-cprod*Saux/kappa3',
  inverseQ4='-AdjReduced/(cprod*Saux)',no_kappa_inverse_or_new_exactdivision=True,
  all16_newAdj_entries_left_lazy=True,fullsource_and_target_inverse_blockproof=True,full57_source_lift_bound_to_frozen181node_operator=True)
 dump(out/'proof.json',proof);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
 result=dict(id=ID,all_passed=True,frozen=True,source_only=True,typed_nodes=len(nodes),actual_samecokernel_shape=[4,11],staged_inverse_shapes=[4,4,3],no_new_scalar_coefficients=True,no_extra_kappa_p_R_h22_inverse=True,one_row_cost_profiles=profiles,no_worker_budget_pass_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
