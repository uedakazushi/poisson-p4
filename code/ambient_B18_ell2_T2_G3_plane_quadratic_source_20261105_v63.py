#!/usr/bin/env python3
"""Factored finite-resolution G3_2 source gate; no signed quadratic values."""
import argparse,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G3_PLANE_QUADRATIC_INPUT_GATE_20261105_v63'
PREFIX='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'factored_known_resolution_source_only_no_signed_quad'})
 hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 for p in ['reports/AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63/report.json',
  'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
  'reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63_REVIEW/report.json']:
  assert read(p)['all_passed']
 rpath='reports/AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63/actual-G3-plane-condition-fixed-rank-and-kernel.json';r=read(rpath)
 t=read(H.PATHS['transport']);pref=read(H.PATHS['active']);active=pref['active_variables'];assert active==[0,1,4,5,7];vars6=active+[8]
 hs[PREFIX]=H.sha(ROOT/PREFIX)
 I3=list(map(H.parse,H.selected(PREFIX,'actual_active_I3_origin_columns')))
 PiY=list(map(H.parse,H.selected(PREFIX,'PiY3_from102_columns')))
 Pi42=list(map(H.parse,t['active_Pi42_columns']));D30=list(map(H.parse,pref['active_d3_columns']))
 fixed=list(map(H.parse,t['active_F2_fixed182_columns']));norm={};reps={}
 for rec in t['active_F2_fixed182_receipts']:
  j=rec['source_column'];p=rec['pivot'];w={j:F(1)}
  for k,c in rec['subtract']:w=add(w,reps[k],-F(c))
  reps[p]=H.guard(scale(w,1/F(rec['divide_by'])));norm[p]=H.parse(rec['normalized']);assert mv(fixed,reps[p])==norm[p]
 # Literal raw differential is affine-linear; only its fixed variable part.
 blocks2=[(2-q,w) for q in range(3) for w in itertools.combinations(active,q)];look={w:k for k,(p,w) in enumerate(blocks2)}
 blocks3=[(3-q,w) for q in range(4) for w in itertools.combinations(active,q)];raw=[]
 for p,w in blocks3:
  for j in range(14):
   v={}
   if p:v[(14*look[w]+j)*6+5]=F(1)
   for s,z in enumerate(w):v=add(v,{(14*look[w[:s]+w[s+1:]]+j)*6+vars6.index(z):F((-1)**(p+s))})
   raw.append(v)
 assert len(raw)==364
 selected=t['active126_source_columns'];P=[]
 for j in selected:
  v={}
  for k,c in raw[j].items():
   row,var=divmod(k,6)
   for target,a in Pi42[row].items():v=add(v,{target*6+var:c*a})
  P.append(H.guard(v))
 kernels=[H.parse(q['original238_source_combination']) for q in r['full_constant_kernel_columns']];assert len(kernels)==r['enlarged_constant_kernel_dimension'];nker=len(kernels)
 A=[];kernel_source=[];first_z=[];linear_support=0
 Efree=[{j:F(1)} for j in t['active_F2_free42']]
 for index,c in enumerate(kernels):
  w=H.guard(mv(I3,c));assert not mv(D30,w)
  v=H.guard(mv(raw,w));corr={};first={}
  for var in range(6):
   q={i//6:a for i,a in v.items() if i%6==var};pi=mv(Pi42,q)
   assert not mv(PiY[:42],pi)
   if var in (4,5):assert not pi
   for i,x in pi.items():first[i*6+var]=x
   original=add(q,mv(Efree,pi),-1);q=dict(original);a={}
   while q:
    p=min(q);assert p in norm;x=q[p];a=H.guard(add(a,reps[p],x));q=H.guard(add(q,norm[p],-x))
   assert mv(fixed,a)==original
   for i,x in a.items():
    if i>=56:corr[(i-56)*6+var]=x
  A.append(H.guard(corr));kernel_source.append(w);first_z.append(first);linear_support+=len(corr)
  live=sum(map(len,A+kernel_source+first_z+P+raw+Pi42+PiY+fixed))+sum(map(len,norm.values()))+sum(map(len,reps.values()))+sum(map(len,kernels))+20000
  assert live<150000,('source_live',index,live)
 # Frozen degree2 normal form is a linear map only; evaluating it on basis
 # monomials does not evaluate any actual G3quadratic relation.
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json')
 nodes=lin['degree2_monomial_rewrite_nodes'];memo={}
 def nf(i):
  if i in memo:return memo[i]
  n=nodes[i]
  if n[0]=='zero_ell2':v={}
  elif n[0]=='standard':v={n[1]:F(1)}
  else:
   v={}
   for j,x in n[2]:v=add(v,nf(j),F(x))
  memo[i]=v;return v
 mons44={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 pairs=[w for w in itertools.combinations_with_replacement(range(6),2) if w!=(5,5)];assert len(pairs)==20
 NF=[]
 for row in range(42):
  for u,v in pairs:
   m=[0]*9;m[vars6[u]]+=1;m[vars6[v]]+=1;mi=mons44[tuple(m)];q={}
   for target,c in PiY[row].items():q=add(q,nf(lin['all264_NF_root_nodes'][target][mi]),c)
   NF.append(q)
 pair_index={w:i for i,w in enumerate(pairs)}
 abs_bounds=[];terms=0;highest=0
 for a in A:
  bound=[F(0) for _ in range(6)];D=1
  for k,x in a.items():
   source,u=divmod(k,6)
   for i,y in P[source].items():
    row,v=divmod(i,6);pair=tuple(sorted((u,v)))
    if pair==(5,5):continue
    for target,z in NF[row*20+pair_index[pair]].items():
     term=abs(x*y*z);D=math.lcm(D,term.denominator);bound[target]+=term;terms+=1
  # A common positive denominator and absolute cleared numerator bound also
  # bound every signed partial sum, without evaluating the signed sum.
  U=max([0]+[x.numerator*(D//x.denominator) for x in bound])
  highest=max(highest,D.bit_length(),U.bit_length());abs_bounds.append([str(D),str(U)])
 sourcecount=sum(map(len,A+P+NF+kernels+kernel_source+first_z))+len(nodes)+sum(map(len,raw+Pi42+PiY+fixed))+sum(map(len,norm.values()))+sum(map(len,reps.values()))
 nextbounds={'factoredA_enlarged':sum(map(len,A)),'selected126_linearP':sum(map(len,P)),
  'PiY_L48_NF840':sum(map(len,NF)),'kernel_and_source_binding':sum(map(len,kernels+kernel_source+first_z)),
  'all_enlarged_projected_quad_output_dense':6*nker,'rank_basis6_by6_and_source_reps6_by_enlarged':36+6*nker,
  'all_enlarged_source_factors6':6*nker,'fixed_input_and_metadata_scratch':18000}
 whole=sum(nextbounds.values());assert whole<60000
 max_input=max([0]+[max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for v in A+P+NF for x in v.values()]);operational=max(highest,3*max_input+2)
 assert operational<2048
 rank_slots=nextbounds['rank_basis6_by6_and_source_reps6_by_enlarged']+nextbounds['all_enlarged_source_factors6']
 raw_bound=(whole-rank_slots)*(2*len(str(2**operational))+24)+rank_slots*(2*len(str(2**2048))+24)+2_000_000;assert raw_bound<32*1024*1024
 data={'source_kernel_enlarged':[srec(v) for v in kernels],'source_original_raw364_combinations':[srec(v) for v in kernel_source],
  'actual_enlarged_selected126_correction_linear_columns':[srec(v) for v in A],
  'actual_selected126_rawD31_Pi42_linear_columns':[srec(v) for v in P],
  'PiY_L48_degree2_NF42_times20_columns':[srec(v) for v in NF],
  'active_variable_global_indices':vars6,'degree2_active_pair_order20':[list(w) for w in pairs],
  'all_enlarged_inner_minus_EfreePi42_vectors_fixed182_span_checked':True,'actual_nonzero_first_z_linear_columns':[srec(v) for v in first_z],
  'typed_actual_quadratic_relation_formula':'NF(PiY0*G3_2*c)=-sum(A126linear[V-EfreePi42V] * Pi42D31selected126linear) then frozenL48NF; ell²termszero',
  'unit_derivative_proof':'P42(t)=last42 inverse[Eselected56,D3selected126(t),Efree42]; changingEfree42 byselected56 multiples leavesP42 unchanged. P42(t)D3selected126(t)=0. ConstantI3 hasD30I3=0, soG3(t)=P42(t)D31I3. P42_1=-Pi42D31selected126*A126. For enlargedkernel firstPi42V is z-linear only andPiY0first0. The unknownfinal99quotient andbottomderivatives times thisz-linear vector vanish modL48 becauseall378 target*z*variable NFzero. ActualG3bottomzero beforeunitclear.',
  'not_a_full_HF2_completeness_claim':True,'signed_actual_quadratic_values_evaluated':False,
  'positive_bound_per_enlarged_common_D_U':abs_bounds,'positive_output_bound_bits':highest,'max_factored_input_bits':max_input,'product_and_output_operational_bound_bits':operational,
  'signed_product_summand_count':terms,'source_phase_owned_slots':sourcecount+20000,
  'next_whole_slot_bound':whole,'next_slot_breakdown':nextbounds,'next_raw_byte_bound':raw_bound,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'QQ_slots':150000,'raw_bytes':32*1024*1024,'guard_bits':2048}}
 ep=out/'actual-enlarged-G3quadratic-factored-source-and-known-NF-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'signed_actual_quadratic_values_evaluated':False,'new_physical_or_implicit_coefficients':False,
  'enlarged_zero_projected_linear_kernel_count':nker,'next_whole_slot_bound':whole,'next_raw_byte_bound':raw_bound,
  'positive_output_bound_bits':highest,'max_factored_input_bits':max_input,'source_phase_owned_slots':sourcecount+20000,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','enlarged_zero_projected_linear_kernel_count','next_whole_slot_bound','next_raw_byte_bound','positive_output_bound_bits','max_factored_input_bits','source_phase_owned_slots','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
