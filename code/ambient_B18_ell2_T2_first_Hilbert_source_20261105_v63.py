#!/usr/bin/env python3
"""Source/known-resolution gate for actual T2 first Hilbert layer.

No X2/M1/implicit coefficient is evaluated, and no rank is searched here.
"""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63'
def parse(v):return H.parse(v)
def quotient(n,basis):
 free=[i for i in range(n) if i not in basis];fi={p:i for i,p in enumerate(free)};red={}
 for p in reversed(range(n)):
  if p in fi:red[p]={fi[p]:F(1)}
  else:
   v={}
   for j,c in basis[p].items():
    if j!=p:v=add(v,red[j],-c)
   red[p]=H.guard(v)
 return free,[red[i] for i in range(n)]
def project_linear(P,cols):
 out=[]
 for v in cols:
  q={}
  for k,c in v.items():
   row,var=divmod(k,9)
   for i,a in P[row].items():q=add(q,{i*9+var:a*c})
  out.append(H.guard(q))
 return out
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'known_resolution_source_profile_no_rank_no_X2'})
 hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 minpath='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json';mini=read(minpath)
 mu=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json')
 read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json')
 t=read(H.PATHS['transport']);pref=read(H.PATHS['active']);res=read(H.PATHS['resolution']);g=read(H.PATHS['G2']);active=pref['active_variables'];zero=pref['zero_actions']
 assert active==[0,1,4,5,7] and zero==[2,3,6]
 Ynorm={p:parse(v) for p,v in mu['normalized_basis_rows']};freeY,PiY=quotient(102,Ynorm);assert freeY==mini['remaining_free3_Y_rows']
 assert all(not mv(PiY,parse(v)) for v in mu['actual_selected_Y0_columns102'])
 inv30=list(map(parse,mini['fixed_bottom30_inverse_columns']));PiBottom=[{} for _ in range(56)]
 for i,row in enumerate(mini['free30_V_indices']):PiBottom[row]={j-27:c for j,c in inv30[i].items() if j>=27}
 # Active I3 is the fixed unit source correction, not a new rank computation.
 D3=list(map(parse,pref['active_d3_columns']));fixed=list(map(parse,t['active_F2_fixed182_columns']));norm={};reprs={}
 for rec in t['active_F2_fixed182_receipts']:
  j=rec['source_column'];p=rec['pivot'];v={j:F(1)}
  for k,c in rec['subtract']:v=add(v,reprs[k],-F(c))
  reprs[p]=scale(v,1/F(rec['divide_by']));norm[p]=parse(rec['normalized']);assert mv(fixed,reprs[p])==norm[p]
 selected=t['active126_source_columns'];complement=[j for j in range(364) if j not in selected];assert len(complement)==238
 I3=[]
 for j in complement:
  v=dict(D3[j]);a={}
  while v:
   p=min(v);assert p in norm;c=v[p];a=add(a,reprs[p],c);v=add(v,norm[p],-c)
  assert all(k>=56 for k in a)
  w=add({j:F(1)},{selected[k-56]:c for k,c in a.items()},-1);assert not mv(D3,w);I3.append(w)
 # Literal canonical d3 linear coefficients, d=dPeriodic+(-1)^p*dKoszul.
 blocks2=[(2-q,w) for q in range(3) for w in itertools.combinations(active,q)];lookup={w:k for k,(p,w) in enumerate(blocks2)}
 blocks3=[(3-q,w) for q in range(4) for w in itertools.combinations(active,q)];rawlin=[]
 for p,w in blocks3:
  for j in range(14):
   v={}
   if p:v[(14*lookup[w]+j)*9+8]=F(1)
   for s,i in enumerate(w):v=add(v,{(14*lookup[w[:s]+w[s+1:]]+j)*9+i:F((-1)**(p+s))})
   rawlin.append(v)
 assert len(rawlin)==364
 Pi42=list(map(parse,t['active_Pi42_columns']));product_bound=0
 for w in I3:
  for j,c in w.items():
   for k in rawlin[j]:product_bound+=len(Pi42[k//9])
 assert product_bound+50000<150000,('active_G3_before_expansion_support_guard',product_bound)
 activeG3=project_linear(Pi42,[mv(rawlin,v) for v in I3]);full=list(activeG3)
 activeG2=[]
 for v in map(parse,g['actual_linear_G2_active102_by42_columns']):activeG2.append({(i//6)*9+(active+[8])[i%6]:c for i,c in v.items()})
 _,G1,_,_=H.full_G2(t,pref,res,g)
 for iz,z in enumerate(zero):
  for j in range(42):
   v={j*9+z:F(1)}
   for k,c in activeG2[j].items():row,var=divmod(k,9);v=add(v,{(42+17*iz+row)*9+var:c})
   full.append(v)
 pairs=list(itertools.combinations(range(3),2))
 for ip,(ia,ib) in enumerate(pairs):
  for j in range(17):
   v={(42+17*ib+j)*9+zero[ia]:F(-1),(42+17*ia+j)*9+zero[ib]:F(1)}
   for k,c in G1[j].items():row,var=divmod(k,9);v=add(v,{(93+3*ip+row)*9+var:c})
   full.append(v)
 for j in range(3):
  full.append({(93+3*pairs.index(tuple(k for k in range(3) if k!=i))+j)*9+zero[i]:F((-1)**i) for i in range(3)})
 assert len(full)==418;projected=project_linear(PiY,full)
 # d2linear*d3linear=0 coefficient-wise (ell²=0), all source identities.
 G2=list(map(parse,read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json')['full_G2_linear234_by102_columns']))
 mon2=[(i,j) for i in range(9) for j in range(i,9) if (i,j)!=(8,8)];mi={m:i for i,m in enumerate(mon2)}
 for col in full:
  image={}
  for k,c in col.items():
   source,var=divmod(k,9)
   for i,a in G2[source].items():
    row,v=divmod(i,9);m=tuple(sorted((var,v)))
    if m!=(8,8):image=add(image,{row*44+mi[m]:a*c})
  assert not image,('literal_G2linear_G3linear_identity',image)
 slots={'active_I3_sparse':sum(map(len,I3)),'canonical_d3linear_sparse':sum(map(len,rawlin)),
  'all_G3linear_sparse':sum(map(len,full)),'projected418_G3linear_sparse':sum(map(len,projected)),
  'PiY_and_PiBottom':sum(map(len,PiY+PiBottom)),
  'all642_original_source_combinations':sum(len(q['original768_combination']) for q in mini['all642_zero_constant_combinations']),
  'old8_actual_X1_inputs':4271,'new640_ell_inputs_upper':7928,'native_effective_ell_inputs_upper':6144,
  'bottom27_by642_dense_upper':27*642,'bottom_rank_basis_and_factors_upper':27*27+27*642,
  'bottom_kernel_original_source_expansion_upper':642*40,'metadata_and_scratch':16000}
 assert sum(slots.values())<150000
 data={'base':'B=(k[ell]/ell²)[[zA0,...,zA6,rA]]','minimal6_target_order':mini['final_minimal6_target_order'],
  'PiY3_from102_columns':[srec(v) for v in PiY],'PiBottom3_from56_columns':[srec(v) for v in PiBottom],
  'active_F3_complement238':complement,'actual_active_I3_origin_columns':[srec(v) for v in I3],
  'actual_active_G3linear42_by238_columns':[srec(v) for v in activeG3],
  'actual_full_G3linear102_by418_columns':[srec(v) for v in full],
  'actual_minimal6_G3linear_top27_by418_columns':[srec(v) for v in projected],
  'all418_G2linear_G3linear_zero_mod_ell2':True,
  'all642_constant_Y0_kernel_source_combinations':mini['all642_zero_constant_combinations'],
  'source_formula':{'active_G3linear':'Pi42*d3linear*I3origin; bothunitderivativesdrop because d3originI3origin=0 andPi42d3origin=0',
   'known_G3_in_minimal6':'(PiY3*G3linear,0); allG3origins0 killsformalunitderivativeandbottom0killsunknownYp',
   'known_bottom3':'PiBottom3*actualXgamma1 for642Y0kernel; allinputsarealreadyfixedfull9X1',
   'remaining_top3':'onconstantkernelofknownbottomlinear, Xgamma1=Xp0*a1; subtractactualXp*a1 soZ=O(m²). G2linear*Y1=M0*(Xgamma2-Xp1*a1); unknownM1drops. Y1ambiguityG3linear isalreadyretainedasactualrelations',
   'next_pair_profile':'afterknownG3linear/bottomrank, deriveadjoints forremainingtoplinearquotient; do notassumeallX2pairsorhonestsecondsavailable'},
  'known_resolution_rank_or_bottom_rank_not_started':True,'new_X2_M1_physical_or_implicit_coefficients':False,
  'active_before_expansion_product_bound':product_bound,'whole_next_known_array_slot_bound':sum(slots.values()),'slot_breakdown':slots,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'coefficients':150000,'bits':384}}
 ep=out/'actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';H.dump(ep,data)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'all418_G2linear_G3linear_source_identities':True,'rank_started':False,
  'new_physical_implicit_X2_or_M1_coefficients':False,'whole_next_known_array_slot_bound':sum(slots.values()),
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_next_known_array_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
