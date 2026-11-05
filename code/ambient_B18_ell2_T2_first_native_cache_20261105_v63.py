#!/usr/bin/env python3
"""Finite C16 native adapter + intrinsicS2rell from accepted caches only."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as U
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_T_native_S2_B11_20261104_v63 import Finite,counts,encode
from ambient_B18_ell2_raw33_quotient_20261105_v63 import mv,srec,add
ROOT=U.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63'
class C16(Finite):
 def __init__(self,nus):
  self.m=2;self.b=4;self.n=16;self.Z=(F(0),)*16;self.ONE=self.c(1);self.nus=nus;self.nilcap=9;self.table={};self.inverse_checks=0
  def base(a,b):
   v={}
   for i,c in a.items():
    for j,d in b.items():
     if not(i&j):v=add(v,{i|j:c*d})
   return v
  for i in range(16):
   ai,ti=divmod(i,4)
   for j in range(16):
    aj,tj=divmod(j,4);v=base({ai:F(1)},{aj:F(1)});carry=ti&tj;mask=ti^tj
    for k in range(2):
     if carry&(1<<k):v=base(v,{a:F(c) for a,c in enumerate(nus[k]) if F(c)})
    if v:self.table[i,j]=[(4*a+mask,c) for a,c in v.items()]
def main(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);U.dump(out/'phase.json',{'phase':'authorizedfiniteC16nativecache_noimplicit'})
 hs={};saved=0;peak=0
 def read(p):hs[p]=U.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def save(name,data,live=()):
  nonlocal saved,peak
  q=counts(data);v=counts(live);peak=max(peak,saved+q+v+15000);assert peak<150000
  def check(x):
   if isinstance(x,F):assert max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=384
   elif isinstance(x,dict):
    for v in x.values():check(v)
   elif isinstance(x,(list,tuple)):
    for v in x:check(v)
  check(data);U.dump(out/name,encode(data));saved+=q;assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
  U.dump(out/'phase.json',{'phase':name,'saved_coefficients':saved,'owned_saved_plus_live_bound':peak})
 gr=read(gate);assert gr['all_passed'] and not gr['actual_native_adapter_or_secondvalues_started'];ip=str(Path(gate).parent/'actual-native-Rell4-cache-source-input-and-adapter-gate.json');D=read(ip);assert hs[ip]==gr['evidence_sha256'][Path(ip).name]
 peer=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63_REVIEW/report.json');assert peer['all_passed']
 for p,h in gr['input_sha256'].items():assert U.sha(ROOT/p)==h;hs[p]=h
 S,I,P,Q,R,K,C,prep,E,Op=[read(D['normalizer_source_references'][k]) for k in ['source','interface','proper','equations','first5','second5','membership','preparation','normalizer','Q_operator']]
 first={i:tuple(map(F,v)) for i,v in D['actual105_Rell4']};nu=list(map(lambda v:tuple(map(F,v)),D['actual_nu_A_B_Rell4']));actual=list(map(lambda v:tuple(map(F,v)),D['actual14_C16']));kappa=list(map(lambda v:tuple(map(F,v)),D['actual_kappa9_C16']));F12=list(map(lambda v:tuple(map(F,v)),D['actual_physicalF12_C16']));ev=list(map(lambda v:tuple(map(F,v)),D['actual_even15_Rell4']));od=list(map(lambda v:tuple(map(F,v)),D['actual_odd15_Rell4']))
 a=C16(nu);a.install();assert [actual[j] for j in I['retained_source9_indices']]==kappa
 prepath='reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json';hs[prepath]=U.sha(ROOT/prepath);GL=U.selected_top_field(ROOT/prepath,'GL18_source_parameter_positions')
 def frame_family(family,tau):
  f=S['families'][family];ss=a.dual(first[f['p_quadratic_relation']['s']]);p=H.add(a.tau(tau),H.scale(ss,F(1,2)));G=H.ident(5)
  for name,(aa,bb) in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma') and name!='gamma0':
    index=int(name[5:]);i,j=divmod(GL[index-1],5);G[i][j]=H.add(G[i][j],a.dual(first[aa]),H.mul(a.dual(first[bb]),p))
  core=H.ident(5)
  if family=='CA':core[0][3]=H.scale(p,F(-2));core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
  else:core[0][3]=H.scale(p,F(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),F(1,4))
  return H.mm(G,core)
 G=H.ident(5)
 for (i,j),x,y in zip(I['actual_group15_positions'],ev,od):G[i][j]=H.add(G[i][j],a.dual(x),H.mul(a.dual(y),a.tau(0)))
 HA=frame_family('CA',0);assert H.frame(H.mm(G,HA))==H.ident(5);FB=H.frame(H.mm(G,frame_family('CB',1)));assert [FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]==actual[9:]
 rows=P['zero_coordinate_rows36'];targets=[18,19,28,29,60,67,72,63,64,68,73,74];transform=H.mm(FB,G);itransform=H.mi(transform);ann=H.constmat(P['CB_annihilator17_by36']);R17=[[a.Z]*12 for _ in range(17)]
 for j,target in enumerate(targets):
  inp=[a.Z]*150;inp[target]=a.ONE
  for k,row in enumerate(E['eliminated_zero_rows4']):inp[row]=H.c(E['eliminated_four_coordinate_reconstruction'][k][rows.index(target)])
  pushed=H.physical_push(inp,transform,itransform);v=H.mm(ann,[[pushed[i]] for i in rows])
  for i in range(17):R17[i][j]=v[i][0]
 assert H.mm(R17,[[p] for p in F12])==[[H.evalpoly(p,actual)] for p in Q['equations14variable_sparse']]
 du=H.add(a.ONE,H.neg(H.mul(kappa[1],kappa[1])));pu=H.add(H.c(5),H.mul(kappa[1],kappa[1]));rr=[]
 for name in ['v0','v1','b','c','e']:
  p=R['five_explicit_rational_roots'][name];rr.append(H.mul(H.evalpoly(p['numerator'],kappa),H.power(du,-p['D_denominator_power'])))
 coeff=[H.telescope(p,actual,[0,1,10,11,13],rr) for p in Q['equations14variable_sparse']];Ai=H.mi(coeff[:5]);BA=H.mm(coeff[5:],Ai)
 R12=[[H.add(R17[i+5][j],H.neg(H.add(*(H.mul(BA[i][k],R17[k][j]) for k in range(5))))) for j in range(12)] for i in range(12)]
 for i,p in enumerate(R['all12_actual_remaining_and_divided_records']):
  unit=H.scale(H.power(du,p['D_denominator_power']),F(p['constant_denominator'],p['primitive_content_unit']));R12[i]=[H.mul(unit,q) for q in R12[i]]
 hpolys=[p['divided_primitive_polynomial'] for p in R['all12_actual_remaining_and_divided_records']];hv=[[H.evalpoly(p,kappa)] for p in hpolys];assert H.mm(R12,[[p] for p in F12])==[[H.mul(kappa[1],p[0])] for p in hv]
 save('phase1-actual-offintersection17-and-R12-C16.json',{'R17':R17,'R12':R12,'full_first5_telescopes':True,'all17_identity':True,'dH12_R12F12':True},live=[coeff,BA,Ai,actual,kappa,F12,hv])
 del R17,coeff,BA,Ai,transform,itransform,inp,pushed,v,HA,FB,G
 core=[kappa[j] for j in [0,1,2,6]];inds=[K['source9_variable_order'].index(name) for name in K['selected5_pivot_order']];rr=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];rr.append(H.mul(H.evalpoly(p['numerator'],core),H.power(du,-p['D_power']),H.power(pu,-p['P_power'])))
 co=[H.telescope(p,kappa,inds,rr) for p in hpolys];perm=K['selected5_source_H_rows']+[i for i in range(12) if i not in K['selected5_source_H_rows']];Ai=H.mi([co[i] for i in perm[:5]]);BA=H.mm([co[i] for i in perm[5:]],Ai);T=[[a.Z]*12 for _ in range(12)]
 for i in range(5):
  for j in range(5):T[i][perm[j]]=Ai[i][j]
 for i in range(7):
  T[i+5][perm[i+5]]=a.ONE
  for j in range(5):T[i+5][perm[j]]=H.neg(BA[i][j])
 for i,p in enumerate(K['remaining7_actual_core_polynomials']):
  lam=H.scale(H.mul(H.power(du,-p['D_denominator_power']),H.power(pu,-p['P_denominator_power'])),F(p['QQ_unit_content_numerator'],p['QQ_unit_content_denominator']));T[i+5]=[H.mul(H.inv(lam),x) for x in T[i+5]]
 E3=H.ident(12);E3[5][5]=H.inv(du);E3[6][6]=H.inv(du);r=H.add(core[2],H.neg(core[0]),H.neg(core[1]));E3[7][5]=H.neg(H.mul(r,H.inv(H.mul(pu,du))));E3[7][7]=H.scale(H.inv(H.mul(pu,du)),F(-2));T=H.mm(E3,T);E4=H.ident(12)
 for j in range(4):
  fac=H.power(du,[2,2,2,3][j])
  for i in range(3):E4[8+j][5+i]=H.neg(H.mul(fac,H.evalpoly(C['cofactor3_by4'][i][j],core)))
 T=H.mm(E4,T);N=H.mm(T,R12);gs=[H.add(kappa[j],H.neg(v)) for j,v in zip(inds,rr)]+[H.evalpoly(p,core) for p in prep['generators3']]
 assert H.mm(T,hv)==[[p] for p in gs]+[[a.Z] for _ in range(4)] and H.mm(N,[[p] for p in F12])==[[H.mul(kappa[1],p)] for p in gs]+[[a.Z] for _ in range(4)]
 npath='reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json';hs[npath]=U.sha(ROOT/npath);oldN=U.selected_top_field(ROOT/npath,'N_C_B44');assert [[list(p[:8]) for p in row] for row in N]==[[[F(c) for c in p[:8]] for p in row] for row in oldN]
 ge=read('reports/AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63/actual30-firstell-and-native-S1ell.json');assert [list(p[8:12]) for p in gs]==[[F(c) for c in v] for v in ge['actual_G8_C4_first_ell']]
 save('phase2-actual-T-N-G8-C16.json',{'T':T,'N':N,'actualG8':gs,'TH12_G8_zero4':True,'NF12_dG8_zero4':True,'all_N0_rfirst_matches':True,'all_G8_firstell_matches':True},live=[co,Ai,BA,R12,hv])
 del co,Ai,BA,R12,hv,E3,E4,T,core,actual,kappa,F12
 Sn=[[a.Z]*32 for _ in range(12)]
 for i,j,g,sign in Op['native_S12_by32']['koszul28_in_first8_rows']:Sn[i][j]=H.scale(gs[g],F(sign))
 for i,j,c in Op['native_S12_by32']['free4_zero_source_columns']:Sn[i][j]=H.c(c)
 S0=[{} for _ in range(128)];S2=[{} for _ in range(128)];Sr=[{} for _ in range(128)];Se=[{} for _ in range(128)]
 for j in range(32):
  for inp in range(4):
   u=[F(0)]*16;u[inp]=F(1)
   for i in range(12):
    p=H.mul(Sn[i][j],tuple(u))
    for mask in range(4):
     for base,cols in [(0,S0),(1,Sr),(2,Se),(3,S2)]:
      if p[4*base+mask]:cols[4*j+inp][4*i+mask]=p[4*base+mask]
 mg=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json');assert S0==list(map(U.parse,mg['native_S0_columns'])) and Se==list(map(U.parse,mg['native_S1ell_columns']))
 L=[]
 for base in range(3):
  columns=[{} for _ in range(48)]
  for i in range(12):
   for k in range(4):
    for j in range(12):
     for mask,c in enumerate(N[i][j][4*base:4*base+4]):
      if c and not(mask&k):columns[4*i+k]=add(columns[4*i+k],{4*j+(mask|k):c})
  L.append(columns)
 assert L[0]==list(map(U.parse,mg['L0_Res4_Ntranspose_columns']))
 fold=read('reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json');Gamma=list(map(U.parse,fold['two_Gamma_on56_actual_X2']));GL0=list(map(U.parse,fold['two_Gamma_L0_on48_intrinsic_native_S2']));values=[]
 def dot(f,v):return sum((c*v.get(i,F(0)) for i,c in f.items()),F(0))
 for row in fold['all619_source_folds']:
  c=U.parse(row['native_intrinsic_S2_source_combination']);assert not mv(S0,c);v=mv(S2,c);omr=U.parse(row['native_Omega_r']);ome=U.parse(row['native_Omega_ell']);vr=mv(L[1],ome);ve=mv(L[2],omr)
  z=[dot(GL0[k],v)+dot(Gamma[k],vr)+dot(Gamma[k],ve) for k in range(2)];U.guard({i:c for i,c in enumerate(z) if c});values.append(z)
 save('phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',{'intrinsic_S2rell_columns':[srec(v) for v in S2],
  'actual_L1r_columns':[srec(v) for v in L[1]],'actual_L1ell_columns':[srec(v) for v in L[2]],
  'all619_two_actual_native_parts':values,'all619_nativeS0_zero_rechecked':True,
  'L1ell_Omega_r_retained':True,'negative_intrinsic_pivot_correction_already_in_old_knownparts_not_added_again':True,
  'actual_combined_HF1_or_T2_length_not_computed':True},live=[N,gs,Sn])
 report={'id':STEM,'all_passed':True,'all619_actual_native_parts_source_bound':True,'all_offintersection_and_cached_source_checks':True,
  'no_new_physical_implicit_or_gauge_coefficients':True,'actual_combined_HF1_or_T2_length_computed':False,
  'input_sha256':hs,'code_sha256':U.sha(Path(__file__)),'note_sha256':U.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{p.name:U.sha(p) for p in out.iterdir() if p.name!='phase.json'},'owned_saved_plus_live_QQ_bound':peak,'seconds':time.monotonic()-start}
 U.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','owned_saved_plus_live_QQ_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);a=p.parse_args();main(a.output.resolve(),a.gate)
