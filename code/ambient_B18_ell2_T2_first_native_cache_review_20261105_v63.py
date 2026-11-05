#!/usr/bin/env python3
"""Independent sparse-C16 native source reconstruction; no source/gauge solve."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_T_native_S2_B11_review_20261104_v63 import B44,V
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63'
P={'primary':f'reports/{A}/report.json','R':f'reports/{A}/phase1-actual-offintersection17-and-R12-C16.json',
'N':f'reports/{A}/phase2-actual-T-N-G8-C16.json','value':f'reports/{A}/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',
'gate':f'reports/{G}/actual-native-Rell4-cache-source-input-and-adapter-gate.json','gate_peer':f'reports/{G}_REVIEW/report.json',
'GL':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
'oldN':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
'firstell':'reports/AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63/actual30-firstell-and-native-S1ell.json',
'mu_gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
'fold':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json',
'fold_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63_REVIEW/report.json',
'oldpeer':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def bmul(a,b):
 out={}
 for i,c in a.items():
  for j,d in b.items():
   if not i&j:out=add(out,{i|j:c*d})
 return out
class C16(B44):
 def mul(self,a,b):
  out={}
  for(i,u),c in a.items():
   for(j,v),d in b.items():
    p=bmul({i:c},{j:d});carry=u&v;mask=u^v
    if carry&1:p=bmul(p,self.nus[0])
    if carry&2:p=bmul(p,self.nus[1])
    out=self.add(out,{(k,mask):x for k,x in p.items()})
  return out

def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'SAVED_NATIVE_C16_SOURCE_REPLAY_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in ['gate_peer','fold_peer','oldpeer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 for k,p in D['gate']['normalizer_source_references'].items():inputs[k]={'path':p,'sha256':sha(ROOT/p)};D[k]=json.loads((ROOT/p).read_text())
 S,I,P0,Q0,H,K,C,prep,E,Op=[D[k] for k in ['source','interface','proper','equations','first5','second5','membership','preparation','normalizer','Q_operator']]
 data=D['gate'];nu=[{i:Q(c) for i,c in enumerate(p) if Q(c)} for p in data['actual_nu_A_B_Rell4']];ring=C16(nu);ring.install()
 first={i:{j:Q(c) for j,c in enumerate(p) if Q(c)} for i,p in data['actual105_Rell4']}
 actual=list(map(ring.flat,data['actual14_C16']));kappa=list(map(ring.flat,data['actual_kappa9_C16']));F12=list(map(ring.flat,data['actual_physicalF12_C16']));assert [actual[j] for j in I['retained_source9_indices']]==kappa
 for i in range(16):
  for j in range(16):
   a,b=divmod(i,4),divmod(j,4);p=ring.mul({a:Q(1)},{b:Q(1)})
   for k in range(16):assert ring.mul(p,{divmod(k,4):Q(1)})==ring.mul({a:Q(1)},ring.mul({b:Q(1)},{divmod(k,4):Q(1)}))
 GL=D['GL']['GL18_source_parameter_positions']
 def family_frame(fam):
  f=S['families'][fam];p=ring.add(ring.tau(0 if fam=='CA' else 1),ring.scale(ring.base(first[f['p_quadratic_relation']['s']]),Q(1,2)));M=V.eye(5)
  for name,(a,b) in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma') and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);M[i][j]=ring.add(M[i][j],ring.base(first[a]),ring.mul(ring.base(first[b]),p))
  F=V.eye(5)
  if fam=='CA':F[0][3]=ring.scale(p,-2);F[1][0]=p;F[1][3]=ring.neg(ring.mul(p,p))
  else:F[0][3]=ring.scale(p,Q(1,2));F[1][0]=p;F[1][3]=ring.scale(ring.mul(p,p),Q(1,4))
  return V.matmul(M,F)
 group=V.eye(5)
 for(i,j),a,b in zip(I['actual_group15_positions'],data['actual_even15_Rell4'],data['actual_odd15_Rell4']):
  aa=ring.base({k:Q(c) for k,c in enumerate(a) if Q(c)});bb=ring.base({k:Q(c) for k,c in enumerate(b) if Q(c)});group[i][j]=ring.add(group[i][j],aa,ring.mul(bb,ring.tau(0)))
 assert V.normalized_flag_frame(V.matmul(group,family_frame('CA')))==V.eye(5);FB=V.normalized_flag_frame(V.matmul(group,family_frame('CB')));assert [FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]==actual[9:]
 rows=P0['zero_coordinate_rows36'];targets=[18,19,28,29,60,67,72,63,64,68,73,74];transform=V.matmul(FB,group);ann=V.constmat(P0['CB_annihilator17_by36']);R17=[[{}]*12 for _ in range(17)]
 for j,target in enumerate(targets):
  vector=[{}]*150;vector[target]=ring.ONE
  for k,row in enumerate(E['eliminated_zero_rows4']):vector[row]=ring.c(E['eliminated_four_coordinate_reconstruction'][k][rows.index(target)])
  moved=V.tensor_push(vector,transform);p=V.matmul(ann,[[moved[i]] for i in rows])
  for i in range(17):R17[i][j]=p[i][0]
 def unpack(mat):return [[ring.flat(p) for p in row] for row in mat]
 assert R17==unpack(D['R']['R17']) and V.matmul(R17,[[p] for p in F12])==[[V.polynomial(p,actual)] for p in Q0['equations14variable_sparse']]
 du=ring.add(ring.ONE,ring.neg(ring.mul(kappa[1],kappa[1])));pu=ring.add(ring.c(5),ring.mul(kappa[1],kappa[1]));roots=[]
 for name in ['v0','v1','b','c','e']:
  rec=H['five_explicit_rational_roots'][name];roots.append(ring.tm(V.polynomial(rec['numerator'],kappa),ring.power(du,-rec['D_denominator_power'])))
 co=[V.synthetic_prefix_coefficients(p,actual,[0,1,10,11,13],roots) for p in Q0['equations14variable_sparse']];Ainv=V.matrix_inverse(co[:5]);B=V.matmul(co[5:],Ainv)
 R12=[[ring.add(R17[i+5][j],ring.neg(ring.add(*(ring.mul(B[i][k],R17[k][j]) for k in range(5))))) for j in range(12)] for i in range(12)]
 for i,p in enumerate(H['all12_actual_remaining_and_divided_records']):
  unit=ring.scale(ring.power(du,p['D_denominator_power']),Q(p['constant_denominator'],p['primitive_content_unit']));R12[i]=[ring.mul(unit,x) for x in R12[i]]
 hpolys=[p['divided_primitive_polynomial'] for p in H['all12_actual_remaining_and_divided_records']];hv=[[V.polynomial(p,kappa)] for p in hpolys]
 assert R12==unpack(D['R']['R12']) and V.matmul(R12,[[p] for p in F12])==[[ring.mul(kappa[1],p[0])] for p in hv]
 dump(o,'phase.json',{'phase':'OWN_FULL17_AND_R12_OFFINTERSECTION_VERIFIED'})
 core=[kappa[j] for j in [0,1,2,6]];axes=[K['source9_variable_order'].index(n) for n in K['selected5_pivot_order']];roots=[]
 for name in K['selected5_pivot_order']:
  rec=K['five_actual_rational_roots'][name];roots.append(ring.tm(V.polynomial(rec['numerator'],core),ring.power(du,-rec['D_power']),ring.power(pu,-rec['P_power'])))
 co=[V.synthetic_prefix_coefficients(p,kappa,axes,roots) for p in hpolys];perm=K['selected5_source_H_rows']+[i for i in range(12) if i not in K['selected5_source_H_rows']];Ainv=V.matrix_inverse([co[i] for i in perm[:5]]);B=V.matmul([co[i] for i in perm[5:]],Ainv);T=[[{}]*12 for _ in range(12)]
 for i in range(5):
  for j in range(5):T[i][perm[j]]=Ainv[i][j]
 for i in range(7):
  T[i+5][perm[i+5]]=ring.ONE
  for j in range(5):T[i+5][perm[j]]=ring.neg(B[i][j])
 for i,p in enumerate(K['remaining7_actual_core_polynomials']):
  unit=ring.scale(ring.tm(ring.power(du,-p['D_denominator_power']),ring.power(pu,-p['P_denominator_power'])),Q(p['QQ_unit_content_numerator'],p['QQ_unit_content_denominator']));T[i+5]=[ring.mul(ring.inv(unit),x) for x in T[i+5]]
 E3=V.eye(12);E3[5][5]=ring.inv(du);E3[6][6]=ring.inv(du);r=ring.add(core[2],ring.neg(core[0]),ring.neg(core[1]));E3[7][5]=ring.neg(ring.tm(r,ring.inv(ring.mul(pu,du))));E3[7][7]=ring.scale(ring.inv(ring.mul(pu,du)),-2);T=V.matmul(E3,T);E4=V.eye(12)
 for j in range(4):
  for i in range(3):E4[8+j][5+i]=ring.neg(ring.tm(ring.power(du,[2,2,2,3][j]),V.polynomial(C['cofactor3_by4'][i][j],core)))
 T=V.matmul(E4,T);N=V.matmul(T,R12);Gs=[ring.add(kappa[j],ring.neg(p)) for j,p in zip(axes,roots)]+[V.polynomial(p,core) for p in prep['generators3']]
 assert T==unpack(D['N']['T']) and N==unpack(D['N']['N']) and Gs==list(map(ring.flat,D['N']['actualG8']))
 assert V.matmul(T,hv)==[[p] for p in Gs]+[[{}] for _ in range(4)] and V.matmul(N,[[p] for p in F12])==[[ring.mul(kappa[1],p)] for p in Gs]+[[{}] for _ in range(4)]
 for row,old in zip(N,D['oldN']['N_C_B44']):
  for p,q in zip(row,old):assert all(p.get((b,k),Q(0))==Q(q[4*b+k]) for b in range(2) for k in range(4))
 assert [[str(p.get((2,k),Q(0))) for k in range(4)] for p in Gs]==D['firstell']['actual_G8_C4_first_ell']
 Sn=[[{}]*32 for _ in range(12)]
 for row,col,g,sign in Op['native_S12_by32']['koszul28_in_first8_rows']:Sn[row][col]=ring.scale(Gs[g],Q(sign))
 for row,col,c in Op['native_S12_by32']['free4_zero_source_columns']:Sn[row][col]=ring.c(c)
 layers=[[{} for _ in range(128)] for _ in range(4)]
 for j in range(32):
  for inp in range(4):
   for i in range(12):
    p=ring.mul(Sn[i][j],{(0,inp):Q(1)})
    for(b,k),c in p.items():layers[b][4*j+inp][4*i+k]=c
 S0,Sr,Se,S2=layers;assert S0==list(map(sp,D['mu_gate']['native_S0_columns'])) and Se==list(map(sp,D['mu_gate']['native_S1ell_columns'])) and list(map(enc,S2))==D['value']['intrinsic_S2rell_columns']
 L=[]
 for b in range(3):
  columns=[{} for _ in range(48)]
  for i in range(12):
   for k in range(4):
    for j in range(12):
     # r/ell monic firsts vanish, so no first-order mask overlap carry.
     p=ring.mul(N[i][j],{(0,k):Q(1)})
     for(base,mask),c in p.items():
      if base==b:columns[4*i+k][4*j+mask]=c
  L.append(columns)
 assert L[0]==list(map(sp,D['mu_gate']['L0_Res4_Ntranspose_columns'])) and list(map(enc,L[1]))==D['value']['actual_L1r_columns'] and list(map(enc,L[2]))==D['value']['actual_L1ell_columns']
 gammas=list(map(sp,D['fold']['two_Gamma_on56_actual_X2']));gl=list(map(sp,D['fold']['two_Gamma_L0_on48_intrinsic_native_S2']));values=[]
 for p in D['fold']['all619_source_folds']:
  c=sp(p['native_intrinsic_S2_source_combination']);assert not mv(S0,c);v=mv(S2,c);rpart=mv(L[1],sp(p['native_Omega_ell']));epart=mv(L[2],sp(p['native_Omega_r']));values.append([str(sum((a*v.get(i,Q(0)) for i,a in gl[k].items()),Q(0))+sum((a*(rpart.get(i,Q(0))+epart.get(i,Q(0))) for i,a in gammas[k].items()),Q(0))) for k in range(2)])
 assert values==D['value']['all619_two_actual_native_parts']
 ep=o/'independent-native-C16-adapter-and-all619-mixed-parts.json';dump(o,ep.name,{'inputs':inputs,'all4096_C16_associativity_identities':True,'same_CA_CB_flag_frame_from_actual105_and_canonical_gauge':True,'all17_and_full_first5_telescopes_R12_offintersection_verified':True,'allT_N_G8_coefficients_and_TH_NF_offintersection_verified':True,'all_N0_rfirst_and_G8_ellfirst_cached_matches':True,'all128_intrinsicS2rell_and48_L1r_L1ell_columns_verified':True,'all619_nativeS0_zero_and_both_native_response_values_verified':True,'native_L1ell_Omega_r_retained_no_individual_effectiveXp_shortcut':True,'no_new_source_gauge_implicit_coefficients':True,'no_complete_combined_HF1_or_T2_length_claim':True})
 report={'id':STEM,'all_passed':True,'all_actual_cached_native_source_and619_parts_accepted':True,'combined_HF1_or_T2_length_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
