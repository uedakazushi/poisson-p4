#!/usr/bin/env python3
"""Kernel-only five cached functionals; no canonicalhead0longedge is evaluated."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'm0peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63_REVIEW/report.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'minimalpeer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'connecting':'reports/AMBIENT_B18_NIL_CONNECTING16_SOURCE_OPERATOR_20261105_v63/actual-formal-connecting16-and-ten-source-obstructions.json',
 'connectingpeer':'reports/AMBIENT_B18_NIL_CONNECTING16_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json',
 'resolution':'reports/AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63/actual-reduced-Cbar-three-sector-Koszul.json',
 'G3':'reports/AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63/actual-Cprime9-actions-and-genuine-G3-target.json',
 'G3peer':'reports/AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63_REVIEW/report.json',
 'first':'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63/phase1-actual-fullnine-first-V-columns.json',
 'quadratic':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase4-all-twenty-actual-quadratic-values.json',
 'dual':'reports/AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_PREFLIGHT_20261105_v63/actual-curve-source-and-one-residual-cubic-functional.json',
 'curve0':'reports/AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_VALUES_20261105_v63/report.json',
 'curvepeer':'reports/AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_VALUES_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def sp(v):return [[i,c]for i,c in enumerate(v)if c]
def prefix(path,key,n):
 text=path.read_text();token=json.dumps(key);assert text.count(token)==1;i=text.index(':',text.index(token))+1
 while text[i].isspace():i+=1
 assert text[i]=='[';i+=1;start=i;out=[];dec=json.JSONDecoder()
 while len(out)<n:
  while text[i].isspace()or text[i]==',':i+=1
  v,i=dec.raw_decode(text,i);out.append(v)
 return out,{'key':key,'retained_prefix_length':n,'byte_interval':[start,i],'full_array_never_decoded':True}
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 ordinary={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k not in['minimal','connecting','resolution']};assert all(ordinary[k]['all_passed']for k in ordinary if k.endswith('peer'));assert ordinary['curve0']['actual_residual_scalar']=='0'
 mn,mc=prefix(ROOT/PATHS['minimal'],'operator_nodes',9);cn,cc=prefix(ROOT/PATHS['connecting'],'matrix_operator_nodes',11);diffs,dc=prefix(ROOT/PATHS['resolution'],'all8_differentials_sparse_columns',2)
 F=ordinary['m0'];m=[[Q(c)for c in row]for row in F['m0_24_by44']];E=F['selected_source_columns24'];iv=[[Q(c)for c in row]for row in F['selected_minor_inverse']];n=[[Q(0)]*24 for _ in range(44)]
 for j,i in enumerate(E):n[i]=iv[j]
 assert all(sum((m[i][k]*n[k][j]for k in range(44)),Q(0))==Q(i==j)for i in range(24)for j in range(24))
 assert mn[5]==['QQ_matrix_rows',[[str(int(i==j))for j in E]for i in range(44)]]
 assert mn[6]==['matrix_product',[0,5]]and mn[7][0:2]==['origin_unit_matrix_inverse',6]and mn[8]==['matrix_product',[5,7]]
 assert cn[0]==['actual_external_matrix','minimal','matrices.m']and cn[1]==['actual_external_matrix','minimal','matrices.direct_unit_section_n']
 assert cn[7]==['matrix_product',[1,2]]and cn[8]==['horizontal_blocks',[6,7]]and cn[10]==['matrix_product',[8,9]]
 retained=[42,43,50,68,70];assert cn[9]==['QQ_matrix_rows',[[str(int(i==j))for j in retained]for i in range(128)]]
 assert all(row[43]==0 for row in m)
 c0=[Q(int(i==42))-sum((n[i][k]*m[k][42]for k in range(24)),Q(0))for i in range(44)]
 pairs=list(itertools.combinations(range(8),2));wedge=[6,24,26];derivs=[]
 for j,source in zip(wedge,[0,4,5]):
  assert pairs[j]==(source,7);col=diffs[1][j];assert all(sum(exp)>=1 for row,exp,c in col);d=[Q(0)]*24
  for row,exp,c in col:
   if exp==[0]*7+[1]:d[row]+=Q(c)
  assert d==[Q(-int(i==source))for i in range(24)];derivs.append(d)
 assert all(sum(exp)>0 for col in diffs[1]for row,exp,c in col)
 L=[rec(v)for v in ordinary['first']['all44_actual_linear_columns27_sparse']];A=[[[col.get(9*t+a,Q(0))for col in L]for t in range(3)]for a in range(9)];B=A[8]
 M=[[[Q(c)for c in col]for col in mat]for mat in ordinary['G3']['actual_R0_actions9_by9_columns']];H=[[-M[a][t][7]for t in range(3)]for a in range(9)]
 phi=rec(ordinary['dual']['residual_actual_cubic_dual495_sparse']);mons=list(itertools.combinations_with_replacement(range(9),3));b=[phi.get(165*t+mons.index((7,7,7)),Q(0))for t in range(3)];a=[phi.get(165*t+mons.index((7,7,8)),Q(0))for t in range(3)];c=[phi.get(165*t+mons.index((7,8,8)),Q(0))for t in range(3)]
 assert a==[Q(-1,75),Q(0),Q(-14,55)]and b==[Q(92,20625),Q(0),Q(3756,75625)]and c==[Q(0),Q(0),Q(1)]
 assert all(sum((a[t]*B[t][j]+b[t]*A[7][t][j]for t in range(3)),Q(0))==0 for j in range(44))
 cz=[Q(2,25)*H[z][0]+H[z][2]/75 for z in range(7)];assert all(H[z][1]==0 for z in range(7))
 assert all(sum((H[z][t]*B[t][j]for t in range(3)),Q(0))==cz[z]*A[7][0][j]for z in range(7)for j in range(44))
 Q2=[[Q(c)for c in col]for col in ordinary['quadratic']['all44_actual_V2_pureplane_nine_values']]
 rho=[sum((b[t]*Q2[j][3*t]+a[t]*Q2[j][3*t+1]+c[t]*Q2[j][3*t+2]-H[7][t]*B[t][j]for t in range(3)),Q(0))+sum((cz[z]*A[z][0][j]for z in range(7)),Q(0))for j in range(44)]
 ga=Q(ordinary['dual']['gA_changed_r_ell_coefficients'][0]);alpha=Q(11,12)*ga;assert alpha==Q(1375,288)
 s1=[[-alpha*n[i][23]for i in range(44)]]+[[-n[i][j]for i in range(44)]for j in[0,4,5]]
 assert all(B[1][43]==0 for _ in[0])
 data={'literal_direct_section_path':'CONNECTING16.nodes1 -> minimal.matrices.direct_unit_section_n -> minimal.nodes8=E*(mE)^-1; NOT olderhonestnpre',
  'minimal_head_inclusion_rows128':retained,'all_origin_wedge_cycles_zero_by_d2_origin':True,'all24_m0_e43_entries_zero':True,'free43_origin_cycle_is_literal_e43':True,
  'source_m0_24_by44':m,'source_n0_44_by24':n,'head42_origin_cycle44_sparse':sp(c0),'lambda0_source_row44_sparse':sp(B[1]),
  'sigma2_replacement_rows_a_b':[a,b],'all44_aV1ell_plus_bV1r_zero':True,'all7_sigma1z_replacement_scalars':cz,'all308_sigma1z_adjoint_identities_zero':True,
  'remaining_pure_r_first_cycle_adjoint44_sparse':sp(rho),'all4_other_head_pure_r_first_cycles44_sparse':list(map(sp,s1)),
  'free43_m1r_e43_only_row23_value':alpha,'wedge_d2_r_sources':[[j,pairs[j]]for j in wedge],
  'actual_five_scalar_recipe':'a=lambda0*c0_head42; bj=rho*s1_otherheadj. DO NOTpairinSOURCEbuilder.',
  'all_missing_nativeV3_terms_absent_for_four_heads':'Wedgesc0=0. Free43c0=e43, higherδV3e43=peeredGAcurve0; nativecurveV2e43=0inALLdirections. ThusH0*c0crossQ2termsalsovanishe. Canonicalhead42lambda1=beta remainsformal/unused.',
  'unknown_beta_retained':'fixedG3imagehead0=(a,beta),rhead0=(0,xi),xi=-11a/16. Othersheads=(0,bj); otherfirst/higherzero. rankminor a*xi independentofbeta.',
  'kernel14_formula':'hj-bj/xi*(rhead0), j1..4; thenfirstbasisindices5,7,8,9; thenhigher10..15. Removeindices0,6 from16basis; retainedfree-tailI14. Allkernelvectorshavehead0coordinate0, sounknownbeta dropsliterally.',
  'genuine_kernel_source_operator':'acceptedCONNECTING16.actual_ell0_cycles44_by16 timesnewQQ16_by14matrix; no assertion offullR0nil lifts orscalar780 expansion',
  'all_needed_values_already_cached':True,'no_new_source_gauge_jet_or_RegDivell_value_job':True,'actual_five_scalar_values_not_computed':True,
  'future_numeric_cap':60000,'future_static_numeric_bound':30000,'future_input_plus_output_bytes':32<<20,'future_seconds':60,'future_memory_bytes':2<<30,'future_bit_guard':2048,
  'selective_array_custody':[mc,cc,dc]}
 ep=out/'five-cache-functionals-and-genuine-fourteen-kernel-source.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n');assert sum((ROOT/p).stat().st_size for p in PATHS.values())+ep.stat().st_size<32<<20
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_five_cached_scalars_for_actual_Nbar14_kernel','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'actual_edge_scalar_or_kernel_matrix_values_computed':False,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':main()
