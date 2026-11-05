#!/usr/bin/env python3
"""Independent whole8psi cache postcheck; no producer/rank search."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,mv,scale
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63';STEM=A+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';P=json.loads(pp.read_text());assert P['all_passed'];inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}}
 for p,h in P['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]={'path':p,'sha256':h}
 assert sha(ROOT/'code/ambient_B18_ell2_T2_active4_cache_bundle_20261105_v63.py')==P['code_sha256'] and sha(ROOT/'research/notes'/f'{A}.md')==P['note_sha256']
 for name,h in P['evidence_sha256'].items():assert sha(pp.parent/name)==h
 gp=ROOT/'reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/report.json';G=json.loads(gp.read_text());assert G['all_passed'] and sha(gp)=='e196d39bab86c6f15073bd655eec3a49f4947bd2ebf6c2dc329fa63336bf3672'
 ip=gp.parent/'whole-cache-eight-psi-four-quadratic-classes-source-input.json';assert sha(ip)==G['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 peer=ROOT/'reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63_REVIEW/report.json';assert sha(peer)=='404d8f2d97ca22cd228fb93e2266d5a26ab2e96082dcf7681a31bf5145fdbfa0';inputs['fresh_corrected_source_peer']={'path':str(peer.relative_to(ROOT)),'sha256':sha(peer)}
 maxbits=0
 def check(x):
  nonlocal maxbits
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length()) if isinstance(x,Q) else abs(x).bit_length();maxbits=max(maxbits,b);assert b<=948
  return x
 def dot(a,b):
  v=Q(0)
  for i,x in a.items():v=check(v+x*b.get(i,Q(0)))
  return v
 table={(i,j):sp(v) for i,j,v in d['raw33_product_table_nonzero']}
 def product(a,b):
  v={}
  for i,x in a.items():
   for j,y in b.items():
    if (i,j) in table:v=add(v,table[i,j],check(x*y))
  for x in v.values():check(x)
  return v
 primitive=d['primitive_V23_source_raw67'];nu=list(map(sp,primitive['nu_A_B']));assert primitive['mask']==3 and primitive['V_column']==23
 phi=[scale(product(nu[1],sp(primitive['bB'])),Q(-1)),product(nu[0],sp(primitive['bA'])),sp(primitive['delta'])]
 inv=[check(sum((dot(sp(bar),v) for bar,v in zip(row,phi,strict=True)),Q(0))) for row in d['two_bar_rho']]
 saved1=json.loads((pp.parent/'phase1-actual-alpha23-beta23-and-primitive23.json').read_text());assert phi==list(map(sp,saved1['actual_primitive_V23_raw33'])) and inv==list(map(Q,saved1['actual_alpha23_beta23']))
 L0,Lr,Le,S0,Sr,Se,S2=[list(map(sp,d[k])) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell']];Gamma=list(map(sp,d['Gamma_two']));lam=list(map(sp,d['lambda_two']));CV={i:Q(x) for i,x in d['literal_C0_second_values14']}
 R=list(map(sp,d['pivot27_r_columns']));E=[sp(v) if mode=='literal' else add(mv(L0,Se[v]),mv(Le,S0[v])) for mode,v in d['pivot27_ell_recipes']]
 Kr=list(map(sp,d['G2_Kr']));Ke=list(map(sp,d['G2_Ke']));kappas={i:Q(x) for i,x in d['M1_kappa_by_source']};psi=[];pieces=[]
 for n,(f,ad,g2) in enumerate(zip(d['eight_source_folds'],d['eight_old_formal_results'],d['G2_eight'],strict=True)):
  assert (f['extra_active_source'],f['contraction'])==(ad['extra_active_source'],ad['contraction'])==(g2['extra_active_source'],g2['contraction'])
  old=Q(ad['known_first_cross_response_only'])+sum((Q(x)*CV[i] for i,x in ad['literal_original_C0_second_rell_sensitivity']),Q(0));native=Q(0);pivot=Q(0);g2sum=Q(0);nchannels=[];gchannels=[]
  for ch in range(2):
   c=sp(f['native_S2rell_two_zero_origin_combinations'][ch]);assert not mv(S0,c)
   v=add(add(mv(L0,mv(S2,c)),mv(Lr,mv(Se,c))),mv(Le,mv(Sr,c)));z=dot(Gamma[ch],v);native=check(native+z);nchannels.append(str(z))
   ar=sp(f['known_Xp1_correction_a_r_two_channels'][ch]);ae=sp(f['known_Xp1_correction_a_ell_two_channels'][ch]);pivot=check(pivot+dot(Gamma[ch],add(mv(R,ae),mv(E,ar))))
   channel=g2['channels'][ch];assert not channel['mixed51_source_coefficients'];v=add(mv(Kr,sp(channel['active_inverse_first67_ell'])),mv(Ke,sp(channel['active_inverse_first67_r'])));z=check(-dot(lam[ch],v));g2sum=check(g2sum+z);gchannels.append(str(z))
  m1=check(kappas[f['extra_active_source']]*(inv[0] if f['contraction']=='A' else -inv[1]));value=check(old+native-pivot+m1-g2sum);psi.append(value)
  pieces.append({'extra_active_source':f['extra_active_source'],'contraction':f['contraction'],'old_cached':str(old),'native_coherent':str(native),'pivot_to_subtract':str(pivot),'M1_gauge_invariant':str(m1),'G2quadratic_to_subtract':str(g2sum),'psi':str(value),'native_two_channels':nchannels,'G2quadratic_two_channels':gchannels,'all_native_S0_zero_verified':True,'canonical_unknown_mu_not_assigned':True})
 saved=json.loads((pp.parent/'phase2-eight-psi-pieces-and-four-actual-quadratic-classes.json').read_text());assert pieces==saved['all_eight_actual_psi_component_identities'] and inv==list(map(Q,saved['actual_alpha23_beta23']))
 P3=list(map(sp,d['G3_P']));Chi=list(map(sp,d['G3_chi2']));pair={tuple(v):i for i,v in enumerate(d['G3_pair20'])};Qcols=list(map(sp,d['Q_two6_columns']));raw=[];coords=[];actual=[]
 for n,r in enumerate(d['G3_four']):
  v={}
  for k,x in sp(r['selected126_inverse_linear_coefficients']).items():
   source,u=divmod(k,6)
   for i,y in P3[source].items():
    row,w=divmod(i,6);p=tuple(sorted((u,w)))
    if p!=(5,5):v=add(v,Chi[20*row+pair[p]],check(-x*y))
  raw.append(v);c=add(v,{0:psi[2*n],1:psi[2*n+1]},Q(-1));coords.append(c);actual.append(mv(Qcols,c))
 assert raw==list(map(sp,saved['raw_known_G3_two_coordinates'])) and coords==list(map(sp,saved['four_actual_G3_quadratic_coordinates2'])) and actual==list(map(sp,saved['four_actual_G3_quadratic_vectors6']))
 rank=json.loads((pp.parent/'actual-active4-quadratic-fixed-integer-rank-and-factorizations.json').read_text());den=int(d['coordinate_global_positive_D']);up=int(d['coordinate_global_abs_integer_U']);I=[]
 for c in coords:
  col=[]
  for i in range(2):x=c.get(i,Q(0));assert den%x.denominator==0;y=x.numerator*(den//x.denominator);assert abs(y)<=up;col.append(y)
  I.append(col)
 assert I==rank['actual_integer_cleared_matrix2_by4_columns'] and den==int(rank['common_denominator'])
 selected=rank['selected_original_columns'];rows=rank['fixed_original_minor_rows'];r=rank['rank'];assert r==len(selected)==len(rows) and 0<=r<=2
 if r==0:assert all(not any(c) for c in I) and not selected and not rows
 elif r==1:assert I[selected[0]][rows[0]]==int(rank['fixed_integer_minor'])!=0
 else:assert I[selected[0]][0]*I[selected[1]][1]-I[selected[0]][1]*I[selected[1]][0]==int(rank['fixed_integer_minor'])!=0
 factors=list(map(sp,rank['all_four_selected_basis_factorizations']));assert len(factors)==4
 for c,a in zip(coords,factors,strict=True):assert mv([coords[j] for j in selected],a)==c
 for j in rank['kernel_free_source_indices']:
  source={j:Q(1)}
  for k,x in factors[j].items():source=add(source,{selected[k]:-x})
  assert not mv(coords,source)
 assert rank['kernel_free_source_indices']==[j for j in range(4) if j not in selected]
 assert P['actual_four_quadratic_initial_relation_rank']==r and P['actual_HF2_upper_from_this_family']==6-r
 assert inv==list(map(Q,P['actual_alpha23_beta23'])) and P['actual_global_HF2_total_length_or_H0_not_computed']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'independent-all8psi-four-classes-and-fixed-rank.json'
 record={'actual_alpha23_beta23':list(map(str,inv)),'all8_old_native_pivot_M1_G2_componentwise_identities':True,'all4_raw_final2_and6_coordinate_vectors':True,'fixed_integer_rank':r,'all4_basis_factorizations_and_kernel':True,'actual_HF2_upper_from_this_family':6-r,
         'source_frame_signs':'CorrectedSUMgate only. psi=old+native-pivot+M1-G2quad; class=rawG3-psi. Allindividualnativepivot carries included;sourceS0zero verified. No canonicalmu assigned.',
         'scope':'Onlyactive4directionsmodaccepted218zero. No nonactive/cycle closure, HF2equality,totalT2length,H0exponent claim. No producer import or new rank search.'}
 path.write_text(json.dumps(record,separators=(',',':'))+'\n');report={'id':STEM,'all_passed':True,'actual_four_quadratic_initial_relation_rank':r,'actual_HF2_upper_from_this_family':6-r,'actual_HF2_equality_or_total_length_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'largest_checked_rational_bits':maxbits,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'rank':r,'HF2_upper':6-r,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
