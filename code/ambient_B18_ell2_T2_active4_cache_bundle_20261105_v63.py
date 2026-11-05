#!/usr/bin/env python3
"""ONE prepared cache-only eight-psi/four-G3-class value bundle."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63'
def main(out,gate,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized_cache_bundle_all_unknown_source_gates_closed'})
 cp=ROOT/custody;cu=json.loads(cp.read_text());assert cu['all_passed'] and cu['safe_help_literal_PASS']
 assert H.sha(Path(__file__))==cu['worker_sha256']
 for p,h in cu['local_worker_and_helper_sha256'].items():assert H.sha(ROOT/p)==h
 gp=ROOT/gate;gr=json.loads(gp.read_text());assert gr['all_passed'] and not gr['value_job_started']
 assert H.sha(gp)==cu['corrected_gate_sha256']
 ip=gp.parent/'whole-cache-eight-psi-four-quadratic-classes-source-input.json';assert H.sha(ip)==gr['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 pp=ROOT/peer;pr=json.loads(pp.read_text());assert pr['all_passed']
 assert H.sha(gp) in json.dumps(pr), 'fresh_source_peer_must_pin_the_corrected_gate_hash'
 hs={gate:H.sha(gp),str(ip.relative_to(ROOT)):H.sha(ip),peer:H.sha(pp),custody:H.sha(cp)}
 for p,h in gr['input_sha256'].items():assert H.sha(ROOT/p)==h
 bound=gr['whole_value_height_bound_bits'];nonrank=max(r['bits'] for r in d['source_positive_height_bounds']);temporary=max(bound,3*nonrank+1)
 assert bound==948 and temporary<=2048
 maxbits=0;temporary_peak=0;cursor={}
 class Stop(Exception):pass
 def check(x,commit=True):
  nonlocal maxbits,temporary_peak
  if isinstance(x,F):b=max(abs(x.numerator).bit_length(),x.denominator.bit_length())
  else:b=abs(x).bit_length()
  temporary_peak=max(temporary_peak,b)
  if b>(bound if commit else temporary):
   cursor['offending_uncommitted_value']=str(x);cursor['bits']=b;raise Stop('source_height_before_commit')
  if commit:maxbits=max(maxbits,b)
  return x
 def parse(v):return {i:check(F(x)) for i,x in v}
 def axpy(a,b,c=F(1)):
  q=dict(a)
  for i,x in b.items():
   cursor['term_index']=i;cursor['partial_accumulator']=srec(q)
   y=check(q.get(i,F(0))+check(c*x,False))
   if y:q[i]=y
   else:q.pop(i,None)
  return q
 def mv(cols,v):
  q={}
  for j,c in v.items():q=axpy(q,cols[j],c)
  return q
 def dot(a,b):
  z=F(0)
  for i,c in a.items():z=check(z+check(c*b.get(i,F(0)),False))
  return z
 table={(i,j):parse(v) for i,j,v in d['raw33_product_table_nonzero']}
 def product(a,b):
  q={}
  for i,x in a.items():
   for j,y in b.items():
    if (i,j) in table:q=axpy(q,table[i,j],check(x*y,False))
  return q
 bounds={r['name']:(int(r['D']),int(r['U'])) for r in d['source_positive_height_bounds']}
 def verify_envelope(x,name):
  D,U=bounds[name];assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U,(name,str(x))
 inv=[];pieces=[];psi=[];rawg3=[];coordinates=[];actual=[]
 try:
  cursor={'stage':'two_primitive_invariants23'};p=d['primitive_V23_source_raw67'];assert p['mask']==3 and p['V_column']==23
  nu=list(map(parse,p['nu_A_B']));phi=[{i:-x for i,x in product(nu[1],parse(p['bB'])).items()},product(nu[0],parse(p['bA'])),parse(p['delta'])]
  for n,row in enumerate(d['two_bar_rho']):
   z=F(0)
   for t,bar in enumerate(row):z=check(z+dot(parse(bar),phi[t]))
   verify_envelope(z,'inv23_'+str(n));inv.append(z)
  H.dump(out/'phase1-actual-alpha23-beta23-and-primitive23.json',{'nu_A_B_source_cached':p['nu_A_B'],
   'actual_primitive_V23_raw33':[srec(v) for v in phi],'actual_alpha23_beta23':list(map(str,inv)),
   'all_higher_source_handled_by_frozen_universal_raw33_descent':True,'no_individual_new_R0_Taylor_Phi_coefficients_queried':True})
  L0,Lr,Le,S0,Sr,Se,S2=[list(map(parse,d[k])) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell']]
  Gamma=list(map(parse,d['Gamma_two']));lam=list(map(parse,d['lambda_two']));CV={i:F(x) for i,x in d['literal_C0_second_values14']}
  R=list(map(parse,d['pivot27_r_columns']));E=[]
  for mode,v in d['pivot27_ell_recipes']:
   if mode=='literal':E.append(parse(v))
   else:
    assert mode=='native';E.append(axpy(mv(L0,Se[v]),mv(Le,S0[v])))
  Kr=list(map(parse,d['G2_Kr']));Ke=list(map(parse,d['G2_Ke']));kappas={i:F(x) for i,x in d['M1_kappa_by_source']}
  assert all(all(row<17 for row in q) for q in lam)
  for n,(fold,ad,g2) in enumerate(zip(d['eight_source_folds'],d['eight_old_formal_results'],d['G2_eight'])):
   cursor={'stage':'eight_psi_cache_components','source8_index':n};extra=fold['extra_active_source'];kind=fold['contraction']
   assert (extra,kind)==(ad['extra_active_source'],ad['contraction'])==(g2['extra_active_source'],g2['contraction'])
   old=check(F(ad['known_first_cross_response_only']))
   for i,x in ad['literal_original_C0_second_rell_sensitivity']:assert i in CV;old=check(old+check(F(x)*CV[i],False))
   verify_envelope(old,'oldWithFirstCross'+str(n));native=F(0);correction=F(0);g2response=F(0);native_channels=[];g2channels=[]
   for ch in range(2):
    c=parse(fold['native_S2rell_two_zero_origin_combinations'][ch]);assert not mv(S0,c)
    v=axpy(axpy(mv(L0,mv(S2,c)),mv(Lr,mv(Se,c))),mv(Le,mv(Sr,c)))
    z=dot(Gamma[ch],v);verify_envelope(z,f'nat_{n}_{ch}');native=check(native+z);native_channels.append(str(z))
    ar=parse(fold['known_Xp1_correction_a_r_two_channels'][ch]);ae=parse(fold['known_Xp1_correction_a_ell_two_channels'][ch])
    z=dot(Gamma[ch],axpy(mv(R,ae),mv(E,ar)));verify_envelope(z,f'pivot_{n}_{ch}');correction=check(correction+z)
    g=g2['channels'][ch];assert not g['mixed51_source_coefficients']
    v=axpy(mv(Kr,parse(g['active_inverse_first67_ell'])),mv(Ke,parse(g['active_inverse_first67_r'])))
    z=check(-dot(lam[ch],v));g2response=check(g2response+z);g2channels.append(str(z))
   m1=check(kappas[extra]*(inv[0] if kind=='A' else -inv[1]));verify_envelope(m1,'M1_'+str(n))
   z=check(old+native-correction+m1-g2response);verify_envelope(z,'psi_'+str(n));psi.append(z)
   pieces.append({'extra_active_source':extra,'contraction':kind,'old_cached':str(old),'native_coherent':str(native),
    'pivot_to_subtract':str(correction),'M1_gauge_invariant':str(m1),'G2quadratic_to_subtract':str(g2response),
    'psi':str(z),'native_two_channels':native_channels,'G2quadratic_two_channels':g2channels,
    'all_native_S0_zero_verified':True,'canonical_unknown_mu_not_assigned':True})
  P=list(map(parse,d['G3_P']));Chi=list(map(parse,d['G3_chi2']));pair={tuple(p):i for i,p in enumerate(d['G3_pair20'])};Q=list(map(parse,d['Q_two6_columns']))
  D=int(d['coordinate_global_positive_D']);U=int(d['coordinate_global_abs_integer_U'])
  for n,r in enumerate(d['G3_four']):
   cursor={'stage':'four_actual_quadratic_classes','extra':r['extra_active_source']};g={}
   for k,x in parse(r['selected126_inverse_linear_coefficients']).items():
    source,u=divmod(k,6)
    for i,y in P[source].items():
     row,v=divmod(i,6);p=tuple(sorted((u,v)))
     if p!=(5,5):g=axpy(g,Chi[row*20+pair[p]],check(-x*y,False))
   c=axpy(g,{0:psi[2*n],1:psi[2*n+1]},F(-1));q=mv(Q,c)
   for x in c.values():assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
   rawg3.append(g);coordinates.append(c);actual.append(q)
  ep=out/'phase2-eight-psi-pieces-and-four-actual-quadratic-classes.json';H.dump(ep,{
   'actual_alpha23_beta23':list(map(str,inv)),'all_eight_actual_psi_component_identities':pieces,
   'raw_known_G3_two_coordinates':[srec(v) for v in rawg3],'four_actual_G3_quadratic_coordinates2':[srec(v) for v in coordinates],
   'four_actual_G3_quadratic_vectors6':[srec(v) for v in actual],'known_bottom_quad_and_all_normalizer_covariance_NF_zero':True,
   'no_new_physical_implicit_gauge_or_jet_coefficients':True,'actual_global_HF2_or_length_not_yet_claimed':True})
  cursor={'stage':'fixed_rank_integer2_by4'};I=[[check(v.get(i,F(0)).numerator*(D//v.get(i,F(0)).denominator)) for i in range(2)] for v in coordinates]
  assert all(abs(x)<=U for col in I for x in col);selected=[];factors=[];minor=0;rows=[]
  nonzero=[j for j,v in enumerate(I) if any(v)]
  if nonzero:
   j=nonzero[0];selected=[j];row=next(i for i,x in enumerate(I[j]) if x);rows=[row];minor=I[j][row]
   for k in nonzero[1:]:
    z=check(I[j][0]*I[k][1]-I[j][1]*I[k][0])
    if z:selected.append(k);rows=[0,1];minor=z;break
  rank=len(selected)
  for col in I:
   if rank==0:assert not any(col);factors.append({})
   elif rank==1:
    a=check(F(col[rows[0]],minor));assert all(check(minor*col[i])==check(I[selected[0]][i]*col[rows[0]]) for i in range(2));factors.append({0:a} if a else {})
   else:
    a,b=[I[j] for j in selected];N=[check(b[1]*col[0]-b[0]*col[1]),check(-a[1]*col[0]+a[0]*col[1])]
    assert all(check(a[i]*N[0]+b[i]*N[1])==check(minor*col[i]) for i in range(2));factors.append({j:check(F(x,minor)) for j,x in enumerate(N) if x})
  free=[j for j in range(4) if j not in selected]
  H.dump(out/'actual-active4-quadratic-fixed-integer-rank-and-factorizations.json',{
   'actual_integer_cleared_matrix2_by4_columns':I,'common_denominator':str(D),'rank':rank,
   'fixed_original_source_extra_labels':[r['extra_active_source'] for r in d['G3_four']],
   'selected_original_columns':selected,'fixed_original_minor_rows':rows,'fixed_integer_minor':str(minor),
   'all_four_selected_basis_factorizations':[srec(v) for v in factors],'kernel_free_source_indices':free,
   'kernel_factored_recipe':'e_freej-sum_k factors[j][k]*e_selectedk; freeI; zeroimage followsintegercrosschecks',
   'six_coordinate_rank_and_factors_proof':'fixedQ2embeddinghasnonzerominor; all840NFfactoridentities andallfourintegercoordinatefactoridentities implyexact6coordinatefactorizations withoutduplicatedrationalmatrixproducts',
   'actual_HF2_upper_from_this_family':6-rank,'actual_HF2_equality_or_total_length_not_claimed':True,
   'scope':'complete4activeG3directionsmodaccepted218zero; remainingnonactiveG3/cyclefamiliesareNOTclaimedclosed'})
 except Stop as e:
  H.dump(out/'guard-frontier.json',{'stop':str(e),'cursor':cursor,'completed_invariants':list(map(str,inv)),
   'completed_psi_pieces':pieces,'completed_rawG3':[srec(v) for v in rawg3],'completed_coordinates':[srec(v) for v in coordinates],
   'input_path':str(ip.relative_to(ROOT)),'input_sha256':H.sha(ip),'no_retry_or_cap_raise':True});print(json.dumps({'all_passed':False,'guard_stop':cursor}));return
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'actual_alpha23_beta23':list(map(str,inv)),'actual_four_quadratic_initial_relation_rank':rank,
  'actual_HF2_upper_from_this_family':6-rank,'actual_global_HF2_total_length_or_H0_not_computed':True,
  'no_new_physical_implicit_gauge_or_jet_coefficients':True,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),
  'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:H.sha(p) for p in out.iterdir() if p.name!='phase.json'},
  'owned_numeric_source_bound':gr['whole_owned_numeric_bound'],'source_committed_and_integer_rank_bound_bits':bound,
  'source_scalar_temporary_bound_bits':temporary,'actual_committed_height_bits':maxbits,'actual_temporary_peak_bits':temporary_peak,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_alpha23_beta23','actual_four_quadratic_initial_relation_rank','actual_HF2_upper_from_this_family','actual_committed_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer,a.custody)
