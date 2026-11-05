#!/usr/bin/env python3
"""Prepared complete nine-class cache bundle; explicit ROOT GO required."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUES_20261105_v63'
def main(out,gate,peer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 cp=ROOT/custody;cu=json.loads(cp.read_text());assert cu['all_passed'] and cu['safe_help_literal_PASS'] and H.sha(Path(__file__))==cu['worker_sha256']
 for p,h in cu['worker_and_helper_sha256'].items():assert H.sha(ROOT/p)==h
 gp=ROOT/gate;g=json.loads(gp.read_text());pp=ROOT/peer;pr=json.loads(pp.read_text());assert g['all_passed'] and pr['all_passed'] and H.sha(gp) in json.dumps(pr);assert H.sha(gp)==cu['source_gate_sha256'] and not g['value_job_started']
 ip=gp.parent/'whole-old-native-pivot-nine-classes-positive-source-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 for p,h in g['input_sha256'].items():assert H.sha(ROOT/p)==h
 cap=g['whole_value_height_bound_bits'];temp=g['temporary_value_height_bound_bits'];assert cap==362<=384 and temp==613<=1153
 assert count([g,pr,cu])<1000 and 2*count(d)+5000<g['whole_value_numeric_bound']<60000
 H.dump(out/'phase.json',{'phase':'authorized_existing_cache_complete9_values','source_gate':H.sha(gp),'source_peer':H.sha(pp),'custody':H.sha(cp),'source_committed_bits':cap,'temporary_bits':temp,'no_new_source_coefficients':True})
 maximum=0;temporary=0;cursor={};call_index=0
 class Stop(Exception):pass
 def ck(x,commit=True):
  nonlocal maximum,temporary
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length()) if isinstance(x,F) else abs(x).bit_length();temporary=max(temporary,b)
  if b>(cap if commit else temp):cursor.update({'offending_uncommitted_value':str(x),'bits':b});raise Stop('source_positive_guard')
  if commit:maximum=max(maximum,b)
  return x
 def parse(v):return {i:ck(F(x)) for i,x in v}
 def add(v,w,c=F(1)):
  q=dict(v)
  for i,x in w.items():
   cursor['partial_accumulator']=srec(q);cursor['term_index']=i;y=ck(q.get(i,F(0))+ck(c*x,False))
   if y:q[i]=y
   else:q.pop(i,None)
  return q
 def image(A,v):
  nonlocal call_index
  call_index+=1;cursor['arithmetic_call_index']=call_index;cursor['arithmetic_kind']='matrix_image'
  q={}
  for j,x in v.items():cursor['source_column_index']=j;cursor['source_coefficient']=str(x);q=add(q,A[j],x)
  return q
 def dot(row,v):
  nonlocal call_index
  call_index+=1;cursor['arithmetic_call_index']=call_index;cursor['arithmetic_kind']='scalar_dot'
  q=F(0)
  for i,x in row.items():cursor['dot_index']=i;cursor['partial_scalar_accumulator']=str(q);q=ck(q+ck(x*v.get(i,F(0)),False))
  return q
 def envelope(x,b):
  D,U=map(int,b);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
 L0,Lr,Le,S0,Sr,Se,S2=[list(map(parse,d[k])) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell']];Gamma=list(map(parse,d['Gamma_two']));Q=list(map(parse,d['Q_two6_columns']));CV={i:F(x) for i,x in d['literal_old104_second_values']};pR=list(map(parse,d['pivot27_r_columns']));pE=[]
 psi=[];pieces=[];coordinates=[];actual=[];carry={}
 try:
  cursor={'stage':'individual_pivot_firsts'}
  for mode,v in d['pivot27_ell_recipes']:
   if mode=='literal':pE.append(parse(v))
   else:assert mode=='native';pE.append(add(image(L0,Se[v]),image(Le,S0[v])))
  f=d['carry_first_inputs'];C={j:parse(v) for j,v in f['selected99_original768']};R={j:parse(v) for j,v in f['original_r_first_columns']};E={j:parse(v) for j,v in f['original640_ell_first_columns']};inv=list(map(parse,f['fixed_B30_inverse_columns']));free=f['free30_V_indices']
  for k in {k for c in C.values() for k in c if k>=640}:E[k]=image(L0,Se[k-640])
  for j,c in C.items():
   cursor={'stage':'cached34_pivot_chart_firsts','source99':j};assert not image(S0,{k-640:x for k,x in c.items() if k>=640});pair=[]
   for raw in [R,E]:
    v=image(raw,c);a=image(inv,{i:v[row] for i,row in enumerate(free) if row in v});pair.append({i:x for i,x in a.items() if i<27})
   carry[j]=pair
  H.dump(out/'phase1-cached34-pivot-chart-firsts.json',{'all34_native_S0_zero_checked':True,'same_original27_free30_B30_chart':True,'all34_a1_r_ell':[[j,[srec(v) for v in pair]] for j,pair in sorted(carry.items())],'no_new_first_coefficients':True})
  for n,(r,wrec,b) in enumerate(zip(d['old18_formal_results'],f['formal_eighteen_weights'],d['all18_component_positive_bounds'])):
   cursor={'stage':'actual18_psi','index':n};assert (r['original418_source'],r['contraction'])==(wrec['original418_source'],wrec['contraction'])==(b['source418'],b['contraction']);aa=[[{},{}],[{},{}]]
   for k,x in parse(wrec['weights']).items():
    j,ch=divmod(k,2)
    for direction in range(2):aa[ch][direction]=add(aa[ch][direction],carry[j][direction],x)
   old=F(0)
   for label,x in r['exact_independent_original_C0_second_sensitivity']:assert label in CV;old=ck(old+ck(F(x)*CV[label],False))
   envelope(old,b['old']);native=F(0);pivot=F(0);nchannels=[];pchannels=[]
   for ch,cnrec in enumerate(r['exact_original_native128_combinations_two_channels']):
    cn=parse(cnrec);assert not image(S0,cn);vec=add(add(image(L0,image(S2,cn)),image(Lr,image(Se,cn))),image(Le,image(Sr,cn)));nv=dot(Gamma[ch],vec);envelope(nv,b['native_channels'][ch]);native=ck(native+nv);nchannels.append(str(nv))
    ar,ae=aa[ch];pv=dot(Gamma[ch],add(image(pR,ae),image(pE,ar)));envelope(pv,b['pivot_channels'][ch]);pivot=ck(pivot+pv);pchannels.append(str(pv))
   z=ck(old+native-pivot);envelope(z,b['psi']);psi.append(z);pieces.append({'source418':r['original418_source'],'contraction':r['contraction'],'old_cached':str(old),'native_coherent':str(native),'pivot_to_subtract':str(pivot),'M1_identically_zero':True,'G2quadratic_identically_zero':True,'psi':str(z),'native_two_channels':nchannels,'pivot_two_channels':pchannels})
  D,U=int(d['coordinate_common_D']),int(d['coordinate_abs_clear_U']);I=[];sources=[]
  for j in range(9):
   assert pieces[2*j]['contraction']=='A' and pieces[2*j+1]['contraction']=='B' and pieces[2*j]['source418']==pieces[2*j+1]['source418'];sources.append(pieces[2*j]['source418']);v={k:-x for k,x in enumerate(psi[2*j:2*j+2]) if x};coordinates.append(v);actual.append(image(Q,v));integer=[ck(x.numerator*(D//x.denominator)) for x in [-psi[2*j],-psi[2*j+1]]];assert all(D%x.denominator==0 for x in psi[2*j:2*j+2]) and max(map(abs,integer))<=U;I.append(integer)
  H.dump(out/'phase2-actual18-psi-and-nine-quadratic-classes.json',{'all18_exact_component_identities':pieces,'nine_source418_labels':sources,'actual_nine_coordinates2':[srec(v) for v in coordinates],'actual_nine_vectors6':[srec(v) for v in actual],'source_M1_G2_rawG3_bottom_zero':True,'no_new_physical_gauge_or_implicit_coefficients':True})
  cursor={'stage':'fixed_integer2_by9_rank'};selected=[];rows=[];minor=0
  nonzero=[j for j,v in enumerate(I) if any(v)]
  if nonzero:
   first=nonzero[0];selected=[first];row=next(i for i,x in enumerate(I[first]) if x);rows=[row];minor=I[first][row]
   for j in nonzero[1:]:
    det=ck(I[first][0]*I[j][1]-I[first][1]*I[j][0])
    if det:selected.append(j);rows=[0,1];minor=det;break
  rank=len(selected);factors=[];numerators=[]
  for v in I:
   if rank==0:assert not any(v);factors.append({});numerators.append([])
   elif rank==1:
    N=v[rows[0]];assert all(ck(minor*v[i])==ck(I[selected[0]][i]*N) for i in range(2));a=ck(F(N,minor));factors.append({0:a} if a else {});numerators.append([N])
   else:
    a,b=[I[j] for j in selected];N=[ck(b[1]*v[0]-b[0]*v[1]),ck(-a[1]*v[0]+a[0]*v[1])];assert all(ck(a[i]*N[0]+b[i]*N[1])==ck(minor*v[i]) for i in range(2));factors.append({k:ck(F(x,minor)) for k,x in enumerate(N) if x});numerators.append(N)
  kernels=[]
  for j in range(9):
   if j in selected:continue
   w={j:minor if rank else 1}
   for k,x in enumerate(numerators[j]):
    if x:w[selected[k]]=w.get(selected[k],0)-x
   assert all(ck(sum(I[k][i]*x for k,x in w.items()))==0 for i in range(2));kernels.append([[k,str(x)] for k,x in sorted(w.items())])
  ep=out/'actual-nine-quadratic-fixed-integer-rank-and-all-factorizations.json';H.dump(ep,{'actual_integer_matrix2_by9':I,'source418_labels':sources,'common_denominator':str(D),'rank':rank,'selected_original9_columns':selected,'fixed_minor_rows':rows,'fixed_integer_minor':str(minor),'all9_selected_basis_factorizations':[srec(v) for v in factors],'cleared_integer_kernel_vectors':kernels,'kernel_free_block':'minor*I onfreecolumns (I ifrank0); normalizedsourcekernelobtainedbydivisionbythenonzerocertifiedminor','all_integer_factor_and_kernel_identities_verified':True,'six_coordinate_source':'fixedQembeddingminorandcoordinateidentitiesgiveall9six-coordinatefactors;nofakeupperzero','complete_G3_projected393_quadratic_rank':rank,'actual_HF2_upper_from_this_family':6-rank,'actual_HF2_equality_or_total_T2_length_not_claimed':True,'remaining619_mixed_cycle_family_open':True})
 except Stop:
  H.dump(out/'guard-frontier.json',{'cursor':cursor,'compiled_carry':[[j,[srec(v) for v in pair]] for j,pair in sorted(carry.items())],'compiled_individual_pivot_ell':[srec(v) for v in pE],'completed_psi':list(map(str,psi)),'completed_pieces':pieces,'completed_coordinates':[srec(v) for v in coordinates],'no_retry_or_cap_raise':True});print(json.dumps({'all_passed':False,'guard_stop':cursor}));return
 assert count([d,carry,pE,psi,pieces,coordinates,actual,I,factors,kernels])+5000<g['whole_value_numeric_bound'];assert sum(p.stat().st_size for p in out.iterdir())<g['whole_RAW_bound']
 report={'id':STEM,'all_passed':True,'actual_nine_quadratic_rank':rank,'complete_G3_projected393_quadratic_rank':rank,'actual_HF2_upper_bound':6-rank,'actual_HF2_equality_or_total_T2_not_computed':True,'no_new_physical_gauge_implicit_or_jet_coefficients':True,'input_sha256':{gate:H.sha(gp),peer:H.sha(pp),custody:H.sha(cp),str(ip.relative_to(ROOT)):H.sha(ip)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:H.sha(p) for p in out.iterdir() if p.name!='report.json'},'source_committed_bit_bound':cap,'source_temporary_bit_bound':temp,'actual_maximum_committed_bits':maximum,'actual_temporary_peak_bits':temporary,'whole_numeric_bound':g['whole_value_numeric_bound'],'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_nine_quadratic_rank','complete_G3_projected393_quadratic_rank','actual_HF2_upper_bound','actual_maximum_committed_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer,a.custody)
