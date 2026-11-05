#!/usr/bin/env python3
"""Exact nonactive171 quotient/source and first-only lower9 input gate."""
import argparse,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={};priorbits=0
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 # Metadata-only pass precedes all new source combinations. Each item is read
 # once at a time and discarded. The upper includes simultaneous raw/parsed
 # source copies, 9*2*99 b1 records, 2*99*3 first outputs, all171 receipts,
 # positive envelopes and serialization scratch; no dense56*768 X0 is read.
 catalog=[
 ('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63/report.json',None),
 ('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63_REVIEW/report.json',None),
 ('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json',None),
 ('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json',['actual_full_G3linear102_by418_columns','PiY3_from102_columns','PiBottom3_from56_columns']),
 (H.PATHS['G2'],None),
 ('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json',None),
 ('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json',None),
 ('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63/actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json',None),
 ('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json',['actual_selected_Y0_columns102','normalized_basis_rows','normalized_basis_original_selected_source_representations','actual_selected_Y0_source_columns']),
 ('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['all741_zero_origin_source_combinations','L0_Res4_Ntranspose_columns','native_S1ell_columns','native_S0_columns']),
 ('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_HEIGHT_BOUND_20261105_v63/report.json',None),
 (H.PATHS['old640'],['source_columns640_rows56_sparse_dual9']),
 (H.PATHS['old128'],['source_columns128_rows48_sparse_dual9']),
 ('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])]
 input_numeric=0;catalog_counts=[]
 for path,keys in catalog:
  n=count(read(path)) if keys is None else sum(count(field(path,k)) for k in keys)
  input_numeric+=n;catalog_counts.append([path,n])
 source_prior_owned=2*input_numeric+24000
 assert source_prior_owned<150000,('whole_source_prior_owned',source_prior_owned)
 source_prior_raw=2_000_000+source_prior_owned*(2*math.ceil(384*math.log10(2))+40)
 assert source_prior_raw<32*1024*1024,('whole_source_prior_RAW',source_prior_raw)
 H.dump(out/'phase.json',{'phase':'SOURCE_only_before_new_combinations','input_numeric':input_numeric,'source_prior_owned':source_prior_owned,'source_prior_RAW':source_prior_raw,'catalog_counts':catalog_counts,'source_code_sha256':H.sha(Path(__file__))})
 def ck(b):
  nonlocal priorbits
  priorbits=max(priorbits,b[0].bit_length(),b[1].bit_length());return b
 def ex(x):x=F(x);return ck((x.denominator,abs(x.numerator)))
 def bp(a,b):
  D=math.lcm(a[0],b[0]);return ck((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def bm(a,b):return ck((a[0]*b[0],a[1]*b[1]))
 def eclass(values,l1=False):
  values=list(map(F,values));D=math.lcm(*(x.denominator for x in values)) if values else 1;u=[abs(x.numerator)*(D//x.denominator) for x in values];return ck((D,sum(u) if l1 else max([0]+u)))
 active=read('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63/report.json');assert active['all_passed'] and active['actual_four_quadratic_initial_relation_rank']==0
 assert read('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63_REVIEW/report.json')['all_passed']
 g=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json')
 pre='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
 G3=list(map(H.parse,field(pre,'actual_full_G3linear102_by418_columns')));PY=list(map(H.parse,field(pre,'PiY3_from102_columns')));PB=list(map(H.parse,field(pre,'PiBottom3_from56_columns')))
 G2=read(H.PATHS['G2']);G2cols=list(map(H.parse,G2['actual_linear_G2_active102_by42_columns']))
 assert sorted({i//6 for v in G2cols for i in v if i%6 in [4,5]})==[12,13,14,15,16]
 support=sorted({j for j,v in enumerate(G2cols) if any(i%6 in [4,5] for i in v)});assert support==[0,1,29,31,36,38]
 assert [i for i,v in enumerate(PY[:42]) if v]==support and not any(PY[42:])
 zero=[2,3,6];assert all(all(i%9 not in zero for i in H.parse(v)) for v in g['actual_projected_G3linear_columns'][:238])
 assert all(all(i%9 in zero for i in H.parse(v)) for v in g['actual_projected_G3linear_columns'][238:])
 selected=g['selected_columns'];assert len([j for j in selected if j<238])==16
 selected_non=[j for j in selected if j>=238];assert selected_non==[238,239,267,280,281,309,322,323,351]
 records=[];nine=[];structural=[];templates={}
 for j in range(238,418):
  if j in selected_non:continue
  a=H.parse(g['all418_original_basis_factors'][j]);assert all(selected[k]>=238 for k in a)
  c=add({j:F(1)},{selected[k]:x for k,x in a.items()},-1);assert all(k>=238 for k in c)
  if j>=364:
   assert not a;reason='G1zero2_or_F0zero3_fullfirstplane0'
  else:
   iz,local=divmod(j-238,42);assert all((k-238)//42==iz for k in c)
   v={(k-238)%42:x for k,x in c.items()}
   if local in templates:assert templates[local]==v
   else:templates[local]=v
   if local not in [31,36,38]:assert not a and local not in support;reason='G2zero1_column_outside_six_plane_sources'
   else:reason=None
  if reason:
   assert all(i%9 not in [7,8] for i in G3[j]) and not g['actual_projected_G3linear_columns'][j]
   structural.append(j)
  else:nine.append({'source_index':j,'original418_combination':srec(c),'zero_variable':zero[(j-238)//42],'local_source':(j-238)%42})
  records.append({'source_index':j,'original418_combination':srec(c),'structural_zero_reason':reason})
 assert len(records)==171 and len(structural)==162 and [r['source_index'] for r in nine]==[269,274,276,311,316,318,353,358,360]
 assert len(templates)==39
 # Four fixed lower-only duals and all378 z*linear zero normal forms.
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json')
 nf=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json');memo={i:H.parse(v) for i,v in nf['all264_linear_NF_columns']};mi={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 for t in range(6):
  for z in range(7):
   for var in range(9):
    m=[0]*9;m[z]+=1;m[var]+=1;assert not memo[lin['all264_NF_root_nodes'][t][mi[tuple(m)]]]
 std={(t,tuple(m)):i for i,(t,m) in enumerate(lin['degree2_linear_standard_module_monomials'])};rr=[0]*9;rr[7]=2;re=[0]*9;re[7]=re[8]=1
 D=[]
 for ell in [[0,1,4],[11,0,270]]:
  for m in [rr,re]:D.append({std[(3+k,tuple(m))]:F(x) for k,x in enumerate(ell) if x})
 psi=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63/actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json');Q=list(map(H.parse,[psi['NF_top2_r_squared'],psi['NF_top2_r_ell']]))
 assert all(not sum((x*q.get(i,F(0)) for i,x in row.items()),F(0)) for row in D for q in Q)
 cols=[std[(3,tuple(rr))],std[(4,tuple(rr))],std[(3,tuple(re))],std[(4,tuple(re))]];det=F(0)
 for p in itertools.permutations(range(4)):
  sign=(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4));det+=sign*math.prod(D[i].get(cols[p[i]],F(0)) for i in range(4))
 assert det==-121
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
 B=list(map(H.parse,field(mu,'actual_selected_Y0_columns102')));norm={p:H.parse(v) for p,v in field(mu,'normalized_basis_rows')};rep={p:H.parse(v) for p,v in field(mu,'normalized_basis_original_selected_source_representations')}
 rE={p:eclass(v.values(),True) for p,v in rep.items()};BE=eclass([x for v in B for x in v.values()])
 GE=eclass([x for v in G3 for x in v.values()])
 # PRIOR point-coefficient bound, before any new b1 is constructed.
 for r in nine:
  c=H.parse(r['original418_combination'])
  bm(eclass(c.values(),True),GE)
  for var in [7,8]:
   qb={}
   for j,x in c.items():
    for k,y in G3[j].items():
     if k%9==var:qb[k//9]=bp(qb.get(k//9,(1,0)),bm(ex(x),ex(y)))
   aE=(1,0)
   for p in sorted(norm):
    x=qb.pop(p,(1,0))
    if not x[1]:continue
    aE=bp(aE,bm(x,rE[p]))
    for k,y in norm[p].items():
     if k!=p:qb[k]=bp(qb.get(k,(1,0)),bm(x,ex(y)))
   bm(aE,BE)
 point_prior=priorbits;assert point_prior<=384,('point_b1_prior_guard',point_prior)
 H.dump(out/'point-prior.json',{'all_passed':True,'point_b1_and_membership_partial_sum_bound_bits':point_prior,'new_b1_not_constructed_yet':True})
 for r in nine:
  c=H.parse(r['original418_combination']);first=H.guard(mv(G3,c));p=[]
  for var in [7,8]:
   q={i//9:x for i,x in first.items() if i%9==var};original=dict(q);a={}
   assert all(42<=i<93 and (i-42)%17 in [12,13,14,15,16] for i in q)
   while q:
    k=min(q);assert k in norm;x=q[k];q=H.guard(add(q,norm[k],-x));a=H.guard(add(a,rep[k],x))
   assert mv(B,a)==original;p.append(a)
  r['b1_r_selected99']=srec(p[0]);r['b1_ell_selected99']=srec(p[1]);r['all_full102_plane_memberships_verified']=True
 # Existing full firsts: prior source-bound before Hstar1 construction.
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 combs=field(mug,'all741_zero_origin_source_combinations');sel=field(mu,'actual_selected_Y0_source_columns');C=[H.parse(combs[j]['original768']) for j in sel]
 height=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_HEIGHT_BOUND_20261105_v63/report.json');XE=next((int(r['common_denominator']),int(r['abs_cleared_numerator_upper'])) for r in height['bounds'] if r['name']=='all_source_firsts_r_ell')
 PE=(1,0)
 for row in range(3):
  e=eclass([v[row] for v in PB if row in v],True);D0=math.lcm(PE[0],e[0]);PE=(D0,max(PE[1]*(D0//PE[0]),e[1]*(D0//e[0])))
 hbounds=[]
 for c in C:hbounds.append(bm(bm(eclass(c.values(),True),XE),PE))
 # Operational L0*Se, PB projection, and the separately checked S0*cN
 # zero identities are covered BEFORE construction, including partial sums.
 Lbound=eclass([x for v in field(mug,'L0_Res4_Ntranspose_columns') for i,x in v]);Sebound=eclass([x for v in field(mug,'native_S1ell_columns') for i,x in v]);S0bound=eclass([x for v in field(mug,'native_S0_columns') for i,x in v])
 native_inner=bm(ex(48),bm(Lbound,Sebound));bm(native_inner,PE)
 for c in C:bm(eclass((x for j,x in c.items() if j>=640),True),S0bound)
 first_prior=priorbits;assert first_prior<=384,('Hstar_first_prior_guard',first_prior)
 H.dump(out/'first-prior.json',{'all_passed':True,'cached_first_and_membership_partial_sum_bound_bits':first_prior,'Hstar1_not_constructed_yet':True})
 rawr=[]
 for path,key in [(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9'),(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9')]:
  for col in field(path,key):rawr.append(H.guard(mv(PB,{row:H.parse(v)[8] for row,v in col if H.parse(v).get(8)})))
 olde=list(map(H.parse,field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json','actual_first640_ell_columns')));olde=[H.guard(mv(PB,v)) for v in olde]
 L=list(map(H.parse,field(mug,'L0_Res4_Ntranspose_columns')));Se=list(map(H.parse,field(mug,'native_S1ell_columns')));S0=list(map(H.parse,field(mug,'native_S0_columns')))
 eff=[H.guard(mv(PB,mv(L,v))) for v in Se];Hr=[];He=[]
 for c in C:
  cn={j-640:x for j,x in c.items() if j>=640};assert not mv(S0,cn);Hr.append(H.guard(mv(rawr,c)));e={}
  for j,x in c.items():e=H.guard(add(e,olde[j] if j<640 else eff[j-640],x))
  He.append(e)
 # Positive lower4x9 envelopes only; no signed lower products evaluated.
 lowerbounds=[]
 for r in nine:
  br=H.parse(r['b1_r_selected99']);be=H.parse(r['b1_ell_selected99'])
  for ell in [[0,1,4],[11,0,270]]:
   for mixed in [False,True]:
    terms=[]
    for j in set(br)|set(be):
     for k,l in enumerate(ell):
      if mixed:
       terms.extend([abs(F(l)*Hr[j].get(k,F(0))*be.get(j,F(0))),abs(F(l)*He[j].get(k,F(0))*br.get(j,F(0)))])
      else:terms.append(abs(F(l)*Hr[j].get(k,F(0))*br.get(j,F(0))))
    D0=math.lcm(*(v.denominator for v in terms)) if terms else 1;U0=sum(v.numerator*(D0//v.denominator) for v in terms);lowerbounds.append([str(D0),str(U0)])
 DG=math.lcm(*(int(v[0]) for v in lowerbounds));UG=max(int(u)*(DG//int(d)) for d,u in lowerbounds);rankbits=5*UG.bit_length()+8;assert rankbits<2048
 data={'negative_completed_active222':{'primary':active['id'],'sha256':hs['reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63/report.json'],'all_actual_quadratic_classes_zero':True},
  'all171_nonactive_projected_kernel_basis':records,'structural_zero162_source_indices':structural,'remaining9_sources':nine,
  'projected_nonactive_selected9':selected_non,'template39_local_combinations':[[j,srec(v)] for j,v in sorted(templates.items())],
  'structural_zero_proof':'rawPiY0G3higher0foreachnonactivebyexacttensor(currentPiYmixed0;onlyactiveblockz*Ilinear). For162fullfirstplane0, everyunknownPi1/bottomderivative/normalizercorrectionisztimeslinear, killedbyall378NFchecks. No coefficient assigned0.',
  'canonical_gauge_carry_for_remaining9':'Pi1/Ystar1/A0upperquadremainformal; fixedlower4annihilatesentireupperQ2, includinganyconstantnormalizerK*w2. Lowerprojectionisexactly-l(Hstar1*b1) onrr/rell.',
  'lower_four_duals6':[srec(v) for v in D],'Q_two_upper_quadratic_columns':[srec(v) for v in Q],'all_eight_dual_Q_annihilations':True,'fixed_lower_dual4_minor_columns':cols,'fixed_lower_dual4_minor_determinant':str(det),
  'Hstar1_r_selected99_columns':[srec(v) for v in Hr],'Hstar1_ell_selected99_columns':[srec(v) for v in He],
  'Hstar_source_binding':'PiBottom3*actualfullX1ofselected99original741combinations; nativeeffectiveL0SeONLYaftereachS0cN0; oldrbarandfullellcacheusedunchanged; noindividualeffectivepivotshortcut',
  'point_b1_prior_bound_bits':point_prior,'cached_first_construction_prior_bound_bits':first_prior,'source_prior_owned_bound':source_prior_owned,'source_prior_RAW_bound':source_prior_raw,'all_input_numeric':input_numeric,
  'all36_lower_positive_D_U':lowerbounds,'common_integer_lower_D':str(DG),'common_integer_lower_abs_U':str(UG),'integer4by9_rank_bound_bits':rankbits,
  'source_actual_lower_values_evaluated':False,'actual_phi_X2_X3_M1_values_not_queried':True,
  'source_scope':'393G3projectedkernel=active222+nonactive171;162nonactivequadraticzero,only9notyetvalued. Remaining619cyclefamilyseparate.',
  'lower_only_rank_scope':'futurelowerimage rank is a LOWERboundontruequadraticrelationrank andanUPPERboundHF2; nofull6relationorHF2equality/total/nilclaim'}
 owned=count(data)+10000;sourceowned=count([g,G3,PY,PB,G2cols,B,norm,rep,C,rawr,olde,L,Se,S0,eff,Hr,He,data])+8000
 assert owned<60000 and sourceowned<=source_prior_owned<150000,('owned',owned,sourceowned,source_prior_owned)
 raw=len(json.dumps(data).encode())+5000*(2*len(str(2**rankbits))+32)+2_000_000;assert raw<32*1024*1024
 data.update({'future_lower_value_owned_bound':owned,'source_owned_bound':sourceowned,'future_lower_value_RAW_bound':raw})
 ep=out/'exact171-decomposition-162-zero-9-lower-first-source-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'structural_zero_count':162,'remaining_source_count':9,'new_rank_or_actual_lower_value_job':False,
  'point_b1_prior_bits':point_prior,'cached_first_prior_bits':first_prior,'future_lower_rank_height_bits':rankbits,
  'future_lower_owned_bound':owned,'future_lower_RAW_bound':raw,'source_owned_bound':sourceowned,'source_prior_owned_bound':source_prior_owned,'source_prior_RAW_bound':source_prior_raw,'input_numeric':input_numeric,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','structural_zero_count','remaining_source_count','point_b1_prior_bits','cached_first_prior_bits','future_lower_rank_height_bits','future_lower_owned_bound','source_owned_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
