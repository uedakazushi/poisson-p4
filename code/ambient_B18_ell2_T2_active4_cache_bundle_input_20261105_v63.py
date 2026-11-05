#!/usr/bin/env python3
"""Whole SOURCE gate and prepared cache-only active4 value bundle."""
import argparse,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;BASE='AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_INPUT_GATE_20261105_v63'
def dot(a,b):return sum((x*b.get(i,F(0)) for i,x in a.items()),F(0))
def gate(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={};boundrows=[]
 H.dump(out/'phase.json',{'phase':'whole_cache_input_and_positive_bounds_only_no_values'})
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 for p in ['reports/AMBIENT_B18_ELL2_T2_G3_SOURCE8_M1_CACHE_REDUCTION_20261105_v63_REVIEW/report.json',
  'reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_QUADRATIC_INPUT_GATE_20261105_v63_REVIEW/report.json',
  'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW/report.json']:
  assert read(p)['all_passed']
 refs=field('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63/actual-combined-invariant-and-A23-cache-source-input.json','existing_cache_paths')
 fold=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json')
 ad=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_AD_20261105_v63/actual-eight-old-formal-rell-sensitivities-with-M1-native-retained.json')
 res=read('reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_QUADRATIC_INPUT_GATE_20261105_v63/actual-factored-G2-rell-eight-and-G3-quadratic-four-source-operators.json')
 red=read('reports/AMBIENT_B18_ELL2_T2_G3_SOURCE8_M1_CACHE_REDUCTION_20261105_v63/actual-M1-alpha23-beta23-gauge-reduction-and14-cache-bindings.json')
 psi=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63/actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json')
 rawbasis=field(refs['raw33'],'actual33_basis');mons=field('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json','monomials660_degree_le4')
 nf=field('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json','all660_actual_raw33_normal_forms');lookup={tuple(m):v for m,v in zip(mons,nf)}
 bm=[]
 for s,r,e in rawbasis:
  m=[0]*9
  if s:m[s-1]=1
  m[7]=r;m[8]=e;bm.append(m)
 table=[]
 for i in range(33):
  for j in range(33):
   m=tuple(a+b for a,b in zip(bm[i],bm[j]));v=[] if m[8]>=2 or sum(m)>=5 else lookup[m]
   if v:table.append([i,j,v])
 nu=field('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63/phase1-actual-two-invariants-from-raw33-primitive.json','nu_A_B_full33')
 physical=field(refs['physical12'],'physical12_coefficients');p67=physical[5];assert p67['physical_row']==67;del physical
 bars=field('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json','actual_two_three_component_bar_rho33')
 oldvals={i:v[3] for i,v,g in field(refs['actual_old_C0'],'all104_original_C0_values_and_opaque_gradients') if not g}
 wanted=ad['unknown_original_C0_second_label_list'];CV={i:oldvals[i] if i in oldvals else '0' for i in wanted};assert set(wanted)-set(oldvals)=={7148,7688,7943};CV[0]='0'
 mug=refs['Mugate'];L0=field(mug,'L0_Res4_Ntranspose_columns');S0=field(mug,'native_S0_columns');Se=field(mug,'native_S1ell_columns')
 S2=field(refs['native'],'intrinsic_S2rell_columns');Lr=field(refs['native'],'actual_L1r_columns');Le=field(refs['native'],'actual_L1ell_columns')
 g8=field('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json','actualG8');sop=field(H.PATHS['S'],'native_S12_by32');Sr=[{} for _ in range(128)]
 for row,col,k,sign in sop['koszul28_in_first8_rows']:
  for inp in range(4):
   for mask,c in enumerate(map(F,g8[k][4:8])):
    if c and not(mask&inp):Sr[4*col+inp]=add(Sr[4*col+inp],{4*row+(mask|inp):F(sign)*c})
 piv=field(refs['mini'],'original27_X_source_columns');pr={}
 for path,key,startindex in [(refs['old640'],'source_columns640_rows56_sparse_dual9',0),(refs['old128'],'source_columns128_rows48_sparse_dual9',640)]:
  for j,col in enumerate(field(path,key)):
   if j+startindex in piv:pr[j+startindex]=[[row,str(H.parse(v)[8])] for row,v in col if H.parse(v).get(8)]
 olde=field(refs['old640ell'],'actual_first640_ell_columns');pivotE=[['literal',olde[j]] if j<640 else ['native',j-640] for j in piv]
 # Fixed source2 basis for every upper quadratic normal form.
 Q=list(map(H.parse,[psi['NF_top2_r_squared'],psi['NF_top2_r_ell']]))
 pair=next((i,j) for i in range(6) for j in range(i+1,6) if Q[0].get(i,F(0))*Q[1].get(j,F(0))-Q[0].get(j,F(0))*Q[1].get(i,F(0)))
 i,j=pair;det=Q[0].get(i,F(0))*Q[1].get(j,F(0))-Q[0].get(j,F(0))*Q[1].get(i,F(0))
 dual=[{i:Q[1].get(j,F(0))/det,j:-Q[1].get(i,F(0))/det},{i:-Q[0].get(j,F(0))/det,j:Q[0].get(i,F(0))/det}]
 chi=[]
 for vs in res['G3_PiY_L48_NF42times20_columns']:
  v=H.parse(vs);a={k:dot(row,v) for k,row in enumerate(dual) if dot(row,v)};assert mv(Q,a)==v;chi.append(srec(a))
 d={'primitive_V23_source_raw67':{'nu_A_B':nu,'bA':p67['bA'],'bB':p67['bB'],'delta':p67['delta'],'mask':3,'V_column':23},
  'raw33_product_table_nonzero':table,'two_bar_rho':bars,'literal_C0_second_values14':[[i,CV[i]] for i in wanted],
  'eight_old_formal_results':ad['all_eight_formal_old_results'],'eight_source_folds':fold['eight_source_folds'],
  'Gamma_two':fold['Gamma56_two_channels'],'lambda_two':fold['lambda26_two_channels'],
  'native_L0':L0,'native_Lr':Lr,'native_Le':Le,'native_S0':S0,'native_Sr':[srec(v) for v in Sr],'native_Se':Se,'native_S2rell':S2,
  'pivot27_source_indices':piv,'pivot27_r_columns':[pr[j] for j in piv],'pivot27_ell_recipes':pivotE,
  'M1_kappa_by_source':[[r['extra_active_source'],r['kappa']] for r in red['all_four_M1_and_gauge_reductions']],
  'G2_Kr':res['Pi17_Kfirst_r_columns67'],'G2_Ke':res['Pi17_Kfirst_ell_columns67'],'G2_eight':res['eight_G2_rell_factored_responses'],
  'G3_four':res['four_G3_quadratic_factored_responses'],'G3_P':res['G3_Pi42_rawD31selected126_linear_columns'],
  'G3_pair20':res['G3_active_pair_order20'],'G3_chi2':chi,'Q_two6_columns':[srec(v) for v in Q],
  'Q_fixed_minor_rows':list(pair),'Q_fixed_minor_determinant':str(det),'Q_dual_two6_rows':[srec(v) for v in dual],
  'all840_upper_NF_factor_Q_chi_checked':True,
  'source_scope':'fouractiveG3quadraticclasses, notall1012initialgenerators. M1normalizerjointlycanceledbyinvariantsandW2plane0; no canonicalmu assigned0.'}
 # Common denominator / absolute numerator envelopes, never signed values.
 def emit(name,D,U,why):
  assert D>0 and U>=0
  if U==0:D=1 # This is an exactly zero envelope, not arbitrary content cancellation.
  boundrows.append({'name':name,'D':str(D),'U':str(U),'bits':max(D.bit_length(),U.bit_length()),'why':why});return D,U
 def exact(name,x):x=F(x);return emit(name,x.denominator,abs(x.numerator),'literal_known_scalar')
 def scalarclass(name,values,l1=False):
  values=list(map(F,values));D=math.lcm(*(x.denominator for x in values)) if values else 1
  us=[abs(x.numerator)*(D//x.denominator) for x in values];return emit(name,D,sum(us) if l1 else max([0]+us),'known_input_LCM_L1' if l1 else 'known_input_LCM_max')
 def mul(name,*bs,times=1):return emit(name,math.prod(x[0] for x in bs),times*math.prod(x[1] for x in bs),'positive_product_sum')
 def summ(name,*bs):
  D=math.lcm(*(x[0] for x in bs));return emit(name,D,sum(x[1]*(D//x[0]) for x in bs),'positive_LCM_sum')
 def union(name,*bs):
  D=math.lcm(*(x[0] for x in bs));return emit(name,D,max(x[1]*(D//x[0]) for x in bs),'positive_union_max')
 def vals(cols):return [x for col in cols for i,x in col]
 tableE=scalarclass('table', [x for i,j,col in table for k,x in col]);nuE=union('nuL1',*(scalarclass('nu'+str(i),[x for k,x in v],True) for i,v in enumerate(nu)))
 bE=union('bL1',scalarclass('bA',[x for i,x in p67['bA']],True),scalarclass('bB',[x for i,x in p67['bB']],True));deltaE=scalarclass('delta',[x for i,x in p67['delta']]);phiE=union('primitive23_component',mul('nu_times_b',nuE,bE,tableE),deltaE)
 invE=[mul('inv23_'+str(i),scalarclass('barL1_'+str(i),[x for c in row for k,x in c],True),phiE) for i,row in enumerate(bars)]
 oldE=scalarclass('old14',CV.values());mats={k:scalarclass(k,vals(d[k])) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell']}
 Gamma=[H.parse(v) for v in d['Gamma_two']];gE=[scalarclass('GammaL1_'+str(i),v.values(),True) for i,v in enumerate(Gamma)]
 pivotR=scalarclass('pivotR',vals(d['pivot27_r_columns']));oldeE=scalarclass('oldpivotE',[x for mode,v in pivotE if mode=='literal' for k,x in v]);newE=union('nativePivotE',mul('L0Se',mats['native_L0'],mats['native_Se'],times=48),mul('LeS0',mats['native_Le'],mats['native_S0'],times=48));pivotEbound=union('allPivotE',oldeE,newE)
 psi_bounds=[]
 for n,(f,a) in enumerate(zip(d['eight_source_folds'],d['eight_old_formal_results'])):
  oldbound=mul('oldpart'+str(n),scalarclass('oldWeightL1'+str(n),[c for i,c in a['literal_original_C0_second_rell_sensitivity']],True),oldE)
  oldbound=summ('oldWithFirstCross'+str(n),oldbound,exact('knownCross'+str(n),a['known_first_cross_response_only']))
  nat=[];cor=[]
  for ch in range(2):
   cE=scalarclass(f'cL1_{n}_{ch}',[x for i,x in f['native_S2rell_two_zero_origin_combinations'][ch]],True)
   nx=summ(f'nativeMatrix_{n}_{ch}',mul(f'L0S2_{n}_{ch}',mats['native_L0'],mats['native_S2rell'],times=48),mul(f'LrSe_{n}_{ch}',mats['native_Lr'],mats['native_Se'],times=48),mul(f'LeSr_{n}_{ch}',mats['native_Le'],mats['native_Sr'],times=48))
   nat.append(mul(f'nat_{n}_{ch}',gE[ch],cE,nx))
   ar=scalarclass(f'aprL1_{n}_{ch}',[x for i,x in f['known_Xp1_correction_a_r_two_channels'][ch]],True);ae=scalarclass(f'apeL1_{n}_{ch}',[x for i,x in f['known_Xp1_correction_a_ell_two_channels'][ch]],True)
   cor.append(summ(f'pivot_{n}_{ch}',mul(f'pr_{n}_{ch}',gE[ch],pivotR,ae),mul(f'pe_{n}_{ch}',gE[ch],pivotEbound,ar)))
  extra=f['extra_active_source'];kap=F(dict(d['M1_kappa_by_source'])[extra]);mi=mul('M1_'+str(n),exact('kappa_'+str(n),kap),invE[0 if f['contraction']=='A' else 1])
  g2bounds=[tuple(map(int,v)) for v in res['positive_envelopes_G2_channel16_then_G3_four'][2*n:2*n+2]]
  psi_bounds.append(summ('psi_'+str(n),oldbound,*nat,*cor,mi,*g2bounds))
 pairindex={tuple(p):i for i,p in enumerate(d['G3_pair20'])};P3=list(map(H.parse,d['G3_P']));Chi=list(map(H.parse,chi));coord_bounds=[]
 for n,r in enumerate(d['G3_four']):
  bs=[]
  for output in range(2):
   terms=[]
   for k,xs in r['selected126_inverse_linear_coefficients']:
    source,u=divmod(k,6);x=F(xs)
    for i,y in P3[source].items():
     row,v=divmod(i,6);p=tuple(sorted((u,v)))
     if p!=(5,5) and Chi[row*20+pairindex[p]].get(output):terms.append(abs(x*y*Chi[row*20+pairindex[p]][output]))
   E=scalarclass(f'G3coordTermL1_{n}_{output}',terms,True);bs.append(summ(f'coord_{n}_{output}',E,psi_bounds[2*n+output]))
  coord_bounds.extend(bs)
 Dglobal=math.lcm(*(b[0] for b in coord_bounds));Uglobal=max(b[1]*(Dglobal//b[0]) for b in coord_bounds);ub=Uglobal.bit_length();rankbound=3*ub+3
 qE=scalarclass('Qbasis',vals(d['Q_two6_columns']));coordE=union('allCoords',*coord_bounds);fullE=mul('all6QuadVectors',qE,coordE,times=2)
 height=max([rankbound]+[r['bits'] for r in boundrows]);assert height<2048,('whole_value_height',height,ub)
 owned=count(d)+count(boundrows)+16000;assert owned<60000,('whole_value_owned',owned)
 source_owner=count([d,fold,ad,res,red,psi,rawbasis,mons,nf,nu,p67,bars,oldvals,L0,S0,Se,S2,Lr,Le,g8,sop,Sr,piv,pr,olde,pivotE,Q,dual,chi,boundrows])+6000
 assert source_owner<60000,('whole_source_builder_owned',source_owner)
 serialized_input_bytes=len((json.dumps(d,separators=(',',':'))+'\n').encode())+len(json.dumps(boundrows).encode())+12000;raw=serialized_input_bytes+5000*(2*len(str(2**height))+32)+2_000_000;assert raw<32*1024*1024
 d.update({'source_positive_height_bounds':boundrows,'coordinate_global_positive_D':str(Dglobal),'coordinate_global_abs_integer_U':str(Uglobal),
  'integer2_by4_rank_all_crosschecks_height_bound_bits':rankbound,'whole_value_height_bound_bits':height,
  'whole_owned_numeric_bound':owned,'whole_raw_byte_bound':raw,
  'source_builder_owned_numeric_upper':source_owner,
  'four_actual_class_formula':'G3raw_coordinates2 - [psi_A,psi_B]; embedwithQ; knownbottomNF0; allsignsfromsource8post27andfiniteinverseproofs',
  'no_new_physical_or_implicit_coefficients':True,'value_job_started':False})
 ep=out/'whole-cache-eight-psi-four-quadratic-classes-source-input.json';H.dump(ep,d)
 report={'id':BASE,'all_passed':True,'value_job_started':False,'new_physical_or_implicit_coefficients':False,
  'whole_owned_numeric_bound':owned,'whole_raw_byte_bound':raw,'whole_value_height_bound_bits':height,
  'coordinate_integer_upper_bits':ub,'all840_Q_embedding_identities':True,'input_sha256':hs,
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{BASE}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_owned_numeric_bound','whole_raw_byte_bound','whole_value_height_bound_bits','coordinate_integer_upper_bits','seconds']}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();gate(a.output.resolve())
