#!/usr/bin/env python3
"""Independent whole8psi source/bounds review, without any new value dot."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
def count(x):
 if isinstance(x,int):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return 0
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';P=json.loads(pp.read_text());assert P['all_passed'];inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}};D={}
 for p,h in P['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]={'path':p,'sha256':h}
 assert sha(ROOT/'code/ambient_B18_ell2_T2_active4_cache_bundle_corrected_input_20261105_v63.py')==P['code_sha256'] and sha(ROOT/'research/notes'/f'{A}.md')==P['note_sha256']
 ep=pp.parent/'whole-cache-eight-psi-four-quadratic-classes-source-input.json';assert sha(ep)==P['evidence_sha256'][ep.name];d=json.loads(ep.read_text());inputs['bundle']={'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)}
 def selected(path,keys):return fields(ROOT/path,keys)
 graph='reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63/actual-combined-invariant-and-A23-cache-source-input.json';refs=selected(graph,['existing_cache_paths'])['existing_cache_paths']
 foldpath='reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json';fold=json.loads((ROOT/foldpath).read_text())
 adpath='reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_AD_20261105_v63/actual-eight-old-formal-rell-sensitivities-with-M1-native-retained.json';ad=json.loads((ROOT/adpath).read_text())
 respath='reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_QUADRATIC_INPUT_GATE_20261105_v63/actual-factored-G2-rell-eight-and-G3-quadratic-four-source-operators.json';res=json.loads((ROOT/respath).read_text())
 for key,want in [('eight_source_folds',fold['eight_source_folds']),('eight_old_formal_results',ad['all_eight_formal_old_results']),('Gamma_two',fold['Gamma56_two_channels']),('lambda_two',fold['lambda26_two_channels']),('G2_Kr',res['Pi17_Kfirst_r_columns67']),('G2_Ke',res['Pi17_Kfirst_ell_columns67']),('G2_eight',res['eight_G2_rell_factored_responses']),('G3_four',res['four_G3_quadratic_factored_responses']),('G3_P',res['G3_Pi42_rawD31selected126_linear_columns']),('G3_pair20',res['G3_active_pair_order20'])]:assert d[key]==want
 raw=selected(refs['raw33'],['actual33_basis','coordinate_actions9_sparse_columns']);basis=raw['actual33_basis'];acts=[[sp(c) for c in M] for M in raw['coordinate_actions9_sparse_columns']];table={(i,j):sp(v) for i,j,v in d['raw33_product_table_nonzero']}
 for i,(sector,r,e) in enumerate(basis):
  for j in range(33):
   v={j:Q(1)}
   for coordinate,power in ([(sector-1,1)] if sector else [])+[(7,r),(8,e)]:
    for k in range(power):v=mv(acts[coordinate],v)
   assert v==table.get((i,j),{})
 assert len(table)==150
 nu=selected('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63/phase1-actual-two-invariants-from-raw33-primitive.json',['nu_A_B_full33'])['nu_A_B_full33'];p67=selected(refs['physical12'],['physical12_coefficients'])['physical12_coefficients'][5];assert p67['physical_row']==67
 assert d['primitive_V23_source_raw67']=={'nu_A_B':nu,'bA':p67['bA'],'bB':p67['bB'],'delta':p67['delta'],'mask':3,'V_column':23}
 bars=selected('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json',['actual_two_three_component_bar_rho33'])['actual_two_three_component_bar_rho33'];assert d['two_bar_rho']==bars
 old=selected(refs['actual_old_C0'],['all104_original_C0_values_and_opaque_gradients'])['all104_original_C0_values_and_opaque_gradients'];cache={i:v[3] for i,v,g in old if not g};wanted=ad['unknown_original_C0_second_label_list'];CV=dict(d['literal_C0_second_values14']);assert set(CV)==set(wanted)
 for i,v in CV.items():assert Q(v)==(Q(cache[i]) if i in cache else Q(0)) and (i in cache or i in (7148,7688,7943))
 assert Q(CV[0])==0
 mug=selected(refs['Mugate'],['L0_Res4_Ntranspose_columns','native_S0_columns','native_S1ell_columns']);nat=selected(refs['native'],['intrinsic_S2rell_columns','actual_L1r_columns','actual_L1ell_columns'])
 for key,value in [('native_L0',mug['L0_Res4_Ntranspose_columns']),('native_S0',mug['native_S0_columns']),('native_Se',mug['native_S1ell_columns']),('native_S2rell',nat['intrinsic_S2rell_columns']),('native_Lr',nat['actual_L1r_columns']),('native_Le',nat['actual_L1ell_columns'])]:assert d[key]==value
 monicpath='reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63/actual-native-Rell4-cache-source-input-and-adapter-gate.json';monics=selected(monicpath,['actual_nu_A_B_Rell4'])['actual_nu_A_B_Rell4'];inputs['native_plane_monic_source']={'path':monicpath,'sha256':sha(ROOT/monicpath)}
 assert all(Q(nu[i])==0 for nu in monics for i in (0,1,2))
 g8=selected('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json',['actualG8'])['actualG8'];Sop=selected('reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',['native_S12_by32'])['native_S12_by32'];Sr=[{} for _ in range(128)]
 for row,col,k,sign in Sop['koszul28_in_first8_rows']:
  for inp in range(4):
   for mask,x in enumerate(g8[k][4:8]):
    q=Q(x)
    if q and not mask&inp:Sr[4*col+inp]=add(Sr[4*col+inp],{4*row+(mask|inp):Q(sign)*q})
 assert Sr==list(map(sp,d['native_Sr']))
 S0=list(map(sp,d['native_S0']))
 for rec in d['eight_source_folds']:
  for cn in rec['native_S2rell_two_zero_origin_combinations']:assert not mv(S0,sp(cn))
 piv=selected(refs['mini'],['original27_X_source_columns'])['original27_X_source_columns'];assert piv==d['pivot27_source_indices'];pr={}
 for path,key,startindex in [(refs['old640'],'source_columns640_rows56_sparse_dual9',0),(refs['old128'],'source_columns128_rows48_sparse_dual9',640)]:
  for j,col in enumerate(selected(path,[key])[key]):
   if j+startindex in piv:pr[j+startindex]={row:x for row,p in col if (x:=sp(p).get(8,Q(0)))}
 assert [pr[j] for j in piv]==list(map(sp,d['pivot27_r_columns']))
 olde=selected(refs['old640ell'],['actual_first640_ell_columns'])['actual_first640_ell_columns']
 assert d['pivot27_ell_recipes']==[['literal',olde[j]] if j<640 else ['native',j-640] for j in piv]
 red=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_SOURCE8_M1_CACHE_REDUCTION_20261105_v63/actual-M1-alpha23-beta23-gauge-reduction-and14-cache-bindings.json').read_text());assert d['M1_kappa_by_source']==[[r['extra_active_source'],r['kappa']] for r in red['all_four_M1_and_gauge_reductions']]
 Qcols=list(map(sp,d['Q_two6_columns']));psi=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63/actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json').read_text());assert Qcols==[sp(psi['NF_top2_r_squared']),sp(psi['NF_top2_r_ell'])]
 i,j=d['Q_fixed_minor_rows'];det=Qcols[0].get(i,Q(0))*Qcols[1].get(j,Q(0))-Qcols[0].get(j,Q(0))*Qcols[1].get(i,Q(0));assert det and det==Q(d['Q_fixed_minor_determinant'])
 dual=list(map(sp,d['Q_dual_two6_rows']));expected_dual=[{i:Qcols[1].get(j,Q(0))/det,j:-Qcols[1].get(i,Q(0))/det},{i:-Qcols[0].get(j,Q(0))/det,j:Qcols[0].get(i,Q(0))/det}]
 assert dual==[{k:v for k,v in row.items() if v} for row in expected_dual]
 def dot(a,b):return sum((q*b.get(k,Q(0)) for k,q in a.items()),Q(0))
 for n,vs in enumerate(res['G3_PiY_L48_NF42times20_columns']):
  v=sp(vs);a={k:x for k,row in enumerate(dual) if (x:=dot(row,v))};assert a==sp(d['G3_chi2'][n]) and mv(Qcols,a)==v
 # Complete independent positive-envelope replay. No signed response, Phi,
 # invariant dot or new source coefficient is evaluated in this phase.
 rows=[]
 def emit(name,D,U,why):
  if not U:D=1
  rows.append({'name':name,'D':str(D),'U':str(U),'bits':max(D.bit_length(),U.bit_length()),'why':why});return D,U
 def exact(name,x):x=Q(x);return emit(name,x.denominator,abs(x.numerator),'literal_known_scalar')
 def cls(name,values,l1=False):
  v=list(map(Q,values));den=math.lcm(*(x.denominator for x in v)) if v else 1;us=[abs(x.numerator)*(den//x.denominator) for x in v];return emit(name,den,sum(us) if l1 else max([0]+us),'known_input_LCM_L1' if l1 else 'known_input_LCM_max')
 def mul(name,*bs,times=1):return emit(name,math.prod(x[0] for x in bs),times*math.prod(x[1] for x in bs),'positive_product_sum')
 def sm(name,*bs):
  den=math.lcm(*(b[0] for b in bs));return emit(name,den,sum(b[1]*(den//b[0]) for b in bs),'positive_LCM_sum')
 def un(name,*bs):
  den=math.lcm(*(b[0] for b in bs));return emit(name,den,max(b[1]*(den//b[0]) for b in bs),'positive_union_max')
 def values(cols):return [x for c in cols for i,x in c]
 tableE=cls('table',[x for i,j,c in d['raw33_product_table_nonzero'] for k,x in c]);nuE=un('nuL1',*(cls('nu'+str(k),[x for i,x in v],True) for k,v in enumerate(nu)));bE=un('bL1',cls('bA',[x for i,x in p67['bA']],True),cls('bB',[x for i,x in p67['bB']],True));deltaE=cls('delta',[x for i,x in p67['delta']]);phiE=un('primitive23_component',mul('nu_times_b',nuE,bE,tableE),deltaE)
 invE=[mul('inv23_'+str(i),cls('barL1_'+str(i),[x for c in row for k,x in c],True),phiE) for i,row in enumerate(bars)];oldE=cls('old14',CV.values());mats={k:cls(k,values(d[k])) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell']};Gamma=list(map(sp,d['Gamma_two']));gE=[cls('GammaL1_'+str(k),v.values(),True) for k,v in enumerate(Gamma)]
 prE=cls('pivotR',values(d['pivot27_r_columns']));oeE=cls('oldpivotE',[x for mode,v in d['pivot27_ell_recipes'] if mode=='literal' for i,x in v]);npE=sm('nativePivotE',mul('L0Se',mats['native_L0'],mats['native_Se'],times=48),mul('LeS0',mats['native_Le'],mats['native_S0'],times=48));peE=un('allPivotE',oeE,npE);psiB=[]
 for n,(f,a) in enumerate(zip(d['eight_source_folds'],d['eight_old_formal_results'],strict=True)):
  ob=mul('oldpart'+str(n),cls('oldWeightL1'+str(n),[c for i,c in a['literal_original_C0_second_rell_sensitivity']],True),oldE);ob=sm('oldWithFirstCross'+str(n),ob,exact('knownCross'+str(n),a['known_first_cross_response_only']));native=[];cor=[]
  for ch in range(2):
   ce=cls(f'cL1_{n}_{ch}',[x for i,x in f['native_S2rell_two_zero_origin_combinations'][ch]],True);nx=sm(f'nativeMatrix_{n}_{ch}',mul(f'L0S2_{n}_{ch}',mats['native_L0'],mats['native_S2rell'],times=48),mul(f'LrSe_{n}_{ch}',mats['native_Lr'],mats['native_Se'],times=48),mul(f'LeSr_{n}_{ch}',mats['native_Le'],mats['native_Sr'],times=48));native.append(mul(f'nat_{n}_{ch}',gE[ch],ce,nx));ar=cls(f'aprL1_{n}_{ch}',[x for i,x in f['known_Xp1_correction_a_r_two_channels'][ch]],True);ae=cls(f'apeL1_{n}_{ch}',[x for i,x in f['known_Xp1_correction_a_ell_two_channels'][ch]],True);cor.append(sm(f'pivot_{n}_{ch}',mul(f'pr_{n}_{ch}',gE[ch],prE,ae),mul(f'pe_{n}_{ch}',gE[ch],peE,ar)))
  kap=Q(dict(d['M1_kappa_by_source'])[f['extra_active_source']]);mi=mul('M1_'+str(n),exact('kappa_'+str(n),kap),invE[0 if f['contraction']=='A' else 1]);psiB.append(sm('psi_'+str(n),ob,*native,*cor,mi,*[tuple(map(int,v)) for v in res['positive_envelopes_G2_channel16_then_G3_four'][2*n:2*n+2]]))
 P3=list(map(sp,d['G3_P']));chi=list(map(sp,d['G3_chi2']));pairs={tuple(p):i for i,p in enumerate(d['G3_pair20'])};coords=[]
 for n,r in enumerate(d['G3_four']):
  for output in range(2):
   terms=[]
   for k,xs in r['selected126_inverse_linear_coefficients']:
    source,u=divmod(k,6);x=Q(xs)
    for i,y in P3[source].items():
     row,v=divmod(i,6);pair=tuple(sorted((u,v)))
     if pair!=(5,5) and chi[20*row+pairs[pair]].get(output):terms.append(abs(x*y*chi[20*row+pairs[pair]][output]))
   coords.append(sm(f'coord_{n}_{output}',cls(f'G3coordTermL1_{n}_{output}',terms,True),psiB[2*n+output]))
 DG=math.lcm(*(b[0] for b in coords));UG=max(b[1]*(DG//b[0]) for b in coords);rankbound=3*UG.bit_length()+3;qE=cls('Qbasis',values(d['Q_two6_columns']));coordE=un('allCoords',*coords);mul('all6QuadVectors',qE,coordE,times=2)
 assert rows==d['source_positive_height_bounds'] and len(rows)==281
 assert str(DG)==d['coordinate_global_positive_D'] and str(UG)==d['coordinate_global_abs_integer_U'] and UG.bit_length()==315
 assert rankbound==948==d['integer2_by4_rank_all_crosschecks_height_bound_bits']==d['whole_value_height_bound_bits']==P['whole_value_height_bound_bits'] and max(r['bits'] for r in rows)<=948<2048
 extras=['source_positive_height_bounds','coordinate_global_positive_D','coordinate_global_abs_integer_U','integer2_by4_rank_all_crosschecks_height_bound_bits','whole_value_height_bound_bits','whole_owned_numeric_bound','whole_raw_byte_bound','source_builder_owned_numeric_upper','four_actual_class_formula','no_new_physical_or_implicit_coefficients','value_job_started'];base={k:v for k,v in d.items() if k not in extras}
 owned=count(base)+count(rows)+16000;assert owned==22891==d['whole_owned_numeric_bound']==P['whole_owned_numeric_bound']<60000
 serialized=len((json.dumps(base,separators=(',',':'))+'\n').encode())+len(json.dumps(rows).encode())+12000;raw=serialized+5000*(2*len(str(2**948))+32)+2_000_000;assert raw==5131577==d['whole_raw_byte_bound']==P['whole_raw_byte_bound']<32*1024*1024
 assert d['source_builder_owned_numeric_upper']==26085<60000 and not d['value_job_started'] and d['no_new_physical_or_implicit_coefficients'] is True
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'independent-whole8psi-custody-frame-sign-height-and-ownership.json'
 proof={'signs':'psi=old+native-pivot+M1-G2quad; class2=rawG3coord2-(psiA,psiB). Native L2 drops only the16 checkedS0zero combinations; allL0S2,LrSe,LeSr kept. IndividualnativeXp_ell=L0Se+LeS0 kept in recipe andheight.',
 'frame':'Prior independent full56kappae23/W2rr/rellzero andall378z*m identities cancel normalizer covariance jointly;alpha23/beta23 invariant. No unknowncanonicalmu/A0setzero.',
 'source':'All14oldcallbacks literalcacheorprovedcoordinatezero;raw33primitive67maskAB usesfullmonicsincludingzfirsts,notplanesubstitution. Srcompiledonlyafterprovingnu_r=nu_ell=0; nonzero nu_rell remains in cachedS2.',
 'bound':'All281positiveenvelopes/globalD/U315and948bit integer-crossbound replayed; table1089/sourcefields/840Qembedding/minor checked.22891ownsprunedinput,parsedscalarsand16kworking/outputreserve;worker mustnotreloadbroaderparents. Freshsignedvaluejob stillrequiresrootGO.',
 'scope':'Fouractiveclassesonly;nonzero addsactualquadraticinitialrelations/uppercover,zeroendsfamily. No completeHF2/total length assertion.'}
 dump={'all_literal_bundle_inputs_and_pivot_carries_verified':True,'all1089_raw33_table_products_verified':True,'all16_nativeS0_combinations_zero':True,'all840_Qfactor_and_fixed2minor_verified':True,'all281_height_rows_replayed':True,'whole_height_bits':948,'native_pivot_bound_uses_actual_SUM':True,'fraction_commit_bits':max(r['bits'] for r in rows),'scalar_temporary_bits':3*max(r['bits'] for r in rows)+1,'whole_numeric':owned,'whole_RAW':raw,'proof_and_scope':proof,'no_new_signed_values_or_phi_dots':True};path.write_text(json.dumps(dump,separators=(',',':'))+'\n')
 erratum=ROOT/'reports/AMBIENT_B18_ELL2_T2_ACTIVE4_WHOLE_HEIGHT_ERRATUM_20261105_v63/report.json';assert sha(erratum)=='be755cc7b5c5ac2b3131231274b4a5f1f58a5170efa50477d331e7b86f6b95f6';inputs['old_gate_and_peer_withdrawal_erratum']={'path':str(erratum.relative_to(ROOT)),'sha256':sha(erratum)}
 report={'id':STEM,'all_passed':True,'whole8psi_cache_input_frame_sign_height_ownership_accepted':True,'actual_values_or_rank_computed':False,'corrected_nativePivotE_SUM_semantics_verified':True,'old_gate_and_old_peer_remain_withdrawn':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'bits':948,'numeric':owned,'RAW':raw}))
if __name__=='__main__':main()
