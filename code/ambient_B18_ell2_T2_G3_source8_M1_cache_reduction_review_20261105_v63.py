#!/usr/bin/env python3
"""Independent SOURCE8 folding, formal AD, gauge and literal-cache review."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,scale,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_G3_SOURCE8_M1_CACHE_REDUCTION_20261105_v63';STEM=A+'_REVIEW'
def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 inputs={};reports={}
 for stem,code in [(A,'ambient_B18_ell2_T2_G3_source8_M1_cache_reduction_20261105_v63.py'),('AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63','ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63.py'),('AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_AD_20261105_v63','ambient_B18_ell2_T2_G3_active4_source8_AD_20261105_v63.py')]:
  pp=ROOT/'reports'/stem/'report.json';r=json.loads(pp.read_text());assert r['all_passed'];reports[stem]=r;inputs[stem]={'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}
  assert sha(ROOT/'code'/code)==r['code_sha256'] and sha(ROOT/'research/notes'/f'{stem}.md')==r['note_sha256']
  for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]={'path':p,'sha256':h}
  for p,h in r['evidence_sha256'].items():assert sha(pp.parent/p)==h
 P=json.loads((ROOT/'reports'/A/'actual-M1-alpha23-beta23-gauge-reduction-and14-cache-bindings.json').read_text())
 foldpath=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json';D=json.loads(foldpath.read_text())
 adpath=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_AD_20261105_v63/actual-eight-old-formal-rell-sensitivities-with-M1-native-retained.json';AD=json.loads(adpath.read_text())
 gp=ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 gate=fields(gp,['Mnew0_columns','full_G2_linear234_by102_columns','all741_zero_origin_source_combinations','L0_Res4_Ntranspose_columns','native_S0_columns','native_S1ell_columns','old640_typed_closure_nodes','old640_source_columns'])
 M=list(map(sp,gate['Mnew0_columns']));G2=list(map(sp,gate['full_G2_linear234_by102_columns']));S0=list(map(sp,gate['native_S0_columns']));L=list(map(sp,gate['L0_Res4_Ntranspose_columns']));Se=list(map(sp,gate['native_S1ell_columns']));nativee=[mv(L,v) for v in Se]
 assert not M[23]
 lam=list(map(sp,D['lambda26_two_channels']))
 assert len(lam)==2 and len(G2)==102
 for j,col in enumerate(G2):
  lr=sum((lam[0].get(i//9,Q(0))*v for i,v in col.items() if i%9==7),Q(0))
  le=sum((lam[0].get(i//9,Q(0))*v for i,v in col.items() if i%9==8),Q(0))
  rr=sum((lam[1].get(i//9,Q(0))*v for i,v in col.items() if i%9==7),Q(0))
  assert lr==0 and le==rr
 Gamma=[{j:c for j,col in enumerate(M) if (c:=sum((x*q.get(i,Q(0)) for i,x in col.items()),Q(0)))} for q in lam]
 assert Gamma==list(map(sp,D['Gamma56_two_channels']))
 mu=fields(ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json',['actual_selected_Y0_source_columns','actual_selected_Y0_columns102'])
 indices=mu['actual_selected_Y0_source_columns'];Y0=list(map(sp,mu['actual_selected_Y0_columns102']));combs=[sp(gate['all741_zero_origin_source_combinations'][i]['original768']) for i in indices]
 mini=fields(ROOT/'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',['original27_X_source_columns','free30_V_indices','fixed_bottom30_inverse_columns'])
 C=mini['free30_V_indices'];inv=list(map(sp,mini['fixed_bottom30_inverse_columns']))
 x0data=json.loads((ROOT/'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json').read_text());X0=[{} for _ in range(768)]
 for row,v in enumerate(x0data['full_X0_56_by768']):
  for j,x in enumerate(v):
   q=Q(x)
   if q:X0[j][row]=q
 del x0data
 Xp=[X0[j] for j in mini['original27_X_source_columns']]
 XR=[]
 for path,key in [('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9'),('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json','source_columns128_rows48_sparse_dual9')]:
  F=fields(ROOT/path,[key])
  for column in F[key]:XR.append({row:v for row,p in column if (v:=sp(p).get(8,Q(0)))})
 ell=fields(ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns']);XE=list(map(sp,ell['actual_first640_ell_columns']))+nativee
 Xeff=[];a1=[]
 for j,c in enumerate(combs):
  assert not mv(X0,c) and not mv(S0,{i-640:x for i,x in c.items() if i>=640})
  v=[mv(XR,c),mv(XE,c)];aa=[];ee=[]
  for var,x in zip((7,8),v):
   a=mv(inv,{k:x[row] for k,row in enumerate(C) if x.get(row,Q(0))});a={k:y for k,y in a.items() if k<27};eff=add(x,mv(Xp,a),Q(-1));aa.append(a);ee.append(eff)
   rhs={}
   for k,y in Y0[j].items():rhs=add(rhs,{i//9:b*y for i,b in G2[k].items() if i%9==var})
   assert mv(M,eff)==rhs
  a1.append(aa);Xeff.append(ee)
 # Original scalar roots are retained as function IDs; no r/g division occurs.
 ops=dict(gate['old640_typed_closure_nodes']);B=[dict(v) for v in gate['old640_source_columns']]
 c0path=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';c0=fields(c0path,['honest_D0_relation_columns','actual_C0_basis8_scalar_nodes']);Pcols=c0['honest_D0_relation_columns']
 weights=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63/actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json').read_text())['all_four_pairs_of_formal_psi_contractions']
 by={};pos=0
 for weight in weights:
  for name in ('A','B'):
   w=sp(weight[f'formal_{name}_weights_on99psi2']);source=[{},{}];vr=[{},{}];ve=[{},{}];apr=[{},{}];ape=[{},{}];yy=[{},{}]
   for slot,x in w.items():
    j,ch=divmod(slot,2);source[ch]=add(source[ch],combs[j],x);vr[ch]=add(vr[ch],Xeff[j][0],x);ve[ch]=add(ve[ch],Xeff[j][1],x);apr[ch]=add(apr[ch],a1[j][0],x);ape[ch]=add(ape[ch],a1[j][1],x);yy[ch]=add(yy[ch],Y0[j],x)
   old={};native=[]
   for ch,c in enumerate(source):
    cn={i-640:x for i,x in c.items() if i>=640};assert not mv(S0,cn);native.append(cn)
    for j,x in c.items():
     if j>=640:continue
     for row,a in Gamma[ch].items():
      if row<48:
       root=B[j].get(row)
       if root is None:continue
       op=ops[root]
       if op[0]=='c':continue
       key=(1,op[1]) if op[0]=='actual_C0_scalar_function' else (0,root)
      else:key=(1,Pcols[j][row-47])
      old=add(old,{key:a*x})
   r=D['eight_source_folds'][pos];pos+=1;assert r['extra_active_source']==weight['source_extra_index'] and r['contraction']==name
   for key,v in [('source_original768_two_channels',source),('native_S2rell_two_zero_origin_combinations',native),('M1_known_Xeff_r_two_channels',vr),('M1_known_Xeff_ell_two_channels',ve),('known_Xp1_correction_a_r_two_channels',apr),('known_Xp1_correction_a_ell_two_channels',ape),('pure_known_G2quadratic_Y0_two_channels',yy)]:assert v==list(map(sp,r[key]))
   assert old=={tuple(k):Q(v) for k,v in r['old_typed_or_literal_C0_root_functional']};by[(weight['source_extra_index'],name)]=(vr,ve)
 kappas=[Q(3,55),Q(0),Q(15,44),Q(3,55)]
 for extra,k in zip([77,101,117,125],kappas):
  ar,ae=by[extra,'A'];br,be=by[extra,'B'];v={23:k} if k else {};nv=scale(v,Q(-1))
  assert ar==[{},{}] and ae==[v,{}] and br==[nv,{}] and be==[{},v]
  assert not ar[0] and not add(ae[0],br[0]) and not be[0]
 # Independent first/second symbolic arithmetic on the pruned 12-node closure.
 atoms={i:list(map(Q,v)) for i,v in D['actual_atom_origin_r_ell_firsts']};vals={};grad={}
 def mul4(a,b):return [a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[2]*b[0],a[0]*b[3]+a[3]*b[0]+a[1]*b[2]+a[2]*b[1]]
 for i,n in D['needed_typed_nodes']:
  if n[0]=='c':v=[Q(n[1],n[2]),Q(0),Q(0),Q(0)];g={}
  elif n[0]=='actual_C0_scalar_function':v=atoms[n[1]];g={n[1]:Q(1)}
  elif n[0]=='sum':
   v=[Q(0)]*4;g={}
   for j,a,b in n[1]:q=Q(a,b);v=[x+q*y for x,y in zip(v,vals[j])];g=add(g,grad[j],q)
  else:
   assert n[0]=='product';v=[Q(n[1],n[2]),Q(0),Q(0),Q(0)];g={}
   for j in n[3]:g=add(scale(g,vals[j][0]),grad[j],v[0]);v=mul4(v,vals[j])
  vals[i]=v;grad[i]=g
 assert [[i,list(map(str,v))] for i,v in sorted(vals.items())]==AD['typed_known4_lower_first_coefficients']
 labels=set()
 for r,saved in zip(D['eight_source_folds'],AD['all_eight_formal_old_results'],strict=True):
  known=Q(0);unknown={}
  for (kind,i),value in r['old_typed_or_literal_C0_root_functional']:
   q=Q(value)
   if kind:unknown=add(unknown,{i:q})
   else:known+=q*vals[i][3];unknown=add(unknown,grad[i],q)
  assert known==Q(saved['known_first_cross_response_only']) and unknown==sp(saved['literal_original_C0_second_rell_sensitivity']);labels|=set(unknown)
 assert sorted(labels)==P['all14_old_second_labels']==AD['unknown_original_C0_second_label_list'] and len(labels)==14
 oldpath=ROOT/'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_INPUT_GATE_20261105_v63/actual104-old-C0-ancestor-cache-bindings-and-Rell4-gate.json';old=json.loads(oldpath.read_text());oldlabels=set(old['formal_original_C0_second_labels104'])
 assert sorted(labels-oldlabels)==[7148,7688,7943] and len(labels&oldlabels)==11 and old['literal_zero_node0']
 assert [c0['actual_C0_basis8_scalar_nodes'][i] for i in (0,4,5)]==[7148,7688,7943]
 actual=fields(ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',['families'])['families']['CA'];gens=actual['actual_zero_odd_I0_generators']+actual['actual_positive_odd_J0_generators'];Y=[gens[i][1] for i in [0,1,2,3,4,5,6,8,9]]
 assert c0['actual_C0_basis8_scalar_nodes'][:7]==Y[:7]
 rho=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json').read_text());assert rho['combined_two_rho_factor_raw33'] and not rho['failures']
 for stem in ['AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW','AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63_REVIEW','AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW']:
  p=ROOT/'reports'/stem/'report.json';assert json.loads(p.read_text())['all_passed'];inputs[stem]={'path':str(p.relative_to(ROOT)),'sha256':sha(p)}
 proof={
  'M1':'A= kappa*lambda0*M1_r*e23; B=-kappa*(lambda0*M1_ell*e23-lambda1*M1_r*e23). All full56 first arrays and source folds independently rebuilt.',
  'gauge':'lambda0G2_r=0 and lambda0G2_ell=lambda1G2_r imply alpha23,beta23 unchanged by M1->M1+G2_1A0 for arbitrary actual admissible A0. No canonical unknown is assigned a value.',
  'W2':'W2rr=A.vr0=0 and W2rell=A.ve0+B.vr0=0. Every remaining term has a z factor or ell². Projection to bottom3 and any constant normalizer preserve z factors, hence known bottom NF and normalizer covariance NF vanish jointly.',
  'same_rho':'For w=e23, M0w=0. Thus F=Mw has F0=0 and obeys the same universal exact G1F=Vw equations used to construct the mixed rho, including nuisance F2. The already independently verified componentwise allorder raw33 factor applies unchanged.',
  'primitive':'Column23=4*5+3 of Q48, namely maskAB of physical67. In the accepted C4 trace frame this is (-nuB*bB67,nuA*bA67,delta67), nu=t+s_trace²/4. No primitive Phi coefficient or dot is evaluated.',
  'coordinates':'Actual C0 source construction basis=Y[:7]+[gA]; literal basis roots and CA source generator IDs checked, not display names. Roots7148/7688/7943 are zA0/zA4/zA5 and their rell coefficient is zero identically.',
  'scope':'All14 old callbacks are eligible, and two new invariant dots are source-defined. Native/Xp1 and contracted G2/G3 finite responses plus whole storage/height source gates remain before values. No HF2/length claim.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'independent-source8-M1-gauge-and-literal-cache-reduction.json'
 ep.write_text(json.dumps({'all99_firsts_fullM0_equations_and_eight_source_folds':True,'all_eight_nativeS0_combinations_zero':True,'kappa':list(map(str,kappas)),
                           'all_lambda102_gauge_identities':True,'all12_formal_old_nodes_and58_sensitivities':True,'old_labels14_literal_cache_binding':True,
                           'source_and_scope_proof':proof,'actual_invariant_or_second_values_not_evaluated':True},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source8_M1_two_invariants_W2_and_old_cache_source_accepted':True,'actual_new_second_or_invariant_values_computed':False,
         'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'kappa':list(map(str,kappas))}))
if __name__=='__main__':main()
