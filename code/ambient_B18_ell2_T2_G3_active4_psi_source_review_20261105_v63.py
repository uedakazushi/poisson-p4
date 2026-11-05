#!/usr/bin/env python3
"""Independent active4/eight-psi source and coherent-frame review."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,mv,scale
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63';STEM=A+'_REVIEW'
PROFILE='AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE_PROFILE_20261105_v63'
def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 inputs={};D={}
 for stem,code in [(A,'ambient_B18_ell2_T2_G3_active4_psi_source_20261105_v63.py'),(PROFILE,'ambient_B18_ell2_T2_G3_active4_source_20261105_v63.py')]:
  path=ROOT/'reports'/stem/'report.json';report=json.loads(path.read_text());assert report['all_passed'];inputs[stem]={'path':str(path.relative_to(ROOT)),'sha256':sha(path)}
  assert sha(ROOT/'code'/code)==report['code_sha256'] and sha(ROOT/'research/notes'/f'{stem}.md')==report['note_sha256']
  for p,h in report['input_sha256'].items():
   assert sha(ROOT/p)==h;inputs[p]={'path':p,'sha256':h}
  for p,h in report['evidence_sha256'].items():
   ep=path.parent/p;assert sha(ep)==h;inputs[p]={'path':str(ep.relative_to(ROOT)),'sha256':h};D[p]=json.loads(ep.read_text())
 # The new packets have no coefficient-query premise. Pin existing complete
 # source/unit and normalized-frame proofs, rather than assuming their values.
 for stem in ['AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63_REVIEW']:
  path=ROOT/'reports'/stem/'report.json';assert json.loads(path.read_text())['all_passed'];inputs[stem]={'path':str(path.relative_to(ROOT)),'sha256':sha(path)}
 prefix=ROOT/'reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
 F=fields(prefix,['actual_full_G3linear102_by418_columns','PiY3_from102_columns']);G3=list(map(sp,F['actual_full_G3linear102_by418_columns']));PiY=list(map(sp,F['PiY3_from102_columns']))
 gp=ROOT/'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json';g=json.loads(gp.read_text())
 cp=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63/actual-G3-plane-condition-fixed-rank-and-kernel.json';c=json.loads(cp.read_text())
 active16=[i for i in g['selected_columns'] if i<238];assert len(active16)==16 and g['selected_columns'][:16]==active16
 assert all(set(sp(v))<=set(range(16)) for v in g['all418_original_basis_factors'][:238])
 assert set(active16)<=set(c['selected_original_source_columns'])
 extras=[j for j in c['selected_original_source_columns'] if j not in active16];assert extras==[77,101,117,125]
 data=D['actual-four-active-source-complements-and-unknown-derivative-profile.json'];factor=D['actual-eight-formal-psi-contractions-and-full-derivative-factor-proof.json']
 assert data['active_projected_rank16_sources']==active16 and data['full_active_projected_kernel_dimension']==222 and data['known_plane_kernel_dimension']==218
 mup=ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
 mu=fields(mup,['actual_selected_Y0_columns102','normalized_basis_rows','normalized_basis_original_selected_source_representations','actual_selected_Y0_source_columns'])
 basis=list(map(sp,mu['actual_selected_Y0_columns102']));norm={p:sp(v) for p,v in mu['normalized_basis_rows']};reps={p:sp(v) for p,v in mu['normalized_basis_original_selected_source_representations']};sources=mu['actual_selected_Y0_source_columns']
 assert len(basis)==len(norm)==len(reps)==len(sources)==99 and len(set(sources))==99
 assert sources==data['selected99_unit_source741_indices']==factor['selected99_source741_indices']
 # All selected99 are among the same741 zero-origin source combinations.
 gatepath=ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 gate=fields(gatepath,['all741_zero_origin_source_combinations','native_S0_columns']);inputs['zero_origin_gate']={'path':str(gatepath.relative_to(ROOT)),'sha256':sha(gatepath)}
 S0=list(map(sp,gate['native_S0_columns']));zero_records=gate['all741_zero_origin_source_combinations']
 assert len(zero_records)==741
 zero_combs=[sp(v['original768']) for v in zero_records]
 x0path=ROOT/'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json'
 x0data=json.loads(x0path.read_text());inputs['original_X0']={'path':str(x0path.relative_to(ROOT)),'sha256':sha(x0path)}
 # The saved source matrix has literal dense row order56-by768; preserve it.
 rows=x0data['full_X0_56_by768'];assert len(rows)==56 and all(len(row)==768 for row in rows)
 X0=[{} for _ in range(768)]
 for i,row in enumerate(rows):
  for j,v in enumerate(row):
   q=Q(v)
   if q:X0[j][i]=q
 assert len(X0)==768
 for j in sources:
  combo=zero_combs[j];assert not mv(X0,combo)
  assert not mv(S0,{i-640:a for i,a in combo.items() if i>=640})
 linpath=ROOT/'reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json';lin=json.loads(linpath.read_text());nodes=lin['degree2_monomial_rewrite_nodes'];memo={}
 def nf(k):
  if k in memo:return memo[k]
  node=nodes[k]
  if node[0]=='standard':v={node[1]:Q(1)}
  elif node[0]=='zero_ell2':v={}
  else:
   v={}
   for j,a in node[2]:v=add(v,nf(j),Q(a))
  memo[k]=v;return v
 monidx={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 def monnf(target,u,v):
  if u==v==8:return {}
  m=[0]*9;m[u]+=1;m[v]+=1;return nf(lin['all264_NF_root_nodes'][target][monidx[tuple(m)]])
 unknown=set();rebuilt_records=[]
 for index,source in enumerate(extras):
  weights=sp(g['all418_original_basis_factors'][source]);assert set(weights)<=set(range(16))
  combo=add({source:Q(1)},{active16[k]:a for k,a in weights.items()},Q(-1));first=mv(G3,combo);b={}
  for variable in range(9):
   q={i//9:a for i,a in first.items() if i%9==variable};original=dict(q);assert not mv(PiY,q);coords={}
   for p in sorted(norm):
    a=q.get(p,Q(0))
    if a:q=add(q,norm[p],-a);coords=add(coords,reps[p],a)
   assert not q and mv(basis,coords)==original
   for j,a in coords.items():b[9*j+variable]=a
  rec=data['remaining_four_source_directions'][index];assert rec['source_extra_index']==source
  assert combo==sp(rec['actual_projected_zero_linear_source238_combination']) and first==sp(rec['full102_by9_first_term']) and b==sp(rec['unique_selected99_linear_correction'])
  sensitivity={};zchecks=0
  for k,a in b.items():
   j,variable=divmod(k,9)
   for target in range(3):
    for u in range(9):
     n=monnf(target,variable,u)
     if variable<7 or u<7:assert not n;zchecks+=1
     for output,v in n.items():sensitivity=add(sensitivity,{output*99*27+j*27+target*9+u:-a*v})
  assert sensitivity==sp(rec['unknown_post27_PiY_Ystar1_sensitivity6_by99by3by9']) and zchecks==rec['literal_z_NF_checks']
  for flat in sensitivity:
   rem=flat%(99*27);j,rem=divmod(rem,27);target,u=divmod(rem,9);unknown.add((j,target,u))
  rebuilt_records.append((b,sensitivity))
 assert [list(w) for w in sorted(unknown)]==data['required_unknown_derivative_labels'] and len(unknown)==145
 Pi=list(map(sp,g['Pi_remaining_top_linear_columns']));assert len(Pi)==27 and g['remaining_top_linear_quotient_coordinates']==[25,26]
 QA=monnf(2,7,7);QB=monnf(2,7,8)
 assert QA==sp(factor['NF_top2_r_squared']) and QB==sp(factor['NF_top2_r_ell'])
 for i,(b,expected) in enumerate(rebuilt_records):
  wa={};wb={}
  for j in range(99):
   r=b.get(9*j+7,Q(0));e=b.get(9*j+8,Q(0))
   if r:wa[2*j]=r;wb[2*j+1]=r
   if e:wb[2*j]=e
  weights=factor['all_four_pairs_of_formal_psi_contractions'][i]
  assert weights['source_extra_index']==extras[i] and wa==sp(weights['formal_A_weights_on99psi2']) and wb==sp(weights['formal_B_weights_on99psi2'])
  rebuilt={}
  for q,w in [(QA,wa),(QB,wb)]:
   for slot,a in w.items():
    j,psi=divmod(slot,2)
    for u,p in enumerate(Pi):
     for output,v in q.items():
      x=p.get(psi,Q(0))
      if x:rebuilt=add(rebuilt,{output*99*27+j*27+u:-a*x*v})
  assert rebuilt==expected
 assert factor['formal_derivative_slots_checked']==4*6*99*3*9==64152
 assert not factor['unknown_values_evaluated'] and not factor['new_rank_search']
 proof={
  'post27_identity':'For selected99, Xj0=0 and a0=0. Expanding G2Yj=MXj and Ystar=Yj-Yp*a gives G2_1 Ystar1=M0(Xj2-Xp1*a1)+M1(Xj1-Xp0*a1)-G2_2Yj0. G2_1Yp0=M0Xp1+M1Xp0 retains the top carry; M0Xp0=0 removes a2. M2Xj0=0.',
  'exact_comparison_gauge':'Let eps_M(u,b)=n_M G2u+Hfree_M b and Mprime=M+G2 A. U_A(u,b)=(u+A eps_M(u,b),b) is a formal unit and eps_Mprime U_A=eps_M. It fixes (G3s,0) exactly. Changing the chosen homotopy lift adds G3-boundaries, never assumes homotopy R-linearity.',
  'normalized_frame_transport':'At the constant level after27+99 unit cancellations, the induced gauge is (top,bottom)->(top+PiY0 A0 Hfree0 E_free3 bottom,bottom). Therefore canonical psi must not be placed into the normalized L48 NF alone. Either use normalized Mprime1=M1+G2_1 A0 and its post27 carries, or transport the entire top/bottom quadratic relation and cover by this typed triangular frame.',
  'lift_ambiguity':'Actual homotopy Yprime may differ from Y+A X by G3C. Its first term differs by G3_1 C0 and is killed by PiTop2. Source unit re-elimination differences are linear multiples of known first relations in quadratic NF.',
  'future_inputs':'Eight newly weighted X2/M1/G2_2 contractions still need a source availability and sensitivity gate. All99 source origins and native S0 combinations vanish, so L2S0 drops, but old honest-second sensitivities must be rechecked under the new weights. Canonical other-two mu and A0 are uncomputed, not set zero.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'independent-active4-eight-psi-and-coherent-gauge-proof.json'
 path.write_text(json.dumps({'active_projected_rank':16,'active_projected_kernel_dimension':222,'plane_kernel_dimension':218,'extra_sources':extras,
                             'all_four_full102_by9_selected99_decompositions':True,'all99_fullX0_and_nativeS0_combinations_zero':True,
                             'all64152_formal_derivative_slots_and_eight_weights':True,'formal_derivative_label_count':145,
                             'post27_identity_and_coherent_frame_proof':proof,'actual_eight_values_or_HF2_not_computed':True},separators=(',',':'))+'\n')
 failure=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_PSI_SOURCE_FACTOR_20261105_v63_REVIEW_FIRST_READER_ATTEMPT/failure.json'
 inputs['preserved_first_reader_attempt']={'path':str(failure.relative_to(ROOT)),'sha256':sha(failure)}
 report={'id':STEM,'all_passed':True,'active4_and_eight_formal_psi_factorization_accepted':True,'post27_and_coherent_gauge_identity_accepted':True,
         'numeric_source_availability_or_sensitivity_for_new_weights_closed':False,'actual_psi_values_or_HF2_computed':False,
         'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
         'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'extra_sources':extras,'formal_slots':64152}))
if __name__=='__main__':main()
