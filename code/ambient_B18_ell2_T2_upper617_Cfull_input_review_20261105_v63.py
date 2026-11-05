#!/usr/bin/env python3
"""Fresh full-coordinate/Psi and Cfull prior/code peer, no signed C construction."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,math,resource,signal,subprocess,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import fields,sp,count,sha
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,d):(o/n).write_text(json.dumps(d,separators=(',',':'))+'\n')
 def manifest(p,expected=None):
  r=rd(p);assert r['all_passed']
  if expected:assert pins[p]==expected
  for k,h in r.get('input_sha256',{}).items():
   if isinstance(h,dict):assert sha(ROOT/h['path'])==h['sha256']
   else:assert sha(ROOT/k)==h
  for k,h in r.get('evidence_sha256',{}).items():assert sha((ROOT/p).parent/k)==h
  return r
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/'
 br=manifest(b+'report.json','b4d00c49577561f7fddb4e928ca905adbc40699703a96dbfd9e204071cec9edb')
 bd=rd(b+'known-Psi-point-correction-and-full-b1-source-annihilation.json')
 g='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/'
 gr=manifest(g+'report.json','7ac51097df48c4416f377fa94380f3abb910f6f10bae77945f96e93ad76b30b4')
 d=rd(g+'pruned-known-Cfull224-input-and-positive-source-envelopes.json')
 put('phase.json',{'phase':'full_Psi_and_Cfull_source_before_any_signed_C','b4d_sha256':pins[b+'report.json'],'gate_sha256':pins[g+'report.json']})
 cache='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json'
 ks=['Gamma_two','pivot27_source_indices','pivot27_r_columns','pivot27_ell_recipes','native_L0','native_Le','native_S0','native_Se','Q_two6_columns','Q_dual_two6_rows']
 raw=fld(cache,ks)
 for k in ks[:8]:assert raw[k]==d[k]
 Qcols=list(map(sp,raw['Q_two6_columns']));L=list(map(sp,raw['Q_dual_two6_rows']))
 assert len(Qcols)==len(L)==2
 def dot(a,b):return sum((x*b.get(i,Q(0)) for i,x in a.items()),Q(0))
 assert [[dot(row,col) for col in Qcols] for row in L]==[[Q(1),Q(0)],[Q(0),Q(1)]]
 assert bd['Q_two6']==raw['Q_two6_columns'] and bd['Q_left_inverse_two6_rows']==raw['Q_dual_two6_rows']
 lin='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json'
 lf=fld(lin,['degree2_monomial_order44','degree2_linear_standard_module_monomials','all264_NF_root_nodes'])
 nf='reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json'
 N={int(i):sp(v) for i,v in fld(nf,['all264_linear_NF_columns'])['all264_linear_NF_columns']}
 rr=[0]*9;rr[7]=2;re=[0]*9;re[7]=re[8]=1
 for j in range(3):
  for k,m in enumerate((rr,re)):
   coord=2*j+k;target=3+j;mi=lf['degree2_monomial_order44'].index(m)
   assert lf['degree2_linear_standard_module_monomials'][coord]==[target,m] and N[lf['all264_NF_root_nodes'][target][mi]]=={coord:Q(1)}
 psi_rr=[{i:row[2*i] for i in range(3) if row.get(2*i)} for row in L];psi_re=[{i:row[2*i+1] for i in range(3) if row.get(2*i+1)} for row in L]
 pre='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
 Pi=list(map(sp,fld(pre,['PiBottom3_from56_columns'])['PiBottom3_from56_columns']));assert len(Pi)==56
 dc=[]
 for ch in range(2):dc.append({'channel':ch,'r_row56':[[j,str(dot(psi_re[ch],col))] for j,col in enumerate(Pi) if dot(psi_re[ch],col)],'ell_row56':[[j,str(dot(psi_rr[ch],col))] for j,col in enumerate(Pi) if dot(psi_rr[ch],col)]})
 assert dc==bd['DeltaC_known_point_rows']==d['DeltaC_known_point_rows']
 assert bd['do_not_use_KU_annihilation_or_B2zero'] and bd['no_actual_CU_K_X2_X3_M1_or_relation_values']
 for p in bd['source_premises'].values():assert rd(p)['all_passed']
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_C_COMPRESSION_SOURCE_20261105_v63/report.json';pure=rd(cp);assert pure['all_passed'] and pins[cp]=='ffe7b0871f93c863490645189549d5959f063a23583bc500e7176d6517ab2caf'
 chart='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json'
 c=fld(chart,['free30_V_indices','fixed_bottom30_inverse_columns','original27_X_source_columns'])
 assert c['original27_X_source_columns']==d['pivot27_source_indices']
 P=[{} for j in range(56)];ownPi=[{} for j in range(56)]
 for j,col in zip(c['free30_V_indices'],map(sp,c['fixed_bottom30_inverse_columns'])):
  P[j]={i:x for i,x in col.items() if i<27};ownPi[j]={i-27:x for i,x in col.items() if i>=27}
 assert P==list(map(sp,d['P27_columns56'])) and ownPi==Pi==list(map(sp,d['Pi3_columns56']))
 invp='reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63/actual-two29-pairings-and-recovered30-kernel-invariants.json'
 iv=fld(invp,['actual_all_three_point_kernel_direction_invariants'])['actual_all_three_point_kernel_direction_invariants']
 alpha=[Q(iv[str(j)][0]) for j in [23,29,30]];beta=[Q(iv[str(j)][1]) for j in [23,29,30]]
 assert alpha==list(map(Q,d['alpha_three'])) and beta==list(map(Q,d['beta_three']))
 # Positive denominator/numerator pairs: all SUM branches retained, never signed C.
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());assert maximum<=147;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def positive_mv(cols,col):
  v={}
  for j,x in col.items():
   for i,y in cols[j].items():v[i]=plus(v.get(i,(1,0)),mul(ex(x),ex(y)))
  return v
 G=list(map(sp,d['Gamma_two']));L0=list(map(sp,d['native_L0']));Le=list(map(sp,d['native_Le']));S0=list(map(sp,d['native_S0']));Se=list(map(sp,d['native_Se']))
 Xr=[{i:ex(x) for i,x in sp(v).items()} for v in d['pivot27_r_columns']];Xe=[];support=0
 for kind,arg in d['pivot27_ell_recipes']:
  if kind=='literal':v={i:ex(x) for i,x in sp(arg).items()};support+=len(arg)
  else:
   assert kind=='native';v=positive_mv(L0,Se[arg]);w=positive_mv(Le,S0[arg]);support+=len(set(v)|set(w))
   for i,x in w.items():v[i]=plus(v.get(i,(1,0)),x)
  Xe.append(v)
 def bdot(a,col):
  v=(1,0)
  for i,x in a.items():v=plus(v,mul(ex(x),col.get(i,(1,0))))
  return v
 gx=[[[bdot(row,v) for v in Xr],[bdot(row,v) for v in Xe]] for row in G]
 delta=[(sp(v['r_row56']),sp(v['ell_row56'])) for v in dc];bounds=[]
 for ch in range(2):
  bch=[]
  for k in range(2):
   bvec=[]
   for j in range(56):
    q=ex(delta[ch][k].get(j,Q(0)))
    if ch==0:
     for i,x in Pi[j].items():q=plus(q,mul(ex((alpha if k==0 else beta)[i]),ex(x)))
    for i,x in P[j].items():q=plus(q,mul(gx[ch][k][i],ex(x)))
    bvec.append([str(q[0]),str(q[1])])
   bch.append(bvec)
  bounds.append(bch)
 assert bounds==d['positive_Cfull_two_by2_by56'] and maximum==gr['source_positive_prior_bits']==147
 assert d['native_pivot_ell_is_SUM_L0Se_LeS0'] and d['full_coordinate_not_pureU'] and count(d)==2401
 custody='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_PREPARED_CUSTODY_20261105_v63/report.json';cu=rd(custody)
 assert cu['all_passed'] and pins[custody]=='18eef8958cb07696c207efe86d6d0561e715389989e6b4cc40a22c5901992ac8'
 for path,h in cu['literal_input_sha256'].items():assert sha(ROOT/path)==h
 worker=ROOT/cu['worker'];assert sha(worker)==cu['worker_sha256']=='6c1c18de55d39ab9fe60d6c8bb9abbd7953a55f414479952e0b7fa07f47636d1'
 helpresult=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,timeout=10);assert helpresult.returncode==0 and hashlib.sha256(helpresult.stdout).hexdigest()==cu['safe_help_stdout_sha256']
 # Recompute literal prospective worker ledger; do not construct signed Xe/C.
 bound=2*count(d)+count([G,L0,Le,S0,Se,P,Pi,list(map(sp,d['pivot27_r_columns'])),delta,alpha,beta])+2*support+108+224+332+600+1000
 assert bound==cu['explicit_prearith_live_bound']==9314<gr['future_Cfull_numeric_bound']==11644
 assert cu['derived_committed_bits']==147 and cu['derived_temporary_bits']==3*147+1==442 and cu['hard_committed_bits']==384 and cu['hard_temporary_bits']==1153
 assert cu['source_helpers']=='NONE:standardlibraryonly'
 # RAW means the complete set of literal files read/hashed by the worker,
 # including its source-gate historical hash closure, plus an entire1MB output.
 closure=set(pins)|set(gr['input_sha256'])|set(br['input_sha256'])|set(cu['literal_input_sha256'])|{cu['worker']}
 closurebytes=sum((ROOT/path).stat().st_size for path in closure)
 own_reserved=1000000 # fresh peer/code/note/custody metadata upper
 raw_upper=closurebytes+own_reserved+1000000;assert raw_upper<33554432
 owned=count([br,bd,gr,d,raw,Qcols,L,lf,N,psi_rr,psi_re,Pi,dc,pure,c,P,ownPi,iv,alpha,beta,G,L0,Le,S0,Se,Xr,Xe,gx,delta,bounds,cu,pins])+20000;assert owned<150000
 proof={'all6_bottom_standard_NF_and_LQ_identity':True,'DeltaC_reconstructed_from_same_Pi3':True,'frozen_ffe7_scope_clarification':'ffe7 compresses f_U only. It never supplies FULL b1 annihilation by itself. Full class Q*f_U+NFbottom(B2) has upper LQcoordinates f_U+Psi(B2), and B2 need not vanish on kerLOWER. Cfull=C_U+DeltaC and Kfull only may be used for full393 source annihilation. Frozen ffe7 bytes preserved.','full393_source_annihilation_proof':'For X=0,Y=G3v first27 coordinates vanish. Fixed25 correction leaves393 source directions. Complete zero quadratic classes imply the full upper coordinate equals -Kfull B99 C G3linear v=0. This is a statement about the COMPLETE coordinate map, and is valid in the coherent source frame; no use of KU-only or B2=0.','conditional_choice_independence':'Only if separately proven kerG2linear918=imG3linear242 may the full covector be descended via an output adjoint. No actual K values or rightsolve are supplied here.','all_pruned_fields_bound_to_original_cache_chart_invariants':True,'individual_native_pivotell_SUM_not_union':True,'all224_positive_envelopes_independently_rebuilt':True,'positive_prior_bits':maximum,'temporary_bound':442,'prospective_worker_live_upper':bound,'worker_cap':11644,'literal_RAW_input_closure_bytes':closurebytes,'fresh_peer_metadata_RAW_reserve':own_reserved,'full_output_RAW_reserve':1000000,'whole_RAW_upper':raw_upper,'no_signed_C_K_or_upper_values_evaluated':True,'no_new_physical_gauge_higher_coefficients':True,'source_peer_owned_bound':owned}
 ep='independent-full-Psi-Cfull-source-code-height-and-whole-bound.json';put(ep,proof)
 report={'id':STEM,'all_passed':True,'SOURCE_prior_code_custody_only_no_signed_C':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'whole_RAW_upper':raw_upper,'seconds':time.monotonic()-t};put('report.json',report);assert sum(x.stat().st_size for x in o.iterdir())<own_reserved;print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned,'RAWupper':raw_upper}))
if __name__=='__main__':main()
