#!/usr/bin/env python3
"""D coefficient gate only: source-frame, ordinary-pair, positive/custody peer."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,itertools,json,math,resource,signal,subprocess,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_D_INPUT_GATE_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_D_INPUT_GATE_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='552d60c970aefd792af4ab4ae2fe622441712e1618be43fbb29d908b8ec9e253'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 d=rd(b+'pruned-D-coefficient-functional-and-ten-pair-source-input.json');assert count(d)==2851
 C=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json',['Cfull_two_by2_by56'])['Cfull_two_by2_by56'];assert d['Cfull_two_by2_by56']==C;C=[[[Q(x) for x in row] for row in ch] for ch in C]
 R=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63/actual-two-rho1144-covectors-and-full918-identities.json',['rho_two_by1144_sparse_rows'])['rho_two_by1144_sparse_rows'];assert d['rho_two_by1144_sparse_rows']==R;R=list(map(sp,R))
 m=fld('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['Mnew0_columns'])['Mnew0_columns'];assert m==d['Mnew0_columns56']
 frame='reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json';assert fld(frame,['actual_Mnew0_columns'])['actual_Mnew0_columns']==m;M=list(map(sp,m));assert len(M)==56 and all(0<=i<26 for c in M for i in c)
 for p,h in [('reports/AMBIENT_B18_ELL2_T2_UPPER617_RHO_20261105_v63_REVIEW/report.json','0e5a9cb0cef8573a4be81c5f2dbb7e1e6771f3a8814b833495fcd76729bd33b8'),('reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63_REVIEW/report.json','76e1aa8842efe0d6d5af78572641c23437ebc902d612c8424876d0058624a48e'),('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json','6428e2178c41cfd7cfcf271402c4d13898879414a0b5129fd409e3c97300d460')]:assert rd(p)['all_passed'] and pins[p]==h
 pairs=[list(v) for v in itertools.combinations_with_replacement(range(9),2) if v!=(8,8)];assert pairs==d['quadratic_pair_order44'];ix={tuple(v):i for i,v in enumerate(pairs)}
 active=sorted({(7,7),(7,8)}|{tuple(pairs[i%44]) for row in R for i in row});assert active==list(map(tuple,d['active_pairs10']))==[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)]
 maximum=0
 def ex(x):
  nonlocal maximum
  x=Q(x);v=(x.denominator,abs(x.numerator)) if x else (1,0);maximum=max(maximum,v[0].bit_length(),v[1].bit_length());assert maximum<=92;return v
 def add(a,b0):
  nonlocal maximum
  if not a[1]:return b0
  if not b0[1]:return a
  D=math.lcm(a[0],b0[0]);v=(D,a[1]*(D//a[0])+b0[1]*(D//b0[0]));maximum=max(maximum,D.bit_length(),v[1].bit_length());assert maximum<=92;return v
 def mul(a,b0):
  nonlocal maximum
  if not a[1] or not b0[1]:return (1,0)
  v=(a[0]*b0[0],a[1]*b0[1]);maximum=max(maximum,v[0].bit_length(),v[1].bit_length());assert maximum<=92;return v
 env=[]
 for t0 in range(2):
  ch=[]
  for pair in active:
   row=[]
   for j,col in enumerate(M):
    q=ex(C[t0][1][j]) if pair==(7,7) else ex(C[t0][0][j]) if pair==(7,8) else (1,0)
    for target,x in col.items():q=add(q,mul(ex(R[t0].get(44*target+ix[pair],Q(0))),ex(x)))
    row.append([str(q[0]),str(q[1])])
   ch.append(row)
  env.append(ch)
 assert env==d['positive_D_two_by10_by56'] and maximum==r['source_positive_prior_bits']==92 and 3*92+1==r['temporary_prior_bits']==277
 assert d['no_OLDpivot_honest2_cache_assumption'] and d['native_individualXp2_full_L2S0_retained']
 cb='reports/AMBIENT_B18_ELL2_T2_UPPER617_D_PREPARED_CUSTODY_20261105_v63/';cu=rd(cb+'report.json');assert cu['all_passed'] and pins[cb+'report.json']=='284d0fee53d8b43857640785bda7bd9c009d6796681b20fdc03a2562ca5132e8'
 for p,h in cu['literal_input_sha256'].items():assert sha(ROOT/p)==h
 frontier='reports/AMBIENT_B18_ELL2_T2_UPPER617_PIVOT4_SECOND_SOURCE_20261105_v63/report.json';fr=rd(frontier);assert fr['all_passed'] and pins[frontier]=='aa578fe00bcfab557c8c727a7a387fe24ba120a43820f94c205e2b5f884fc0d9' and fr['unknown_boundaries']==4
 worker=ROOT/cu['worker'];assert sha(worker)==cu['worker_sha256']=='435a9f4c184cb309acdc19792587e4500f5c536be49d012cdd224d0516b49e9e'
 helper=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,timeout=10);assert helper.returncode==0 and hashlib.sha256(helper.stdout).hexdigest()==cu['safe_help_stdout_sha256']
 code=worker.read_text();assert 'q=C[t][1][j] if pair==(7,7) else C[t][0][j] if pair==(7,8)' in code and 'term=mul(rho[t].get(44*target+mi,F(0)),x);q=sub(q,term)' in code and "if b>92:stop('source92_before_commit',x)" in code
 assert all(s in code for s in ["'pending_operation':dict(pending)","'accumulator_before_operation':str(acc)","'completed_D_prefix':D","'input_sha256':pins","'no_retry_or_cap_raise':True"])
 ledger=rd(cb+'allstage-live-RAW-prior-ledger.json');assert pins[cb+'allstage-live-RAW-prior-ledger.json']==cu['ownership_ledger_sha256'];parsed=count([C,R,M,ix,active]);assert parsed==652
 assert 2*2851+parsed+2240+1000==ledger['success_live_upper']==cu['success_live_upper']==9594<14334
 assert ledger['guard_state_prefix1120_plus40meta']==1160 and 2*2851+parsed+3*1160+1000==ledger['guard_live_upper']==cu['guard_live_upper']==10834<14334
 assert cu['source_positive_bits']==92 and cu['temporary_bits']==277 and cu['hard_guards']==[384,1153] and cu['helpers']=='NONE:standardlibraryonly'
 rawset=set(cu['literal_input_sha256'])|{cu['worker'],cb+'report.json',cb+'allstage-live-RAW-prior-ledger.json'};rawbytes=sum((ROOT/p).stat().st_size for p in rawset);assert rawbytes+2100000<33554432
 proof={'all_Cfull_Psi_rho_and_M0_fields_literal_original':True,'M0_matches_actual_active_transport_and_is_unchanged_by_free3_gauge':'MtoM+G2A0changesnoM0becauseminimalG2origin0;no canonicalM1zeroassumed.','ordinary_pair_formula':'At r²ell, C_r multipliesX2_rell and C_ell multipliesX2_rr. Dpair=Ccontribution-rho_pair*M0. No factor2;allother44pairsoutsideactive10aresourcezero.','all1120_positive_envelopes_reconstructed92_temp277':True,'outer_source_frontier':'OLD604/613/622/631outsideK99old80;thisgateneverqueriesXc2orXp2. Fourg2frontiersandindividualnativeL2S0remainseparateuntilwhole617sensitivity.','success9594_guard10834_cap14334':True,'whole_RAW_independent_upper':rawbytes+2100000,'custody_whole_RAW2137396':cu['whole_RAW_upper'],'worker_literal_guard_help_and_custody':True,'no_signed_D_Xc2_RHS_honestgradient_firstcross_b1_upper_or_higher_values':True}
 ep='independent-D-source-frame-ordinarypairs-positive-and-custody.json';put(ep,proof)
 owned=count([r,d,C,R,m,M,pairs,ix,active,env,cu,fr,ledger,proof,pins])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'SOURCE_coefficient_gate_only_no_signed_D':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()
