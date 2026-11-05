#!/usr/bin/env python3
"""Fresh full918 alpha clearing SOURCE/control/positive/custody peer only."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,math,resource,signal,subprocess,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_INPUT_GATE_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_INPUT_GATE_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='6e68e307280750ea713beebb73fa1c3979cc7d122b6389c899e11484c9c804c9'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 d=rd(b+'actual-four-row-to-full918-alpha-clearing-pruned-input.json');assert count(d)==2911
 rows=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63/actual-four-specific-B99-point-dual-rows-STOP-only384.json',['actual_four_rows_channel_pair'])['actual_four_rows_channel_pair'];assert d['actual_four_point_rows']==rows
 rp='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63_REVIEW/report.json';assert rd(rp)['all_passed'] and pins[rp]=='9468089e150e31cd367fb115f08654cd807ff4fdaca6635fe6a64623250c0e3f'
 s={(v['channel'],v['pair_rr_rell']):sp(v['point_row102']) for v in rows};pre=[]
 for t0 in range(2):pre.append({**{9*i+7:x for i,x in s[t0,1].items()},**{9*i+8:x for i,x in s[t0,0].items()}})
 assert pre==list(map(sp,d['pre_two_full918_rows']))
 h='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';z=fld(h,['PiY3_from102_columns','actual_full_G3linear102_by418_columns'])
 assert z['PiY3_from102_columns']==d['PiY3_from102_columns'] and z['actual_full_G3linear102_by418_columns']==d['all418_G3_full918_columns']
 jsrc=fld('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json',['selected_columns'])['selected_columns'];assert len(jsrc)==25 and jsrc==d['selected25_G3_source_indices'] and [z['actual_full_G3linear102_by418_columns'][j] for j in jsrc]==d['J25_full918_columns']
 fpath='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_FIXED_ROWS_20261105_v63/actual-fixed-specific-dual-rows-and-F25-with-next-clearing-prior.json';assert fld(fpath,['actual_F25_columns27'])['actual_F25_columns27']==d['F25_columns27']
 # The existing fixedF25 source proof is reused, never a new row inversion.
 fpr='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63_REVIEW/report.json';assert rd(fpr)['all_passed']
 Pi=list(map(sp,d['PiY3_from102_columns']));J=list(map(sp,d['J25_full918_columns']));F=list(map(sp,d['F25_columns27']));G=list(map(sp,d['all418_G3_full918_columns']));assert [len(Pi),len(J),len(F),len(G)]==[102,25,27,418]
 assert all(0<=i<3 for c in Pi for i in c) and all(0<=i<25 for c in F for i in c) and all(0<=i<918 for c in J+G for i in c)
 maxbits=0
 def mark(a):
  nonlocal maxbits
  maxbits=max(maxbits,a[0].bit_length(),a[1].bit_length());assert maxbits<=137;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def plus(a,b0):
  if not a[1]:return b0
  if not b0[1]:return a
  D=math.lcm(a[0],b0[0]);return mark((D,a[1]*(D//a[0])+b0[1]*(D//b0[0])))
 def mul(a,b0):return mark((a[0]*b0[0],a[1]*b0[1])) if a[1] and b0[1] else (1,0)
 Z=(1,0);lj=[];bounds=[]
 for t0 in range(2):
  qj=[]
  for c in J:
   q=Z
   for i,x in c.items():q=plus(q,mul(ex(pre[t0].get(i,Q(0))),ex(x)))
   qj.append(q)
  bvec={}
  for i in range(102):
   for v in range(9):
    q=ex(pre[t0].get(9*i+v,Q(0)))
    for target,pi in Pi[i].items():
     for j,x in F[9*target+v].items():q=plus(q,mul(qj[j],mul(ex(pi),ex(x))))
    if q[1]:bvec[9*i+v]=q
  lj.append(qj);bounds.append(bvec)
 clearing=maxbits;assert clearing==110==r['clearing_positive_prior_bits']
 for row in bounds:
  for col in G:
   q=Z
   for i,x in col.items():q=plus(q,mul(row.get(i,Z),ex(x)))
 assert maxbits==137==r['all836_G3_check_positive_prior_bits'] and 3*137+1==r['temporary_prior_bits']==412
 assert [[[str(D),str(U)] for D,U in row] for row in lj]==d['positive25_LJ_envelopes']
 assert [[[i,str(D),str(U)] for i,(D,U) in sorted(row.items())] for row in bounds]==d['positive_two_full918_alpha_envelopes']
 cb='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_PREPARED_CUSTODY_20261105_v63/';cu=rd(cb+'report.json');assert cu['all_passed'] and pins[cb+'report.json']=='043cdd12e0aa824fe9d1b81b5098fd23537a4fa4c5d6b0131509691d33374fb5'
 for p,h0 in cu['literal_input_sha256'].items():assert sha(ROOT/p)==h0
 worker=ROOT/cu['worker'];assert sha(worker)==cu['worker_sha256']=='68c7271ec7a5fe739c69ca7e6b3b0662edfb6d6a8a18aac973b0d197ccc7adfc'
 helper=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,timeout=10);assert helper.returncode==0 and hashlib.sha256(helper.stdout).hexdigest()==cu['safe_help_stdout_sha256']
 code=worker.read_text();assert "if b>137:stop('source137_before_commit',x)" in code and 'pending.update' in code and 'return ck(a*b)' in code and 'return ck(a+b)' in code
 assert 'term=mul(pi,x);term=mul(lj[b],term);q=add(q,-term)' in code and "if q:stop('actual_Kfull_G3_annihilation_identity_failed',q)" in code
 assert all(v in code for v in ["'pending_operation':dict(pending)","'partial_accumulator_before_operation':str(accumulator)","'completed_alpha_prefix':[rec(v) for v in alpha]","'completed_LJ25'","'completed_G3_literal_dots'","'input_sha256':pins","'no_retry_or_cap_raise':True"])
 parsed=count([Pi,J,F,G,pre,[{i:(int(a),int(b0)) for i,a,b0 in row} for row in d['positive_two_full918_alpha_envelopes']]])
 live=2*count(d)+parsed+2*3672+2*836+50+2000;assert live==cu['success_live_upper']==19554<r['future_value_numeric_bound']==24230
 ledger=rd(cb+'explicit-clearing-live-RAW-and-guard-state-ledger.json') if (ROOT/cb/'explicit-clearing-live-RAW-and-guard-state-ledger.json').exists() else None
 if ledger is None:
  options=[p for p in (ROOT/cb).glob('*.json') if p.name!='report.json'];assert len(options)==1;lp=options[0];pins[str(lp.relative_to(ROOT))]=sha(lp);ledger=json.loads(lp.read_text())
 assert any(h0==cu['ownership_ledger_sha256'] for h0 in pins.values())
 assert sum(map(len,bounds))==ledger['positive_alpha_support_bound']==86
 assert ledger['full_guard_state_numeric_upper']==2*86+50+836+30==1088
 assert ledger['guard_state_three_copies']==3*1088==3264
 assert ledger['guard_path_live_upper']==2*2911+parsed+3264+1000==12752==cu['guard_live_upper']<24230
 assert ledger['success_prearith_live_upper']==live and ledger['source_committed_prior']==137 and ledger['temporary_integer_prior']==412
 assert max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for A in [Pi,J,F,G,pre] for c in A for x in c.values())<=137
 rawset=set(cu['literal_input_sha256'])|{cu['worker'],cb+'report.json'}|{p for p in pins if p.startswith(cb)};rawbytes=sum((ROOT/p).stat().st_size for p in rawset);assert rawbytes+1000000+100000+2000000<33554432
 assert cu['helpers']=='NONE:standardlibraryonly' and cu['source_committed_bits']==137 and cu['temporary_bits']==412 and cu['hard_guards']==[384,1153] and cu['whole_numeric_cap']==24230
 proof={'Krell_to_r_Krr_to_ell_ordinaryslot918_source_bound':True,'actual4rows_original9468_and_fixed25_PiY_J_F_all_literal':True,'all1836_clearing_positive_envelopes_and50_LJ_rebuilt':True,'all836_G3_positive_check_bounds':True,'clearing110_G3_137_temp412':True,'worker_multiply_and_add_guard_before_use_and_all_G3_checks_required':True,'full_state_cursor_operands_accumulator_partialalpha_LJ_Gdots_hashes_preserved':True,'success_live19554_cap24230':True,'custody_guard_live12752':cu['guard_live_upper'],'whole_literal_RAW_plus_fresh_metadata_and_output':rawbytes+3100000,'custody_whole_RAW3146494':cu['whole_RAW_upper'],'source_scope':'Onlyfull918alphaand836knownG3dots. No transposeR/rho/b1/sourceZ2/upper/HF. No canonicalcomplementM1coefficientisassignedzero.','no_signed_alpha_G3dots_or_new_response_evaluated':True}
 ep='independent-alpha918-clearing-source-prior-custody-and-whole-bound.json';put(ep,proof)
 owned=count([r,d,rows,z,jsrc,s,pre,Pi,J,F,G,lj,bounds,cu,ledger,proof,pins])+12000;assert owned<60000
 report={'id':STEM,'all_passed':True,'SOURCE_code_custody_positive_only_no_signedalpha':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()
