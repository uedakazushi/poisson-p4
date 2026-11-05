#!/usr/bin/env python3
"""SOURCE only: fixed row58 col56 Q/R/E integer pseudo-division budget."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,ast,gc,hashlib,json,os,re,resource,signal,subprocess,sys,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode,encode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_STRIP_PREFLIGHT_20261105_v64'
WORKER='code/s61_sparse_generic5_constrained83_row58_existing_QRE_strip_20261105_v64.py'
TARGET='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_EXISTING_QRE_STRIP_20261105_v64'
VALUE='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64'
VALUE_REVIEW='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_QQ_VALUE_REVIEW_20261105_v64'
CODEC_GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STREAM_CODEC_SOURCE_GATE_20261105_v64'
EXPECTED={VALUE:'20524cd69e98bc7f5334a07544035fc04ea4fc1486ac64ce885670c0ce70b173',VALUE_REVIEW:'788daf79ec359ed9819656bff6c7a624bfdc4ce96363475a574c519f7c1c7e06',CODEC_GATE:'c9921d8f8db77f5a4cf2739b26d88a4728cb871c4e58b4ed9aa748a6d09e2d6f'}
def dump(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for z in iter(lambda:f.read(65536),b''):h.update(z)
 return h.hexdigest()
def scalar_count(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(scalar_count,v))
 if isinstance(v,dict):return sum(map(scalar_count,v.values()))
 return 0
def main(a):
 assert a.freeze
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in()).throw(TimeoutError('SOURCE60s')));signal.alarm(60)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={};sizes={}
 def add(path):
  p=ROOT/path;inputs[path]=sha(p);sizes[path]=p.stat().st_size;return p
 try:
  source_meta=0
  for stem,h in EXPECTED.items():
   path=f'reports/{stem}/report.json';p=add(path);assert inputs[path]==h;v=json.loads(p.read_text());assert v['all_passed'];source_meta=max(source_meta,scalar_count(v));del v
  api_path=f'reports/{VALUE}/actual15-col56-and-once-shared-row35.json'
  api=json.loads(add(api_path).read_text());rec=next(z for z in api['actual15_ninth_col56_numerators']if z['row']==58)
  assert rec['terms']==951 and rec['content']=='1/353894400000'
  path=f'reports/{VALUE}/'+rec['file'];ip=add(path);assert inputs[path]==rec['sha256']=='5774ebc24c84576c4d7b21f5e245644a4e1a6ad224447d9ab82ec0d1bf2a035b'
  binary=ip.read_bytes();assert binary[:9]==b'Q5CELL1\n\x01';p=decode(binary[9:]);assert b'Q5CELL1\n\x01'+encode(p)==binary
  q=Q(p['content']);assert str(q)=='1/353894400000'and len(p['terms'])==951
  valuation=[min(p['shift'][j]+t[j+1]for t in p['terms'])for j in range(5)];assert valuation==[-1,0,0,0,-1]
  P={};g=0;N0=0;groups={}
  for t in p['terms']:
   c=int(t[0]);e=tuple(p['shift'][j]+t[j+1]-valuation[j]for j in range(5));assert c and e not in P and min(e)>=0;P[e]=c;g=gcd(g,abs(c));N0+=abs(c);f=(e[0],e[1],e[4]);groups[f]=max(groups.get(f,0),e[2]+e[3])
  assert g==1 and len(P)==951 and len(groups)==17
  U=sum((d+1)*(d+2)//2 for d in groups.values());D=max(groups.values());assert(U,D)==(1739,16)
  box=[max(e[j]for e in P)for j in range(5)];assert box==[6,2,8,13,6]
  source_numeric_upper=25000+8*len(p['terms'])+7*len(P)+5*U+source_meta+2048
  assert source_numeric_upper<150000
  # No polynomial division is performed anywhere in this SOURCE program.
  NP=N0*(D+1)**D*2**(D*(D+1)//2)
  Mdiv=32**D*NP*(1+50*(D+1)*2**D)
  Mrebuild=NP*82**D;Mglobal=max(Mdiv,Mrebuild);Midentity=83*Mdiv;h=Mglobal.bit_length()
  trials=D+3;cells=2*trials+2;meta=25000;runtime_reserve=5000;root_reserve=2000
  meta_proof=5*U+8*951+64*trials+runtime_reserve+root_reserve
  ordinary=meta+7*(3*U)+2048;codec=meta+7*(3*U)+3*U+2048
  assert meta_proof==24519<=meta and ordinary==63567 and codec==68784<150000 and 2*h<2048
  # Every exponent is <=16, including each fixed H/L/t fibre; five one-byte
  # unsigned exponent integers. sint length-tag and all headers are overbounded.
  assert max(box)<=16 and D<128
  one_cell=2048+U*((h+7)//8+5+16)
  for f in [WORKER,'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py',str(Path(__file__).relative_to(ROOT)),f'research/notes/{ID}.md',f'research/notes/{TARGET}.md']:add(f)
  wp=ROOT/WORKER;tree=ast.parse(wp.read_text());assert not(ROOT/'reports'/TARGET).exists()
  assert sum(isinstance(n,ast.ImportFrom)and n.module=='s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64'for n in ast.walk(tree))==1
  help_result=subprocess.run([sys.executable,str(wp),'--help'],capture_output=True,text=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=10)
  assert help_result.returncode==0 and '--execute' in help_result.stdout and not(ROOT/'reports'/TARGET).exists()
  actual_read=sum(sizes.values());future_headers=3*65536;trace_json_reserve=2*1024*1024
  raw=actual_read+future_headers+cells*one_cell+trace_json_reserve
  assert raw<32*1024*1024
  proof={
   'id':ID,'source_only':True,'no_polynomial_trial_division_or_new_matrix_coefficient':True,
   'worker':WORKER,'worker_sha256':inputs[WORKER],
   'original_binary_path':path,'original_binary_sha256':inputs[path],
   'original_QQ_content':str(q),'actual_Laurent_monomial_valuation5':valuation,
   'source_terms':951,'primitive_integer_coefficient_gcd':1,
   'fixed_H_L_t_fibres':[[list(k),d]for k,d in sorted(groups.items())],
   'support_universe_terms':U,'normalized_exponent_box5':box,'max_s_q_total_degree':D,
   'input_integer_l1_norm':str(N0),'all_normalized_quotient_l1_bound':str(NP),
   'all_pseudo_maps_integer_l1_bound':str(Mdiv),'all_rebuild_maps_integer_l1_bound':str(Mrebuild),
   'all_identity_partial_integer_bound':str(Midentity),
   'integer_height_bits':{'initial_l1':N0.bit_length(),'normalized_quotient':NP.bit_length(),'pseudo_map_and_partial':Mdiv.bit_length(),'reverse_product':Mrebuild.bit_length(),'identity_partial':Midentity.bit_length(),'codec_qq_or_gcd':h,'codec_naive_coefficient_times_G':2*h,'codec_hard_integer_bit_guard':2048},
   'height_proof':[
    'For a*x+B, norm(B)/a<=2; a true quotient has l1 norm <= norm(P)*(d+1)*2^d.',
    'Each success lowers every nonzero fixed-fibre total(s,q) degree by one; <=16 successes.',
    'NP=N0*17^16*2^(16*17/2) bounds all primitive exact quotients.',
    'Scale each trial by a^d (a<=32,d<=16). Leading-variable order gives integer divisions by a; processed degree k residual coefficients are multiples of a^k.',
    'The leading term is removed before tail updates, so every residual partial has l1 <=32^16*NP*(1+50*17*2^16). Every quotient and multiplication register fits that or83times that bound.',
    'If remainder=0, Gauss lemma for primitive integer F and P gives V/a^d primitive integer; no parameter inverse or rational-content inflation.',
    'Reverse products have l1 <=NP*82^16. Constant terms of Q,R,E prevent any coordinate monomial valuation change.',
    'For codec q=1, G and primitive c each <=Mglobal; naive c*G has <=2*h=928 bits, below2048, including signed integer normalization and valuehash registers.'
   ],
   'fixed_factor_order':['Q=5*q0+1','R=q0+1','E=32*s0+15*q0+35'],
   'not_divided':['F9','any_unknown_or_normal_coordinate'],'no_general_GCD_or_new_row_search':True,
   'max_successes':D,'max_trials':trials,'max_integer_trace_cells':cells,
   'metadata_upper':meta,'metadata_derived_upper':meta_proof,
   'metadata_budget_terms':{'universe_tuple_scalars':5*U,'decoded_original_numeric_scalars':8*951,'max_trial_index':64*trials,'runtime_pr_proof_misc':runtime_reserve,'root_review_max_scalars':root_reserve},
   'whole_owned_upper':codec,'ordinary_three_map_upper':ordinary,
   'ownership_proof':[
    'At most three P/residual/quotient integer maps; seven numeric cells charged per nonzero term, five key coordinates plus integer and one slot.',
    'Only one existing-key alias-sort array and codec state; additional3*n+2048, never another copied coefficient map.',
    'Reverse reconstruction retains only B and new T maps; local P alias is explicitly deleted.',
    'The original decoded951 term list is charged inside metadata25000, and released before trials.',
    'Root review must be <=65536 RAW bytes and <=2000 numeric metadata scalars.'
   ],
   'one_integer_cell_RAW_upper':one_cell,'actual_READ_input_RAW_bytes':actual_read,
   'actual_READ_files_bytes':sizes,'future_source_proof_report_rootreview_RAW_reserve':future_headers,
   'all_trace_report_receipts_RAW_reserve':trace_json_reserve,'fresh_RAW_upper':raw,
   'RAW_proof':'Distinct actual READ closure + <=40 uncompressed canonical integer cells + all trial/result/failure/report JSON and source/root metadata reserves. No gzip, deletion, or omitted versions.',
   'source_numeric_upper':source_numeric_upper,'whole_limits':{'seconds':60,'AS_bytes':2000000000,'owned_numeric':150000,'fresh_RAW_bytes':32*1024*1024},
   'matrix_scope':{'source_full_rows':168,'source_full_columns':15,'other_columns_unchanged':14,'extra_top_RHS_rows':87,'constraint_zero_rows':66,'P290_applied_once':True,'F9_complement_and_all_old_complements_open':True,'no_tenth_matrix_update_or_rank_claim':True},
   'prepared_worker_help_passed':True,'prepared_worker_help_sha256':hashlib.sha256(help_result.stdout.encode()).hexdigest()
  }
  bp=out/'integer-pseudo-division-budget-and-proof.json';bp.write_bytes(dump(proof));assert bp.stat().st_size<=65536
  report={'id':ID,'all_passed':True,'source_only':True,'prepared_fixed_QRE_worker_safe':True,'no_trial_division_executed':True,'no_new_matrix_coefficient':True,'source_terms':951,'fixed_fibres':17,'support_universe_upper':U,'max_trial_divisions':trials,'future_owned_upper':codec,'future_RAW_upper':raw,'all_codec_temporary_bits_upper':2*h,'source_numeric_upper':source_numeric_upper,'source_metadata_report_max':source_meta,'worker_sha256':inputs[WORKER],'input_sha256':inputs,'evidence_sha256':{bp.name:sha(bp)},'seconds':time.monotonic()-start}
  rp=out/'report.json';rp.write_bytes(dump(report));assert rp.stat().st_size<=65536
  proof_scalars=scalar_count(proof);report_scalars=scalar_count(report)
  assert proof_scalars+report_scalars+3000<=runtime_reserve
  assert sum(p.stat().st_size for p in out.iterdir())+actual_read<=32*1024*1024
  signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp));print('proof_sha256',sha(bp))
 except Exception as exc:
  (out/'failure.json').write_bytes(dump({'all_passed':False,'source_only':True,'error':repr(exc),'no_trial_division_executed':True,'seconds':time.monotonic()-start}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
