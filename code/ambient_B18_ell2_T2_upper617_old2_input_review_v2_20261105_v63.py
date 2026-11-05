#!/usr/bin/env python3
"""Independent source/positive/code/custody audit; no signed OLD2 fold."""
import argparse,ast,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_UPPER617_OLD2_INPUT_GATE_20261105_v63'
STEM=A+'_REVIEW_V2'
G='reports/'+A+'/'
C='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_PREPARED_CUSTODY_V2_20261105_v63/'
W='code/ambient_B18_ell2_T2_upper617_old2_formal_20261105_v63.py'
S='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def fields(p,names):
 text=p.read_text();dec=json.JSONDecoder();i=text.index('{')+1;z={}
 while True:
  while text[i].isspace() or text[i]==',':i+=1
  if text[i]=='}':break
  k,i=dec.raw_decode(text,i)
  while text[i].isspace() or text[i]==':':i+=1
  v,i=dec.raw_decode(text,i)
  if k in names:z[k]=v
  del v
 assert set(z)==set(names),(p,names,set(z))
 return z
def main(out):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;maximum=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=q.stat().st_size
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,x):(out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 r=read(G+'report.json');assert pins[G+'report.json']=='e193894d22d0da56a2420c3441c7aecdfd523405c84149316104713cde48c0ce'
 d=read(G+'617-combined-OLD2-independent-symbol-positive-input.json');assert pins[G+'617-combined-OLD2-independent-symbol-positive-input.json']==r['evidence_sha256']['617-combined-OLD2-independent-symbol-positive-input.json']
 c=read(C+'report.json');assert pins[C+'report.json']=='1db52feced9adb8ae2f4853efcb080b08e96e64295ef9612981fad0365442fac'
 ledger=read(C+'explicit-whole-ownership-and-RAW-ledger.json');assert pins[C+'explicit-whole-ownership-and-RAW-ledger.json']==c['ownership_ledger_sha256']
 assert r['all_passed'] and c['all_passed'] and sha(ROOT/W)==c['worker_sha256']=='62c9ea29d111ae2e9d4eb711351443c12638fb80b29c160b33466bd1cc7e2aa9'
 for p,h in {**r['input_sha256'],**c['literal_input_sha256']}.items():assert sha(ROOT/p)==h
 sr=read(S+'report.json');s=read(S+'combined-OLD2-D-and-four-pivot-source-factor-inventory.json')
 assert pins[S+'report.json']=='addc65c2eb868114a773f8dbed9e7c513cc5f21b7e9e46d7604f23498d9d8e16'
 assert pins[S+'combined-OLD2-D-and-four-pivot-source-factor-inventory.json']==sr['evidence_sha256']['combined-OLD2-D-and-four-pivot-source-factor-inventory.json']
 assert count([r,d,c,ledger,s,sr,pins])+50000<150000
 dump('phase.json',{'phase':'literal_source_inputs_custody_and_preoperation_owner_PASS','no_signed_gradient_or_response':True})
 # Direct literal source and current617 subtractions; no corrected coefficients
 # are evaluated with their signs in this peer.
 for x in ['D_two_by10_by56','pair_order10','original619_c0_source_combinations','kernel617_source_recipes']:
  assert d[x]==s[x],x
 assert d['literal_outer_four_pivot_source']==s['outer_four_pivot_source']
 assert d['raw619_old4_plane_a1_nonempty28coefficients']==s['raw619_old4_plane_a1_nonempty28coefficients']
 assert sum(len(v) for j,v in d['raw619_old4_plane_a1_nonempty28coefficients'])==28
 assert d['pair_order10']==[[0,7],[0,8],[1,7],[1,8],[4,7],[4,8],[5,7],[5,8],[7,7],[7,8]]
 # Independently check a1 slots directly against original actual known firsts.
 firstp='reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json'
 first=read(firstp,['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations']
 a=[]
 for j,v in enumerate(first):
  q=[[i,x] for i,x in v['actual_a1_27_times9'] if i//9<4 and i%9 in [7,8] and F(x)]
  if q:a.append([j,q])
 assert a==d['raw619_old4_plane_a1_nonempty28coefficients'];del first,a
 # Every leaf origin is read literally from already accepted caches. No second
 # value is inferred for a known-zero point, including all36 only79 functions.
 p79='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json'
 v79=read(p79,['first640'])['first640']['all79_atom_values_full10'];point79={i:dict(v).get(0,'0') for i,v in v79};del v79
 p164='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json'
 v164=read(p164,['all_requested164_root_values'])['all_requested164_root_values'];point164={i:dict(v).get(0,'0') for i,v in v164};del v164
 expected=[]
 for i,k in s['typed_actual_C0_function_classification']:
  xs=[F(v[i]) for v in [point79,point164] if i in v];assert xs and all(x==0 for x in xs)
  assert k==('source164' if i in point164 else 'C0_79_only');expected.append([i,k,'0'])
 assert sum(k=='C0_79_only' for i,k,x in expected)==36
 gp='reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63_REVIEW/'
 gr=read(gp+'report.json');proof=read(gp+'independent-two57-fifth-input-and-source-order-audit.json')
 assert pins[gp+'report.json']=='917f2492740394b7ebbdeca99886cfa4cbbe18018c2171eaa5788af0c463fdd5'
 for rec in proof['both_families_proof_records']:
  assert rec['global_bbar_order_at_least4'] and rec['cofactor_origin_zero'] and rec['all_signed_original_g_and_bbar_roots_match_frozen_quotients']
 expected.extend([i,'outer_C0_regular_g_original_source917f','0'] for i in range(566400,566404))
 assert expected==d['literal_function_origins71'] and len(expected)==71
 # All typed origins vanish. Thus the dependence of a second coefficient on
 # independent leaf seconds is its point differential: multiproducts of arity
 # >=2 contribute first-cross terms, which remain outside this formal job.
 Z=(1,0)
 def guard(n):
  nonlocal maximum
  maximum=max(maximum,abs(n).bit_length());assert abs(n).bit_length()<=120
  return n
 def ex(x):
  q=F(x);return (guard(q.denominator),guard(abs(q.numerator))) if q else Z
 def add(x,y):
  if not x[1]:return y
  if not y[1]:return x
  den=guard(math.lcm(x[0],y[0]));return den,guard(guard(x[1]*(den//x[0]))+guard(y[1]*(den//y[0])))
 def mul(x,y):return (guard(x[0]*y[0]),guard(x[1]*y[1])) if x[1] and y[1] else Z
 def axpy(dst,src,scale):
  if not scale[1]:return
  for key,x in src.items():dst[key]=add(dst.get(key,Z),mul(x,scale))
 grad={};linear=[];zeros=[]
 for i,n in s['typed_program']:
  if n[0]=='actual_C0_scalar_function':assert n[1] in {x[0] for x in expected};z={n[1]:ex(1)}
  elif n[0]=='sum':
   z={}
   for j,a,b in n[1]:assert j<i;axpy(z,grad[j],ex(F(a,b)))
  else:
   assert n[0]=='product' and n[3] and all(j<i for j in n[3]);z={}
   if len(n[3])==1:axpy(z,grad[n[3][0]],ex(F(n[1],n[2])))
  grad[i]=z
  if z:
   nn=['sum',[[j,a,b] for j,a,b in n[1] if grad[j] and a]] if n[0]=='sum' else n
   linear.append([i,nn])
  else:zeros.append(i)
 assert len(grad)==757 and len(linear)==104 and sum(map(len,grad.values()))==112
 assert linear==d['linear_typed_program104'] and zeros==d['all_zero_gradient_typed_node_ids653']
 ownnodes=[[i,[[f,str(a),str(b)] for f,(a,b) in sorted(z.items())]] for i,z in grad.items() if z]
 assert ownnodes==d['positive_linear_node_gradients']
 cols=[[j,[[row,node] for row,node in v if grad[node]]] for j,v in s['D_old_root_columns']]
 cols=[v for v in cols if v[1]];assert cols==d['linear_D_old_root_columns45'] and len(cols)==45 and sum(len(v) for j,v in cols)==62
 # Absolute source combinations are assembled first, then D/Gamma are applied;
 # this is a different grouping from the producer's per-case root route.
 cg={i:{i:ex(1)} for i,k,v in expected};outer=d['literal_outer_four_pivot_source']
 for i,n in outer['literal_C0_program']:
  if n[0]=='g':assert i in cg
  elif n[0]=='n':cg[i]=dict(cg[n[1]])
  else:assert n[0]=='a';cg[i]=dict(cg[n[1]]);axpy(cg[i],cg[n[2]],ex(1))
 c0=[{j:ex(x) for j,x in v if j<640} for v in d['original619_c0_source_combinations']]
 a0={j:{i:ex(x) for i,x in v} for j,v in d['raw619_old4_plane_a1_nonempty28coefficients']}
 D=d['D_two_by10_by56'];Gm=[dict(v) for v in outer['Gamma_two']];piv=dict(outer['four_source_root_columns']);old4=outer['original4_old_pivots'];cols=dict(cols)
 saved={(j,t):{(p,f):(int(a),int(b)) for p,f,a,b in vv} for j,t,vv in d['positive617_combined_symbol_coefficients']}
 checks=0;channels=0;rootmax=0
 for rec in d['kernel617_source_recipes']:
  j=rec['free_case619'];cc=dict(c0[j]);aa=dict(a0.get(j,{}))
  for b,q in zip([563,609],rec['subtract_actual_basis563609'],strict=True):axpy(cc,c0[b],ex(q));axpy(aa,a0.get(b,{}),ex(q))
  for t0 in range(2):
   roots={}
   for col,w in cc.items():
    for row,node in cols.get(col,[]):
     for pair in range(10):
      if F(D[t0][pair][row]):key=('typed',node,pair);roots[key]=add(roots.get(key,Z),mul(w,ex(D[t0][pair][row])))
   for slot,w in aa.items():
    ip,var=divmod(slot,9);assert ip<4 and var in [7,8];pair=8 if var==8 else 9
    for row,(kind,node) in piv[old4[ip]]:
     if row in Gm[t0]:key=(kind,node,pair);roots[key]=add(roots.get(key,Z),mul(w,ex(Gm[t0][row])))
   rootmax=max(rootmax,len(roots));z={}
   for (kind,node,pair),w in roots.items():
    for f,v in (grad[node] if kind=='typed' else cg[node]).items():key=(pair,f);z[key]=add(z.get(key,Z),mul(w,v))
   env=saved.get((j,t0),{});assert set(z)==set(env)
   for key,(den,u) in z.items():
    a,b=env[key];assert F(u,den)==F(b,a),(j,t0,key,'positive envelope');checks+=1
   channels+=bool(z)
 dump('positive-source-checkpoint.json',{'checks':checks,'channels':channels,'root_support_upper':rootmax,'own_positive_maximum_bits':maximum,'producer_prior_upper':120,'no_signed_gradient_or_response':True})
 assert checks==1926 and channels==106 and rootmax<=72 and maximum<=120
 assert ledger['prior_input_numeric_count_including9_final_metadata']==count(d)==14866
 assert ledger['success_and_guard_upper']==2*count(d)+d['future_scratch_upper']+ledger['numeric_output_and_JSON_codec_upper']+2000==48788
 assert ledger['success_and_guard_upper']<ledger['unchanged_future_cap']==49698
 assert ledger['numeric_output_and_JSON_codec_upper']==2*(6395+237)==13264
 assert ledger['whole_RAW_upper']==4610546<r['future_RAW_cap']<32*1024*1024
 # Code/custody only: no producer imports or invocation beyond safe argparse help.
 code=(ROOT/W).read_text();tree=ast.parse(code)
 imports=[n for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom))]
 assert all(not isinstance(n,ast.ImportFrom) or n.module in ['fractions','pathlib'] for n in imports)
 assert "pair=8 if var==8 else 9" in code and 'times(-x,G[ch][row])' in code
 assert 'def plus(a,b):return ck(a+b)' in code and 'def times(a,b):return ck(a*b)' in code
 assert code.index('if b>120:')<code.index('return q',code.index('def ck('))
 h=subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/W),'--help'],capture_output=True,timeout=10)
 assert h.returncode==0 and hashlib.sha256(h.stdout).hexdigest()==c['safe_help_stdout_sha256']
 owned=count([r,d,c,ledger,s,sr,point79,point164,proof,gr,expected,grad,linear,zeros,ownnodes,cols,cg,c0,a0,D,Gm,piv,saved,pins])+12000
 assert owned<150000 and raw+1_000_000<32*1024*1024
 ep='independent-OLD2-source-zero-origins-positive617-signs-and-custody.json'
 ev={'all71_literal_ids_origin0_bound_no_seconds_zero_filled':True,'all757_origins0_no_empty_products':True,'all104_positive_linear_routes112_supports':True,
  'all1926_positive_coefficient_envelopes106_channels_checked':True,'all28_original_a1_plane_slots_bound':True,'outer_sign_and_rr_rell_ordinary_orientation_checked':True,
  'no_signed_OLD2_gradient_or_second_response_values':True,'native23_L2S0_and_nonlinear_firstcross_terms_separate':True,
  'all_literal_worker_custody_help_and_metadata_counts_verified':True,'source_prior120_temporary361':True,
  'worker_execution_whole48788_below49698':True,'worker_RAW4610546_below_declared_cap':True,'independent_owned_numeric_upper':owned}
 dump(ep,ev);rep={'id':STEM,'all_passed':True,'source_code_custody_only':True,'input_sha256':pins,
  'worker_sha256':sha(ROOT/W),'gate_sha256':pins[G+'report.json'],'custody_sha256':pins[C+'report.json'],
  'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
  'whole_numeric_upper':owned,'no_signed_sensitivity_or_response_evaluated':True,'seconds':time.monotonic()-t}
 dump('report.json',rep);dump('phase.json',{'phase':'FINAL_SOURCE_ONLY','report_sha256':sha(out/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owned,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
