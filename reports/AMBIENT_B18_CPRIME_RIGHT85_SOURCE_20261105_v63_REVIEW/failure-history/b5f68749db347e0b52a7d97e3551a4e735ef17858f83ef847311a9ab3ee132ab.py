#!/usr/bin/env python3
"""Independent full kernel witnesses and cached-only adjoint SOURCE audit."""
import hashlib,itertools,json,math,re,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_RIGHT85_SOURCE_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def plus(v,w,c=F(1)):
 z=dict(v)
 for i,a in w.items():
  z[i]=z.get(i,F(0))+c*a
  if not z[i]:del z[i]
 return z
def numeric(v):
 if isinstance(v,(F,int)):return 1
 if isinstance(v,str):return int(re.fullmatch(r'-?\d+(?:/\d+)?',v)is not None)
 if isinstance(v,(list,tuple)):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_RIGHT85_SOURCE_PREFLIGHT_20261105_v63';r=load(primary,h='242319f5b5fffdd0902820d2a0928786b9f2d44906477a6c6356c3e8be55478f')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 e=load(primary,'complete-right85-and-cached-plane-source-adjoints.json',h='26aaa12be22f0aa761a7afc91a43a33f5252673a2d586c7dc6abbd4479d3600f')
 first=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase1-actual-fullnine-first-V-columns.json');cp=first['fixed_pivot_columns'];L=[pairs(c)for c in first['all44_actual_linear_columns27_sparse']];assert len(cp)==24 and cp==e['actual_V1_fixed24_source_columns']
 assert load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW',h='5faf10001ed74c2b9e107d7fa1f4fe873600423affa75305a75a8f1d3fad5bb3')['all_passed']
 qv=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63','phase4-all-twenty-actual-quadratic-values.json');raw=qv['all44_actual_V2_pureplane_nine_values'];m2=list(itertools.combinations_with_replacement(range(9),2));i2={m:i for i,m in enumerate(m2)};m3=list(itertools.combinations_with_replacement(range(9),3));i3={m:i for i,m in enumerate(m3)};columns=[]
 for j in cp:
  for x in range(9):
   v={}
   for i,c in L[j].items():t,y=divmod(i,9);v=plus(v,{45*t+i2[tuple(sorted((x,y)))]:c})
   columns.append(v)
 assert columns==[pairs(c)for c in e['all216_first_product_columns135_sparse']]
 pc=e['rank131_pivot_original_columns'];pr=e['rank131_pivot_rows_creation_order'];K=[pairs(k)for k in e['all85_syzygy_columns216_sparse']];free=e['free_source_rows85'];assert len(set(pc))==len(pc)==len(set(pr))==len(pr)==131 and free==[j for j in range(216)if j not in pc]and len(K)==len(free)==85
 for k,v in enumerate(K):
  out={}
  for j,c in v.items():out=plus(out,columns[j],c)
  assert not out and [v.get(j,F(0))for j in free]==[F(i==k)for i in range(85)]
 seen=set()
 for row,w,a in e['rank131_sparse_normalized_pivot_and_source_witnesses']:
  w,a=pairs(w),pairs(a);out={}
  for j,c in a.items():out=plus(out,columns[j],c)
  assert out==w and min(w)==row and w[row]==1 and row not in seen;seen.add(row);k=pr.index(row);assert set(a).issubset(pc[:k+1])and a.get(pc[k],F(0))
 assert seen==set(pr)
 digest=hashlib.sha256()
 for k in K:digest.update(json.dumps(enc([[i,c]for i,c in sorted(k.items())]),separators=(',',':')).encode()+b'\n')
 assert digest.hexdigest()==e['kernel_sha256']
 old=load('AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW','independent-known-cubic-fixed-certificate-and-open-dual.json');dual={int(i):F(c)for i,c in old['known_quotient_dual495'].items()};assert all(m3[i%165][0]>=7 for i in dual)
 folded=[];zcount=0
 for k in K:
  adj={}
  for source,c in k.items():
   j,x=divmod(source,9)
   if x<7:zcount+=1;continue
   for target in range(3):
    for p,mon in enumerate([(7,7),(7,8),(8,8)]):
     row=165*target+i3[tuple(sorted(mon+(x,)))];a=c*dual.get(row,F(0))
     if a:adj=plus(adj,{9*cp[j]+3*target+p:a})
  folded.append(adj)
 assert folded==[pairs(a)for a in e['all85_cached_plane_Q2_adjoints396_sparse']]and zcount==e['z_factor_callback_terms_proved_zero']
 leaves=sorted(set(i for a in folded for i in a));assert leaves==e['all_required_actual_leaves']and len(leaves)==115 and sum(map(len,folded))==463 and sum(bool(a)for a in folded)==40 and max(map(len,folded))==25
 # Check cache availability/height only; do not form any adjoint-value dot.
 actual={i:F(raw[i//9][i%9])for i in leaves};assert max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())for c in actual.values())<=14 and max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())for a in folded for c in a.values())<=23
 prep=load('AMBIENT_B18_CPRIME_RIGHT85_PREPARED_WORKER_20261105_v63',h='31d88529263efefceff2a103859c5d61bfb1481c0462277baee0d8558c2089d3e');worker=ROOT/prep['prepared_code'];assert sha(worker)==prep['prepared_code_sha256']=='67c0a19e9c79cb0a67d2a726b0f0f1d5a6f69808a95ff2060a4cd616057182bc';snap=ROOT/'reports/AMBIENT_B18_CPRIME_RIGHT85_PREPARED_WORKER_20261105_v63/source-history'/f'{sha(worker)}.py';assert sha(snap)==sha(worker);pins[str(snap.relative_to(ROOT))]=sha(snap)
 assert prep['static_numeric_envelope']==22000<prep['numeric_cap']==60000 and 25*(23+14)+math.ceil(math.log2(25))==930<2048
 conservative=numeric(e)+numeric(qv)+numeric(r)+1000+396+85*(1+2*12)+2000+8000;assert conservative<22000
 valueout=ROOT/'reports/AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63';assert not valueout.exists();run=subprocess.run(['python',str(worker),'--help'],capture_output=True,text=True,timeout=15,cwd=ROOT);assert run.returncode==0 and not run.stderr and not valueout.exists();(OUT/'independent-help-stdout.txt').write_text(run.stdout)
 save('independent-complete85-source-kernel-and-cached-adjoint-gate.json',dict(all216_literal_firstproducts=True,all85_full135_multiplybacks=True,all85_free_tail_identity_rows=True,all131_sparse_pivot_source_multiplybacks=True,rank131_selected_column_transform_triangular=True,complete_kernel_dimension85=True,all463_folded_coefficients=True,active_adjoint_rows=40,available_cached_leaf_count=115,z_factor_terms_vanish_by_literal_pureplane_support=zcount,actual_values_NOT_paired=True,initial_source_proof='Each kernel vector is a homogeneouslinear coefficient vector on24literal originalVcolumns. V0=0 andV1u1=0 inALL135positions, soV*u1 beginsin degree3 with V2u1. No higher source-column change,M3,constant-kernel correction or gauge3 enters.',nonplane_proof='Residual495dual is supported onlypure r,ell cubics. A z factor or a nonplane quadratic cannot hit it; all retained leaf values are the accepted originalGamma-frame Q20 pureplane coefficients.',scope='Complete85 first-linear-syzygy block only. The19 constant-kernel cubic sources remain outside this test; allzero is onlypartialnegative.',conservative_numeric_bound=conservative,accepted_envelope=22000,hard_numeric_cap=60000,max_adjoint_bits=23,max_cached_leaf_bits=14,max_dot_terms=25,rational_dot_bit_bound=930,worker_import_only_pass=True))
 for p in[Path(__file__),worker,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rawbytes=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir()if p.is_file());assert rawbytes<2_000_000 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,complete_right85_kernel_source_and_adjoint_gate_passes=True,all85_values_NOT_evaluated=True,one_cached_value_job_ready=True,whole_numeric_envelope=22000,hard_numeric_cap=60000,derived_dot_bits=930,no_new_physical_source_gauge_or_jet=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=rawbytes,seconds=time.monotonic()-st))
 print(json.dumps(dict(all_passed=True,complete_kernel85=True,values_not_evaluated=True,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
