#!/usr/bin/env python3
"""Independent literal/source/tensor/custody peer. No residual rank query."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,math,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple,set)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def fields(path,wanted):
 s=path.read_text();dec=json.JSONDecoder();i=s.index('{')+1;out={}
 def ws(i):
  while s[i].isspace():i+=1
  return i
 while True:
  i=ws(i)
  if s[i]=='}':break
  if s[i]==',':i=ws(i+1)
  k,i=dec.raw_decode(s,i);i=ws(i);assert s[i]==':';i=ws(i+1)
  if k in wanted:out[k],i=dec.raw_decode(s,i)
  else:
   quoted=escape=False;depth=0
   while i<len(s):
    c=s[i]
    if quoted:
     if escape:escape=False
     elif c=='\\':escape=True
     elif c=='"':quoted=False
    elif c=='"':quoted=True
    elif c in '[{':depth+=1
    elif c in ']}':
     if not depth:break
     depth-=1
    elif c==',' and not depth:break
    i+=1
 assert set(out)==set(wanted);return out
def sp(a):return {int(i):Q(v) for i,v in a if Q(v)}
def add(a,b):
 c=a.copy()
 for i,v in b.items():
  c[i]=c.get(i,Q(0))+v
  if not c[i]:del c[i]
 return c
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,d):(out/n).write_text(json.dumps(d,separators=(',',':'))+'\n')
 def manifest(p,expected=None):
  d=rd(p);assert d['all_passed']
  if expected:assert pins[p]==expected
  for k,h in d.get('input_sha256',{}).items():
   if isinstance(h,dict):assert sha(ROOT/h['path'])==h['sha256']
   else:assert sha(ROOT/k)==h
  for k,h in d.get('evidence_sha256',{}).items():assert sha((ROOT/p).parent/k)==h
  if 'code_sha256' in d:
   stem=d['id'];paths=list((ROOT/'code').glob('*'))
   assert any(sha(q)==d['code_sha256'] for q in paths if q.is_file())
  return d
 pb='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_SOURCE_PROFILE_20261105_v63/'
 gb='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63/'
 pr=manifest(pb+'report.json','effebaaa9b673f6a19bdde737ea9c06819b96db755f9cbfac136787a22ec9c6a')
 gr=manifest(gb+'report.json','9f1627a8f293869e96ab9854befc9e5e9c320f1a822103416d37fd63cbda2772')
 d=rd(pb+'literal-active-G2-linear-strand-matrix-singleton-units-and-residual-source.json')
 p=rd(gb+'pruned46-by41-integer-rank-source-and-live-bound.json');ad=rd(gb+'source-profile-count-addendum.json')
 put('phase.json',{'phase':'independent_source_before_rank','gate_sha256':pins[gb+'report.json'],'no_residual_rank_query':True})
 g2p='reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json'
 z=fld(g2p,['variable_order6','actual_linear_G2_active102_by42_columns'])
 assert z['variable_order6']==['zA0','zA1','zA4','zA5','rA','ell']
 G2=list(map(sp,z['actual_linear_G2_active102_by42_columns']));assert len(G2)==42
 mon=list(itertools.combinations_with_replacement(range(6),2));mon.remove((5,5));mi={m:i for i,m in enumerate(mon)}
 A=[]
 for col in G2:
  for v in range(6):
   q={}
   for k,x in col.items():
    row,u=divmod(k,6);m=tuple(sorted((u,v)))
    if m==(5,5):continue
    k=20*row+mi[m];assert k not in q;q[k]=x
   A.append(q)
 assert A==list(map(sp,d['literal_matrix_columns'])) and sum(map(len,A))==708
 assert d['ordinary_quadratic_pairs20']==list(map(list,mon))
 rows=[{j:col[i] for j,col in enumerate(A) if i in col} for i in range(340)]
 lr=set(range(340));lc=set(range(252));piv=[]
 while True:
  hit=next(((i,next(j for j in rows[i] if j in lc),'row_singleton') for i in sorted(lr) if sum(j in lc for j in rows[i])==1),None)
  if hit is None:hit=next(((next(i for i in A[j] if i in lr),j,'column_singleton') for j in sorted(lc) if sum(i in lr for i in A[j])==1),None)
  if hit is None:break
  i,j,k=hit;piv.append((i,j,k,A[j][i]));lr.remove(i);lc.remove(j)
 assert [(r['row340'],r['column252'],r['singleton_kind'],Q(r['literal_nonzero_unit'])) for r in d['fixed_singleton_unit_pivots']]==piv and len(piv)==162
 R=sorted(i for i in lr if any(j in lc for j in rows[i]));C=sorted(j for j in lc if any(i in lr for i in A[j]));assert [len(R),len(C)]==[46,41]
 assert R==d['remaining_nonzero_row_indices']==p['original340_row_indices46'] and C==d['remaining_nonzero_column_indices']==p['original252_column_indices41']
 assert [[j,[[i,str(A[j][i])] for i in sorted(A[j]) if i in lr]] for j in C]==d['remaining_literal_matrix_columns']
 clear=[];M=[]
 for i in R:
  D=math.lcm(*(x.denominator for j,x in rows[i].items() if j in lc));row=[int(A[j].get(i,Q(0))*D) for j in C];U=sum(map(abs,row));clear.append([i,str(D),str(U),U.bit_length()]);M.append(row)
 assert clear==d['remaining_row_integer_D_L1_bounds'] and M==p['residual_integer_matrix46_by41']
 assert sum(bool(x) for row in M for x in row)==102 and p['row_integer_scaling_D46']==[v[1] for v in clear]
 zeros=sorted(lc-set(C));assert len(zeros)==49 and zeros==p['additional_zero_source_columns49']
 bits=sorted([v[3] for v in clear],reverse=True);assert sum(bits[:28])==304 and sum(bits[:29])==306 and 2*306+1==613
 # Fixed known composition, bounded by independent common-denominator absolute sums BEFORE signed sums.
 g3p='reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_INPUT_GATE_20261105_v63/actual-full-G3-linear-zero3-split-and-active-kernel-input.json'
 gg=fld(g3p,['active_coordinates6','zero_variable_coordinates3','actual_active_G3linear252_by238_columns'])
 G3=list(map(sp,gg['actual_active_G3linear252_by238_columns']));assert len(G3)==238
 envelopes=[];prior=0
 for col in G3:
  terms={}
  for j,c in col.items():
   for i,x in A[j].items():terms.setdefault(i,[]).append((c,x))
  env={}
  for i,ts in terms.items():
   den=math.lcm(*(c.denominator*x.denominator for c,x in ts));U=sum(abs(c.numerator*x.numerator)*(den//(c.denominator*x.denominator)) for c,x in ts);prior=max(prior,den.bit_length(),U.bit_length());env[i]=(den,U)
  envelopes.append(env)
 assert prior<=gr['positive_product_prior_bits']<=384
 actual=0
 for col,env in zip(G3,envelopes):
  q={}
  for j,c in col.items():
   for i,x in A[j].items():q=add(q,{i:c*x});actual=max(actual,abs((c*x).numerator).bit_length(),(c*x).denominator.bit_length())
  assert not q
 # Reconstruct the SAME full tensor coefficients, without a rank computation.
 va=gg['active_coordinates6'];vz=gg['zero_variable_coordinates3'];assert va==[0,1,4,5,7,8] and vz==[2,3,6]
 gp='reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_INPUT_GATE_20261105_v63/actual-full-G1-and186-mixed-adjoint-source-input.json'
 g1=fld(gp,['exact_full_G1_polynomial_columns'])['exact_full_G1_polynomial_columns'];L1=[]
 for col in g1[:17]:
  q={}
  for target,ex,x in col:
   if sum(ex)==1:q[target*9+ex.index(1)]=Q(x)
  L1.append(q)
 full2=[{(k//6)*9+va[k%6]:x for k,x in c.items()} for c in G2]
 for iz,zv in enumerate(vz):
  for j in range(17):full2.append(add({j*9+zv:Q(-1)},{(17+3*iz+k//9)*9+k%9:x for k,x in L1[j].items()}))
 pairs=list(itertools.combinations(range(3),2))
 for a,b in pairs:
  for j in range(3):full2.append({(17+3*b+j)*9+vz[a]:Q(1),(17+3*a+j)*9+vz[b]:Q(-1)})
 fp='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 assert full2==list(map(sp,fld(fp,['full_G2_linear234_by102_columns'])['full_G2_linear234_by102_columns']))
 full3=[{(k//6)*9+va[k%6]:x for k,x in c.items()} for c in G3]
 for iz,zv in enumerate(vz):
  for j in range(42):full3.append(add({j*9+zv:Q(1)},{(42+17*iz+k//6)*9+va[k%6]:x for k,x in G2[j].items()}))
 for ip,(a,b) in enumerate(pairs):
  for j in range(17):full3.append(add({(42+17*b+j)*9+vz[a]:Q(-1),(42+17*a+j)*9+vz[b]:Q(1)},{(93+3*ip+k//9)*9+k%9:x for k,x in L1[j].items()}))
 for a in range(3):full3.append({(93+3*2+a)*9+vz[0]:Q(1),(93+3*1+a)*9+vz[1]:Q(-1),(93+a)*9+vz[2]:Q(1)})
 fp='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'
 assert len(full3)==418 and full3==list(map(sp,fld(fp,['actual_full_G3linear102_by418_columns'])['actual_full_G3linear102_by418_columns']))
 rp='reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63_REVIEW/report.json';rr=rd(rp);assert rr['all_passed'] and pins[rp]=='8c4bd870deb749165c9afe69d7a06bd7d32823fcb2cedcfc1a8474922e517864'
 rank=fld('reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63/actual-full-G3linear-active-fixed-rank-and-kernel.json',['active_rank','full_G3linear_rank']);assert rank=={'active_rank':62,'full_G3linear_rank':242}
 assert ad['corrected_source_phase_upper']==16784+1664+2008==20456 and ad['rank_worker_cap_unchanged']==16784
 c=rd('reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_PREPARED_CUSTODY_20261105_v63/report.json')
 assert c['all_passed'] and pins['reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_PREPARED_CUSTODY_20261105_v63/report.json']=='3fc1f4e7aa9e7ff66d039a34f6e11b29bb9e12c29268dd539bfc15dcf274bff1'
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 worker=ROOT/c['worker'];assert sha(worker)==c['worker_sha256']=='ff2546259e773dccf58c3dabbd10666fdfacba8a959ae98d422d3a098107bacc'
 helpresult=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,timeout=10);assert helpresult.returncode==0 and hashlib.sha256(helpresult.stdout).hexdigest()==c['safe_help_stdout_sha256']
 assert p['failure_or_success_live_bound']['sum']==14886==7000+1886+500+4500+1000<16784 and count(p)<3500
 assert p['hard_committed384']==384 and p['hard_temporary1153']==1153 and p['all_minors_order_le29_bits306']==306 and p['integer_product_difference_bits613']==613
 owned=count([pr,gr,d,p,ad,z,G2,A,rows,piv,R,C,clear,M,gg,G3,envelopes,g1,L1,full2,full3,rank,c,pins])+12000;assert owned<150000
 proof={'all252_literal_source_columns_reconstructed':True,'all162_singleton_units_replayed_without_fill':True,'residual_integer_matrix_exact_no_rank_query':True,'all238_A_G3_zero_checked':True,'positive_known_product_prior_bits':prior,'all102_fullG2_and418_fullG3_same_tensor_frame':True,'conditional_tensor_theorem':'linear complex over active graded ring k[z0,z1,z4,z5,r,ell]/ell², extended by three independent zero variables, is total tensor with their Koszul resolution. The augmentation is a quasi-isomorphism for this bounded-below free complex. Therefore active H2 internal degree3 zero (IF residual rank28) implies full918 coefficient-degree1 kernel equals image fullG3linear242. No general associated-graded exactness.','Bareiss_proof':'Each current Schur entry is an ordered selected minor (arbitrary consistent pivot row/column orders). Integer rowL1 Hadamard bounds all minors through29 by306bits; p*a-b*c uses at most613bits, exact previous-pivot division yields next minor. 162 singleton Laplace units plus residual fixed minor prove rank. Zero Schur verifies full span; reverse source substitutions give complete kernel.','worker_manual_audit':'only pruned46x41 integer arithmetic; no residual rank arithmetic invoked by peer; exact division/frontier/STOP and source reverse recipes reviewed','numeric_worker_bound':14886,'source_count_addendum':20456,'source_peer_owned_bound':owned,'full_b1_functional_annihilation_separate':True,'no_C_K_upper_or_higher_response_values':True}
 ep='independent-literal-tensor-custody-and-rank-gate-proof.json';put(ep,proof)
 report={'id':STEM,'all_passed':True,'SOURCE_code_custody_only_no_residual_rank':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(out/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned}))
if __name__=='__main__':main()
