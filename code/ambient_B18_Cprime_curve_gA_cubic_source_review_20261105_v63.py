#!/usr/bin/env python3
"""Independent curve-column SOURCE/frame and unrun one-dot audit."""
import hashlib,itertools,json,math,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def selected(path,names):
 text=path.read_text();d={};dec=json.JSONDecoder()
 for name in names:
  token=json.dumps(name);assert text.count(token)==1;pos=text.index(':',text.index(token))+1
  while text[pos].isspace():pos+=1
  d[name]=dec.raw_decode(text,pos)[0]
 return d
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_PREFLIGHT_20261105_v63';r=load(primary,h='1119ab28d3dcf8f2b4175adcc59394fb6e62516ae7ed62ab6b47efca3b10a21d')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 e=load(primary,'actual-curve-source-and-one-residual-cubic-functional.json',h='d3005ffbd5e82a0f37227884e13c6e9ed67b66d54a9137181814656638d86f7b')
 prior='AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW';assert load(prior,h='481862051b9e604ac92163b452f88e0cb46682d043393668b6648a89610e3559')['known_cubic_rank']==26
 evidence=load(prior,'independent-known-cubic-fixed-certificate-and-open-dual.json');dual=pairs(e['residual_actual_cubic_dual495_sparse']);q27=pairs(e['residual_compressed_dual27_sparse']);assert dual=={int(i):F(c)for i,c in evidence['known_quotient_dual495'].items()}and q27=={int(i):F(c)for i,c in evidence['known_quotient_dual27'].items()}
 cpath=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';C=selected(cpath,['actual_C0_basis8_scalar_nodes','all35_regular_coordinate_divisions','cached_actual_CA_gA_all_degree3_coefficients']);assert load('AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW',h='9e3f2de68224ba1696888b21a51de61f1cded2dfbca258cd856bdc1d99bd2882')['all_passed']
 tail=C['all35_regular_coordinate_divisions'][0];records=tail['all7_regular_coordinate_divisions'];assert tail['original_pure_CA_source_function']==7790 and tail['restriction_r_s_node']==505938==C['actual_C0_basis8_scalar_nodes'][-1]and len(records)==7
 subset=load('AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE2_SUBSET_20261105_v63','minimal-two57-source2-and-labels.json',h='9c3cd8d7a91b9f079bf36b348daabaa5348b91d0ce5f4674cd9d9cbfbcaeaf1a');assert subset['both57_source2']['CA']['F9_nodes'][:7]==C['actual_C0_basis8_scalar_nodes'][:7]and subset['both57_source2']['CA']['F9_nodes'][7:]==[17493,19106]
 ids={rec[k]for rec in records for k in ['node','numerator','coordinate','current','previous']}|{rec['numerator']-1 for rec in records};text=cpath.read_text();pos=text.index(':',text.index('"DAG_scalar_nodes"'))+1
 while text[pos].isspace():pos+=1
 assert text[pos]=='[';pos+=1;dec=json.JSONDecoder();nodes={}
 for j in range(max(ids)+1):
  while text[pos].isspace()or text[pos]==',':pos+=1
  node,pos=dec.raw_decode(text,pos)
  if j in ids:nodes[j]=node
 del text
 previous=505938
 for j,rec in enumerate(records):
  assert rec['previous']==previous and rec['coordinate']==C['actual_C0_basis8_scalar_nodes'][j]and F(rec['origin'])==0
  assert nodes[rec['node']]==['r',rec['numerator'],rec['coordinate'],f'C0_source7790:zero{j}']and nodes[rec['numerator']]==['a',rec['current'],rec['numerator']-1]and nodes[rec['numerator']-1]==['n',previous];previous=rec['current']
 assert previous==7790 and records==e['all7_ordered_tail_division_source_records']
 joint=load('AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63','joint-degree3-certificate.json');assert joint['monomials']==[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]]and joint['g_CA_all10_coefficients']==C['cached_actual_CA_gA_all_degree3_coefficients']
 old=list(map(F,C['cached_actual_CA_gA_all_degree3_coefficients']));assert old[:6]==[F(0)]*6;changed=[sum((c*F(math.comb(s,k),5**k)for s,c in enumerate(old[6:])if s>=k),F(0))for k in range(4)];assert changed==list(map(F,e['gA_changed_r_ell_coefficients']))==[F(125,24),F(5,3),F(1,6),F(2,375)]
 # Constant native -> gamma basis inverse, applied symbolically to -gA.
 P=[[F(1),0,0],[F(11,12),0,1],[0,1,0]];Pi=[[1,0,0],[0,0,1],[F(-11,12),1,0]];assert all(sum(P[i][k]*Pi[k][j]for k in range(3))==int(i==j)for i in range(3)for j in range(3));assert [-row[0]for row in Pi]==[F(-1),0,F(11,12)]
 m3=list(itertools.combinations_with_replacement(range(9),3));curve={}
 for k,c in enumerate(changed):j=m3.index((7,)*(3-k)+(8,)*k);curve[j]=-c;curve[330+j]=F(11,12)*c
 assert curve==pairs(e['actual_curve_relation495_sparse'])and len(curve)==8
 assert all(m3[i%165][0]>=7 for i in dual)and all(dual.get(165*t+m3.index((8,8,8)),F(0))==0 for t in range(3))
 prep=load('AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_PREPARED_WORKER_20261105_v63',h='bee6d11169ed6a1803a6a69c5a2129477bf0d1558b7d59d958ec716e890ce05f');worker=ROOT/prep['prepared_code'];assert sha(worker)==prep['prepared_code_sha256']=='c3107dce2c780651b94ebd4d8e2910ad8c346cb9c8f778ba6bcda15b483e7c01'
 assert max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())for c in list(dual.values())+list(curve.values()))<=17 and 8*34+3==275<2048 and 275+15<300
 assert prep['numeric_cap']==8000 and prep['full_hashed_input_raw_bound_bytes']==15000000<32<<20
 out=ROOT/'reports/AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_VALUES_20261105_v63';assert not out.exists();run=subprocess.run(['python',str(worker),'--help'],capture_output=True,text=True,timeout=15,cwd=ROOT);assert run.returncode==0 and not run.stderr and not out.exists();(OUT/'independent-help-stdout.txt').write_text(run.stdout)
 save('independent-curve-source-unit-column-change-and-one-dot-gate.json',dict(actual_coordinate_source_nodes=C['actual_C0_basis8_scalar_nodes'],telescoping_nodes=nodes,all7_literal_difference_and_regular_division_nodes=True,inherited_divisibility='AcceptedhonestC0review proves eachorderedIFTrestriction difference vanishes at itsregularF9coordinate. The r nodes denote these certified Hadamardquotients, notunknownfactor inverses.',identity='f63=gA+sum zi*ei; curve depends only CA F9(r,s).',source_column_change='Replacecol43 bycol43−sum ei*col(36+i). Itsinverse adds those terms; offdiagonalupdate square0 anddet1. It preserves actualCprimeimage without evaluating ei.',old_curve_coefficients=old,changed_curve_coefficients=changed,Gamma_source_vector=curve,native_gamma_basis_matrix=P,native_gamma_inverse=Pi,residual_dual27=q27,residual_dual495=dual,source_pairing_NOT_evaluated=True,known26_plus_nonzero_pairing_proof='In rows oldfixed26 then25 and columns oldfixed26 thencurve, det27=det26*(q27*curve). Nonzero killsall27knowncubicclasses; combinedwithactualHF333/lower9 proveslength9/e2. Zeroleavesonlythiscandidate negative.',unrun_worker_import_only_pass=True,QQ_cap=8000,input_bit_bound=17,scalar_bit_bound=275,source_or_gauge_or_jet_evaluation=False))
 for p in[Path(__file__),worker,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir()if p.is_file());assert raw<32<<20 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,actual_curve_source_unipotent_column_change_and_frame_accepted=True,single_scalar_worker_gate_ready=True,source_pairing_NOT_evaluated=True,numeric_cap=8000,no_new_source_gauge_or_physical_jet=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=raw,seconds=time.monotonic()-st))
 print(json.dumps(dict(all_passed=True,pairing_not_evaluated=True,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
