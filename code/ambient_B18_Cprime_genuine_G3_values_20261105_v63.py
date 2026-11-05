#!/usr/bin/env python3
"""Prepared cached9 actions and genuineG3 using fixed unit block witnesses."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];GATE='AMBIENT_B18_CPRIME_GENUINE_G3_SOURCE_PREFLIGHT_20261105_v63';STEM='AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def count(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def guard(v):
 if isinstance(v,Q):assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=2048
 elif isinstance(v,dict):
  for x in v.values():guard(x)
 elif isinstance(v,(tuple,list)):
  for x in v:guard(x)
 return v
def main(peer_path):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={}
 def read(path):p=ROOT/path;pins[path]=sha(p);return json.loads(p.read_text())
 gp=f'reports/{GATE}/report.json';r=read(gp);peer=read(peer_path);assert r['all_passed']and peer['all_passed']and peer['input_sha256'][gp]==sha(ROOT/gp)
 exact_peer=read('reports/AMBIENT_B18_CPRIME_EXACT_TORSION_12_20261105_v63_REVIEW/report.json');assert exact_peer['all_passed']
 for path,h in r['input_sha256'].items():assert sha(ROOT/path)==h; pins[path]=h
 ep=f'reports/{GATE}/actual-nine-action-and-genuine-G3-source-interface.json';e=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 qp=next(p for p in r['input_sha256']if p.endswith('phase4-all-twenty-actual-quadratic-values.json'));raw=read(qp)['all44_actual_V2_pureplane_nine_values'];L=[rec(v)for v in e['source_L44_columns27_sparse']];cp=e['source_V1_pivot24_columns'];rows=e['source_V1_pivot24_rows'];iv=[[Q(c)for c in row]for row in e['source_V1_fixed24_inverse']];rowpos={v:i for i,v in enumerate(rows)}
 q1=[rec(v)for v in e['first_NF3_rows27_sparse']];q3=[rec(v)for v in e['quadratic_NF3_rows135_sparse']];E1=e['first_basis_embedding_rows27'];E2=e['quadratic_basis_embedding_rows135'];m2=list(itertools.combinations_with_replacement(range(9),2));mi={m:i for i,m in enumerate(m2)}
 Qnf=[]
 for col in raw:
  v={45*target+mi[mon]:Q(col[3*target+k])for target in range(3)for k,mon in enumerate([(7,7),(7,8),(8,8)])if Q(col[3*target+k])};Qnf.append([sum((c*v.get(i,Q(0))for i,c in q.items()),Q(0))for q in q3])
 def unit(i):return [Q(j==i)for j in range(9)]
 zero=[Q(0)]*9
 def linear(row):
  v=zero.copy();v[3:6]=[q.get(row,Q(0))for q in q1]
  if row in rowpos:
   k=rowpos[row];v[6:9]=[-sum((Qnf[cp[j]][f]*iv[j][k]for j in range(24)),Q(0))for f in range(3)]
  return v
 def quadratic(row):return [Q(0)]*6+[q.get(row,Q(0))for q in q3]
 actions=[]
 for var in range(9):
  columns=[linear(9*target+var)for target in range(3)]
  for row in E1:target,x=divmod(row,9);columns.append(quadratic(45*target+mi[tuple(sorted((x,var)))]))
  columns+= [zero.copy()for _ in range(3)];actions.append(guard(columns))
 def apply(A,v):return guard([sum((c*A[j][i]for j,c in enumerate(v)if c),Q(0))for i in range(9)])
 relationhash=hashlib.sha256()
 for j,col in enumerate(L):
  result=[Q(0)]*6+Qnf[j]
  for row,c in col.items():result=[a+c*b for a,b in zip(result,linear(row))]
  assert result==zero;relationhash.update(json.dumps(enc(result)).encode())
 commhash=hashlib.sha256()
 for a in range(9):
  for b in range(9):
   for j in range(9):
    v=apply(actions[a],actions[b][j]);assert v==apply(actions[b],actions[a][j]);commhash.update(json.dumps(enc(v),separators=(',',':')).encode())
 # Thisstrictblocksupport provesall729tripleproducts identicallyzero.
 assert all(not any(col[:3])for A in actions for col in A)
 assert all(not any(actions[v][j][:6])for v in range(9)for j in range(3,6))
 assert all(actions[v][j]==zero for v in range(9)for j in range(6,9))
 cyclic=[unit(j)for j in range(3)]
 for row in E1:target,var=divmod(row,9);cyclic.append(actions[var][target])
 for row in E2:target,k=divmod(row,45);u,v=m2[k];cyclic.append(apply(actions[v],actions[u][target]))
 assert cyclic==[unit(j)for j in range(9)]
 ell=actions[8];S=[[ell[j][6+i]for j in range(3)]for i in range(3)]
 assert [col[3:6]for col in ell[:3]]==[[Q(1),Q(0),Q(0)],[Q(0)]*3,[Q(0),Q(0),Q(1)]]
 assert ell[3:6]==[unit(6+i)for i in range(3)]
 ell2=[apply(ell,col)for col in ell];assert ell2==[unit(6),zero,unit(8)]+[zero]*6
 K=unit(1)
 for i in range(3):K[3+i]=-S[i][1]
 assert apply(ell,K)==zero;kernel=[K,unit(6),unit(7),unit(8)];assert all(apply(ell,c)==zero for c in kernel)
 quotient=[[Q(i==1)for i in range(9)],[Q(i==7)for i in range(9)]];basis=[K,unit(7)]
 def project(v):return [v[1],v[7]]
 assert [project(c)for c in basis]==[[Q(1),Q(0)],[Q(0),Q(1)]]and all(project(c)==[Q(0),Q(0)]for c in ell2)
 target=[];alphas=[]
 for var,A in enumerate(actions):
  AK=apply(A,K);assert not any(AK[:6])and apply(ell,AK)==zero and apply(A,basis[1])==zero
  col=project(AK);assert col[0]==0;target.append([col,[Q(0),Q(0)]]);alphas.append(col[1])
 assert alphas[8]==0
 data={'basis9_order':['H0','H1','H2','F0=gamma0ell','F1=gamma2r','F2=gamma2ell','Q0=gamma0ell2','Q1=gamma2rell','Q2=gamma2ell2'],
  'actual_R0_actions9_by9_columns':actions,'all44_actual_relation_reduction_sha256':relationhash.hexdigest(),'all81_commutation_checks_sha256':commhash.hexdigest(),'strict3layer_support_proves_m3zero':True,'all9_literal_cyclic_basis_vectors_I9':True,
  'ell_head_quadratic_block_S3_by3':S,'ell_rank5_unit_minor_rows':[3,5,6,7,8],'ell_rank5_unit_minor_columns':[0,2,3,4,5],'ell_rank5_determinant':Q(1),
  'ell_square_columns9':ell2,'ell_square_rank2_I2_rows':[6,8],'ell_square_rank2_I2_columns':[0,2],
  'all4_actual_kernel_ell_basis_vectors9':kernel,'image_ell_square_basis_vectors9':[unit(6),unit(8)],'genuine_Cprime_ell_torsion_basis_vectors9':basis,
  'genuine_G3_projection_on_kernel_rows2_by9':quotient,'actual_genuine_G3_R0_actions2_by2_columns':target,'actual_genuine_G3_action_scalars9':alphas,
  'target_actions_K_to_alpha_Q1_and_Q1_to0':True,'target_m2zero_by_source_block_form':True,'genuine_target_dimension':2,
  'Gamma_from_native_C4drop1':'(a,b,c) native(tauA,tauB,tauAB) -> (a,c,b-11*a/12); apply BEFORE thisfinite9projection',
  'new_source_gauge_jet_or_connecting_values':False,'no_original_nil_generators_or_cube_claim':True}
 owned=count([r,peer,exact_peer,e,raw,L,iv,q1,q3,Qnf,actions,ell2,S,K,kernel,basis,target,data]);assert owned<=24000,owned;guard(data)
 blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert sum((ROOT/p).stat().st_size for p in pins)+len(blob.encode())<32<<20;ev=out/'actual-Cprime9-actions-and-genuine-G3-target.json';ev.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_cached_actual_Cprime9_actions_and_genuine_G3','actual_Cprime_ell3_length':9,'genuine_Cprime_ell_torsion_dimension':2,'genuine_G3_action_scalars9':list(map(str,alphas)),'new_source_gauge_jet_or_connecting_values':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'owned_numeric_slots':owned,'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source-peer-report',required=True);main(p.parse_args().source_peer_report)
