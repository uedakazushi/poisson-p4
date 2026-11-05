#!/usr/bin/env python3
"""Independent cached finite9/genuine target SOURCE theorem and custody."""
import hashlib,itertools,json,math,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_GENUINE_G3_SOURCE_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_GENUINE_G3_SOURCE_PREFLIGHT_20261105_v63';r=load(primary,h='dcf9eb35182a9ed73a7f281e5f0678037a54b20f15ce440c6de2a26d197ec813')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 e=load(primary,'actual-nine-action-and-genuine-G3-source-interface.json',h='286f72eae2164fda084913036c4416097def4b772ba8ce9e06d1318ff0132baf')
 assert load('AMBIENT_B18_CPRIME_EXACT_TORSION_12_20261105_v63_REVIEW',h='d7ec8b85b6622793fd4d5c06dfe96b1fe4b02d7ddea9ef9be3ebcfd837587f0b')['genuine_target_kerell_mod_imell2_canonical']
 first=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase1-actual-fullnine-first-V-columns.json');L=[pairs(v)for v in first['all44_actual_linear_columns27_sparse']];cp=first['fixed_pivot_columns'];rr=first['fixed_pivot_rows'];iv=[[F(c)for c in row]for row in first['fixed_minor_inverse']];free=[i for i in range(27)if i not in rr];assert free==[8,25,26]
 assert e['source_L44_columns27_sparse']==first['all44_actual_linear_columns27_sparse']and e['source_V1_pivot24_columns']==cp and e['source_V1_pivot24_rows']==rr and e['source_V1_fixed24_inverse']==first['fixed_minor_inverse']
 A=[[L[j].get(i,F(0))for j in cp]for i in rr];assert all(sum((A[i][k]*iv[k][j]for k in range(24)),F(0))==F(i==j)and sum((iv[i][k]*A[k][j]for k in range(24)),F(0))==F(i==j)for i in range(24)for j in range(24))
 q1=[]
 for f in free:
  row={f:F(1)}
  for j,k in enumerate(rr):
   c=-sum((L[cp[i]].get(f,F(0))*iv[i][j]for i in range(24)),F(0))
   if c:row[k]=c
  q1.append(row)
 assert q1==[pairs(v)for v in e['first_NF3_rows27_sparse']]and all(sum((c*v.get(i,F(0))for i,c in row.items()),F(0))==0 for row in q1 for v in L)
 q2=[pairs(v)for v in e['quadratic_NF3_rows135_sparse']];E2=e['quadratic_basis_embedding_rows135'];assert E2==[44,133,134]and [[q.get(i,F(0))for i in E2]for q in q2]==[[F(i==j)for j in range(3)]for i in range(3)]
 assert [[q.get(9*t+8,F(0))for q in q1]for t in range(3)]==[[1,0,0],[0,0,0],[0,0,1]]and all(q.get(9+v,F(0))==0 for q in q1 for v in range(9))
 # No action values are formed: this is a bound on their known input terms.
 Q=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63','phase4-all-twenty-actual-quadratic-values.json');raw=[list(map(F,v))for v in Q['all44_actual_V2_pureplane_nine_values']];m2=list(itertools.combinations_with_replacement(range(9),2));plane={(7,7):0,(7,8):1,(8,8):2}
 prep=load('AMBIENT_B18_CPRIME_GENUINE_G3_PREPARED_WORKER_20261105_v63',h='e5e1292215cd6d15d4d28016e9ee351e1a1bb9333353bae782bc5b51d1fe5c6c');den=int(prep['source_action_common_denominator']);bound=int(prep['source_integer_action_numerator_bound']);maxterm=F(0);product_count=0
 for q in q2:
  for i,c in q.items():
   t,mi=divmod(i,45);idx=3*t+plane[m2[mi]]
   for j,source in enumerate(cp):
    for a in iv[j]:
     x=c*raw[source][idx]*a
     if x:assert (x*den).denominator==1;maxterm=max(maxterm,abs(x));product_count+=1
 assert all((c*den).denominator==1 for q in q1+q2 for c in q.values())
 coarse=math.ceil(den*max(F(1),216*maxterm));assert coarse<=bound and (9*bound**2).bit_length()+32<=181<2048
 assert prep['static_numeric_envelope']==18000<prep['hard_numeric_cap']==e['future_numeric_cap']==24000 and prep['raw_plus_output_envelope_bytes']==500000<32<<20
 worker=ROOT/prep['prepared_code'];assert sha(worker)==prep['prepared_code_sha256']=='5830ab8de63620ba279c0b7dcf9e302383eb631a02a2bb3ab356735b8c070f36';snap=ROOT/'reports/AMBIENT_B18_CPRIME_GENUINE_G3_PREPARED_WORKER_20261105_v63/source-history'/f'{sha(worker)}.py';assert sha(snap)==sha(worker);pins[str(snap.relative_to(ROOT))]=sha(snap)
 text=worker.read_text()
 for token in ['assert result==zero','for a in range(9):','for b in range(9):','assert cyclic==[unit(j)for j in range(9)]','assert ell2==[unit(6),zero,unit(8)]+[zero]*6','assert owned<=24000','assert alphas[8]==0','signal.alarm(60)']:assert token in text
 valueout=ROOT/'reports/AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63';assert not valueout.exists();run=subprocess.run(['python',str(worker),'--help'],capture_output=True,text=True,timeout=15,cwd=ROOT);assert run.returncode==0 and not run.stderr and not valueout.exists();(OUT/'independent-help-stdout.txt').write_text(run.stdout)
 proof=dict(exact_source_basis9=['H0','H1','H2','F0=gamma0ell','F1=gamma2r','F2=gamma2ell','Q0=gamma0ell²','Q1=gamma2rell','Q2=gamma2ell²'],all24_first_inverse_products=True,all132_first_quotient_annihilations=True,quadratic_knownK2_exact_source='InM/m³ each V_j=L_j+Q_j, each variable multiple is its linearpartproduct, andtruequadraticinitials becomeexactrelations. Thuskerpi3 is anactualrelation space; removing Qnonplane meanssubtracting suchrelations, neverzero-fillingunknowns.',head_normalform='For linearinpute, subtractactualpivotrelations witht=I24*e_R. Output(0,q1e,−pi3Q_piv*t). Forquadraticf output(0,0,pi3f). ActualHF333 makes thecyclic9-dimensional quotientan isomorphism.',constant_kernel_compatibility='WholeconstantkernelK20 has rank1NF4image killedbyNF3. Thereforepi3(Q_j−Q_piv*I24*L_j,R)=0 forALL44; worker mustrechecktheliteral44relations.',mandatory_module_value_checks=['All44 source relations','All81 actioncommutators onall9basisvectors','Strict3-layer support andcyclicI9','Fixedunitrank5ellminor andell²I2','All4kerellvectors andimell²span','All9inducedG3actions'],fixed_ell_block_proof='E mapsFi→Qi, Q→0, Hj→delta_j0 F0+delta_j2 F2+sumSijQi. Rank5minor(rows3,5,6,7,8;cols0,2,3,4,5) islowertriangular withdiag1. E²=[Q0,0,Q2,0...],rank2. KerE=(K=H1−sumSi1Fi,Q0,Q1,Q2), withunit4projection coords1,6,7,8. ImageE²=(Q0,Q2). Genuinebasis[K,Q1],quotientcoords1,7.',genuine_canonical_source='AcceptedM/m³=M/ell³ andexponents1,2 identifykerE/imE² canonically withfullM[ell]. Neitherrawkernel4 norfree-blocktopisidentifiedwithgenuineG3.',target_action_shape='Sinceallq1gamma1-columnsare0, everyA_vK ispurequadratic andA_vQ1=0. CommutationmakesA_vK lieinkerE; quotientgivesK→alpha_vQ1,Q1→0,m²G3=0. Noalpha evaluated.',source_term_count=product_count,source_action_common_denominator=den,independent_integer_bound=coarse,accepted_integer_bound=bound,derived_arithmetic_bits_upper=181,static_numeric_envelope=18000,hard_numeric_cap=24000,worker_import_only_pass=True,actual_actions_kernel_or_G3_NOT_computed=True,no_new_source_gauge_jet_or_connecting_values=True)
 save('independent-cached-nine-module-and-genuine-G3-source-proof.json',proof)
 for p in[Path(__file__),worker,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 total=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir()if p.is_file());assert total<500000 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,cached9_action_source_and_genuine_G3_constructor_accepted=True,prepared_worker_and_fixed_unit_block_gate_passes=True,actual_action_or_G3_values_NOT_computed=True,static_numeric_envelope=18000,hard_numeric_cap=24000,derived_bit_bound=181,no_new_source_gauge_jet_or_connecting_values=True,parent_value_GO_required=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=total,seconds=time.monotonic()-st))
 print(json.dumps(dict(all_passed=True,actions_not_evaluated=True,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
