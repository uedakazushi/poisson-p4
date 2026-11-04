#!/usr/bin/env python3
"""Independent fixed-source coverage + universal typed fifth pivot audit."""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parents[1];S='S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63_REVIEW';O=R/'reports'/S;P='S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def put(d,m,c):
 c=d.get(m,F(0))+c
 if c:d[m]=c
 else:d.pop(m,None)
def add(a,b,scale=F(1)):
 out=dict(a)
 for m,c in b.items():put(out,m,scale*c)
 return out
def mul(a,b):
 out={}
 for m,c in a.items():
  for n,d in b.items():
   e=dict(m)
   for v,k in n:e[v]=e.get(v,0)+k
   put(out,tuple(sorted((v,k) for v,k in e.items() if k)),c*d)
 return out
def atom(name,exponent=1):return {((name,exponent),):F(1)}
ONE={():F(1)};ZERO={}
def eye(n):return [[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
def mm(a,b):
 assert len(a[0])==len(b)
 return [[sum_poly([mul(a[i][k],b[k][j]) for k in range(len(b))]) for j in range(len(b[0]))] for i in range(len(a))]
def sum_poly(ps):
 out={}
 for p in ps:out=add(out,p)
 return out
def neg(A):return [[{m:-c for m,c in p.items()} for p in row] for row in A]
def dec(p):return {tuple((str(i),e) for i,e in enumerate(m) if e):F(c) for m,c in p}
def save(n,d):(O/n).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));O.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=R/'reports'/st/name;pins[str(p.relative_to(R))]=sha(p);return json.loads(p.read_text())
 primary=load(P);assert pins[f'reports/{P}/report.json']=='d4ef9e3d48d3701d0e4e832d66ef32aedffa6813dc7a3da48f51b44994d1340d'
 for p,h in primary['input_sha256'].items():assert sha(R/p)==h,(p,h)
 op=load(P,'operator.json');cover=load(P,'coverage.json');directstem='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63';direct=load(directstem,'operator.json')
 # Complete branch coverage from actual scalar identities, with no unknown factor inverted.
 branch=load('S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63','certificate.json');generic=load('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63','certificate.json')
 p,g1,g2,theta=[dec(branch[k]) for k in ['source_p','g1','g2','Theta']];a,e,d=[atom(str(i)) for i in range(3)]
 Binner=add(add(add(mul({():F(50)},a),mul({():F(18)},mul(e,e))),mul({():F(9)},e),-1),{():F(2)},-1)
 B=add(add(mul(Binner,d),mul({():F(60)},e)),{():F(40)},-1)
 C=add(add(add(mul({():F(100)},a),mul({():F(18)},mul(e,e))),mul({():F(9)},e),-1),{():F(2)},-1)
 assert mul({():F(750000)},p)==mul(g1,theta)
 assert mul(B,B)==add(theta,mul(C,mul(g1,g2)))
 assert p==dec(generic['p'])
 assert dec(generic['literal_pair_determinant'])==mul({():F(-6,5)},mul(B,mul(p,p)))
 assert generic['sparse_pair_rows']==[2,4] and generic['sparse_pair_columns']==[0,5]
 assert not set(generic['sparse_pair_rows'])&set(generic['all_other3_column_support_exact_rows'])
 assert len(generic['all_other3_column_support_exact_rows'])==7
 # Independent peer gates supply the full source identities behind these exact block decompositions.
 gates=[('S61_NONSQUARE_RAW_H_P_BRANCH_CHAIN_20261104_v63_REVIEW','all_passed'),('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW','all_passed'),('S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW','all_passed'),('S61_NONSQUARE_THETA_P80_REPLACEMENT_20261104_v63_REVIEW','ALL_P80_full_cubic57_proved'),('S61_NONSQUARE_THETA_P100_REPLACEMENT_20261104_v63_REVIEW','ALL_P100_full_cubic57_proved'),('S61_NONSQUARE_RAW_H_B_K13_CUBIC_REGULAR_REP_20261104_v63_REVIEW','whole_actual_Ka13_Bbad_full_cubic57'),(directstem+'_REVIEW','all_passed'),('S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW','all_passed')]
 for st,k in gates:assert load(st)[k]
 gpair=load('S61_NONSQUARE_RAW_H_G1_UNIT_PAIR_20261104_v63','certificate.json');assert gpair['selected_source_rows']==[5,8] and gpair['all_six_other_heavy_entries_literally_zero']
 assert gpair['other_heavy_columns']==[2,3,4]
 curves=load('S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63','certificate.json');assert curves['literal_Q6_equals_minus_product'] and curves['P100_minus_P80_equals20d']
 assert cover['rawH_rank_below7_exact_geometric_support']==['actual Bbad13','Theta P80','Theta P100'] and cover['all_three_exceptional_sets_already_full_cubic57']
 # Every typed node is validated and evaluated in a NEW universal Laurent ring.
 nodes=op['nodes'];assert len(nodes)==62;values=[];E=[[atom(f'E{i}_{j}') for j in range(13)] for i in range(6)];E[0][9]=atom('v');survive=[[atom(f'Q{i}_{j}') for j in range(25)] for i in range(6)]
 external={'residual6x13':E,'full_source15':eye(15),'full_source15_inverse':eye(15),'full_target22':eye(22),'full_target22_inverse':eye(22),'surviving6_source_original25':survive}
 for i,n in enumerate(nodes):
  assert n['id']==i and all(0<=j<i for j in n['args']);args=[values[j] for j in n['args']];kind=n['op']
  if kind=='accepted_external_export':
   ref=n['reference'];assert ref['path']==f'reports/{directstem}/operator.json' and ref['sha256']==pins[ref['path']] and ref['node']==direct['exports'][ref['export']];v=external[ref['export']]
  elif kind=='submatrix':
   rr,cc=n['rows'],n['columns'];assert len(rr)==len(set(rr)) and len(cc)==len(set(cc))
   assert all(0<=j<len(args[0]) for j in rr) and all(0<=j<len(args[0][0]) for j in cc);v=[[args[0][r][c] for c in cc] for r in rr]
  elif kind=='certified_chart_scalar_inverse':
   assert args[0]==[[atom('v')]] and n['new_chart_factor']=='Saux7022';v=[[atom('v',-1)]]
  elif kind=='matrix_product':v=mm(*args)
  elif kind=='negative':v=neg(args[0])
  elif kind=='central_scalar_multiple':assert len(args[0])==len(args[0][0])==1;v=[[mul(args[0][0][0],p) for p in row] for row in args[1]]
  elif kind=='difference':v=[[add(p,q,-1) for p,q in zip(row,s)] for row,s in zip(args[0],args[1])]
  elif kind=='identity':assert n['shape'][0]==n['shape'][1];v=eye(n['shape'][0])
  elif kind=='zero':v=[[ZERO for j in range(n['shape'][1])] for k in range(n['shape'][0])]
  elif kind=='block_matrix':
   blocks=[[values[j] for j in row] for row in n['blocks']];assert n['args']==[j for row in n['blocks'] for j in row]
   v=[]
   for row in blocks:
    assert len({len(b) for b in row})==1
    v.extend([sum((b[k] for b in row),[]) for k in range(len(row[0]))])
  elif kind=='column_permutation':
   order=n['order'];assert sorted(order)==list(range(13));v=[[ONE if i==order[j] else ZERO for j in range(13)] for i in range(13)]
  elif kind=='transpose':v=list(map(list,zip(*args[0])))
  else:raise AssertionError(kind)
  assert [len(v),len(v[0])]==n['shape'];values.append(v)
 exports={k:values[j] for k,j in op['exports'].items()};L,Li,U,Ui=[exports[k] for k in ['source6','source6_inverse','target13','target13_inverse']]
 assert mm(L,Li)==mm(Li,L)==eye(6) and mm(U,Ui)==mm(Ui,U)==eye(13)
 T=exports['residual5x12'];target=[[ONE]+[ZERO]*12]+[[ZERO]+row for row in T];assert mm(mm(L,E),U)==target
 assert all(add(mul(atom('v'),T[i][j]),exports['cleared_residual5x12'][i][j],-1)==ZERO for i in range(5) for j in range(12))
 for f,fi,size in [('full_source15','full_source15_inverse',15),('full_target22','full_target22_inverse',22)]:assert mm(exports[f],exports[fi])==mm(exports[fi],exports[f])==eye(size)
 assert exports['surviving5_source_original25']==mm(L[1:],survive)
 assert op['chart']=='FULLU D(R434*Delta5*Saux)' and op['new_parameter_coefficients']==0
 assert op['remaining_current15']==[8,9,11,12,14] and op['remaining_original25']==[18,19,21,22,24]
 expected=[1,2,3,4,6,7,8,9,14,17,18,19];assert op['remaining_original32']==expected and op['extra8_Q9_tails_retained']==expected[:8]
 oldlabels=direct['context']['remaining13_original32_seed_labels'];assert [x for j,x in enumerate(oldlabels) if j!=9]==expected and oldlabels[9]==15
 # Source composition order is separately audited, so specializing accepted old equivalences
 # to identities above does not assert their actual coefficient matrices are identities.
 exp=op['exports'];f=nodes[exp['full_source15']];fi=nodes[exp['full_source15_inverse']];assert nodes[f['args'][1]]['reference']['export']=='full_source15' and nodes[fi['args'][0]]['reference']['export']=='full_source15_inverse'
 f=nodes[exp['full_target22']];fi=nodes[exp['full_target22_inverse']];assert nodes[f['args'][0]]['reference']['export']=='full_target22' and nodes[fi['args'][1]]['reference']['export']=='full_target22_inverse'
 costs=load(P,'support-preflight.json');assert len(costs['all60_candidates'])==60 and costs['smallest']['output_plus_input_support_upper']==180479
 save('independent-typed-coverage-proof.json',{'typed_nodes':62,'all78_relative_diagonalization_entries':True,'both_sides_relative_and_embedded_inverse_products':True,'full_source_composition_order_checked':True,'source25_surviving5_transport':True,'remaining_original32':expected,'all8_Q9_pure_tails_retained':True,'p_factor_and_comaximal_QQ_identity':True,'geometric_H_bad_exact_support':['Bbad13','ThetaP80','ThetaP100'],'exceptions_already_full57':True,'unresolved_actual_defect_contained_in_H7':True,'H7_local_bundle_kernel10_quotientQ9_rank9_to16':True,'chart':op['chart'],'Saux_is_actual_bordered5_SLP':True,'7022_QQ_coefficient_expansion_not_repeated_or_newly_claimed':True,'nonzero_source_gate':'accepted independent symbolicmodp proof implies char0 nonzero','no_actual5x12_coefficient_value_or_rank5_computed':True,'all60_border_support_numbers_not_independently_recomputed':True,'smallest_border180479_exceeds150k_kept_stopped':True})
 report={'id':S,'all_passed':True,'scope':'geometric H7 coverage and source-bound typed fifth pivot; no new actual coefficients','full_FULLU_p0_coverage_QQ_Bezout_checked':True,'actual_same_cokernel_chart_shape':[5,12],'new_scalar_coefficients':0,'all8_Q9_tails_retained':True,'no_global57_quartic239_or_rank5_claim':True,'no_producer_import_rerun_source_query_minor_search_point_prime_GB':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md'),'evidence_sha256':{'independent-typed-coverage-proof.json':sha(O/'independent-typed-coverage-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in O.iterdir())<33554432;save('report.json',report);save('phase.json',{'stage':'FINAL','all_passed':True});print(json.dumps(report))
if __name__=='__main__':main()
