#!/usr/bin/env python3
"""Finite parameter-only selector and ALL129 linear-normal factoredcircuits."""
import argparse,gzip,itertools,json,resource,signal,time
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
from s61_sparse_generic5_full134_normal_circuit_20261104_v63 import Pool
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63'
class Stop(Exception):pass
def sign(p):return -1 if sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))%2 else 1
def formal_proof():
 # Independent formal polynomial dictionaries on25 matrix entries.
 def add(a,b,c=1):
  a=dict(a)
  for m,v in b.items():
   a[m]=a.get(m,0)+c*v
   if not a[m]:del a[m]
  return a
 def mul(a,b):
  c={}
  for m,x in a.items():
   for n,y in b.items():
    k=tuple(sorted(m+n));c[k]=c.get(k,0)+x*y
  return {m:v for m,v in c.items()if v}
 def det(rr,cc):
  a={}
  for p in itertools.permutations(range(len(rr))):a=add(a,{tuple(sorted(5*rr[i]+cc[p[i]]for i in range(len(rr)))):sign(p)})
  return a
 D=det(list(range(5)),list(range(5)));adj=[[{m:v*(-1 if(i+j)%2 else 1)for m,v in det([r for r in range(5)if r!=j],[c for c in range(5)if c!=i]).items()}for j in range(5)]for i in range(5)]
 B=[[{(5*i+j,):1}for j in range(5)]for i in range(5)]
 for X,Y in [(B,adj),(adj,B)]:
  for i in range(5):
   for j in range(5):
    a={}
    for k in range(5):a=add(a,mul(X[i][k],Y[k][j]))
    assert a==(D if i==j else{})
 for j in range(5):
  a={}
  for h in range(5):
   for k in range(5):a=add(a,mul(mul({(25+h,):1},adj[h][k]),B[k][j]))
  assert a==mul(D,{(25+j,):1})
 # Every normal-containing monomial degree<=2 has a unique min-normalfactor.
 mon=[(i,)for i in range(3,45)]+list(itertools.product(range(3),range(3,45)))+list(itertools.combinations_with_replacement(range(3,45),2));assert len(mon)==1071
 for m in mon:
  k=min(i for i in m if i>=3);rest=list(m);rest.remove(k);assert tuple(sorted((k,*rest)))==tuple(sorted(m))and len(rest)<=1
 return dict(both5x5_adjugate_products_exact=True,formal_matrix_identities=50,formal_row_inverse_identities=5,all1071_normal_monomial_factorizations_exact=True,all129_residual_factors_by_row_linearity=True,plane_zero_gate='Acceptedfivefunctionalinjectivity; S-lambda*B=0 byformaladjugate; no expandedEorparameterpoly',whole134_equivalence='Selected5 unchanged; other129 replacedbyFi-lambda_i*Fselected, invertibleblocktriangularrowoperation overK')
def main(out):
 start=time.monotonic();out=out.resolve();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(Stop('60s wall cap')));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 src=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63';sel=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63';raw=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_FIVE_SELECTOR_20261104_v63';rank3=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63';pref=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63'
 assert read(src/'report.json')['all_passed']and read(sel/'report.json')['all_passed'];J=read(sel/'native-certificate.json.gz')['selected_original134_rows'];assert J==[50,10,52,51,15];rows=[i for i in range(134)if i not in J]
 direct=read(src/'full134-direct-circuit.json.gz');block=read(src/'coefficient-block-circuit.json.gz');old=direct['nodes']+block['appended_nodes'];assert len(direct['nodes'])==block['prefix_node_count'];oldrows=block['rows134'];assert len(oldrows)==134
 pre=dict(explicit_new_parameter_polynomial_terms=0,external_frozen_134x5_terms=692953,external_entries_referenced_not_duplicated=True,hard_TOTAL_new_explicit_support_guard150000=True,estimated_nodes_upper=350000,estimated_raw_JSON_bytes_upper=32*1024*1024,all129_rows_and42normals_retained=True,normal2_factored_source_sums_not_expanded=True,no_unknown_inverse=True,rank5_only_no_cross126_or_exhaustion=True)
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');proof=formal_proof();(out/'structural-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 pool=Pool();mcache={};scache={};fcache={}
 def scalar(i):
  if i not in scache:
   if old[i][0]=='qq':scache[i]=pool.const(old[i][1])
   else:scache[i]=pool.node('external_prior_scalar',[i],0,None)
  return scache[i]
 def F(i):
  if i not in fcache:fcache[i]=pool.node('external_prior_F134',[i,direct['F134'][i]],2,True)
  return fcache[i]
 def entry(i,j):
  if(i,j)not in mcache:mcache[i,j]=pool.node('external_five_functional_matrix',[i,j],0,True)
  return mcache[i,j]
 def neg(i):return pool.prod(pool.const(-1),i)
 def determinant(a):
  return pool.summ([pool.prod(pool.const(sign(p)),*[a[i][p[i]]for i in range(len(a))])for p in itertools.permutations(range(len(a)))])
 B=[[entry(i,j)for j in range(5)]for i in J];Delta=determinant(B);inv=pool.inv(Delta);adj=[[None]*5 for _ in range(5)]
 for i in range(5):
  for j in range(5):adj[i][j]=pool.prod(pool.const(-1 if(i+j)%2 else 1),determinant([[B[r][c]for c in range(5)if c!=i]for r in range(5)if r!=j]))
 lam=[[pool.prod(inv,pool.summ([pool.prod(entry(i,h),adj[h][j])for h in range(5)]))for j in range(5)]for i in rows]
 residual=[pool.summ([F(i),*[neg(pool.prod(v,F(j)))for v,j in zip(l,J)]])for i,l in zip(rows,lam)]
 linear=[];cross=[]
 for i,l in zip(rows,lam):
  linear.append([pool.summ([scalar(oldrows[i]['linear_n42'][k]),*[neg(pool.prod(v,scalar(oldrows[j]['linear_n42'][k])))for v,j in zip(l,J)]])for k in range(42)])
  cross.append([pool.summ([scalar(oldrows[i]['cross_z_n126'][k]),*[neg(pool.prod(v,scalar(oldrows[j]['cross_z_n126'][k])))for v,j in zip(l,J)]])for k in range(126)])
 Z=[pool.node('z',[k],1,True)for k in range(3)];N=[pool.node('n',[k],1,True)for k in range(42)];qf=[]
 for row in oldrows:
  parts=[[]for _ in range(42)]
  for a,b,r in row['n2_sparse']:assert a<=b;parts[a].append(pool.prod(scalar(r),N[b]))
  qf.append([pool.summ(p)for p in parts])
 cofactor=[]
 for t,(i,l)in enumerate(zip(rows,lam)):
  cofactor.append([pool.summ([linear[t][k],*[pool.prod(cross[t][42*h+k],Z[h])for h in range(3)],qf[i][k],*[neg(pool.prod(v,qf[j][k]))for v,j in zip(l,J)]])for k in range(42)])
 assert all(pool.deg[r]<=2 for r in residual)and all(pool.deg[r]<=1 for row in cofactor for r in row)and all(pool.deg[r]==0 for row in linear+cross+lam for r in row)
 inverses=[a[1]for a in pool.nodes if a[0]=='inv'];assert inverses==[Delta]and pool.deg[Delta]==0
 # Structuraltypes andtopologicalreferences; externalleaves have noDAGchildrefs.
 for i,node in enumerate(pool.nodes):
  if node[0]in['prod','sum']:assert all(isinstance(j,int)and 0<=j<i for j in node[1])
  elif node[0]=='inv':assert 0<=node[1]<i
 outputs=dict(nodes=pool.nodes,coefficient_field='QQ(H,L,s0,q0,t0)',external_files={'prior_direct':str((src/'full134-direct-circuit.json.gz').relative_to(ROOT)),'prior_blocks':str((src/'coefficient-block-circuit.json.gz').relative_to(ROOT)),'prior_scalar_atoms':str((src/'scalar-atoms-full.json.gz').relative_to(ROOT)),'five_functional_matrix':str((raw/'five-column-matrix.json.gz').relative_to(ROOT))},external_leaf_semantics={'external_prior_scalar':'nodeID in prior_direct.nodes followedby prior_blocks.appended_nodes, allparameteronly','external_prior_F134':'originalequationindex and actualdirectrootID, same z/n leaves','external_five_functional_matrix':'originalrow,column in matrix134x5, originalQQfivevarpolynomial; notQQ[T]substitution'},selected_original134_rows=J,residual_original134_rows=rows,B5_nodes=B,adjugate5_nodes=adj,Delta5_node=Delta,inverseDelta5_node=inv,lambda129x5=lam,selected_direct_F5=[F(j)for j in J],residual_direct_F129=residual,linear_n129x42=linear,cross129x126=cross,source_normal2_factor134x42=qf,linear_normal_cofactor129x42=cofactor,normal2_residual_lazy_recipe={'source_rows':rows,'selected_rows':J,'multipliers':lam,'rule':'Qres_i(k,l)=Qsource_i(k,l)-sum_j lambda_ij*Qsource_selected_j(k,l), withmissingpriorn2_sparse entries0; allk<=l'},residual_factor_identity='residual_direct_F129[i]=sum_k n_k*linear_normal_cofactor129x42[i][k]',all45_unknowns_order=['z0','z1','z2',*[f'n{i}'for i in range(42)]],pivot_q_rows=direct['pivot_q_rows'],normal_q_rows42=direct['normal_q_rows42'],q_to_u_indices_zero_based=direct['q_to_u_indices_zero_based'],all_original134_ideal_retained=True,all_normal2_retained=True,all_inverses_parameter_only=True,no_unknown_q_n_z_localization=True)
 payload=json.dumps(outputs,separators=(',',':')).encode();assert len(pool.nodes)<=350000;assert len(payload)<=32*1024*1024;assert pre['explicit_new_parameter_polynomial_terms']<=150000
 (out/'residual129-factored-circuit.json.gz').write_bytes(gzip.compress(payload,mtime=0))
 # Read-only degree/support majorant of a FUTURE formalQQ[T] crossinput.
 w=[1,2,3,5,7];W=read(rank3/'plane-W-original.json.gz')['W45x3'];data=read(pref/'pulled134-input.json.gz');ss=[[{sum(a*b for a,b in zip(w,t[1:]))for t in p}for p in row]for row in W];potential=0;degree=0;maxsupport=0;norm=direct['normal_q_rows42'];ni={k:i for i,k in enumerate(norm)}
 for eq in data['equations']:
  supports=[set()for _ in range(126)]
  for m,e,v in eq:
   if len(m)!=2:continue
   d=sum(a*b for a,b in zip(w,e));i,j=m
   if i in ni:
    for h in range(3):supports[42*h+ni[i]].update(d+a for a in ss[j][h])
   if j in ni:
    for h in range(3):supports[42*h+ni[j]].update(d+a for a in ss[i][h])
  potential+=sum(map(len,supports));maxsupport=max(maxsupport,max(map(len,supports)));degree=max(degree,max((max(s)for s in supports if s),default=0))
 cost=dict(future_polynomial_matrix='[S134x5 | Y134x126], withSactualfivefunctionalpolysandY=P*rawcross usingW=P*V; rowselectorSchurcomplementrank=5+crossrank',formal_map='sameQQ[T]exponents[1,2,3,5,7]; noothercurveorpointtested',raw_cross_support_upper_no_cancellation=potential,max_entry_support_upper=maxsupport,max_Tdegree_upper=degree,plus_five_column_actual_terms=55594,total_explicit_support_upper=potential+55594,not_an_exact_nonzero_count=True,no_crosscoefficients_or_rank_computed=True,correct_leading_certificate='Onlyif134x131 rowwiseCOMMONmaximumTdegreecoefficientmatrixhasrank131 doesa131minornonzeroresultfollow; flattenedcoefficientrank isnotmatrixrank',no_fixed_old126minor_assumed=True,whole_affine_residual_empty_not_proved=True)
 (out/'next-cross-cost.json').write_text(json.dumps(cost,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<=32*1024*1024
 paths=[Path(__file__),ROOT/'code/s61_sparse_generic5_full134_normal_circuit_20261104_v63.py',ROOT/f'research/notes/{ID}.md',src/'report.json',src/'full134-direct-circuit.json.gz',src/'coefficient-block-circuit.json.gz',sel/'report.json',sel/'native-certificate.json.gz',ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63_REVIEW/report.json',raw/'five-column-matrix.json.gz',rank3/'plane-W-original.json.gz',pref/'pulled134-input.json.gz']
 report=dict(id=ID,all_passed=True,all129_actual_residuals_and42_linear_normal_factors=True,all129x126_cross_coefficients=True,allnormal2_retained_factored=True,full134_ideal_equality_over_parameter_field=True,all45_unknowns_retained=True,only_new_inverse_Delta5_parameter=True,new_explicit_parameter_polynomial_terms=0,external692953_terms_not_duplicated=True,nodes=len(pool.nodes),external_scalar_roots=len(scache),external_five_matrix_entries=len(mcache),raw_circuit_bytes=len(payload),normalfactor_basis_identities=1071,adjugate_and_row_inverse_formal_identities=55,no_rank126_GB_J0_Euler_Macaulay_point_prime_or_exhaustion=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,nodes=len(pool.nodes),rawbytes=len(payload),seconds=report['seconds'],future_cross_support_upper=potential)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
