#!/usr/bin/env python3
"""Independent anchor reconstruction and GIVEN rank certificate for mu(T)."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,resource,time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63';STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def sparse(raw):return{i:F(q)for i,q in raw if F(q)}
def plus(a,b,c=1):
 a=dict(a)
 for i,q in b.items():
  x=a.get(i,F(0))+c*q
  if x:a[i]=x
  else:a.pop(i,None)
 return a

def main():
 t=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 rp='reports/'+AUTHOR+'/report.json';R=read(rp);assert R['all_passed'];hashes={rp:sha(rp)}
 for p,h in R['input_sha256'].items():assert sha(p)==h,p;hashes[p]=h
 for name,h in R['evidence_sha256'].items():
  p='reports/'+AUTHOR+'/'+name;assert sha(p)==h,p;hashes[p]=h
 for p,h in[('code/ambient_B18_comparison_T_constant_rank_20261104_v63.py',R['code_sha256']),('research/notes/'+AUTHOR+'.md',R['note_sha256'])]:assert sha(p)==h;hashes[p]=h
 for st in['AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63_REVIEW','AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW','AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63_REVIEW','AMBIENT_B18_NATIVE128_FIRST_20261104_v63_REVIEW']:
  p='reports/'+st+'/report.json';q=read(p);assert q['all_passed'];hashes[p]=sha(p)
  for key in['input_sha256','checked_dependency_sha256']:
   for p,h in q.get(key,{}).items():
    if p in R['input_sha256']:assert R['input_sha256'][p]==h
 def get(st,name):return read('reports/'+st+'_20261104_v63/'+name)
 X0=get('AMBIENT_B18_COMPARISON_X_ORIGIN','actual-X0-and-unit-cancellation.json')
 A=get('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR','actual-first640-X-mod-m2.json')
 N=get('AMBIENT_B18_NATIVE128_FIRST','actual-native128-X-mod-m2.json')
 C=get('AMBIENT_B18_CBAR_LITERAL_RESOLUTION','actual-Cbar-minimal-resolution-and-three-targets.json')
 V=get('AMBIENT_B18_FULLV56_COMPARISON_UNIT','actual-fullV56-m0-and-typed-comparison.json')
 E=get('AMBIENT_B18_COMPARISON_T_CONSTANT_RANK','actual-linear-cycle-muT-certificate.json')
 anchor_saved=get('AMBIENT_B18_COMPARISON_T_CONSTANT_RANK','literal-d2-linear-anchors.json')
 columns=A['source_columns640_rows56_sparse_dual9']+N['source_columns128_rows48_sparse_dual9'];assert len(columns)==768
 X=[]
 for j,column in enumerate(columns):
  col={row:sparse(p)for row,p in column};assert len(col)==len(column)
  assert all(row<(56 if j<640 else 48)for row in col)
  for i in range(56):assert col.get(i,{}).get(0,F(0))==F(X0['full_X0_56_by768'][i][j])
  X.append(col)
 relations=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations']if r.get('zero')]
 assert len(relations)==741 and X0['rank_X0']==27
 assert[r['source_column']for r in relations]==E['actual_remaining_source_column_labels']
 assert[r['constant_original_source_combination']for r in relations]==E['actual_constant_source_combinations']
 basis=C['resolution_basis_sector_and_wedge'][2];assert len(basis)==84
 linear=[]
 for col in C['differentials8_column_sparse'][1]:
  result={}
  for row,p in col:
   for mon,q in p:
    if sum(mon)==1:result[8*row+mon.index(1)]=F(q)
  linear.append(result)
 anchors=[8*(8*s+j)+i for s,(i,j)in basis]
 assert len(set(anchors))==84
 for i,k in enumerate(anchors):assert[poly.get(k,F(0))for poly in linear]==[F(i==j)for j in range(84)]
 assert anchors==anchor_saved['selected_derivative_coordinates84']
 assert linear==[sparse(x)for x in E['literal_d2_linear_columns84']]
 assert basis==E['literal_wedge_basis84']
 (out/'preflight.json').write_text(json.dumps({'source_columns':768,'constant_null_source_combinations':741,'literal_d2_anchors':84,
  'complete_linear_constraint_slots':741*24*8,'no_rank_search_or_source_replay':True,'input_sha256':hashes},indent=2)+'\n')
 m0=[[F(q)for q in row]for row in V['all24_by56_m0']]
 column_m0=[[(row,m0[row][i])for row in range(24)if m0[row][i]]for i in range(56)]
 wcols=[];constraints=0;w_nonzero=0
 for j,relation in enumerate(relations):
  xp={}
  for k,c in sparse(relation['constant_original_source_combination']).items():
   for row,dual in X[k].items():
    for axis,q in dual.items():
     index=(row,axis);value=xp.get(index,F(0))+c*q
     if value:xp[index]=value
     else:xp.pop(index,None)
  assert all(axis for row,axis in xp)
  cycle={}
  for(i,axis),q in xp.items():
   for row,c in column_m0[i]:
    key=8*row+axis-1;value=cycle.get(key,F(0))+c*q
    if value:cycle[key]=value
    else:cycle.pop(key,None)
  assert cycle==sparse(E['actual24_by8_linear_cycles_sparse'][j])
  # Independent extraction by the literal identity minor, rather than the
  # author's skew-difference/2 formula. Then check ALL192 coefficients.
  w={i:cycle[k]for i,k in enumerate(anchors)if cycle.get(k,F(0))}
  rebuilt={}
  for i,q in w.items():rebuilt=plus(rebuilt,linear[i],q)
  assert rebuilt==cycle
  assert w==sparse(E['actual84_by741_w0_columns_sparse'][j])
  wcols.append(w);w_nonzero+=len(w);constraints+=24*8
 assert constraints==142272 and w_nonzero==715
 rank=E['complete_rank_certificate'];normalized=[sparse(p)for p in rank['basis_normalized_columns']]
 factors=[sparse(p)for p in rank['all_column_factorization_in_normalized_basis']]
 assert len(normalized)==84 and len(factors)==741 and rank['rank']==84
 leads=[min(p)for p in normalized];assert len(set(leads))==84 and set(leads)==set(range(84))
 assert all(p[k]==1 for p,k in zip(normalized,leads))
 chosen=rank['selected_original_columns'];minor_rows=rank['selected_minor_rows'];assert minor_rows==leads and len(chosen)==84
 for col,factor in zip(wcols,factors):
  value={}
  for p,q in factor.items():value=plus(value,normalized[p],q)
  assert value==col
 receipts=rank['fixed_elimination_receipts'];atrow={};basis_count=0;selected=[];det=F(1)
 for j,(col,record)in enumerate(zip(wcols,receipts)):
  assert record['column']==j;v=dict(col)
  for p,c in record['subtract']:
   c=F(c);assert p<basis_count and v.get(leads[p],F(0))==c
   v=plus(v,normalized[p],-c)
  if record.get('zero'):assert not v
  else:
   row=record['pivot_row'];q=F(record['divide']);assert row==min(v)and q==v[row]and row not in atrow
   assert{i:x/q for i,x in v.items()}==normalized[basis_count]
   assert sparse(record['normalized'])==normalized[basis_count]
   assert all(i<=basis_count for i in factors[j])and factors[j][basis_count]==q
   det*=q;selected.append(j);atrow[row]=basis_count;basis_count+=1
 assert basis_count==84 and selected==chosen and det==F(rank['selected_minor_determinant'])and det
 # Lower minor proof with no new determinant computation: normalized columns
 # sorted by distinct leading row give a triangular unit diagonal matrix.
 # Sorting rows and columns together has determinant sign+1. Selected-column
 # factor matrix is triangular with the recorded divides on its diagonal.
 for i,col in enumerate(normalized):assert all(row>=leads[i]for row in col)
 assert E['mu_Rbar_T_exact']==89-84==5
 proof={'all_passed':True,'all768_source_origins_and27_constant_cancel_combinations_checked':True,
  'all741_Xprime_zero_constants_and_m0_Xprime_linear_checked':True,
  'all142272_actual_d2_linear_cycle_constraints_checked':constraints,
  'all84_literal_identity_anchor_rows_checked':True,
  'independent_wedge_extraction':'literal d2 anchor identity minor, not skew formula',
  'actual_W_shape':[84,741],'actual_W_nonzeros':w_nonzero,
  'all_given_rank_receipts_and_all_column_factorizations_verified':True,
  'fixed84_minor_nonzero_determinant':str(det),'rank_W_exact':84,'mu_Rbar_T_exact':5,
  'unknown_m1_cancels_by_Xprime0_zero':True,
  'unknown_first27_top_wedge_constants_do_not_change_remaining_wedge_constants':True,
  'length_or_depth_from_mu5_claim':False,
  'next_Schur_theory':{
   'Y0_zero_removes_W1_from_linear_R':True,
   'linear_formula':'R1=Y_rest1-Y_S1*W_S0^-1*W_rest0',
   'd3_quadratic_formula':'-Y_S1*W_S0^-1*d3_linear',
   'unknown_first27_top_constants_irrelevant':'all remaining Xprime0 vanish and the bottom5 projection annihilates selected27 constants; corrections start too high for Y1/W0',
   'rank_linear_at_most36_implies_length_at_least9':True,
   'if_rank37_and_known_quadratics_span180':'m²T=m³T; finite-module Nakayama gives m²T=0 and length5+3=8',
   'partial_quadratic_span_does_not_prove_nonzero_higher_layer':True,
   'no_actual_next_rank_or_quadratic_computation_here':True}}
 (out/'independent-constant-cancel-wedge-rank-and-Schur-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,'actual_mu_Rbar_T_exact':5,'actual_W_rank':84,
  'all741_cycles142272_constraints_and_fixed_rank_factorization_verified':True,
  'next_linear_and_d3_quadratic_sufficiency_theory_accepted':True,
  'actual_T_length_K_torsion_or_depth_computed':False,'input_sha256':hashes,
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
