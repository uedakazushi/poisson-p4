#!/usr/bin/env python3
"""Independent finite Schur coefficient and GIVEN rank witness replay."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
import hashlib,json,resource,time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63';STEM=AUTHOR+'_REVIEW'
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
def check_given_rank(columns,certificate):
 basis=[sparse(p)for p in certificate['basis_normalized_columns']]
 factors=[sparse(p)for p in certificate['all_column_factorization_in_normalized_basis']]
 leads=[min(p)for p in basis];assert len(set(leads))==len(basis)
 assert leads==certificate['selected_minor_rows']
 assert all(p[k]==1 for p,k in zip(basis,leads))
 assert len(columns)==len(factors)==len(certificate['fixed_elimination_receipts'])
 for col,factor in zip(columns,factors):
  value={}
  for i,q in factor.items():value=plus(value,basis[i],q)
  assert value==col
 selected=[];count=0;det=F(1)
 for j,(col,r)in enumerate(zip(columns,certificate['fixed_elimination_receipts'])):
  assert r['column']==j;v=dict(col)
  for p,q in r['subtract']:
   q=F(q);assert p<count and v.get(leads[p],F(0))==q;v=plus(v,basis[p],-q)
  if r.get('zero'):assert not v
  else:
   k=r['pivot_row'];q=F(r['divide']);assert k==min(v)and v[k]==q and q
   assert{i:c/q for i,c in v.items()}==basis[count]==sparse(r['normalized'])
   assert all(i<=count for i in factors[j])and factors[j][count]==q
   selected.append(j);count+=1;det*=q
 assert selected==certificate['selected_original_columns']and count==certificate['rank']
 assert det==F(certificate['selected_minor_determinant'])and det
 return count

def main():
 t=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 rp='reports/'+AUTHOR+'/report.json';R=read(rp);assert R['all_passed'];hashes={rp:sha(rp)}
 for p,h in R['input_sha256'].items():assert sha(p)==h,p;hashes[p]=h
 for name,h in R['evidence_sha256'].items():
  p='reports/'+AUTHOR+'/'+name;assert sha(p)==h,p;hashes[p]=h
 for p,h in[('code/ambient_B18_comparison_T_linear_quadratic_20261104_v63.py',R['code_sha256']),('research/notes/'+AUTHOR+'.md',R['note_sha256'])]:assert sha(p)==h;hashes[p]=h
 peer='reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63_REVIEW/report.json';P=read(peer);assert P['all_passed']and P['actual_mu_Rbar_T_exact']==5;hashes[peer]=sha(peer)
 def get(st,f):return read('reports/'+st+'_20261104_v63/'+f)
 X0=get('AMBIENT_B18_COMPARISON_X_ORIGIN','actual-X0-and-unit-cancellation.json')
 A=get('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR','actual-first640-X-mod-m2.json')
 N=get('AMBIENT_B18_NATIVE128_FIRST','actual-native128-X-mod-m2.json')
 Mu=get('AMBIENT_B18_COMPARISON_T_CONSTANT_RANK','actual-linear-cycle-muT-certificate.json')
 C=get('AMBIENT_B18_CBAR_LITERAL_RESOLUTION','actual-Cbar-minimal-resolution-and-three-targets.json')
 S1=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC','phase1-actual-five-target-linear-relations.json')
 S2=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC','phase2-certified-quadratic-input.json')
 E=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC','actual-linear-quadratic-T-certificate.json')
 columns=[{i:sparse(p)for i,p in c}for c in A['source_columns640_rows56_sparse_dual9']+N['source_columns128_rows48_sparse_dual9']]
 relations=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations']if r.get('zero')]
 free32=X0['free_kernel_rows32'];basis27={r['pivot']:sparse(r['normalized'])for r in X0['all768_fixed_elimination_receipts_and_constant_combinations']if 'pivot'in r}
 tail=[i for i in range(32)if i not in basis27];assert len(tail)==5
 assert tail==S1['constant_projection_tail5']and[free32[i]for i in tail]==S1['five_remaining_free_source_rows']
 # Build the five LEFT annihilators by descending back-substitution, rather
 # than the producer's forward subtraction on every source column.
 projection=[]
 for row in tail:
  p={row:F(1)}
  for k in sorted(basis27,reverse=True):
   q=-sum((c*p.get(i,F(0))for i,c in basis27[k].items()if i!=k),F(0))
   if q:p[k]=q
  assert all(sum((c*p.get(i,F(0))for i,c in col.items()),F(0))==0 for col in basis27.values())
  assert[p.get(i,F(0))for i in tail]==[F(i==row)for i in tail]
  projection.append(p)
 Q1=[]
 for r in relations:
  xp={}
  for j,c in sparse(r['constant_original_source_combination']).items():
   for i,row in enumerate(free32):
    for k,a in columns[j].get(row,{}).items():
     if k:
      index=(i,k-1);v=xp.get(index,F(0))+c*a
      if v:xp[index]=v
      else:xp.pop(index,None)
  q={}
  for row,p in enumerate(projection):
   for(i,k),a in xp.items():
    c=p.get(i,F(0));index=8*row+k;v=q.get(index,F(0))+c*a
    if v:q[index]=v
    else:q.pop(index,None)
  Q1.append(q)
 assert Q1==[sparse(p)for p in S1['all741_Q1_columns_sparse']]
 wc=Mu['complete_rank_certificate'];wbasis=[sparse(p)for p in wc['basis_normalized_columns']];wfactors=[sparse(p)for p in wc['all_column_factorization_in_normalized_basis']]
 selected=wc['selected_original_columns'];norm=[]
 for k,j in enumerate(selected):
  value=dict(Q1[j])
  for i,c in wfactors[j].items():
   if i<k:value=plus(value,norm[i],-c)
   else:assert i==k
  diag=wfactors[j][k];assert diag;norm.append({i:c/diag for i,c in value.items()})
 assert norm==[sparse(p)for p in S1['normalized84_unitcycle_Q1_columns_sparse']]
 linear=[];labels=[]
 for j,q in enumerate(Q1):
  if j in selected:continue
  value=dict(q)
  for p,c in wfactors[j].items():value=plus(value,norm[p],-c)
  linear.append(value);labels.append(relations[j]['source_column'])
 assert linear==[sparse(p)for p in S1['actual40_by657_linear_relation_columns']]
 assert labels==S1['actual657_minimal_relation_source_labels']and len(linear)==657
 assert check_given_rank(linear,E['complete_linear_rank_certificate'])==37
 mons=list(combinations_with_replacement(range(8),2));assert[list(m)for m in mons]==S2['quadratic_monomials36'];index={m:i for i,m in enumerate(mons)}
 quadratic=[];label=[]
 for j in E['complete_linear_rank_certificate']['selected_original_columns']:
  for k in range(8):
   outcol={}
   for row,c in linear[j].items():
    target,var=divmod(row,8);i=36*target+index[tuple(sorted((var,k)))];outcol[i]=outcol.get(i,F(0))+c
   quadratic.append({i:c for i,c in outcol.items()if c});label.append({'kind':'variable_times_actual_linear_relation','source_relation':labels[j],'variable':k})
 beta_checks=0
 for d3j,col in enumerate(C['differentials8_column_sparse'][2]):
  beta=sparse(S2['all168_d3_linear_coordinates_in_unit_wedge_basis'][d3j])
  literal=[{}for _ in range(8)]
  for row,p in col:
   for m,c in p:
    if sum(m)==1:literal[m.index(1)][row]=F(c)
  reconstructed=[{}for _ in range(8)];outcol={}
  for i,c in beta.items():
   p,k=divmod(i,8);reconstructed[k]=plus(reconstructed[k],wbasis[p],c)
   for row,a in norm[p].items():
    target,var=divmod(row,8);q=36*target+index[tuple(sorted((var,k)))];outcol[q]=outcol.get(q,F(0))-c*a
  assert literal==reconstructed;beta_checks+=84*8
  quadratic.append({i:c for i,c in outcol.items()if c});label.append({'kind':'certified_actual_d3_quadratic','literal_d3_column':d3j})
 assert len(quadratic)==464 and label==S2['source_columns464']
 assert quadratic==[sparse(p)for p in S2['all464_certified_quadratic_columns_sparse']]
 assert check_given_rank(quadratic,E['quadratic_certificate'])==177
 assert E['actual_dim_mT_mod_m2T']==3 and E['actual_length_T_mod_m2T']==8
 assert not E['certified_quadratic_full180']and E['actual_total_length_T_proved']is None and E['actual_origin_cone_depth_proved']is None
 proof={'all_passed':True,'actual_constant_five_row_projection_independently_left_annihilator_rebuilt':True,
  'all741_Q1_and84_normalized_Q1_and657_linear_relations_checked':True,
  'independent_Schur_reconstruction':'paired triangular factor coefficients, not replay of producer Q1 elimination',
  'complete_given_linear_rank37_factorization_and_minor_verified':True,
  'actual_length_T_mod_m2_exact':8,'actual_HF_degrees0_1':[5,3],
  'all296_linear_variable_multiples_and168_d3_quadratics_source_bound':True,
  'd3_beta_identity_coefficient_slots_checked':beta_checks,
  'all464_known_quadratic_columns_match':True,
  'complete_given_known_quadratic_rank177_factorization_and_minor_verified':True,
  'known_quadratic_span_codimension3_NOT_actual_HF2':True,
  'no_m2T_zero_total_length_depth_conclusion':True,
  'no_new_source_X2_m1_jet_GB_348_or_rank_search':True}
 (out/'independent-five-row-Schur-and-fixed-rank-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,
  'actual_linear_relation_rank':37,'actual_length_T_mod_m2T':8,'certified_known_quadratic_rank':177,
  'actual_total_length_T_or_depth_claim':False,'input_sha256':hashes,
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
