#!/usr/bin/env python3
"""Exact QQ-constant combinations divisible by one positive local variable.

All135 columns and all parameter coefficients are retained. No regeneration,
point, field, Groebner basis, normal graph or pivot is computed.
"""
import argparse,collections,hashlib,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_BIANCHI_CONSTANT_COMBINATION_DIVISION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def padd(a,b):
 z=dict(a)
 for m,c in b.items():z[m]=z.get(m,Q())+c
 return {m:c for m,c in z.items() if c}
def pmul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():
   k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,Q())+c*d
 return {m:c for m,c in z.items() if c}
def decode(ts):
 z={}
 for m,c in ts:
  assert len(m)==3 and min(m)>=0;z=padd(z,{tuple(m):Q(c)})
 return z
def enc(p):return [[list(m),str(c)] for m,c in sorted(p.items())]
def vector_add(target,p,scale=Q(1)):
 for key,c in p.items():
  z=target.get(key,Q())+scale*c
  if z:target[key]=z
  else:target.pop(key,None)
def kernel(constraints,n):
 """Selected exact constraint rows plus a full free-coordinate kernel basis."""
 basis={};selected=[]
 for label,row in sorted(constraints.items()):
  r=dict(row)
  while r:
   pivot=min(r)
   if pivot not in basis:
    lead=r[pivot];basis[pivot]={j:c/lead for j,c in r.items()};selected.append(label);break
   vector_add(r,basis[pivot],-r[pivot])
  if len(basis)==n:break
 free=[j for j in range(n) if j not in basis];vectors=[]
 for f in free:
  v={f:Q(1)}
  for pivot,row in sorted(basis.items(),reverse=True):
   value=-sum((c*v.get(j,Q()) for j,c in row.items() if j!=pivot),Q())
   if value:v[pivot]=value
  assert all(not sum((c*v.get(j,Q())for j,c in row.items()),Q())for row in constraints.values())
  vectors.append(v)
 assert len(vectors)==n-len(basis)
 return {'rank':len(basis),'selected_constraint_keys':[list(k)for k in selected],'free_row_indices':free,'kernel_vectors':[[str(v.get(j,Q()))for j in range(n)]for v in vectors]},vectors

def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s preflight cap')));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False)
 source=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63/operator.json';parent=source.parent/'report.json';pr=json.loads(parent.read_text());assert pr['all_full_raw_map_and_literal_identity_checks_passed']and sha(source)==pr['operator_sha256']
 for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 data=json.loads(source.read_text());weights=data['local_coordinate_weights'];ew=data['equation_weights'];fixed=data['local_fixed_q_index'];columns=data['local_other_q_order'];cp=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=json.loads(cp.read_text());assert [i for i,w in enumerate(weights)if w==0]==[19,23]
 for i in [19,23]:
  assert all(Q(chart['T34'][i][j])==0 for j in chart['normal34']);assert all(Q(x)==0 for x in chart['U34'][i])
 F=[]
 for g,poly in enumerate(data['full_homogeneous135_equations']):
  row={}
  for mon,ts in poly:
   local=tuple(sorted(44 if v==45 else columns.index(v)for v in mon if v!=fixed));assert sum(weights[v]for v in local)==ew[g]
   if any(v in [19,23]for v in local):continue
   row[local]=padd(row.get(local,{}),decode(ts))
  F.append({m:c for m,c in row.items()if c})
 assert len(F)==135
 groups={w:[]for w in [5,6,7,8]}
 for index,rec in enumerate(data['local_weight_pieces']):
  if rec['weight'] not in groups:continue
  vector={}
  for g,var,ts in rec['coefficients']:
   assert 0<=g<135 and -1<=var<45 and ew[g]+(0 if var==-1 else weights[var])==rec['weight']
   if var in [19,23]:continue
   for ex,c in decode(ts).items():vector_add(vector,{(g,var,*ex):c})
  groups[rec['weight']].append((index,vector))
 assert [len(groups[w])for w in [5,6,7,8]]==[20,22,18,11]
 candidates=[];results=[];tested=0
 for w,rows in groups.items():
  n=len(rows);full={}
  for j,(_,vector)in enumerate(rows):
   for key,c in vector.items():full.setdefault(key,{})[j]=c
  fullcert,fullker=kernel(full,n);systems=[]
  for variable,wt in enumerate(weights):
   if wt!=w-4:continue
   bad={key:row for key,row in full.items()if key[1]!=variable};cert,ks=kernel(bad,n);expected=fullcert['rank']-cert['rank'];assert expected>=0
   image_basis={};chosen=[]
   for alpha in ks:
    combined={}
    for j,coef in alpha.items():vector_add(combined,rows[j][1],coef)
    assert all(key[1]==variable for key in combined)
    reduced=dict(combined)
    while reduced:
     pivot=min(reduced)
     if pivot not in image_basis:
      lead=reduced[pivot];image_basis[pivot]={k:c/lead for k,c in reduced.items()};chosen.append((alpha,combined));break
     vector_add(reduced,image_basis[pivot],-reduced[pivot])
   assert len(chosen)==expected
   indices=[]
   for alpha,combined in chosen:
    divided={}
    for(g,var,a,e,d),c in combined.items():assert var==variable and ew[g]==4;divided[g]=padd(divided.get(g,{}),{(a,e,d):c})
    assert divided
    check={}
    for g,c in divided.items():
     for mon,p in F[g].items():check[mon]=padd(check.get(mon,{}),pmul(c,p))
    assert not {m:c for m,c in check.items()if c}
    indices.append(len(candidates));candidates.append({'source_weight':w,'factor_local_variable':variable,'factor_weight':wt,'quotient_weight':4,'source_weight_piece_indices':[i for i,_ in rows],'constant_QQ_combination':[str(alpha.get(j,Q()))for j in range(n)],'all135_columns_retained_in_divisibility_constraints':True,'divided_coefficients':[[g,enc(c)]for g,c in sorted(divided.items())],'literal_divided_full135_identity_verified':True,'normal_graph_or319_transport_not_performed':True})
   systems.append({'factor_local_variable':variable,'factor_weight':wt,'constraint_count':len(bad),'exact_constraint_kernel':cert,'identically_zero_full_combination_dimension':len(fullker),'quotient_nonzero_combination_dimension':expected,'candidate_indices':indices});tested+=1
  results.append({'total_weight':w,'source_piece_indices':[i for i,_ in rows],'number_of_unknown_constant_QQ_row_coefficients':n,'full_coefficient_map_constraint_count':len(full),'full_map_exact_kernel':fullcert,'variable_divisibility_systems':systems})
  (out/'progress.json').write_text(json.dumps({'completed_weights':[r['total_weight']for r in results],'completed_variable_systems':tested,'nonzero_candidates':len(candidates),'seconds':time.monotonic()-start},indent=2)+'\n');print('GROUP',w,'FULLRANK',fullcert['rank'],'SYSTEMS',len(systems),'CANDIDATES',sum(s['quotient_nonzero_combination_dimension']for s in systems),'SEC',time.monotonic()-start,flush=True)
 evidence={'source_operator_sha256':sha(source),'verified_zero_nuisance_coordinates':[19,23],'coefficient_keys':'(equation_index, local_multiplier_variable_or_minus1, a_exponent,e_exponent,d_exponent)','groups':results,'candidates':candidates,'total_variable_systems':tested,'all_parameter_coefficients_split_exactly_over_QQ':True,'zero_combinations_quotiented_using_full_coefficient_map':True,'all135_equations_and_F5_F6_columns_retained':True,'scope':'QQ-constant combinations of same-weight frozen identities followed by one actual variable factor of weight(total-4); no polynomial parameter combinations, arbitrary gcd or saturation'}
 text=json.dumps(evidence,indent=2)+'\n';assert len(text)<32*1024*1024;(out/'preflight.json').write_text(text)
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),source,parent,cp,ROOT/'reports/S61_NONSQUARE_BIANCHI_MONOMIAL_DIVISION_PREFLIGHT_20261104_v63/report.json',ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW/report.json'];report={'id':ID,'all_exact_constant_QQ_kernel_and_division_checks_passed':True,'weight_group_sizes':[len(groups[w])for w in [5,6,7,8]],'variable_systems':tested,'nonzero_divided_weight4_candidate_count':len(candidates),'all135_columns_and_parameter_coefficients_retained':True,'identically_zero_combinations_quotiented':True,'all_nonzero_candidates_literal135_verified':True,'no_normal_graph_transport_pivot_or_rank_claim':True,'no_regeneration_GB_field_point_or_parameter_polynomial_combination':True,'precise_class_only_constant_QQ_combinations_then_single_variable_factor':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'preflight_sha256':sha(out/'preflight.json'),'artifact_bytes':len(text),'seconds':time.monotonic()-start,'frozen':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);args=p.parse_args()
 try:main(args.output.resolve())
 except (TimeoutError,MemoryError)as error:
  args.output.mkdir(parents=True,exist_ok=True);(args.output/'BLOCKED.json').write_text(json.dumps({'completed':False,'constraint':type(error).__name__,'message':str(error),'no_auto_retry':True},indent=2)+'\n');raise
