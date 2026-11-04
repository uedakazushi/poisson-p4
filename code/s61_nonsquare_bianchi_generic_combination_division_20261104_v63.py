#!/usr/bin/env python3
"""One old-point exact QQ rank certificate for unsplit divisibility maps."""
import argparse, hashlib, importlib.util, json, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_BIANCHI_GENERIC_COMBINATION_DIVISION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def determinant(a):
 a=[list(r)for r in a];ans=Q(1);n=len(a)
 for j in range(n):
  k=next((i for i in range(j,n)if a[i][j]),None)
  if k is None:return Q()
  if k!=j:a[k],a[j]=a[j],a[k];ans=-ans
  p=a[j][j];ans*=p
  for i in range(j+1,n):
   t=a[i][j]/p
   for k in range(j+1,n):a[i][k]-=t*a[j][k]
   a[i][j]=Q()
 return ans
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(500000000,500000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('30s cap')));signal.alarm(30)
 out.mkdir(parents=True,exist_ok=False)
 old=ROOT/'code/s61_nonsquare_bianchi_constant_combination_division_20261104_v63.py'
 spec=importlib.util.spec_from_file_location('old_kernel',old);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
 op=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63/operator.json'
 prior=ROOT/'reports/S61_NONSQUARE_BIANCHI_CONSTANT_COMBINATION_DIVISION_20261104_v63/report.json'
 prev=json.loads(prior.read_text())
 for rel,h in prev['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 assert sha(op)==prev['input_sha256'][str(op.relative_to(ROOT))]
 generic=ROOT/'reports/S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.json'
 gp=json.loads(generic.read_text());assert gp['all_passed']and gp['point']==[1,2,3]
 generic_source=ROOT/'code/s61_binary_hyperplane_generic_ranks_20261102.py'
 assert sha(generic_source)==gp['input_sha256'][str(generic_source.relative_to(ROOT))]
 assert "a,e,d=s.symbols('a e d');point={a:1,e:2,d:3}"in generic_source.read_text()
 point=tuple(map(Q,gp['point']));data=json.loads(op.read_text());weights=data['local_coordinate_weights'];ew=data['equation_weights']
 assert data['matrix_entry_coefficient_core_order']==['a','e','d']
 cp=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=json.loads(cp.read_text())
 assert [i for i,w in enumerate(weights)if w==0]==[19,23]
 for i in [19,23]:
  assert all(Q(chart['T34'][i][j])==0 for j in chart['normal34'])
  assert all(Q(x)==0 for x in chart['U34'][i])
 groups={w:[]for w in [5,6,7,8]}
 for pi,rec in enumerate(data['local_weight_pieces']):
  if rec['weight']not in groups:continue
  vector={}
  for g,var,terms in rec['coefficients']:
   assert 0<=g<135 and -1<=var<45
   assert ew[g]+(0 if var==-1 else weights[var])==rec['weight']
   if var in [19,23]:continue
   c=sum((Q(c)*point[0]**ex[0]*point[1]**ex[1]*point[2]**ex[2]for ex,c in terms),Q())
   helper.vector_add(vector,{(g,var):c})
  groups[rec['weight']].append((pi,vector))
 assert [len(groups[w])for w in groups]==[20,22,18,11]
 systems=[]
 for w,rows in groups.items():
  n=len(rows);full={}
  for j,(_,v)in enumerate(rows):
   for key,c in v.items():full.setdefault(key,{})[j]=c
  for variable,wt in enumerate(weights):
   if wt!=w-4:continue
   bad={key:row for key,row in full.items()if key[1]!=variable}
   cert,k=helper.kernel(bad,n)
   square=[];det=None
   if cert['rank']==n:
    square=[[bad[tuple(key)].get(j,Q())for j in range(n)]for key in cert['selected_constraint_keys']]
    det=determinant(square);assert det
   systems.append({'total_weight':w,'factor_local_variable':variable,'factor_weight':wt,'quotient_weight':4,
      'source_weight_piece_indices':[i for i,_ in rows],'unknown_count':n,'nonzero_unsplit_constraint_rows':len(bad),
      'exact_QQ_rank_at_old_point':cert['rank'],'selected_constraint_keys':cert['selected_constraint_keys'],
      'square_minor_column_order':list(range(n))if square else[],
      'square_minor_at_old_point':[[str(c)for c in row]for row in square],
      'exact_QQ_minor_determinant':str(det)if det is not None else None,
      'generic_QQ_a_e_d_map_injective_certified':cert['rank']==n,
      'rank_drop_implies_no_generic_drop_claim':cert['rank']<n})
 allfull=all(s['generic_QQ_a_e_d_map_injective_certified']for s in systems);assert len(systems)==43
 evidence={'old_point':{'report':str(generic.relative_to(ROOT)),'parameter_order':['a','e','d'],'point':[str(v)for v in point]},
 'verified_zero_nuisance_coordinates':[19,23],
 'unsplit_constraint_keys':'(all135_equation_index, local_multiplier_variable_or_minus1); parameter coefficients evaluated together, not split',
 'systems':systems,'all43_generic_maps_injective':allfull,
 'scope':'Same-weight QQ(a,e,d)-linear combinations of frozen rows; division by one actual local variable of weight(w-4); source-level exact weight4 only',
 'no_new_point_field_source_regeneration_GB_normalgraph_or_pivot':True}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 assert (out/'certificate.json').stat().st_size<32*1024*1024
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),old,op,prior,generic,generic_source,cp]
 report={'id':ID,'all_exact_rank_minor_checks_passed':True,'old_accepted_point_only':True,'parameter_order':['a','e','d'],
 'point':[1,2,3],'systems':len(systems),'full_generic_rank_count':sum(s['generic_QQ_a_e_d_map_injective_certified']for s in systems),
 'all43_QQ_a_e_d_variable_division_maps_injective':allfull,'no_global_pointwise_or_general_saturation_claim':True,
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),
 'artifact_bytes':(out/'certificate.json').stat().st_size,'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);args=p.parse_args()
 try:main(args.output.resolve())
 except(TimeoutError,MemoryError)as err:
  args.output.mkdir(parents=True,exist_ok=True);(args.output/'BLOCKED.json').write_text(json.dumps({'completed':False,'constraint':type(err).__name__,'message':str(err),'no_auto_retry':True},indent=2)+'\n');raise
