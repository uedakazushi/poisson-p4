#!/usr/bin/env python3
"""ONE bounded residual5 DAG rank profile; execute ONLY under explicitROOTbudget."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_RANK_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 gate='reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_FINITE_SOURCE_20261105_v63';gr=src(gate+'/report.json');assert gr['all_passed']and pins[gate+'/report.json']=='3d2cd38f949f5bc29507a542409638c17031e90bddac59e4e3b7280af96d3dbb';dag=src(gate+'/finite-residual5-DAG.json');budget=src(gate+'/resource-gate.json');skip=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_SUPPORT_ZERO_20261105_v63/report.json');assert skip['Newton_sharpened_order_lower_bound']==26705
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('lin',hp);lin=importlib.util.module_from_spec(sp);sp.loader.exec_module(lin)
 t2=src('reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63/typed-next-residual-Schur.json');R2=[[Q(v)for v in rr]for rr in t2['row_transform_R2']];Q2=[[Q(v)for v in rr]for rr in t2['column_transform_Q2']];t3=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_SOURCE_STRUCTURE_20261105_v63/typed-residual5-and-next-source.json');R3=[[Q(v)for v in rr]for rr in t3['R3']];Q3=[[Q(v)for v in rr]for rr in t3['Q3']]
 requests={};literal_values={};receipts=[];filecache={}
 def cachedfile(rel):
  if rel not in filecache:filecache[rel]=src(rel)
  return filecache[rel]
 for n in dag['nodes']:
  if n['op']!='literal':continue
  s=n['source'];value=Q(n['value']);actual=None
  if 'constant'in s:actual=Q(s['constant'])
  elif s.get('kind')=='saved_formal_coefficient':requests.setdefault(s['path'],{}).setdefault(s['line'],[]).append(n);continue
  elif s.get('kind')in['R2Q2','R3Q3']:
   u,i,j,v=s['indices'];actual=(R2[u][i]*Q2[j][v])if s['kind']=='R2Q2'else(R3[u][i]*Q3[j][v])
  elif 'saved_scalar_inverse_G1'in s:actual=Q(cachedfile(s['saved_scalar_inverse_G1'])['G_ONLY0and1_QQ']['1'])
  elif 'source_A11_gap69'in s:
   p=s['source_A11_gap69'];requests.setdefault(p['path'],{}).setdefault(p['line'],[]).append(dict(n,source=dict(p,kind='saved_formal_coefficient',gap=69,block='A',row=1,column=1)));continue
  elif s.get('kind')in['completed_H','completed_E']:
   f=cachedfile(s['path']);key=tuple(s['key'])
   if s['kind']=='completed_H':m={ (e['Hrow'],e['Hcolumn'],e['relative_gap']):Q(e['exact_QQ_coefficient'])for e in f['all31']}
   else:m={(e['row'],e['column'],e['relative_gap']):Q(e['exact_QQ_coefficient'])for e in f['all25']}
   actual=m[key]
  else:raise AssertionError(s)
  assert actual==value;literal_values[n['id']]=value;receipts.append(dict(node=n['id'],source=s,verified_value=str(value)))
 for rel,byline in requests.items():
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    if ln not in byline:continue
    e=json.loads(line);c=Q(e['QQ_content']);co={int(k):int(v)for k,v in e['sparse_integer_polynomial']}
    for n in byline[ln]:
     s=n['source'];assert e['entry_index']==s['entry_index']and c==Q(s['QQ_content']);actual=c*co.get(s['gap'],0);assert actual==Q(n['value']);literal_values[n['id']]=actual;receipts.append(dict(node=n['id'],source=s,source_integer=str(co.get(s['gap'],0)),verified_value=str(actual)))
 for r in dag['cached_completed_coefficients_used']:
  f=cachedfile(r['source']);key=tuple(r['key']);m={(e['Hrow'],e['Hcolumn'],e['relative_gap']):Q(e['exact_QQ_coefficient'])for e in f['all31']};assert m[key]==Q(r['value'])
 values=[];den=[];nodeproof=[];power_exponents=[]
 for n in dag['nodes']:
  assert n['id']==len(values)and all(i<n['id']for i in n['args']);D=int(n['denominator_upper']);N=int(n['integer_abs_upper']);args=n['args']
  if n['op']=='literal':q=literal_values[n['id']];assert D==q.denominator;z=q.numerator
  elif n['op']=='sum':
   assert D==math.lcm(*(den[i]for i in args));z=0
   for w,i in zip(n['weights'],args):z+=w*values[i]*(D//den[i]);assert abs(z)<=N
  elif n['op']=='product':
   assert D==math.prod(den[i]for i in args);z=math.prod(values[i]for i in args)
  elif n['op']=='power':
   exponent=n['exponent'];power_exponents.append(exponent);assert exponent==2 and D==den[args[0]]**exponent;z=values[args[0]]**exponent
  else:raise AssertionError(n['op'])
  assert abs(z)<=N;values.append(z);den.append(D);nodeproof.append(dict(node=n['id'],op=n['op'],integer_numerator=str(z),known_denominator=str(D),exact_QQ_value=str(Q(z,D))))
  if time.monotonic()-start>45:raise RuntimeError('45s_guard')
 assert power_exponents==[2]
 RD=int(dag['root_integer_common_denominator']);Z=[[0]*5 for i in range(5)]
 for r in dag['residual5_roots']:
  k=r['node'];assert RD%den[k]==0;z=values[k]*(RD//den[k]);assert abs(z).bit_length()<=budget['root_integer_bit_upper'];Z[r['row']][r['column']]=z
 rf,P,pcols,prows,ops=lin.rank_ops(Z);rank=len(pcols);assert rank<=3;free=[j for j in range(5)if j not in pcols];kernel=[[Q(0)for j in free]for i in range(5)]
 for k,j in enumerate(free):
  kernel[j][k]=Q(1)
  for i,c in enumerate(pcols):kernel[c][k]=-rf[i][j]
 assert lin.mm(Z,kernel)==[[Q(0)for j in free]for i in range(5)]
 minor=[[Q(Z[i][j])for j in pcols]for i in prows];mdet=lin.determinant(minor)if rank else Q(1);minv=lin.inverse(minor)if rank else[];order=pcols+free;Q4=[[Q(0)for j in range(5)]for i in range(5)]
 for j,c in enumerate(order):Q4[c][j]=Q(1)
 U=lin.eye(5)
 for i in range(rank):
  for j,c in enumerate(free):U[i][rank+j]=-rf[i][c]
 Q4=lin.mm(Q4,U);R4=[[RD*v for v in rr]for rr in P];ZQ=[[Q(v,RD)for v in rr]for rr in Z];assert lin.mm(lin.mm(R4,ZQ),Q4)==[[Q(i==j and i<rank)for j in range(5)]for i in range(5)]
 dump(out/'all-literal-source-checks.json',dict(receipts=receipts,cached_source_zero_receipts=dag['cached_completed_coefficients_used'],all116_literal_nodes_bound=True,ONLY_new_closedpower_exponent=2,no_inverse_series=True))
 dump(out/'all213-integer-node-identities.json',dict(nodes=nodeproof,all_metadata_denominator_and_integer_bounds_verified=True,new_G2_finite_power_only=True))
 dump(out/'rank-kernel-unitminor.json',dict(integer_Z1=Z,common_denominator=str(RD),actual_QQ_Z1=ZQ,exact_rank=rank,right_kernel=kernel,row_transform_integerZ=P,row_reduced=rf,row_operations=ops,pivot_rows=prows,pivot_columns=pcols,
  unit_minor=minor,unit_minor_inverse=minv,unit_minor_determinant=str(mdet),R4_for_actualQQ_Z1=R4,Q4_for_actualQQ_Z1=Q4,R4_Z1_Q4='diag(Irank,0)',known_determinant_zero_from_rank_bound=True))
 n=5-rank
 dump(out/'next-source-residual.json',dict(new_residual_size=n,rank=rank,typed_source='Z5=t*L5; R4*L5*Q4=[[U,V],[W,J]],U0=Irank,V0=W0=J0=0; residual=J-W*U^-1*V',R4=R4,Q4=Q4,
  outer_order_before_new_residual=26653+5,known_whole_Newton_order=26705,next_source_goal='next possible whole coefficient a26e91d19 => degree47 determinant of newresidual (unitfactors constant1 since all previousweights excluded); derive bounded finite operator/cost beforeany coefficient job',
  no_new_coefficient_or_det_window_weight_source572=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);actual=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert actual<budget['new_RAW_upper']and actual+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,ONE_finite213node_rank_profile=True,all116literal_sourcechecks=True,all213_integer_DAG_identities=True,ONLY_new_scalar_power_G2=True,exact_firstpositive5_rank=rank,remaining_typed_residual_size=n,kernel_and_nonzero_unitminor_verified=True,
  no_5det_nonzero_attempt=True,no_inverse_series_otherweight_precisionwindow_or_source572=True,known_Newton_order26705_unchanged=True,no_full_det11_zero_nonzero_or_global57_claim=True,
  actual_new_RAW_before_report=actual,actual_external_read_RAW=external,actual_new_plus_read_RAW=actual+external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
