#!/usr/bin/env python3
"""Read-only audit of the proposed114 univariate fifth coefficients.

Evaluates saved coefficients only throughorder4. No fifth unknown is solved.
Exact semantic DAG signatures bind all signed numerators to original sources.
"""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63_REVIEW'
PATHS={
 'author_report':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63/report.json',
 'input':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63/exact-two57-univariate5-input.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'fourth':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'fourth_peer':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json',
 'first':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
 'copy_peer':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json',
 'full':'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(x):return {int(i):Q(c) for i,c in x if Q(c)}
def jm(a,b):return tuple(sum(a[i]*b[n-i] for i in range(n+1)) for n in range(5))
def ji(a):
 assert a[0];v=[1/a[0]]
 for n in range(1,5):v.append(-sum(a[i]*v[n-i] for i in range(1,n+1))/a[0])
 assert jm(a,v)==(1,0,0,0,0)
 return tuple(v)
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'independent_source_order_audit_inputs_saved','inputs':inputs},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['author_report','fourth_peer','copy_peer','C0_peer'])
 for p,h in D['author_report']['input_sha256'].items():assert sha(ROOT/p)==h
 src=D['source']['DAG_scalar_nodes'];mon=D['fourth']['monomials_bivariate_through4']
 y=D['copy']['actual_YB_last2_from_saved105'];yv=[rec(p) for i,p in y]
 assert [v.get(1,Q(0)) for v in yv]==[Q(17,12),Q(7,12)] and all(not v.get(2,Q(0)) for v in yv)
 assert [i for i,p in y]==D['first']['both57_source2']['CB']['F9_nodes'][7:9]
 # Structural semantics keep variable NODE IDs, never display names. Shared
 # canonical expression hashes handle duplicate negations in later DAGs.
 intern={};forms=[]
 def put(k):
  if k not in intern:intern[k]=len(forms);forms.append(k)
  return intern[k]
 def cc(q):return put(('C',Q(q)))
 z=cc(0);one=cc(1)
 def ng(a):
  n=forms[a]
  if n[0]=='C':return cc(-n[1])
  if n[0]=='N':return n[1]
  return put(('N',a))
 def aa(a,b):
  if a==z:return b
  if b==z:return a
  if a==ng(b) or b==ng(a):return z
  if forms[a][0]==forms[b][0]=='C':return cc(forms[a][1]+forms[b][1])
  return put(('A',min(a,b),max(a,b)))
 def mm(a,b):
  if a==z or b==z:return z
  if a==one:return b
  if b==one:return a
  if forms[a][0]==forms[b][0]=='C':return cc(forms[a][1]*forms[b][1])
  return put(('M',min(a,b),max(a,b)))
 def inv(a):
  if forms[a][0]=='C':assert forms[a][1];return cc(1/forms[a][1])
  return put(('I',a))
 def signatures(nodes,roots,bindings=None):
  need=set();todo=list(roots)
  while todo:
   i=todo.pop()
   if i in need:continue
   need.add(i);n=nodes[i]
   if n[0] in ['a','m']:todo.extend(n[1:3])
   elif n[0] in ['n','i']:todo.append(n[1])
   else:assert n[0] in ['c','v']
  vals={}
  for i in sorted(need):
   n=nodes[i];op=n[0]
   if op=='c':v=cc(Q(n[1],n[2]))
   elif op=='v':v=bindings[i] if bindings is not None else put(('V',i))
   elif op=='a':v=aa(vals[n[1]],vals[n[2]])
   elif op=='m':v=mm(vals[n[1]],vals[n[2]])
   elif op=='n':v=ng(vals[n[1]])
   elif op=='i':v=inv(vals[n[1]])
   vals[i]=v
  return {i:vals[i] for i in roots}
 records=[];resource_bound=0
 for f in D['input']['families']:
  name=f['family'];old=D['fourth']['families'][name];frame=D['first']['both57_source2'][name]
  assert f['original_unknown57_order']==frame['unknown57']==old['actual_unknown57_order']
  assert f['original_equations57']==frame['equations57']==old['actual_equations57']
  assert len(frame['equations57'])==57
  slope=list(map(Q,f['last2_target_slopes']));assert slope==([Q(1),Q(1)] if name=='CA' else [Q(17,12),Q(7,12)])
  def line(p):
   v=[Q(0)]*5
   for (a,b),c in zip(mon,p):v[a+b]+=Q(c)*slope[0]**a*slope[1]**b
   return tuple(v)
  sol={i:line(p) for i,p in old['all57_source4_jets']};assert [[i,list(map(str,sol[i]))] for i in frame['unknown57']]==f['cached_univariate_source_through4']
  needed=f['source_needed_original_DAG_nodes'];roots=frame['equations57']+[f['actual_g_original_node']]+f['actual_bbar4_original_nodes'];uses=collections.Counter()
  for i in needed:
   if src[i][0] in ['a','m']:uses.update(src[i][1:3])
   elif src[i][0] in ['n','i']:uses.update(src[i][1:])
  uses.update(roots);vals={};peak=0;unit_checks=0
  for i in needed:
   n=src[i];op=n[0]
   if op=='c':v=(Q(n[1],n[2]),0,0,0,0)
   elif op=='v':v=sol[i]
   elif op=='a':v=tuple(a+b for a,b in zip(vals[n[1]],vals[n[2]]))
   elif op=='m':v=jm(vals[n[1]],vals[n[2]])
   elif op=='n':v=tuple(-a for a in vals[n[1]])
   elif op=='i':v=ji(vals[n[1]]);unit_checks+=1
   else:raise AssertionError(n)
   vals[i]=v;peak=max(peak,len(vals))
   for j in n[1:3] if op in ['a','m'] else n[1:] if op in ['n','i'] else []:
    uses[j]-=1
    if not uses[j]:vals.pop(j,None)
  for k,i in enumerate(frame['equations57']):assert vals[i]==(0,slope[k-55] if k>=55 else 0,0,0,0)
  gp=vals[f['actual_g_original_node']];nums=[vals[i] for i in f['actual_bbar4_original_nodes']]
  assert list(map(str,gp))==f['cached_g_line_through4'] and gp[3] and gp[:3]==(0,0,0) and gp[4]==0
  assert all(all(Q(c)==0 for c in p[:10]) for p in old['four_numerator_all15_coefficients'])
  assert all(p[:4]==(0,0,0,0) for p in nums)
  h1=list(map(Q,f['cached_four_quotient_first_line']));assert all(b[4]==gp[3]*h for b,h in zip(nums,h1))
  # Verify the fixed inverse without doing any new inversion or solving u5.
  jac=[list(map(Q,r)) for r in old['fixed_J57']];iv=[{int(j):Q(c) for j,c in row} for row in f['fixed_J57_inverse_sparse_rows']]
  for i in range(57):
   for j in range(57):assert sum(jac[i][k]*iv[k].get(j,Q(0)) for k in range(57))==Q(i==j)
  # Exact signed source-function identities to the frozen old/new quotients.
  if name=='CB':
   target=D['full']['DAG_scalar_nodes'];inst=D['full']['all16_pointed_unit57_inverse_instances'][14];assert inst['retained_Y_indices']==[7,8]
   bind=dict(inst['bindings']);quot=D['full']['all4_known_regular_g_divisions'];tn=[quot[0]['denominator_g']]+[q['numerator'] for q in quot]
  else:
   target=D['C0']['DAG_scalar_nodes'];inst=D['C0']['all7_pure_CA_tail_prefix57_systems'][0];assert inst['retained_zero_prefix_length']==0
   bind=dict(inst['bindings']);quot=D['C0']['all4_known_regular_gA_divisions'];tn=[quot[0]['denominator']]+[q['numerator'] for q in quot]
  assert set(bind)==set(frame['unknown57'])
  target_values=signatures(target,tn+list(bind.values()))
  source_values=signatures(src,[f['actual_g_original_node']]+f['actual_bbar4_original_nodes'],{i:target_values[j] for i,j in bind.items()})
  original_roots=[f['actual_g_original_node']]+f['actual_bbar4_original_nodes']
  assert [source_values[i] for i in original_roots]==[target_values[i] for i in tn]
  # The static cost bound anticipates a sequential6-coefficient evaluation;
  # it is not a claim that missing fifth heights are already known.
  nnz=sum(len(r) for r in iv);bound=6*(peak+57+len(roots))+nnz+57*6+len(needed)+5000
  assert bound<=f['conservative_live_QQ_slots_and_records_bound']<150000;resource_bound=max(resource_bound,bound)
  records.append({'family':name,'slope':list(map(str,slope)),'all57_saved_source_equations_zero_through4':True,
   'all_signed_original_g_and_bbar_roots_match_frozen_quotients':True,
   'g_line_through4':list(map(str,gp)),'global_bbar_order_at_least4':True,'cofactor_origin_zero':True,
   'scalar_line_h1':list(map(str,h1)),'future_h2_formula':'bbar5/g3',
   'unit_inverse_source_checks':unit_checks,'source_peak_live_nodes':peak,'fixed_J57_times_inverse_identity':True,
   'prospective_source5_unknown_count':57,'prospective_live_bound':bound})
 ev={'source_manifest':inputs,'both_families_proof_records':records,
  'actual_CB_coordinate_curvature':'YB7,YB8 = (17/12,7/12)*r+0*r²+O(r³), bound to ownCBF9 node IDs',
  'source_order_proof':'regular h has globalh0=0 because full own-plane bbar degree3 iszero while g3 nonzero. Along either slope b5=g3*h2+g4*h1+g5*h0; g4=h0=0 implies h2=b5/g3. Pullback second differs from own-slope second only by Dh0*Y2, andY2=0.',
  'no_new_g5_common105_4_or5_needed':True,'ordinary_coefficient_not_factorial':True,
  'old_minus_new_sign':'gdiff = oldCB quotient - newCA quotient in nodeIDs; Xp55_2 is negative gdiff2 in the separate Schur gate',
  'unit_source5_formula':'with u5=0 evaluate degree5 residual of the SAME57 original equations; u5=-fixed_J0^-1 residual. J(u)-J0 timesu5*r5 startsdegree6; aftercorrection verify ALL57 equations through5.',
  'strict_guard_note':'height and elapsedtime of missing fifth values are notknown; any authorized futureworker must retain 60s/2GB/32MiB/150k guards and failedprefix if exceeded.',
  'new_source5_or_g2_unknowns_solved':False,'prospective_new_unknowns114':True,
  'prospective_live_bound_max':resource_bound,'this_preflight_does_not_authorize_worker':True}
 ep=out/'independent-two57-fifth-input-and-source-order-audit.json';blob=json.dumps(ev,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'actual_CA_CB_slopes_and_zero_CB_curvature_accepted':True,
  'both57_equations_saved0_through4_and_fixed_inverse_accepted':True,
  'all_signed_old_new_numerator_denominator_source_bindings_accepted':True,
  'global_h0_zero_and_g4_zero_imply_h2_b5_over_g3_accepted':True,
  'one_order114_univariate_coefficients_suffice_if_authorized':True,
  'new_source5_g2_or114_unknown_values_computed':False,'max_prospective_live_bound':resource_bound,
  'input_hashes':inputs,'code_sha256':sha(Path(__file__)),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
