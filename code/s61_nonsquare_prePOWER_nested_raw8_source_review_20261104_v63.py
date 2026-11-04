#!/usr/bin/env python3
"""Independent fixed raw8 source/unit, integer-bound and nested-lifetime review.
No evaluation of the parameter determinant or worker duplication.
"""
import collections,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parents[1];S='S61_NONSQUARE_PREPOWER_NESTED_RAW8_SOURCE_20261104_v63_REVIEW';O=R/'reports'/S
NEST='S61_NONSQUARE_PREPOWER_NESTED_RAW8_SUPPORT_PREFLIGHT_20261104_v63';FUSED='S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_FUSED_PREFLIGHT_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2)+'\n')
def lattice(b):
 lo,hi=b['degree'];box=b['box'];return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1) for a in range(box[0][0],box[0][1]+1) for e in range(box[1][0],box[1][1]+1))
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));O.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=R/'reports'/st/name;pins[str(p.relative_to(R))]=sha(p);return json.loads(p.read_text())
 nr=load(NEST);assert pins[f'reports/{NEST}/report.json']=='a3f480f9f1904cade9ffc42af2e746ad07c19dd72cfb3bed669f6c40eeaec230'
 fr=load(FUSED);assert fr['all_passed'] and not fr['resource_gate_passes'];plan=load(NEST,'nested-plan.json');op=load(FUSED,'operator-and-cost.json');fproof=load(FUSED,'fused-accumulation-proof.json')
 for rep in [nr,fr]:
  for p,h in rep['input_sha256'].items():assert sha(R/p)==h
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json');h9=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');fx=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');st=load('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json')
 assert load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63_REVIEW')['all_passed']
 assert load('S61_NONSQUARE_PREPOWER_DIRECT_P5_H5_SCHUR_20261104_v63_REVIEW')['all_passed']
 assert op['source_current15_rows']==[0,1,3,6,2,4,5,8] and op['source_original32_columns']==[22,11,13,25,15,19,27,31]
 for i in [1,3,4,7]:assert hc['known_row_units_rho14'][i]['a_d_F_B_C_J_powers']==[1,0,1,0,2,0] and all(not hc['raw_compact14x5'][i][j] for j in [1,2,3])
 # Reconstruct every literal atom from the immutable source, including all QQ/a normalization.
 N=op['scalar_nodes'];cache={v['old_DAG_root']:v['Laurent_polynomial'] for v in fx['node_values']};old=q['scalar_DAG']['nodes'];polys={}
 for i,n in enumerate(N):
  assert n['id']==i and all(j<i for j in n['args'])
  if n['op']!='atom':continue
  ref=n['reference'];path=ref.get('path')
  if 'constant' in ref:p=[] if not F(ref['constant']) else [[[0,0,0],str(F(ref['constant']))]]
  elif path=='accepted_Q2L_source_atom':
   j=ref['old_DAG_root'];d=old[j];assert ref['definition']==d
   a,b=cp['matrix25x32'][d[1]][d[2]] if d[0]=='checkpoint' else pd['Lframe_dual6x2'][d[1]][1] if d[0]=='accepted_pair_second' else q['scalar_DAG']['new_unit_atoms_num_den'][d[1]]
   assert len(b)==1 and b[0][0][1:]==[0,0] and n['QQ_unit']==str(F(b[0][1])) and n['a_shift']==-b[0][0][0];p=a
  elif path=='reports/S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63/native.json':p=cache[ref['old_DAG_root']]
  elif path=='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json':
   r,c=ref['index'];assert r in [1,3,4,7] and c in [0,4] and ref['divided_heavy_columns_by_only_FULLU_units'];p=hc['raw_compact14x5'][r][c];assert F(n['QQ_unit'])==1/F(hc['known_row_units_rho14'][r]['constant'])
  elif path=='reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json':p=h9['source_h9'][h9['original_target_labels'].index(ref['h9_target'])]
  else:assert path=='reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63/native.json' and ref['key']=='exact_unscaled_pure4_minor';p=st['exact_unscaled_pure4_minor']
  p={tuple([m[0]+n['a_shift'],m[1],m[2]]):F(c)/F(n['QQ_unit']) for m,c in p if F(c)};polys[i]=p
 # Exact source-wire fingerprints, with source polynomial atoms canonicalized but not multiplied.
 fingerprints={};descriptor=[]
 def key(v):
  if v not in fingerprints:fingerprints[v]=len(descriptor);descriptor.append(v)
  return fingerprints[v]
 def atomkey(poly):return key(('poly',tuple(sorted((m,str(v)) for m,v in poly.items() if v))))
 zero=atomkey({});one=atomkey({(0,0,0):F(1)});minus=atomkey({(0,0,0):F(-1)})
 def em(x,y):
  if zero in [x,y]:return zero
  if x==one:return y
  if y==one:return x
  return key(('mul',*sorted((x,y))))
 def ea(x,y,sign=1):
  if x==zero:return y if sign==1 else em(minus,y)
  if y==zero:return x
  if sign==-1 and x==y:return zero
  return key(('add' if sign==1 else 'sub',x,y))
 def es(x,k):
  if x==zero or not k:return x
  d=descriptor[x]
  if d[0]=='poly':return atomkey({(m[0]+k,m[1],m[2]):F(c) for m,c in d[1]})
  return key(('a_shift',x,k))
 expr=[]
 for i,n in enumerate(N):
  if n['op']=='atom':v=atomkey(polys[i])
  elif n['op']=='known_a_shift':v=es(expr[n['args'][0]],n['a_exponent_shift'])
  else:
   a,b=[expr[j] for j in n['args']];v=em(a,b) if n['op']=='mul' else ea(a,b,-1 if n['op']=='sub' else 1)
  expr.append(v)
 qcache={}
 def oldexpr(i):
  if i in qcache:return qcache[i]
  n=old[i];kind=n[0]
  if i in cache:v=atomkey({tuple(m):F(c) for m,c in cache[i]})
  elif kind=='z':v=zero
  elif kind=='one':v=one
  elif kind in ['a','s','m']:
   a,b=oldexpr(n[1]),oldexpr(n[2]);v=em(a,b) if kind=='m' else ea(a,b,-1 if kind=='s' else 1)
  else:
   a,b=cp['matrix25x32'][n[1]][n[2]] if kind=='checkpoint' else pd['Lframe_dual6x2'][n[1]][1] if kind=='accepted_pair_second' else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
   assert len(b)==1 and b[0][0][1:]==[0,0]
   v=atomkey({(m[0]-b[0][0][0],m[1],m[2]):F(c)/F(b[0][1]) for m,c in a if F(c)})
  qcache[i]=v;return v
 rows=op['source_current15_rows'];labels=op['source_original32_columns'];cols=q['new22_target_original32_columns'];mat=op['matrix8x8_after_row_a_normalization']
 assert op['raw_input_roots']==[[q['remaining15x22_scalar_roots'][r][cols.index(j)] for j in labels] for r in rows]
 for ri,r in enumerate(rows):
  for ci,c in enumerate(labels):
   if ri<4 and c in [27,31]:v=zero
   elif ri>=4 and c in [27,31]:
    rr=[1,3,4,7][ri-4];k=0 if c==27 else 4;constant=F(hc['known_row_units_rho14'][rr]['constant']);v=atomkey({tuple(m):F(x)*constant for m,x in hc['raw_compact14x5'][rr][k]})
   elif r==0 and c in h9['original_target_labels']:v=atomkey({tuple(m):F(x) for m,x in h9['source_h9'][h9['original_target_labels'].index(c)]})
   else:v=oldexpr(q['remaining15x22_scalar_roots'][r][cols.index(c)])
   assert expr[mat[ri][ci]]==es(v,-op['row_known_a_factors'][ri]),(ri,ci)
 # Fixed arithmetic support/denominator/height recursion. No parameter products evaluated.
 bnds=[];dens=[];norms=[]
 for i,n in enumerate(N):
  kind=n['op'];args=n['args']
  if kind=='atom':
   p=polys[i];d=math.lcm(*(v.denominator for v in p.values())) if p else 1;nn=sum(abs(v.numerator)*(d//v.denominator) for v in p.values());b=None if not p else {'support':len(p),'degree':[min(map(sum,p)),max(map(sum,p))],'box':[[min(m[j] for m in p),max(m[j] for m in p)] for j in range(3)]}
  elif kind=='known_a_shift':
   j=args[0];d,nn=dens[j],norms[j];bb=bnds[j];k=n['a_exponent_shift'];b={'support':bb['support'],'degree':[v+k for v in bb['degree']],'box':[[v+k for v in bb['box'][0]],*bb['box'][1:]]}
  else:
   j,k=args;A,B=bnds[j],bnds[k];assert A and B
   if kind=='mul':
    d=dens[j]*dens[k];nn=norms[j]*norms[k];lo,hi=[A['degree'][a]+B['degree'][a] for a in [0,1]];box=[[A['box'][a][v]+B['box'][a][v] for v in [0,1]] for a in range(3)];limit=A['support']*B['support']
   else:
    assert kind in ['add','sub'];d=math.lcm(dens[j],dens[k]);nn=norms[j]*(d//dens[j])+norms[k]*(d//dens[k]);lo=min(A['degree'][0],B['degree'][0]);hi=max(A['degree'][1],B['degree'][1]);box=[[min(A['box'][a][0],B['box'][a][0]),max(A['box'][a][1],B['box'][a][1])] for a in range(3)];limit=A['support']+B['support']
   b={'degree':[lo,hi],'box':box};b['support']=min(limit,lattice(b),nn)
  assert b==n['bound'] and d==int(n['QQ_denominator']) and nn==int(n['integer_l1_upper']),(i,kind)
  bnds.append(b);dens.append(d);norms.append(nn)
 size=lambda i:bnds[i]['support'] if bnds[i] else 0
 # Independently check all90 signed inner terms and15 outer scales/events.
 terms=op['laplace15'];heavy={x['root'] for x in op['heavy2_roots']};assert len(terms)==15 and len(heavy)==6
 def linearize(i,c=1):
  n=N[i];kind=n['op'];a=n['args']
  if bnds[i] is None:return []
  if kind=='mul' and any(j in heavy for j in a):
   x=a[0] if a[0] in heavy else a[1];y=a[1] if a[0] in heavy else a[0];return [(c,x,y)]
  if kind in ['add','sub']:return linearize(a[0],c)+linearize(a[1],c if kind=='add' else -c)
  if kind=='mul':
   for j,k in [a,a[::-1]]:
    p=polys.get(j)
    if p is not None and len(p)==1 and (0,0,0) in p and p[(0,0,0)] in [1,-1]:return linearize(k,c*int(p[(0,0,0)]))
  raise AssertionError(('bottom grammar',i))
 inner_events=collections.defaultdict(list);outer_events=collections.defaultdict(list);needed_roots=[]
 for k,(t,sup) in enumerate(zip(terms,plan['inner15'])):
  assert t['pure_top4_columns']==list(list(itertools.combinations(range(6),4))[k]) and t['pure_bottom2_columns']==[i for i in range(6) if i not in t['pure_top4_columns']]
  b=t['bottom4_root'];top=t['top4_root'];pieces=linearize(b);assert len(pieces)==6 and {h for c,h,p in pieces}==heavy
  expected=[]
  for c,h,p in pieces:
   assert dens[b]%(dens[h]*dens[p])==0;v={'sign':c,'H2':h,'P2':p,'integer_QQ_scale':str(dens[b]//(dens[h]*dens[p]))};expected.append(v);inner_events[max(h,p)].append({'bottom_index':k,**v});needed_roots.extend([h,p])
  assert expected==sup['six_terms'] and sup['dense_slots']==lattice(bnds[b])
  assert dens[op['output']]%(dens[top]*dens[b])==0;v={'bottom_index':k,'top4':top,'bottom4':b,'integer_QQ_scale':str(dens[op['output']]//(dens[top]*dens[b])),'sign':(-1)**(sum(t['pure_top4_columns']))};outer_events[max(top,b)].append(v);needed_roots.append(top)
 assert {str(k):v for k,v in inner_events.items()}==plan['inner_events'] and {str(k):v for k,v in outer_events.items()}==plan['outer_events']
 need=set();todo=needed_roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);todo.extend(N[i]['args'])
 assert sorted(need)==plan['ordinary_nodes'] and len(need)==1898 and all(t['bottom4_root'] not in need for t in terms)
 uses=collections.Counter()
 for i in need:uses.update(N[i]['args'])
 for es in inner_events.values():
  for e in es:uses.update([e['H2'],e['P2']])
 for es in outer_events.values():
  for e in es:uses[e['top4']]+=1
 outer=lattice(bnds[op['output']]);atom=sum(size(i) for i in need if N[i]['op']=='atom');live=peak=0;active={};completed=[];records=[];peakid=None
 def drop(i):
  nonlocal live
  uses[i]-=1;assert uses[i]>=0
  if not uses[i]:live-=size(i)
 for i in sorted(need):
  live+=size(i);scratch=size(i) if N[i]['op'] in ['add','sub','mul'] else 0;tot=live+sum(active.values())+scratch+outer+atom+15000
  if tot>peak:peak=tot;peakid=i
  for j in N[i]['args']:drop(j)
  for e in inner_events[i]:
   k=e['bottom_index']
   if k not in active:assert not active;active[k]=plan['inner15'][k]['dense_slots']
   tot=live+sum(active.values())+outer+atom+15000
   if tot>peak:peak=tot;peakid=i
   drop(e['H2']);drop(e['P2'])
  for e in outer_events[i]:drop(e['top4']);del active[e['bottom_index']];completed.append(e['bottom_index'])
  records.append({'node':i,'ordinary_live_upper':live,'active_bottom_buffers':{str(k):v for k,v in active.items()},'output_outer_slots':outer,'inner_event_count':len(inner_events[i]),'outer_event_count':len(outer_events[i])})
 assert completed==list(range(15)) and not active and live==0 and records==plan['last_use_records']
 assert (peak,peakid,outer,atom)==(131590,305,58259,8048) and peak<150000 and nr['unchanged_work_reserve']==15000
 # Conditional source11 block: all newly appended upper-right3 entries are literal zero.
 assert op['row_known_a_factors']==[0,1,1,0,0,0,0,0]
 future=op['future_pivot11'];assert future['rows']==[0,1,3,6,2,4,5,8,7,13,10] and future['original32_columns']==[22,11,13,25,15,19,27,31,28,30,29]
 assert future['remaining_current15']==[9,11,12,14] and future['remaining_original32']==[1,2,3,4,6,7,8,9,14,17,18]
 for i in [1,3,4,7]:assert all(not hc['raw_compact14x5'][i][j] for j in [1,2,3])
 # Universal Jacobi complementary-minor identity is a polynomial identity: prove on D(delta2)
 # by block Schur, then descend in the integral universal polynomial ring. No localization used
 # in the resulting formula/raw8 operator.
 evidence={'all2094_fixed_scalar_bounds_integer_denominators_l1':True,'all_literal_source_atoms_reconstructed':len(polys),'all15_same_raw8_source_terms_and90_inner105_scales':True,'nested_ordinary_nodes':1898,'last_use_peak':peak,'outer_dense_slots':outer,'max_one_bottom_dense_slots':max(x['dense_slots'] for x in plan['inner15']),'source_atom_reserve':atom,'unchanged_work_reserve':15000,'artifact_bound_bytes':nr['raw_output_and_atom_artifact_bytes_upper'],'no_runtime_completion_guarantee':True,'known_column_units_only':['L=a*F*C^2','L*gamma'],'source_row_a_factors':op['row_known_a_factors'],'raw8_to_normalized_factor':'L^2*gamma*a^2','no_p_R_h22_inverse_or_guessed_pure_C_factor':True,'upper8_by3_zero_from_actual_source70_identity':True,'future11_determinant_conditional':'knownFULLUunit * normalizedD8 * Delta407','future4x11_requires_actual_D8_nonzero_and_chart_unit':True,'no_actual_raw8_value_or_nonzero_claim':True,'old445860_and158134_failed_gates_not_overwritten':True}
 save('independent-source-lifetime-and-conditional11-proof.json',evidence)
 report={'id':S,'all_passed':True,'source_only':True,'all_actual_source_units_orders_and_nested_plan_verified':True,'all2094_bounds_and105_integer_scales':True,'total_live_upper':peak,'artifact_upper':nr['raw_output_and_atom_artifact_bytes_upper'],'prospective_11pivot_only_conditional':True,'actual_D8_value_not_evaluated_or_nonzero_assumed':True,'no_producer_import_rerun_coefficient_sourcequery_point_prime_rank_or_GB':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md'),'evidence_sha256':{'independent-source-lifetime-and-conditional11-proof.json':sha(O/'independent-source-lifetime-and-conditional11-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};save('report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
