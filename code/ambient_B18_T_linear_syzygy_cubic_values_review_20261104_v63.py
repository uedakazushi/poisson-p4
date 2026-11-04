#!/usr/bin/env python3
"""Independent actual cubic witness71: GB quotient, adjoint, and forward B11 AD.

Only one nonzero witness is asserted. Cached108 source values are an explicit
already independent dependency. Old typed roots are reevaluated forwards with
all other atom second coefficients symbolic, never assigned actual zero.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_T_known_initial_module_length_upper24_review_20261104_v63 import divide,lead
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW'
PATHS={
 'report':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63/report.json',
 'values':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63/actual-linear-syzygy-cubic-source-values.json',
 'adjoint':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase1-cubic-Pi2-and-source-adjoint.json',
 'fold':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase2-linear-syzygy-source-root-contractions.json',
 'GB':'reports/AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63/known-initial-module-and-length-bound.json',
 'GB_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63_REVIEW/report.json',
 'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'w0':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
 'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'S':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'a2':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'theorem':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/report.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'atom1':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'cache':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
 'cache_peer':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json',
 'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
 'native_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(v):return {int(i):Q(c) for i,c in v if Q(c)}
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(c)] for k,c in sorted(v.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,x in b.items():
  z=a.get(i,Q(0))+c*x
  if z:a[i]=z
  else:a.pop(i,None)
 return a
def put(a,k,c):
 z=a.get(k,Q(0))+c
 if z:a[k]=z
 else:a.pop(k,None)
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'independent_source_witness_inputs_saved','inputs':inputs},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['report','GB_peer','theorem','cache_peer','native_peer'])
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 assert D['cache_peer']['all48_copy_leaf_labels_and_all108_composite_values_accepted'] if 'all48_copy_leaf_labels_and_all108_composite_values_accepted' in D['cache_peer'] else D['cache_peer']['all_passed']
 # Independent Pi2 from the full39 homogeneous GB; no producer's cubic
 # rewrite/rank is used to construct these quotient maps.
 gs=[{tuple(e)+(int(t),):Q(c) for t,e,c in p} for p in D['GB']['homogeneous_generators39']]
 ms2=list(itertools.combinations_with_replacement(range(8),2));i2={m:i for i,m in enumerate(ms2)}
 ms3=list(itertools.combinations_with_replacement(range(8),3));i3={m:i for i,m in enumerate(ms3)}
 def ex(m):return tuple(m.count(i) for i in range(8))
 standard=[ex((5,7,7))+(4,),ex((7,7,7))+(4,)];Pi=[{},{}]
 for t in range(5):
  for j,m in enumerate(ms3):
   r,tr=divide({ex(m)+(t,):Q(1)},gs);assert set(r)<=set(standard)
   for k,s in enumerate(standard):
    if r.get(s):Pi[k][120*t+j]=r[s]
 assert [enc(p) for p in Pi]==D['adjoint']['two_quotient_dual_rows600']
 # Check every known1424 cubic relation in its quotient5 standard frame.
 # The same39 GB membership independently covers its actual source quads.
 # The scalar48 origin native source cannot alter this quotient selection.
 qnorm=list(map(rec,D['low']['normalized84_unitcycle_Q1_columns_sparse']))
 wb=list(map(rec,D['w0']['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):j for j,v in enumerate(wb)}
 J=[]
 for w in range(84):
  v={w:Q(1)};coeff={}
  while v:
   k=min(v);p=wl[k];a=v[k];put(coeff,p,a);v=add(v,wb[p],-a)
  y={}
  for p,a in coeff.items():y=add(y,qnorm[p],a)
  J.append(y)
 for p,col in enumerate(wb):
  check={}
  for w,c in col.items():check=add(check,J[w],c)
  assert check==qnorm[p]
 # Forward sparse derivative/Euler scattering. Qnorm*U0^-1 is applied first,
 # then each input cubic monomial is differentiated with OUTPUT Euler shifts.
 vadj=[{},{}];wedge=D['w0']['literal_wedge_basis84']
 for w,(s,pair) in enumerate(wedge):
  i,j=pair;terms=[(8*s+7,i,1,1)] if j==7 else [(8*s+j,i,1,2),(8*s+i,j,-1,2)]
  for row,var,sign,shift in terms:
   for mindex,m in enumerate(ms3):
    if var not in m:continue
    small=list(m);small.remove(var);den=sum(k<7 for k in small)+shift;dc=Q(sign*m.count(var),den)
    for tv,c in J[w].items():
     t,v=divmod(tv,8);output=120*t+i3[tuple(sorted(small+[v]))]
     for channel,d in enumerate(Pi):
      if output in d:put(vadj[channel],120*row+mindex,-dc*c*d[output])
 m0=[list(map(Q,r)) for r in D['V']['all24_by56_m0']];U=D['X0']['free_kernel_rows32'];pi=[{U[j]:c for j,c in rec(r).items()} for r in D['S']['constant_free32_quotient_projection5_rows']]
 kap=[{},{}]
 for ch in range(2):
  for index,c in Pi[ch].items():
   t,m=divmod(index,120)
   for src,a in pi[t].items():put(kap[ch],120*src+m,c*a)
  for index,c in vadj[ch].items():
   row,m=divmod(index,120)
   for src,a in enumerate(m0[row]):
    if a:put(kap[ch],120*src+m,c*a)
 assert [enc(k) for k in kap]==D['adjoint']['two_Xsource3_adjoint_rows56_times120']
 A=D['a2'];rnum=71;syz=rec(A['all119_linear_syzygies_source296'][rnum]['coefficients296_sparse']);src37=list(map(rec,A['actual_selected37_original768_constant_combinations']))
 rows=[{} for _ in range(8)]
 for jv,c in syz.items():j,v=divmod(jv,8);rows[v]=add(rows[v],src37[j],c)
 assert [enc(r) for r in rows]==D['fold']['all119_linear_coefficient_original768_combinations'][rnum]
 X0=D['X0']['full_X0_56_by768']
 for coeff in rows:
  for r in range(56):
   assert sum(Q(X0[r][j])*c for j,c in coeff.items() if j<640)==0
   assert sum(Q(X0[r][j])*c for j,c in coeff.items() if j>=640)==0
 bcols=[dict(c) for c in D['chain']['B48_by640_sparse_columns']];ops=D['chain']['operator_scalar_nodes'];cn=D['C0']['DAG_scalar_nodes'];pcols=D['C0']['honest_D0_relation_columns']
 def label(j,row,pair):
  p=i2[pair]
  if j>=640:return (2,j-640,row,p) if row<48 else None
  if row>=48:
   root=pcols[j][row-48+1];return None if cn[root][0]=='c' else (1,root,p)
  root=bcols[j].get(row)
  if root is None or ops[root][0]=='c':return None
  return (1,ops[root][1],p) if ops[root][0]=='actual_C0_scalar_function' else (0,root,p)
 xp=list(map(rec,A['actual_Xp1_sparse27_columns']));a2=rec(A['all119_actual_a2_sparse27_times36'][rnum]);func=[];cor=[]
 for d in kap:
  f={};known=Q(0)
  for index,c in d.items():
   src,m=divmod(index,120);mon=ms3[m]
   for v in set(mon):
    pair=list(mon);pair.remove(v);pair=tuple(pair)
    for j,b in rows[v].items():
     k=label(j,src,pair)
     if k is not None:put(f,k,c*b)
    for p,first in enumerate(xp):known-=c*first.get(8*src+v,Q(0))*a2.get(36*p+i2[pair],Q(0))
  func.append(f);cor.append(known)
 authorrow=D['fold']['all119_two_cubic_source_functionals'][rnum]
 assert [enc(f) for f in func]==authorrow['two_scalar_root_X2_functionals']
 assert list(map(str,cor))==authorrow['known_Xp1_a2_correction']
 allowed={(0,7):7,(1,7):8,(4,7):9,(5,7):10,(7,7):2}
 assert all(ms2[k[-1]] in allowed for f in func for k in f)
 # Independent FORWARD B11 algebra on all required typed roots. Every atom
 # not in the already peer-accepted108 cache keeps five symbolic seconds.
 bm=[(0,)*8,ex((7,)),ex((7,7)),ex((0,)),ex((1,)),ex((4,)),ex((5,)),ex((0,7)),ex((1,7)),ex((4,7)),ex((5,7))];bi={m:i for i,m in enumerate(bm)}
 def times(a,b):
  z={}
  for i,c in a.items():
   for j,d in b.items():
    k=bi.get(tuple(x+y for x,y in zip(bm[i],bm[j])))
    if k is not None:put(z,k,c*d)
  return z
 cv={int(i):rec(v) for i,v,g in D['cache']['all108_composite_values_and_opaque_second_gradients']};assert all(not g for i,v,g in D['cache']['all108_composite_values_and_opaque_second_gradients'])
 atom={int(i):rec(v) for i,v in D['atom1']['all79_atom_values']};linmap={0:0,8:1,1:3,2:4,5:5,6:6}
 needed=set();todo=[k[1] for f in func for k in f if k[0]==0]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
 val={};unk={}
 for i in sorted(needed):
  n=ops[i]
  if n[0]=='c':v={0:Q(n[1],n[2])} if n[1] else {};u={}
  elif n[0]=='actual_C0_scalar_function':
   k=n[1]
   if k in cv:
    v=cv[k];u={};assert all(v.get(b,Q(0))==atom[k].get(j,Q(0)) for j,b in linmap.items())
   else:
    v={b:atom[k][j] for j,b in linmap.items() if atom[k].get(j,Q(0))};u={(b,k):Q(1) for b in [2,7,8,9,10]}
  elif n[0]=='sum':
   v={};u={}
   for j,a,b in n[1]:v=add(v,val[j],Q(a,b));u=add(u,unk[j],Q(a,b))
  elif n[0]=='product':
   v={0:Q(n[1],n[2])};u={}
   for j in n[3]:
    u=add({k:c*val[j].get(0,Q(0)) for k,c in u.items()},unk[j],v.get(0,Q(0)));v=times(v,val[j])
  else:raise AssertionError(n)
  val[i]=v;unk[i]=u
 native=[{r:rec(v) for r,v in col} for col in D['native']['source_columns128_rows48_sparse_B11_effective_second']]
 assert D['native']['valid_for_actual_zero_origin_source_combinations_only_by_L2_cancellation']
 results=[];ugrad=[];contributions=[]
 for ch,f in enumerate(func):
  z=cor[ch];unknown={};terms=[]
  for k,c in f.items():
   b=allowed[ms2[k[-1]]]
   if k[0]==2:q=native[k[1]].get(k[2],{}).get(b,Q(0));uu={}
   elif k[0]==1:
    assert k[1] in cv;q=cv[k[1]].get(b,Q(0));uu={}
   else:q=val[k[1]].get(b,Q(0));uu={a:v for (bb,a),v in unk[k[1]].items() if bb==b}
   z+=c*q;unknown=add(unknown,uu,c);terms.append([list(k),str(c),str(q),str(c*q)])
  assert not unknown;results.append(z);ugrad.append(enc(unknown));contributions.append(terms)
 assert results==[Q(-251099375,1051776),Q(0)]
 assert list(map(str,results))==D['values']['all119_two_actual_cubic_relation_values'][rnum]['two_actual_cubic_relation_values']
 support=sum(map(len,Pi))+sum(map(len,kap))+sum(map(len,func))+sum(len(v)+len(unk[i]) for i,v in val.items())+sum(map(len,contributions))
 assert support<150000
 ev={'source_manifest':inputs,'independent_source_method':'GB39 reverse-polynomial quotient + direct derivative/Euler scattering + forward typed-root B11 arithmetic with symbolic unknown atom2',
  'canonical_Pi2_600_rows_independent':list(map(enc,Pi)),'source56_Sym3_adjoint_rows_independent':list(map(enc,kap)),
  'all_actual_source_combinations_and_Xp1_a2_signs_matched':True,'witness_index':71,
  'witness_two_actual_cubic_values':list(map(str,results)),'witness_known_pivot_correction':list(map(str,cor)),
  'witness_root_functionals_independent':list(map(enc,func)),'witness_each_scalar_root_contribution':contributions,
  'all_symbolic_uncached_atom_second_coefficients_cancel':True,'division_second_values_never_assigned_actual_zero':True,
  'actual_cached108_source_values_dependency':'independent EXISTING_CACHE_TRANSPORT_REVIEW, same node IDs/full57 inverse/C0 frame; cached source units already verified, no nonunit redivision',
  'native_effective_values_used_after_separate_old_and_native_zero_origin_checks':True,
  'new_true_cubic_initial_relation_nonzero_mod_known_quadratic_times_variables':True,
  'actual_gr3_T_dimension_upper_bound':1,'actual_HF3_or_total_length_computed':False,
  'all119_value_rank_claim_repeated':False,'only_one_nonzero_witness_needed':True,
  'no_new_source_jet_IFTsolve_rank_search_or_348_matrix':True,'expanded_QQ_support':support}
 ep=out/'independent-cubic-witness71-source-and-channel.json';blob=json.dumps(ev,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'own_GB39_Pi2_and_forward_Euler_source_adjoint_accepted':True,
  'own_actual_source_combination_and_Xp1_a2_correction_accepted':True,
  'actual_nonzero_cubic_witness71_accepted':['-251099375/1051776','0'],
  'actual_gr3_T_dimension_upper_bound_accepted':1,'actual_HF3_or_total_length_not_claimed':True,
  'full119_values_rank_not_recomputed':True,'input_hashes':inputs,'code_sha256':sha(Path(__file__)),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
