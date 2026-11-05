#!/usr/bin/env python3
"""Independent source-fourth adjoint replay; no source coefficient solve."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,json,resource,signal,time
from ambient_B18_T_gdiff_second_univariate5_review_20261104_v63 import numeric_eval
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63'
STEM=A+'_REVIEW'
PATHS={
 'primary':f'reports/{A}/report.json','answer':f'reports/{A}/actual-nine-fourth-source-functionals.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'cache3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
 'cache2':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
 'cache2_peer':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW/report.json',
 'CA6':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63/CA-source6-and-regular-quotient-third.json',
 'own_CA6_peer':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63_REVIEW/report.json',
 'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
 'highfold':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/phase4-all619-high-old-source-sensitivities.json',
 'highchain':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/phase1-source-closure-and-bindings.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'atom':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,n,x):
 (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
def plus(a,b,c=Q(1)):
 d=dict(a)
 for k,v in b.items():
  z=d.get(k,Q(0))+c*v
  if z:d[k]=z
  else:d.pop(k,None)
 return d
def scaled(a,c):return {k:c*v for k,v in a.items() if c*v}
def forward(nodes,roots,bound,callback):
 """Independent sparse value/derivative dictionaries with last-use release."""
 seen=set();todo=list(roots)
 def children(i):
  if i in bound:return []
  if i in callback:return callback[i][1]
  n=nodes[i]
  assert n[0] in ['c','a','m','n','i'],('unbound_source',i,n)
  return n[1:] if n[0]!='c' else []
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i);todo.extend(children(i))
 uses=collections.Counter(j for i in seen for j in children(i));uses.update(roots)
 vals={};grads={};live=0;peak=0;invchecks=0
 for i in sorted(seen):
  n=nodes[i];o=n[0]
  if i in bound:v,g=bound[i]
  elif i in callback:
   v,ys,ps=callback[i];g={}
   for y,p in zip(ys,ps):g=plus(g,grads[y],p)
  elif o=='c':v=Q(n[1],n[2]);g={}
  elif o=='a':v=vals[n[1]]+vals[n[2]];g=plus(grads[n[1]],grads[n[2]])
  elif o=='n':v=-vals[n[1]];g=scaled(grads[n[1]],Q(-1))
  elif o=='m':v=vals[n[1]]*vals[n[2]];g=plus(scaled(grads[n[1]],vals[n[2]]),grads[n[2]],vals[n[1]])
  elif o=='i':
   assert vals[n[1]];v=1/vals[n[1]];assert v*vals[n[1]]==1
   g=scaled(grads[n[1]],-v*v);invchecks+=1
  else:raise AssertionError((i,n))
  vals[i]=v;grads[i]=g;live+=1+len(g);peak=max(peak,live);assert peak<125000
  for j in children(i):
   uses[j]-=1
   if not uses[j]:live-=1+len(grads[j]);del vals[j];del grads[j]
 return {i:(vals[i],grads[i]) for i in roots},len(seen),peak,invchecks
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'INDEPENDENT_SOURCE_ADJOINT_POSTVERIFICATION_INPUTS','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 for k in ['primary','joint_peer','cache2_peer','own_CA6_peer']:assert D[k]['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 J=D['joint'];P=D['answer'];unknown=J['all105_unknown_original_node_order'];eq=J['all105_equation_original_node_order'];ui={i:k for k,i in enumerate(unknown)}
 assert unknown==P['source_unknown_order105'] and eq==P['source_equation_order105']
 # Reconstruct the r=s=T restriction from the original bivariate source3.
 mon=J['monomials'];bound={}
 for i,coeffs in J['all105_unknown_jets']:
  a=[Q(0)]*5
  for m,c in zip(mon,coeffs):a[sum(m)]+=Q(c)
  bound[i]=tuple(a)
 assert list(bound)==unknown
 assert [[i,list(map(str,bound[i][:4]))] for i in unknown]==D['cache3']['actual_common105_rline_through3']
 mons=[tuple(m) for m in D['cache2']['base20_monomials_original8']];ix={m:i for i,m in enumerate(mons)}
 for i,coeffs in D['cache2']['all105_actual_source2_values']:
  a={j:Q(c) for j,c in coeffs};assert bound[i][:3]==tuple(a.get(k,Q(0)) for k in [0,ix[(7,)],ix[(7,7)]])
 # Independent full source substitution, with every unknown fourth zero.
 src=D['source']['DAG_scalar_nodes'];ev,units=numeric_eval(src,eq,bound,5)
 expected=[(Q(0),Q(k in [55,56]),Q(0),Q(0)) for k in range(105)]
 assert all(ev[i][:4]==a for i,a in zip(eq,expected))
 residual=[ev[i][4] for i in eq];assert list(map(str,residual))==P['known_F4_residual_with_unknown_fourths_zero']
 # Reconstruct the fixed Jacobian by sparse source AD, not its saved inverse.
 base={i:(bound[i][0],{('physical',i):Q(1)}) for i in unknown}
 ad,closure,peak,iv=forward(src,eq,base,{})
 jac=[[ad[e][1].get(('physical',i),Q(0)) for i in unknown] for e in eq]
 assert jac==[list(map(Q,p)) for p in J['fixed105_Jacobian']]
 assert all(ad[e][0]==0 for e in eq)
 targets=P['target_nodes9'];assert targets==[29,20538,20555,20556,20571,20577,20579,20581,20583]
 rho=[list(map(Q,p)) for p in P['rho9_by105']]
 for t,row in zip(targets,rho):
  assert row==list(map(Q,J['fixed105_Jacobian_inverse'][ui[t]]))
  assert all(sum(row[k]*jac[k][j] for k in range(105))==Q(j==ui[t]) for j in range(105))
 result=[-sum(a*b for a,b in zip(row,residual)) for row in rho]
 assert [[i,str(a)] for i,a in zip(targets,result)]==P['actual_nine_fourth_values']
 own={i:list(map(Q,p)) for i,p in D['CA6']['all57_source_jets_through6']}
 cross=[i for i in targets if i in own];assert cross==[29] and own[29][4]==result[0]
 dump(out,'phase.json',{'phase':'ALL105_SOURCE_EQUATIONS_AND945_ADJOINT_IDENTITIES_PASS','targets9':targets})
 # Bind every fourth-order sensitivity to physical original105 leaves.
 # Recompute all619 fourth sensitivities with a fresh sparse source AD;
 # copiedCB57 has no independent high variable, only psi_1*dYB.
 cnode=D['C0']['DAG_scalar_nodes'];assert cnode[:len(src)]==src
 C=D['C0'];fj=D['FJ'];opaque={}
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in fj[key].items():opaque[int(i)]=Q(p[0])
 for row in D['atom']['all35_ordered_CA_prefix_firsts']:
  opaque[row['node']]=dict((int(j),Q(c)) for j,c in row['first9_bar']).get(0,Q(0))
 for i,p in D['atom']['all8_g_quotient_overrides']:opaque[i]=dict((int(j),Q(c)) for j,c in p).get(0,Q(0))
 for i,a in opaque.items():
  if i not in base:base[i]=(a,{('honest',i):Q(1)})
 for i in C['actual_C0_basis8_scalar_nodes']:
  if i not in bound:base[i]=(Q(0),{})
 cb=fj['both57_source2']['CB'];om=[tuple(m) for m in fj['source2_monomials']];oi={m:i for i,m in enumerate(om)}
 psi={i:list(map(Q,p)) for i,p in cb['all57_source2jets']};ys=cb['F9_nodes'][7:9];callback={}
 copy_pairs=D['copy']['explicit57_source_copy_bindings'];assert [i for i,j in copy_pairs]==cb['unknown57']
 for old,new in copy_pairs:callback[new]=(psi[old][0],ys,[psi[old][oi[(7,)]],psi[old][oi[(8,)]]])
 assert D['copy']['all57_copied_equations_literal_match'] and D['highchain']['missing_source_leaves']==[]
 proots=sorted({i for row in D['high']['all619'] for i,c in row['rawP55_pure_r4_roots'] if Q(c)})
 gradients,nclosure,npeak,nunit=forward(cnode,proots,base,callback)
 folds=[];seen=set()
 for high,claimed in zip(D['high']['all619'],D['highfold']['all619']):
  g={}
  for root,c in high['rawP55_pure_r4_roots']:g=plus(g,gradients[root][1],Q(c))
  actual=[[list(k),str(v)] for k,v in sorted(g.items())]
  assert high['kernel619_index']==claimed['kernel619_index'] and actual==claimed['old_high_source4_sensitivity']
  assert all(k[0]=='physical' and k[1] in ui for k in g)
  seen.update(k[1] for k in g);folds.append({'kernel619_index':high['kernel619_index'],'physical_source4_sensitivity':actual})
 assert sorted(seen)==targets
 ep=out/'independent-nine-adjoint-source-and-copy-checks.json'
 dump(out,ep.name,{'input_manifest':inputs,'all105_source_rline0_through3_rebuilt_from_joint_bivariate3':True,
  'all105_source_equations_through3_checked_with_fourths_zero':True,'independent_source_F4_residual':list(map(str,residual)),
  'independent_source_J105_nonzero_entries':[[i,j,str(a)] for i,p in enumerate(jac) for j,a in enumerate(p) if a],
  'all945_adjoint_selector_identities_verified':True,'actual_nine_fourth_source_functionals':[[i,str(a)] for i,a in zip(targets,result)],
  'independent_CA6_overlap':cross,'all619_source4_sensitivity_folds':folds,
  'copied57_rule':'psi(0^7,YB7,YB8), derivative at0 is psi7*dYB7+psi8*dYB8; source57 unit uniqueness pinned by accepted nonlinear source2/copy peer',
  'all57_literal_copy_node_bindings':copy_pairs,'source_equation_closure':closure,'high_AD_closure':nclosure,
  'peak_live_sparse_QQ':max(peak,npeak),'unit_inverse_checks':units+iv+nunit,
  'no_source_unknown_solve_or_full105fourth_vector':True,
  'scope':'NEW higher source functionals primary independently accepted; no source3/gauge3, full619quartics, HF4, or exactlength assertion'})
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'actual_nine_source_fourth_functionals_independently_accepted':True,
  'all105_equations_and945_adjoint_identities_and_coherent_copy619_folds_accepted':True,
  'new_source_solve_by_reviewer':False,'full105_fourth_vector_constructed':False,'actual619_quartics_or_HF4':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
