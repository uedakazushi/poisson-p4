#!/usr/bin/env python3
"""Independent one-witness580 source/unknown-polynomial/native audit."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time
from ambient_B18_T_joint135_selected_functionals_review_20261104_v63 import Cover
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_T_COMBINED619_QUARTIC_CLOSURE_20261105_v63'
STEM='AMBIENT_B18_T_COMBINED619_QUARTIC_WITNESS_20261105_v63_REVIEW'
P={
 'primary':f'reports/{A}/report.json','witness':f'reports/{A}/actual-nonzero-quartic-witness.json',
 'cases':f'reports/{A}/phase3-combined619-evaluated-cases.json','record':f'reports/{A}/phase1-C0-symbolic-quartic-roots.json',
 'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
 'high_peer':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63_REVIEW/report.json',
 'bridge':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/actual619-quartic-source-order-and-unknown-gates.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'source2':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
 'source2_peer':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW/report.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'own4':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'atom':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
 'g3':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63/actual-eight-g-third-values.json',
 'g3_peer':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63_REVIEW/report.json',
 'CA6':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63/CA-source6-and-regular-quotient-third.json',
 'phi4':'reports/AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63/actual-nine-fourth-source-functionals.json',
 'phi4_peer':'reports/AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63_REVIEW/report.json',
 'N2':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
 'N2_peer':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW/report.json',
 'response':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase2-known-third-RHS-and-proper-response.json',
 'implicit':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase4-selected-joint135-third-adjoint-values.json',
 'implicit_peer':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63_REVIEW/report.json',
 'Op':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
 'lower':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/all619-lower-contribution-values.json',
 'lower_peer':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63_REVIEW/report.json',
 'theta':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json',
 'c4_peer':'reports/AMBIENT_B18_T_CORRECTED168_QUARTIC_VALUES_20261104_v63_REVIEW/report.json',
 'tail':'reports/AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(out,n,x):
 (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
def put(out,k,c):
 a=out.get(k,Q(0))+c
 if a:out[k]=a
 else:out.pop(k,None)
class Jets:
 """Flat exact polynomial ring in base monomials and named h2 symbols."""
 def __init__(self,mons):self.mons=set(mons)
 def c(self,c):return {((),()):Q(c)} if c else {}
 def add(self,a,b,c=Q(1)):
  x=dict(a)
  for k,v in b.items():put(x,k,c*v)
  return x
 def scale(self,a,c):return {k:v*c for k,v in a.items() if v*c}
 def mul(self,a,b):
  out={}
  for (m,p),c in a.items():
   for (n,q),d in b.items():
    mon=tuple(sorted(m+n))
    if mon not in self.mons:continue
    pars=tuple(sorted(p+q));assert len(pars)<=2
    put(out,(mon,pars),c*d)
  return out
 def inv(self,a):
  c=a.get(((),()),Q(0));assert c and all(m or not p for m,p in a)
  h=self.scale(self.add(a,self.c(-c)),-1/c);v=t=self.c(1)
  for _ in range(4):t=self.mul(t,h);v=self.add(v,t)
  v=self.scale(v,1/c);assert self.mul(a,v)==self.c(1);return v
 def piece(self,a,m):return {p:c for (n,p),c in a.items() if n==m}
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(out,'phase.json',{'phase':'ONE580_LITERAL_SOURCE_INPUTS','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in ['primary','high_peer','source2_peer','g3_peer','phi4_peer','N2_peer','implicit_peer','lower_peer','c4_peer','tail'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 high=D['high']['all619'][580];w=D['witness'];bridge=D['bridge']['all619'][580]
 c0=[[406,'-8/3'],[647,'60'],[710,'1']]
 assert high['kernel619_index']==w['kernel619_index']==bridge['constant619_index']==580
 assert high['actual_original768_constant_combination']==bridge['actual_original768_constant_combination_c0']==w['actual_original768_constant_combination']==c0
 assert bridge['full56_Xc1_equals_Xp0_a1'] and bridge['old_and_native_origin_parts_separately_zero']
 assert D['bridge']['Theta4_source']['c4']==str(Q(4565,68832)) if 'c4' in D['bridge']['Theta4_source'] else True
 nodes=D['C0']['DAG_scalar_nodes'];ops=D['chain']['operator_scalar_nodes'];mons=[tuple(m) for m in D['record']['R0_monomials']];J=Jets(mons)
 oldmons=[tuple(m) for m in D['source2']['base20_monomials_original8']]
 bind={i:{(oldmons[int(k)],()):Q(c) for k,c in p if Q(c)} for i,p in D['source2']['all105_actual_source2_values']}
 for i,c in D['phi4']['actual_nine_fourth_values']:
  if Q(c):bind[i][((7,7,7,7),())]=Q(c)
 g2={i:Q(c) for i,c in D['g2']['actual_eight_regular_g_second_values']};g3={i:Q(c) for i,c in D['g3']['actual_eight_common_regular_g_third_values']}
 opaque={}
 def first(p):
  p=list(map(Q,p));x=J.c(p[0])
  for k in [0,1,4,5,7]:
   c=p[8]+p[9] if k==7 else p[k+1]
   if c:x[((k,),())]=c
  return x
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in D['FJ'][key].items():opaque[int(i)]=first(p)
 for row in D['atom']['all35_ordered_CA_prefix_firsts']:
  p={int(i):Q(c) for i,c in row['first9_bar']};x=J.c(p.get(0,Q(0)))
  for k in [0,1,4,5,7]:
   if p.get(k+1,Q(0)):x[((k,),())]=p[k+1]
  opaque[row['node']]=x
 for i,p0 in D['atom']['all8_g_quotient_overrides']:
  p={int(i):Q(c) for i,c in p0};x=J.c(p.get(0,Q(0)))
  for k in [0,1,4,5,7]:
   if p.get(k+1,Q(0)):x[((k,),())]=p[k+1]
  opaque[i]=x
 coords={}
 for k,i in enumerate(D['C0']['actual_C0_basis8_scalar_nodes']):
  if k<7:coords[i]={((k,),()):Q(1)} if (k,) in J.mons else {}
  else:
   g=list(map(Q,D['CA6']['g_through6']));assert g[3]==Q(125,24) and not g[4];coords[i]={((7,7,7),()):g[3]}
 cb=D['FJ']['both57_source2']['CB'];ys=cb['F9_nodes'][7:9];callbacks={j:i for i,j in D['copy']['explicit57_source_copy_bindings']}
 own={i:list(map(Q,p)) for i,p in D['own4']['families']['CB']['all57_source4_jets']};om=D['own4']['monomials_bivariate_through4']
 cr=set(i for field in ['rawP55_pure_r4_roots','Xp55_pure_r3_times_a1r_roots'] for i,c in high[field]);tr=set()
 for field in ['actual_X3_high_order_source_functional','actual_Xp2_times_a1_source_functional']:
  for key,c in high[field]:
   if key[0]==0:tr.add(key[1])
   elif key[0]==1:cr.add(key[1])
 tc=set();todo=list(tr)
 while todo:
  i=todo.pop()
  if i in tc:continue
  tc.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  elif n[0]=='actual_C0_scalar_function':cr.add(n[1])
  else:assert n[0]=='c'
 def children(i):
  if i in bind or i in coords or i in opaque:return []
  if i in callbacks:return ys
  n=nodes[i];assert n[0] in ['c','a','m','n','i'],('unbound_source',i,n)
  return n[1:] if n[0]!='c' else []
 cc=set();todo=list(cr)
 while todo:
  i=todo.pop()
  if i in cc:continue
  cc.add(i);todo.extend(children(i))
 uses=collections.Counter(j for i in cc for j in children(i));uses.update(cr);cv={};peak=0;live=0;parameters=set()
 for i in sorted(cc):
  n=nodes[i];op=n[0]
  if i in bind:v=bind[i]
  elif i in coords:v=coords[i]
  elif i in opaque:
   v=dict(opaque[i])
   for m in mons:
    if len(m)!=2:continue
    if m==(7,7) and i in g2:
     if g2[i]:v[(m,())]=g2[i]
    else:
     param=(i,m);parameters.add(param);v[(m,(param,))]=Q(1)
   if i in g3 and g3[i]:v[((7,7,7),())]=g3[i]
  elif i in callbacks:
   powers=[[J.c(1)] for _ in range(2)]
   for k in range(2):
    for _ in range(4):powers[k].append(J.mul(powers[k][-1],cv[ys[k]]))
   v={}
   for m,c in zip(om,own[callbacks[i]]):
    if c:v=J.add(v,J.mul(powers[0][m[0]],powers[1][m[1]]),c)
  elif op=='c':v=J.c(Q(n[1],n[2]))
  elif op=='a':v=J.add(cv[n[1]],cv[n[2]])
  elif op=='n':v=J.scale(cv[n[1]],-1)
  elif op=='m':v=J.mul(cv[n[1]],cv[n[2]])
  elif op=='i':v=J.inv(cv[n[1]])
  else:raise AssertionError(n)
  cv[i]=v;live+=len(v);peak=max(peak,live);assert peak<100000
  for j in children(i):
   uses[j]-=1
   if not uses[j]:live-=len(cv[j]);del cv[j]
 usesT=collections.Counter()
 for i in tc:
  n=ops[i]
  if n[0]=='sum':usesT.update(j for j,a,b in n[1])
  elif n[0]=='product':usesT.update(n[3])
 usesT.update(tr);tv={}
 for i in sorted(tc):
  n=ops[i];op=n[0]
  if op=='c':v=J.c(Q(n[1],n[2]))
  elif op=='actual_C0_scalar_function':v=cv[n[1]]
  elif op=='sum':
   v={}
   for j,a,b in n[1]:v=J.add(v,tv[j],Q(a,b))
  elif op=='product':
   v=J.c(Q(n[1],n[2]))
   for j in n[3]:v=J.mul(v,tv[j])
  else:raise AssertionError(n)
  tv[i]=v
  for j in [j for j,a,b in n[1]] if op=='sum' else n[3] if op=='product' else []:
   usesT[j]-=1
   if not usesT[j]:del tv[j]
 polynomial={}
 def addpoly(a,b,c):
  z=dict(a)
  for k,v in b.items():put(z,k,c*v)
  return z
 for root,c in high['rawP55_pure_r4_roots']:polynomial=addpoly(polynomial,J.piece(cv[root],(7,7,7,7)),Q(c))
 for root,c in high['Xp55_pure_r3_times_a1r_roots']:polynomial=addpoly(polynomial,J.piece(cv[root],(7,7,7)),Q(c))
 native_keys=[];m2=list(itertools.combinations_with_replacement(range(8),2));m3=list(itertools.combinations_with_replacement(range(8),3))
 for field in ['actual_X3_high_order_source_functional','actual_Xp2_times_a1_source_functional']:
  for key,c in high[field]:
   if key[0]==2:native_keys.append((tuple(key),Q(c)));continue
   mon=(m3 if key[-2]==3 else m2)[key[-1]]
   polynomial=addpoly(polynomial,J.piece((tv if key[0]==0 else cv)[key[1]],mon),Q(c))
 assert polynomial=={():Q(45300359575,158588928)}
 assert not any(key for key in polynomial)
 dump(out,'phase1-independent-honest-LIN-QUAD-cancellation.json',{'case580_source_closure':len(cc),'case580_typed_closure':len(tc),
  'actual_honest_second_parameters':[[i,list(m)] for i,m in sorted(parameters)],'all_linear_AND_quadratic_parameter_coefficients_zero':True,
  'old_known_piece':str(polynomial[()]),'peak_flat_QQ':peak,'coherent_copy57_own4_composed':True})
 # Independent native known response in the same monic algebra, with N2
 # limited to precisely the sufficient B11 coefficients (at most one z).
 mon25=mons[:25];ix={m:i for i,m in enumerate(mon25)};nu=[{i:Q(c) for i,c in enumerate(p) if Q(c)} for p in D['response']['known_nu2_base25']]
 ring=Cover(mon25,nu);gg=[{(i,b):Q(p[4*i+b]) for i in range(25) for b in range(4) if Q(p[4*i+b])} for p in D['response']['known_properG8_C100']]
 b11=[tuple(i for i,n in enumerate(m) for _ in range(n)) for m in D['B11']['B11_basis_monomials8']]
 N=[[{(ix[m],b):Q(p[4*j+b]) for j,m in enumerate(b11) for b in range(4) if Q(p[4*j+b])} for p in row] for row in D['N2']['N_C_B44']]
 Sn=[[{} for j in range(32)] for i in range(12)];OP=D['Op']['native_S12_by32']
 for i,j,k,s in OP['koszul28_in_first8_rows']:Sn[i][j]=ring.scale(gg[k],Q(s))
 for i,j,c in OP['free4_zero_source_columns']:Sn[i][j]=ring.c(Q(c))
 ncols={};native=Q(0);native_evidence=[]
 for key,c in native_keys:
  _,j,outrow,d,k=key
  if j not in ncols:
   jn,b=divmod(j,4);vec=[]
   for o in range(12):
    v={}
    for a in range(12):v=ring.add(v,ring.mul(N[a][o],Sn[a][jn]))
    vec.append(ring.mul(v,{(0,b):Q(1)}))
   ncols[j]=vec
  m=(m3 if d==3 else m2)[k];assert sum(v!=7 for v in m)<=1
  row,b=divmod(outrow,4);value=ncols[j][row].get((ix[m],b),Q(0));native+=c*value
  native_evidence.append([list(key),str(value),str(c)])
 assert native==Q(-491134351145455,17127604224)
 implicit=sum(map(Q,D['implicit']['all619_actual_implicit_third_functional_values_by5'][580]));lower=Q(D['lower']['all619_actual_lower_contribution_values'][580])
 assert implicit==Q(97525955,9911808) and lower==Q(489966482601115,17127604224)
 total=polynomial[()]+native+implicit+lower;assert total==Q(40553074375,178412544) and total
 assert [str(polynomial[()]),str(native),str(implicit),str(lower),str(total)]==[w[k] if k!='old_known_part_with_honest_second_polynomial' else w[k][0][1] for k in ['old_known_part_with_honest_second_polynomial','native_known_response','selected135_implicit_contribution','actual_common_D_lower','actual_corrected_pure_quartic_if_bound']]
 ep=out/'independent-one580-complete-source-and-nonzero-quartic.json'
 dump(out,ep.name,{'source_manifest':inputs,'actual_original768_c0':c0,'case580':580,
  'honest_second_linear_and_quadratic_symbols_all_cancel':True,'independent_native_known_values':native_evidence,
  'old_known':str(polynomial[()]),'native_q_known':str(native),'selected_implicit_minus_rho_f':str(implicit),'common_D_lower':str(lower),
  'actual_corrected_pure_e4_r4':str(total),'nonzero':True,
  'source_inputs':'literal68 Hadamards/oldf63 coefficients, true opaque germs, actual ownCB57 zero7 section and all provided higher values retained',
  'normalization':'accepted exact619 source bridge, Theta4 with b3=-42215/3456,c4=4565/68832, normalized pivots retained',
  'proof':'sole gr4 candidate is e4*r4; nonzero actual relation kills it, m4T=m5T, finite generated Nakayama gives m4T=0',
  'actual_HF':[5,3,2,1,0],'actual_T_total_length':11,'all619_evaluated_by_reviewer':False,
  'new_source_or_implicit_coefficient_solve_by_reviewer':False,'nilcube_or_H0_exponent_claim':False})
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'actual_case580_nonzero_quartic_independently_accepted':True,
  'actual_T_total_length':11,'actual_T_HF':[5,3,2,1,0],'no_all619_or_newsource_solve_claim':True,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
