#!/usr/bin/env python3
"""Independent typed-chart and ONE actual fixed Schur-entry audit."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import argparse,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
CHART='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63'
ENTRY='S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63'
STEM='S61_NONSQUARE_PREPOWER_DIRECT_SCHUR_ENTRY_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def poly(p):return {tuple(m):Q(c) for m,c in p if Q(c)}
def main(chart_only=False):
 global STEM
 if chart_only:STEM=CHART+'_REVIEW'
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 work=ROOT/'work/prePOWER-direct-one-entry-independent-review';work.mkdir(parents=True,exist_ok=True)
 hashes={};cache={}
 def load(p):
  hashes[p]=sha(p)
  if p not in cache:cache[p]=read(p)
  return cache[p]
 def get(stem,name='certificate.json'):return load('reports/'+stem+'/'+name)
 def save(name,d):(out/name).write_text(json.dumps(d,indent=2)+'\n')
 for st,want in [(CHART,'264049acfbbc7abeae10dfd5994b03967b599d5a0ea18f0a10c07efb4eb57b11'),(ENTRY,'5c4c200fef2e029a4a5a3992895f327591801b2475c737b3e670f920cc885f9e')]:
  d=get(st,'report.json');assert sha('reports/'+st+'/report.json')==want
  assert d.get('all_passed',d.get('all_structural_source_and_type_checks_passed'))
  for group in ['input_sha256','source_sha256']:
   for p,h in d.get(group,{}).items():assert sha(p)==h,p;hashes[p]=h
  for n,h in d['evidence_sha256'].items():p='reports/'+st+'/'+n;assert sha(p)==h,p;hashes[p]=h
 for st in ['S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION','S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION','S61_NONSQUARE_CUBIC_ZW25_UNITS','S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT','S61_NONSQUARE_RAW_H_GENERIC_RANK','S61_NONSQUARE_PREPOWER_H22_FIXED_MINORS']:
  d=get(st+'_20261104_v63_REVIEW','report.json');assert d['all_passed']
 assert hashes['reports/S61_NONSQUARE_PREPOWER_H22_FIXED_MINORS_20261104_v63_REVIEW/report.json']=='f2d235534778c7bc868c6e013ee85c9bd1f81e45c02575fdaf4968178c9b4037'
 q=get('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63')
 pair=get('S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63')
 cp=get('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json')
 pd=get('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63')
 h=get('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63')
 rank=get('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63')
 h9=get('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json')
 old=get('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json')
 strip=get('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json')
 chart=get(CHART,'operator.json');N=chart['nodes'];E=chart['exports'];ctx=chart['context']
 roots=q['remaining15x22_scalar_roots'];cols=q['new22_target_original32_columns'];labels=q['new15_source_original25_labels']
 assert len(N)==151 and labels==[9]+h['final14_source_labels'] and cols[17:]==[27,28,29,30,31]
 pr=[0,1,3,6];pc=[22,11,13,25];hr=[2,4,7,13,10];hi=[1,3,6,12,9];hcols=[27,31,28,30,29]
 piv=pr+hr;remaining=[i for i in range(15) if i not in piv];targetpiv=[cols.index(c) for c in pc+hcols];remainingc=[j for j in range(22) if j not in targetpiv]
 assert ctx['P4_current15_rows']==pr and ctx['P4_original32_targets']==pc
 assert ctx['H5_current15_rows']==hr and ctx['H5_compact14_rows']==hi and ctx['H5_original32_targets']==hcols
 assert ctx['remaining6_current15']==remaining==[5,8,9,11,12,14]
 assert ctx['remaining13_original32_seed_labels']==[cols[j] for j in remainingc]==[1,2,3,4,6,7,8,9,14,15,17,18,19]
 assert ctx['pivot9_source_order']==piv and ctx['pivot9_target_current22_order']==targetpiv
 # Independent complete typed-node reader, including all foreign pointers.
 for i,n in enumerate(N):
  assert n['id']==i and len(n['shape'])==2 and all(x>0 for x in n['shape'])
  assert all(0<=j<i for j in n['args']);op=n['op'];s=n['shape'];a=n['args']
  if op=='external_frozen_operator':
   r=n['reference'];d=load(r['path']);assert hashes[r['path']]==r['sha256']
   for k in r['json_pointer'].strip('/').split('/'):d=d[int(k)] if isinstance(d,list) else d[k]
  elif op=='matrix_product':assert N[a[0]]['shape'][1]==N[a[1]]['shape'][0] and s==[N[a[0]]['shape'][0],N[a[1]]['shape'][1]]
  elif op in ('negative','transpose'):
   assert s==(N[a[0]]['shape'] if op=='negative' else N[a[0]]['shape'][::-1])
  elif op=='difference':assert s==N[a[0]]['shape']==N[a[1]]['shape']
  elif op=='central_scalar_multiple':assert N[a[0]]['shape']==[1,1] and s==N[a[1]]['shape']
  elif op=='submatrix':assert s==[len(n['rows']),len(n['columns'])] and all(0<=j<N[a[0]]['shape'][0] for j in n['rows']) and all(0<=j<N[a[0]]['shape'][1] for j in n['columns'])
  elif op=='universal_polynomial_cofactor_adjugate':assert s==N[a[0]]['shape'] and s[0]==s[1] and n['expanded']==False
  elif op=='certified_chart_unit_inverse':assert s==N[a[0]]['shape']==[1,1] and n['expanded']==False
  elif op=='block_matrix':
   br=n['block_rows'];heights=[N[row[0]]['shape'][0] for row in br];widths=[N[k]['shape'][1] for k in br[0]]
   assert s==[sum(heights),sum(widths)] and a==[j for row in br for j in row]
   assert all(N[j]['shape']==[heights[r],widths[c]] for r,row in enumerate(br) for c,j in enumerate(row))
  elif op=='diagonal':assert s[0]==s[1]==len(n['diagonal']) and a==n['diagonal'] and all(N[j]['shape']==[1,1] for j in a)
  elif op in ('row_permutation','column_permutation'):assert s[0]==s[1] and sorted(n['order'])==list(range(s[0]))
  elif op=='coordinate_row_selector':assert s[0]==len(n['rows']) and all(0<=j<s[1] for j in n['rows'])
  elif op=='coordinate_row_embedding':assert s[0]==len(n['columns']) and all(0<=j<s[1] for j in n['columns'])
  elif op=='accepted_forward_Gaussian_source_step':assert s==[25,25] and len(a)==2 and n['pivot_row'] not in n['targets']
  elif op=='identity':assert s[0]==s[1]
  else:assert op in ('zero','one'),op
 assert N[E['row_permutation']]['order']==piv+remaining
 assert N[E['column_permutation']]['order']==targetpiv+remainingc
 diag=N[E['hybrid_row_unit']]['diagonal'];assert len(set(diag[i] for i in range(15) if i not in hr))==1
 assert all(N[diag[i]]['op']=='certified_chart_unit_inverse' for i in hr)
 td=N[E['heavy_column_unit']]['diagonal'];assert len(set(td[:21]))==1 and N[td[21]]['op']=='certified_chart_unit_inverse'
 inverse_refs=[]
 for n in N:
  if n['op']=='certified_chart_unit_inverse':
   r=N[n['args'][0]]['reference'];p=r['json_pointer']
   assert (r['path'].endswith('certificate.json') and (p.startswith('/known_row_units_rho14/') or p=='/known_column_units/4' or p=='/nonzero_minor_of14x5')) or (r['path'].endswith('checkpoint.json') and p.startswith('/all_completed_unit_steps/')) or (r['path'].endswith('native.json') and p=='/exact_unscaled_pure4_minor')
   inverse_refs.append([r['path'],p])
 assert not any('h22' in p.lower() for _,p in inverse_refs)
 assert N[E['residual6x13']]['op']=='difference' and N[E['residual6x13']]['shape']==[6,13]
 assert N[E['surviving6_source_original25']]['op']=='matrix_product' and N[E['surviving6_source_original25']]['args']==[E['surviving6_source_current15'],E['source_current15_to_original25']]
 # Source structural identities and every literal heavy-root correspondence.
 QN=q['scalar_DAG']['nodes'];resolved={}
 def raw(n):
  if n[0]=='checkpoint':return cp['matrix25x32'][n[1]][n[2]]
  if n[0]=='accepted_pair_second':return pd['Lframe_dual6x2'][n[1]][1]
  assert n[0]=='new_unit_atom';return q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
 def resolve(i):
  if i in resolved:return resolved[i]
  n=QN[i];o=n[0]
  if o=='z':v=0
  elif o=='one':v=1
  elif o in ('checkpoint','accepted_pair_second','new_unit_atom'):v=tuple(n) if raw(n)[0] else 0
  else:
   a,b=resolve(n[1]),resolve(n[2])
   if o=='m':v=0 if not a or not b else b if a==1 else a if b==1 else ('m',a,b)
   elif o=='s':v=0 if a==b else a if not b else ('s',a,b)
   else:assert o=='a';v=b if not a else a if not b else ('a',a,b)
  resolved[i]=v;return v
 assert all(resolve(roots[i][j])==0 for i in pr for j in range(17,22))
 for lab,col,trace in h['all70_literal_external_root_traces']:
  assert resolve(roots[labels.index(lab)][cols.index(col)])==(0 if trace==[0] else tuple(trace[-1]))
 compact=h['raw_compact14x5'];assert all(not compact[i][j] for i in hi[:2] for j in [1,3,2])
 assert [rank['input14_remaining_rows'][i] for i in hi]==rank['nonzero_minor_of14x5_rows']
 assert [rank['input5_columns'][j] for j in [0,4,1,3,2]]==rank['nonzero_minor_of14x5_cols']
 for i,br in zip(hi[2:],rank['selected_rows7']):
  for j,bc in zip([1,3,2],rank['selected_cols3']):assert poly(compact[i][j])==poly(rank['actual_block7x3'][br][bc])
 assert strip['named_unit_powers_a_d_F_B_C_J']==[4,0,0,0,3,0]
 assert ctx['Q9_eight_other_pure_tails']==[{'original32':j,'actual_root':roots[0][cols.index(j)],'remaining13_position':[cols[k] for k in remainingc].index(j)} for j in [1,2,3,4,6,7,8,9]]
 assert labels==q['old16_source_original25_rows'][1:] and q['old16_source_original25_rows']==pair['old17_source_original25_rows'][1:]==cp['remaining_rows'][1:]
 op=get('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63','operator-and-cost.json');D=op['nodes'];native=get(ENTRY,'native.json')
 assert len(D)==768 and op['output']==767 and op['current15_row']==5 and op['target_original32']==15
 assert op['source_current15_rows']==[5,2,4] and op['source_coefficients3']==[33,35,37]
 assert native['numerator_terms']==7022 and native['known_a_clearing_power']==0 and native['B_denominator_power']==0
 units=[[[[1,0,0],'1']],[[[0,0,1],'1']],[[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']],[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']],[[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']],[[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']]]
 names=['a','d','F','B','C','J'];oldvals={v['old_DAG_root']:v['Laurent_polynomial'] for v in old['node_values']}
 assert set(h9['source_prePOWER_DAG_roots'])==set(roots[0][cols.index(j)] for j in h9['original_target_labels'])
 # Canonical expression fingerprints bind all imported old Q2L roots. They
 # keep actual cached roots atomic, already coefficient-verified by f2d235.
 fcache={};oldfp={};localfp=[];atompolys={};atommanifest=[]
 def const(c):return ('QQ',Q(c))
 def combine(o,a,b):
  if o=='mul':
   if a==const(0) or b==const(0):return const(0)
   if a==const(1):return b
   if b==const(1):return a
  elif o=='sub':
   if a==b:return const(0)
   if b==const(0):return a
  else:
   if a==const(0):return b
   if b==const(0):return a
  return (o,a,b)
 h9byold=dict(zip(h9['source_prePOWER_DAG_roots'],h9['source_h9']))
 def fOld(i):
  if i in oldfp:return oldfp[i]
  n=QN[i];o=n[0]
  if i in oldvals:v=('Qroot',i) if oldvals[i] else const(0)
  elif i in h9byold:v=('Qroot',i) if h9byold[i] else const(0)
  elif o=='z':v=const(0)
  elif o=='one':v=const(1)
  elif o in ('checkpoint','accepted_pair_second','new_unit_atom'):v=('Qroot',i) if raw(n)[0] else const(0)
  else:v=combine({'a':'add','s':'sub','m':'mul'}[o],fOld(n[1]),fOld(n[2]))
  oldfp[i]=v;return v
 for i,n in enumerate(D):
  assert n['id']==i and all(j<i for j in n['args']) and n['bound'].get('known_B_denominator_power',0)==0
  o=n['op']
  if o!='atom':
   assert o in ('add','sub','mul') and len(n['args'])==2
   localfp.append(combine(o,localfp[n['args'][0]],localfp[n['args'][1]]));continue
  r=n['reference'];qq=Q(n['QQ_denominator_unit']);ashift=n['a_shift'];assert qq and n['B_denominator_power']==0
  if 'constant' in r:p=[[[0,0,0],r['constant']]];fp=const(r['constant']);assert qq==1 and ashift==0
  elif 'namedFULLUunit' in r:p=units[names.index(r['namedFULLUunit'])];fp=('unit',r['namedFULLUunit']);assert qq==1 and ashift==0
  elif r['path']=='accepted_Q2L_source_atom':
   oldid=r['old_DAG_root'];assert QN[oldid]==r['definition'];num,den=raw(QN[oldid]);assert len(den)==1 and den[0][0][1:]==[0,0]
   assert ashift==-den[0][0][0] and qq==Q(den[0][1]);p=num;fp=fOld(oldid)
  else:
   x=load(r['path']);k=r['key'];assert qq==1 and ashift==0
   if k=='node_values':p=oldvals[r['old_DAG_root']];fp=fOld(r['old_DAG_root'])
   elif k=='source_h9':
    j=h9['original_target_labels'].index(r['target']);p=h9['source_h9'][j];fp=fOld(h9['source_prePOWER_DAG_roots'][j])
   else:
    p=x[k]
    for j in r.get('index',[]):p=p[j]
    fp=('bound_scalar',r['path'],k,tuple(r.get('index',[])))
  atompolys[i]=[[[m[0]+ashift,m[1],m[2]],str(Q(c)/qq)] for m,c in p if Q(c)]
  localfp.append(fp);atommanifest.append({'node':i,'reference':r,'literal_terms':len(p),'a_shift':ashift,'QQ_divisor':str(qq)})
 lookup={x:i for i,x in enumerate(localfp)}
 def P(i,j):return lookup[fOld(roots[i][cols.index(j)])]
 A=[[P(i,j) for j in pc] for i in pr];target=[P(i,15) for i in pr];Ps=[[P(i,j) for j in pc+[15]] for i in [5,2,4]]
 keep=set(sum(A,[])+target+sum(Ps,[])+op['source_coefficients3']+op['Bhat_pivot_and_target_roots']+[767])
 # All polynomial values and input data are streamed; no duplicate 7022-term
 # result is committed. The independent native guard counts real live terms.
 text=[];literal=sum(len(p) for p in atompolys.values())
 def emit(p):
  text.append(str(len(p)))
  for m,c in p:assert len(m)==3 and m[1]>=0 and m[2]>=0;text.append(' '.join(map(str,m))+' '+str(c))
 text.append(f'{literal} 768')
 for i,n in enumerate(D):
  if n['op']=='atom':text.append('p');emit(atompolys[i])
  else:text.append({'add':'a','sub':'s','mul':'m'}[n['op']]+' '+' '.join(map(str,n['args'])))
 text.append(str(len(keep))+' '+' '.join(map(str,sorted(keep))));text.append('3')
 for i in op['source_coefficients3']:text.append(str(i));emit(get(ENTRY,f'source-cofactor-node{i}.json')['polynomial'])
 text.append('767');emit(native['exact_Saux_Laurent_polynomial']);emit(native['cleared_numerator_polynomial'])
 for row in A:text.append(' '.join(map(str,row)))
 text.append(' '.join(map(str,target)))
 for row in Ps:text.append(' '.join(map(str,row)))
 text.append(' '.join(map(str,op['source_coefficients3'])));text.append(' '.join(map(str,op['Bhat_pivot_and_target_roots'])))
 emit(rank['literal_pair_determinant'])
 for i in [4,1,3]:
  for j in [0,4]:emit(compact[i][j])
 for p in units:emit(p)
 rows=[h['known_row_units_rho14'][i] for i in [4,1,3]];powers=[r['a_d_F_B_C_J_powers'] for r in rows]
 lcm=[max(p[j] for p in powers) for j in range(6)];assert op['known_rho_LCM_powers']==lcm
 for p in powers:text.append(' '.join(map(str,p)))
 text.extend(r['constant'] for r in rows);text.append(' '.join(map(str,lcm)))
 ip=work/'input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<32<<20
 save('preflight.json',{'typed_nodes':151,'fixed_source_operations':768,'literal_atom_terms':literal,'retained_nodes':sorted(keep),'polynomial_output_terms':7022,'independent_fixed_border_det_size':5,'no_author_import_solver_or_new_entry':True,'input_bytes':ip.stat().st_size,'live_coefficient_guard':150000})
 if chart_only:
  proof={'all151_typed_nodes_source_bound':True,'all70_actual_heavy_root_traces_checked':True,
   'P4_source_rows':pr,'P4_original32_targets':pc,'H5_source_rows':hr,'H5_original32_targets':hcols,
   'remaining_source_rows':remaining,'remaining13_original32_seed_columns':[cols[j] for j in remainingc],
   'hybrid_only_selected_H5_rows_divided_by_rho':True,'only_last_heavy31_divided_by_gamma':True,
   'H5_inverse':'Hnum=[[delta3 Adj2,0],[-Adj3 C32 Adj2,delta2 Adj3]]; both products with H5 equal Delta5 I5.',
   'P4_inverse':'Direct actual Adj4(P4)/D4, D4=a^4 C^3 R434. No h22 inverse.',
   'K9_inverse':'[[P4inv,0],[-H5inv C54 P4inv,H5inv]], inverse of actual lower-block K9.',
   'complete_equivalence':'For M=[[K,B],[C,D]], L=[[Kinv,0],[-C Kinv,I]], U=[[I,-Kinv B],[0,I]] give diag(I9,D-C Kinv B). Their inverse blocks are [[K,0],[C,I]] and [[I,Kinv B],[0,I]].',
   'source25_transport':'Actual forward Gaussian step product, remaining17 selector, then actual S17/deletion and S16/deletion; remaining source covectors multiply this forward operator.',
   'all8_Q9_extra_pure_tails_retained':True,'inverse_references':inverse_refs,
   'chart':'FULLU D(R434*Delta5), including h22=0 points when these chart units are nonzero.',
   'one_entry_arithmetic_peer':'PENDING: separate bounded native arithmetic stopped at48s; no numeric entry identity accepted here.'}
  save('source-and-block-proof.json',proof)
  note='research/notes/'+STEM+'.md'
  for p in [str(Path(__file__).relative_to(ROOT)),note]:hashes[p]=sha(p)
  result={'id':STEM,'all_passed':True,'scope':'Typed actual direct P4/H5 chart and source/cokernel proof only; ONEentry arithmetic remains pending.',
   'actual_samecoker_shape':[6,13],'all151_types_and_actual_source_indices_verified':True,
   'hybrid_scaling_and_forward_source25_covector_transport_accepted':True,
   'no_h22_inverse_all8_extra_Q9_tails_retained':True,'chart':'FULLU D(R434*Delta5)',
   'entry7022_coefficients_or_residual_rank6_or_full57_claimed':False,
   'input_sha256':hashes,'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},
   'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
  save('report.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['input_sha256','evidence_sha256']}),flush=True);return
 if (work/'native.log').exists() and '48 second guard' in (work/'native.log').read_text():raise RuntimeError('Prior bounded arithmetic stopped; a new worker requires fresh authorization. Use --chart-only for the completed structural review.')
 (work/'phase.json').write_text(json.dumps({'phase':'independent source/type/input checks complete; before arithmetic','input_sha256':sha(str(ip.relative_to(ROOT)))},indent=2)+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=work/'checker';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-started)))
 with (work/'native.log').open('w') as log:subprocess.run([str(exe),str(ip),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,57-(time.monotonic()-started)))
 arithmetic=read('reports/'+STEM+'/native.json');assert arithmetic['all_passed']
 proof={'chart':'FULLU intersect D(R434*Delta5); no h22 inversion and no chart complement removed.',
  'typed_H5_inverse':'H5=[P2,0;C32,B3]; Hnum=[delta3 Adj2,0;-Adj3 C32 Adj2,delta2 Adj3], both products=Delta5 I5.',
  'typed_P4_inverse':'Direct Adj4(P4)/D4 with actual source-bound D4=a^4 C^3 R434; not the h22 Schur3 chart.',
  'complete_current15_equivalence':'K9=[P4,0;C54,H5]; L=[Kinv,0;-C Kinv,I6], U=[I9,-Kinv B;0,I13], giving diag(I9,D-C Kinv B), with the displayed inverses.',
  'source_transport':'Forward source Gaussian steps in their saved chronological order, then actual S17/S16 and deletion selectors. Covectors multiply this forward operator; no transpose or unrelated inverse.',
  'one_entry_source':'u=compact4/current5, original32 target15, residual(0,9). N[u,remaining3]=0, so only actual Pair2 cofactors enter. Source coefficients and all5 Bhat values independently checked.',
  'one_entry_formula':'hybrid(0,9)=rho_u/(delta2*rhoLCM*D4)*Saux; hybrid rescales ONLY selected5 H source rows, so unselected u row carries rho_u.',
  'independent_arithmetic':'768 source operations in a separate P/a^k DEGREVLEX engine with no min-a normalization; independent fixed bordered5 determinant by120 permutations equals all7022 QQ coefficients.',
  'eight_Q9_other_pure_tails_retained':True,'not_claimed':['Any other entry/minor/factor/point/prime','Residual rank6 or global fullcubic57/quartic239','Removal of R434=0 or Delta5=0','Free basis of an unknown projective kernel']}
 save('source-and-block-proof.json',{'proof':proof,'source_atoms':atommanifest,'P4_actual_source_node_ids':A,'P4_target15_source_nodes':target,'actual_pair_only_source_nodes':Ps,'all151_node_shapes_and_source_bindings':True})
 for p in [str(Path(__file__).relative_to(ROOT)),str(cpp.relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'typed_DIRECT_P4_H5_actual_samecoker_6_by13_accepted':True,
  'all151_types_source_indices_hybrid_units_and_forward_covector_transport_accepted':True,
  'actual_one_u4_target15_entry_all768_operations_and7022_QQcoefficients_verified':True,
  'independent_bordered5_formula_and_three_actual_source_cofactors_verified':True,
  'no_h22_inverse_all8_Q9_extra_tails_retained':True,'chart':'FULLU D(R434*Delta5)',
  'residual_rank6_global57_other_entry_or_complement_removal_claimed':False,
  'input_sha256':hashes,'ignored_reproducible_input':{'path':str(ip.relative_to(ROOT)),'sha256':sha(str(ip.relative_to(ROOT)))},
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']}),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--chart-only',action='store_true');main(p.parse_args().chart_only)
