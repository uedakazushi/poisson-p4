#!/usr/bin/env python3
"""Actual honest ell² comparison: source gate, then cached X1/minimal origin.

No new physical implicit coefficients, normalizer jets, or N1 are requested.
The native ell column is used only after separate intrinsic S0 annihilation.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
from ambient_B18_ell2_native_gauge_first_ell_20261105_v63 import selected_top_field
ROOT=Path(__file__).resolve().parents[1]
GATE='AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63'
STEM='AMBIENT_B18_ELL2_T2_MU_20261105_v63'
R='reports/'
PATHS={
 'chain':R+'AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'C0':R+'AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'firstjets':R+'AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'quartic':R+'AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'B56':R+'AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json',
 'old640':R+'AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'old128':R+'AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'origin':R+'AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'transport':R+'AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
 'active':R+'AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
 'resolution':R+'AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
 'G2':R+'AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json',
 'N':R+'AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
 'G0':R+'AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'ell':R+'AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63/actual30-firstell-and-native-S1ell.json',
 'S':R+'AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json'}
PEERS=[
 R+'AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
 R+'AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 R+'AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63_REVIEW/report.json',
 R+'AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63_REVIEW/report.json',
 R+'AMBIENT_B18_NATIVE128_FIRST_20261104_v63_REVIEW/report.json',
 R+'AMBIENT_B18_ELL2_T2_ORIGIN_BRIDGE_20261105_v63_REVIEW/report.json',
 R+'AMBIENT_B18_ELL2_HONEST_T2_DVR_TYPE_20261105_v63_REVIEW/report.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {int(i):F(c) for i,c in v if F(c)}
def guard(v):
 assert all(max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=384 for c in v.values()),'384bit guard'
 return v
def dump(p,d):p.write_text(json.dumps(d,separators=(',',':'))+'\n')
def fullpair(p):assert len(p)==10;return {k:F(p[i]) for k,i in enumerate([0,9]) if F(p[i])}
def pairmul(a,b):
 a0=a.get(0,F(0));b0=b.get(0,F(0));v={}
 if a0*b0:v[0]=a0*b0
 c=a0*b.get(1,F(0))+b0*a.get(1,F(0))
 if c:v[1]=c
 return guard(v)
def selected(path,key):return selected_top_field(ROOT/path,key)
def closure(ops,roots):
 seen=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0] in ['c','actual_C0_scalar_function']
 return seen
def original_nodes(ids):
 # Only retain a bounded arithmetic window, never expand the 567k-node graph.
 s=(ROOT/PATHS['C0']).read_text();key='"DAG_scalar_nodes"';p=s.index('[',s.index(key))+1;decoder=json.JSONDecoder();nodes={}
 for i in range(max(ids)+1):
  while s[p].isspace() or s[p]==',':p+=1
  n,p=decoder.raw_decode(s,p)
  if 330000<=i<=max(ids):nodes[i]=n
 return nodes
def full_G2(t,pref,res,g):
 active=pref['active_variables'];zero=pref['zero_actions'];assert active==[0,1,4,5,7] and zero==[2,3,6]
 cols=[];vmap=active+[8]
 for v in map(parse,g['actual_linear_G2_active102_by42_columns']):
  cols.append({(i//6)*9+vmap[i%6]:c for i,c in v.items()})
 D1=list(map(parse,pref['active_d1_columns']));a0=list(map(parse,t['active_F1_fixed67_columns'][:11]));eta=res['physical_target3_basis_positions'];ei={i:j for j,i in enumerate(eta)}
 nonphys=[i for i in range(14) if i not in eta];I1=[];linear=[]
 for j in t['active_F1_free17']:
  v={j:F(1)}
  for k,i in enumerate(nonphys):v=add(v,a0[k],-D1[j].get(i,F(0)))
  assert not mv(D1,v);I1.append(v);q={}
  for row,c in v.items():
   block,i=divmod(row,14)
   if i in ei:q=add(q,{ei[i]*9+(8 if block==0 else active[block-1]):c})
  linear.append(q)
 for iz,z in enumerate(zero):
  for j in range(17):
   v={j*9+z:F(-1)}
   for k,c in linear[j].items():
    i,var=divmod(k,9);v=add(v,{(17+3*iz+i)*9+var:c})
   cols.append(v)
 for ia,ib in itertools.combinations(range(3),2):
  for j in range(3):cols.append({(17+3*ib+j)*9+zero[ia]:F(1),(17+3*ia+j)*9+zero[ib]:F(-1)})
 assert len(cols)==102
 # Convert fixed active rank receipts to original42-source representations.
 norms={};reps={}
 for rec in g['fixed_rank_receipts']:
  j=rec['source_column'];p=rec['pivot'];v={j:F(1)}
  for k,c in rec['subtract']:v=add(v,reps[k],-F(c))
  reps[p]=guard(scale(v,1/F(rec['divide_by'])));norms[p]=parse(rec['normalized'])
  assert mv([parse(v) for v in g['actual_linear_G2_active102_by42_columns']],reps[p])==norms[p]
 return cols,linear,norms,reps
def solve_G2(rhs,cols,norms,reps):
 zero=[2,3,6];active=[0,1,4,5,7,8];y={}
 for iz,z in enumerate(zero):
  for j in range(17):
   c=-rhs.get(j*9+z,F(0))
   if c:y[42+17*iz+j]=c
 rest=add(rhs,mv(cols,y),-1);v={}
 for j in range(17):
  for k,var in enumerate(active):
   c=rest.get(j*9+var,F(0))
   if c:v[j*6+k]=c
 a={}
 while v:
  p=min(v);assert p in norms,('outside active linearimage',p);c=v[p];a=add(a,reps[p],c);v=add(v,norms[p],-c)
 y=add(y,a);rest=add(rhs,mv(cols,y),-1)
 for k,(ia,ib) in enumerate(itertools.combinations(range(3),2)):
  for j in range(3):
   c=rest.get((17+3*ib+j)*9+zero[ia],F(0))
   if c:y[93+3*k+j]=c
 assert mv(cols,y)==rhs,('all234G2linearequations',add(mv(cols,y),rhs,-1))
 return guard(y)
def preflight(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out.mkdir(parents=True,exist_ok=False);dump(out/'phase.json',{'phase':'source_gate_started_no741rank'})
 hs={}
 def readkey(name,key=None):
  path=PATHS[name];hs[path]=sha(ROOT/path);return selected(path,key) if key else json.loads((ROOT/path).read_text())
 for path in PEERS:
  peer=json.loads((ROOT/path).read_text());assert peer['all_passed'];hs[path]=sha(ROOT/path)
 b=readkey('chain');ops=b['operator_scalar_nodes'];roots=[j for col in b['B48_by640_sparse_columns'] for i,j in col];need=closure(ops,roots)
 atoms=sorted({ops[i][1] for i in need if ops[i][0]=='actual_C0_scalar_function'});assert (len(roots),len(need),len(atoms))==(2808,1302,79)
 # Literal bindings of existing scalar firsts. Keep only the needed source leaves.
 h={}
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  h.update({int(i):(PATHS['firstjets'],key,i,fullpair(p)) for i,p in readkey('firstjets',key).items()})
 b56=readkey('B56');basis=[tuple(m) for m in b56['basis56']];ell=[0]*9;ell[8]=1;ellidx=basis.index(tuple(ell))
 ordinary={i:(PATHS['B56'],'all164_original_source_roots',i,{k:c for k,c in [(0,parse(p).get(0,F(0))),(1,5*parse(p).get(ellidx,F(0)))] if c}) for i,p in b56['all164_original_source_roots']}
 Cbasis=readkey('C0','actual_C0_basis8_scalar_nodes');divs=readkey('C0','all35_regular_coordinate_divisions')
 vals=dict(h);vals.update(ordinary);prefix=[]
 for j,i in enumerate(Cbasis):vals[i]=(PATHS['C0'],'actual_C0_basis8_scalar_nodes',j,{})
 words=readkey('firstjets','source2_monomials');mi={tuple(m):i for i,m in enumerate(words)};src=readkey('firstjets','both57_source2')['CA']['pure_function2jets']
 for record in divs:
  f=record['original_pure_CA_source_function'];p=src[str(f)]
  for j,q in enumerate(record['all7_regular_coordinate_divisions']):
   v={k:F(c) for k,c in [(0,p[mi[(j,)]]),(1,p[mi[(j,8)]])] if F(c)}
   assert v.get(0,F(0))==F(q['origin']);vals[q['node']]=(PATHS['firstjets'],'ordered_CA_prefix',q['node'],v)
   prefix.append({'node':q['node'],'source_function':f,'prefix_index':j,'origin_and_ds':srec(v),'no_square_factor2':True})
 ids=readkey('quartic','original_regular_function_nodes_unmodified');grads=readkey('quartic','all8_regular_function_first_gradients_in_FULL9_CA_F9');gnodes=ids['old_CB']+ids['new_CA']
 for i,g in zip(gnodes,grads):vals[i]=(PATHS['quartic'],'all8_regular_function_first_gradients_in_FULL9_CA_F9',i,{1:F(g[8])} if F(g[8]) else {})
 gb=next(x['gB'] for x in b['accepted_actual_regular_function_identities'] if x['kind']=='accepted_actual_CB_f63_tail');vals[gb]=(PATHS['chain'],'accepted_actual_CB_f63_tail',gb,{})
 missing=[i for i in atoms if i not in vals];ast=original_nodes(missing);arithmetic={};leaves=set();todo=list(missing)
 while todo:
  i=todo.pop()
  if i in vals:leaves.add(i);continue
  if i in arithmetic:continue
  assert i in ast,('uncached_actual_original_function',i)
  n=ast[i];assert n[0] in ['a','m','n','c'],('unbound_original_function',i,n);arithmetic[i]=n
  if n[0] in ['a','m']:todo.extend(n[1:3])
  elif n[0]=='n':todo.append(n[1])
 assert len(missing)==28 and len(arithmetic)==84
 leaves.update(i for i in atoms if i in vals)
 bindings=[{'original_node':i,'path':vals[i][0],'field':vals[i][1],'literal_key_or_position':vals[i][2],'origin_and_ds':srec(vals[i][3])} for i in sorted(leaves)]
 del h,ordinary,src,ast,b56,vals
 old640=readkey('old640');old128=readkey('old128');krows=readkey('origin','all768_fixed_elimination_receipts_and_constant_combinations');kernels=[r for r in krows if r.get('zero',False)]
 assert len(kernels)==741 and sum(len(r['constant_original_source_combination']) for r in kernels)==759
 X0=readkey('origin','full_X0_56_by768');Xcols=[{i:F(X0[i][j]) for i in range(56) if F(X0[i][j])} for j in range(768)]
 assert sum(map(len,Xcols))==97
 for r in kernels:assert not mv(Xcols,parse(r['constant_original_source_combination']))
 t=readkey('transport');pref=readkey('active');res=readkey('resolution');g=readkey('G2');G2,linear,norms,reps=full_G2(t,pref,res,g)
 for j,v in enumerate(G2):assert solve_G2(v,G2,norms,reps)=={j:F(1)}
 M=list(map(parse,t['actual_Mnew0_columns']));assert all(not mv(M,v) for v in Xcols)
 N=readkey('N','N_C_B44');L=[{} for _ in range(48)]
 for i in range(12):
  for k in range(4):
   for j in range(12):
    for mask,c in enumerate(map(F,N[i][j][:4])):
     if c and not(mask&k):L[4*i+k]=add(L[4*i+k],{4*j+(mask|k):c})
 del N
 G=readkey('G0','G8');S=readkey('S','native_S12_by32');S0=[{} for _ in range(128)]
 for row,col,k,sign in S['koszul28_in_first8_rows']:
  for inp in range(4):
   for mask,c in enumerate(map(F,G[k][:4])):
    if c and not(mask&inp):S0[4*col+inp]=add(S0[4*col+inp],{4*row+(mask|inp):F(sign)*c})
 for row,col,c in S['free4_zero_source_columns']:
  for mask in range(4):S0[4*col+mask][4*row+mask]=F(c)
 assert all(mv(L,v)==Xcols[640+j] for j,v in enumerate(S0))
 for r in kernels:
  c={j-640:a for j,a in parse(r['constant_original_source_combination']).items() if j>=640};assert not mv(S0,c)
 ell=readkey('ell');S1=list(map(parse,ell['actual_native_S1_ell_48_by128_sparse_columns']));assert sum(map(len,S1))==126
 oldnz=sum(len(p) for col in old640['source_columns640_rows56_sparse_dual9']+old128['source_columns128_rows48_sparse_dual9'] for row,p in col)
 # Worst case one live rank basis, its original-source representation, original
 # selected columns and all741 factorizations. No dense Y102x741 is retained.
 slots={'three_rank_bases':3*102*102,'all741_original_basis_factors':102*741,
  'old_source_first_nonzeros':oldnz,'new640_ell_upper':56*640 if False else 2808+8*640,
  'native_effective_ell_upper':48*128,'full_L0_sparse':sum(map(len,L)),
  'all741_constant_kernel_coefficients':759,'X0_nonzeros':97,'Mnew0_sparse':sum(map(len,M)),
  'G2_full_linear_sparse':sum(map(len,G2)),'G2_active_rank_representations':sum(map(len,reps.values())),
  'G2_active_normalized_basis':sum(map(len,norms.values())),
  'old_ell_compiler_nodes_and_leaves':2*(len(need)+84+len(leaves)+len(atoms)),
  'current_rhs_y_and_rank_scratch':234+3*102+384,'metadata_and_exact_cache_selected_reserve':6000}
 assert sum(slots.values())<150000,slots
 ep=out/'actual-T2-mu-input-source-and-cache-gate.json'
 data={'coordinate_order9':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],
  'old640_typed_closure_nodes':[[i,ops[i]] for i in sorted(need)],'old640_source_columns':b['B48_by640_sparse_columns'],
  'old79_atom_ids':atoms,'old_exact84_original_arithmetic_nodes':[[i,n] for i,n in sorted(arithmetic.items())],
  'all_exact_cached_leaf_bindings_origin_and_ds':bindings,'all35_ordered_prefix_bindings':prefix,
  'C0_original_source_hash':hs[PATHS['C0']],'full9_CB_h_is_already_in_YA':True,
  'all741_zero_origin_source_combinations':[{'source_column':r['source_column'],'original768':r['constant_original_source_combination']} for r in kernels],
  'all741_fullX0_and_intrinsic_nativeS0_zero_verified':True,
  'Mnew0_columns':[srec(v) for v in M],'full_G2_linear234_by102_columns':[srec(v) for v in G2],
  'active42_normalized_basis':[[p,srec(v)] for p,v in sorted(norms.items())],
  'active42_original_source_representations':[[p,srec(v)] for p,v in sorted(reps.items())],
  'L0_Res4_Ntranspose_columns':[srec(v) for v in L],'native_S0_columns':[srec(v) for v in S0],
  'native_S1ell_columns':[srec(v) for v in S1],
  'all102_G2_unique_leftsolve_fixtures_exact':True,'all6144_L0S0_actual_origin_entries_match':True,
  'source_proof':{'source_presentation':'complete honestT2 132x1186 via active3/26/102/418 resolution; all768 original fullR0 source columns',
   'unit27':'X0 source27 pivots survive in30 free coordinates since M0X0=0 and splitV; cancel them overlocalB',
   'remaining741':'constant source elimination gives actualXeff0=0; applying M gives(G2Y)1=M0Xeff1, no M1',
   'unique_Y0':'fullG2linear102 is injective andG3origin0; all234 equations checked for each suppliedY0',
   'mu':'132-27-rank(Y0)=105-rank(Y0)',
   'native_ell':'L1*S0cN vanishes by separatelycheckedS0cN=0; useL0*S1ell onlyaftercombination',
   'old_ell':'existingactualFULL9sourcefirsts; ds/5, CBh alreadyinYA; orderedHadamardretainsr/s andz<=prefix; oldf63distinctfromnewplane505938',
   'DVR_check':'peered pureDVR type givesmu_D(T2)=11; thereforemu_B<=11 andrankY>=94; not a length inference'},
  'whole_slot_bound':sum(slots.values()),'slot_breakdown':slots,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'coefficient_bits':384},
  'new_physical_or_implicit_coefficients':False,'actual741_Y_or_rank_started':False}
 dump(ep,data);assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':GATE,'all_passed':True,'actual741_Y_or_rank_started':False,'all741_zero_origin_and_nativeS0_verified':True,
  'fullG2linear_unique_solve_fixtures':102,'original_full9_source_cache_closed':True,'whole_slot_bound':sum(slots.values()),
  'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{GATE}.md'),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
 dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();preflight(a.output.resolve())
