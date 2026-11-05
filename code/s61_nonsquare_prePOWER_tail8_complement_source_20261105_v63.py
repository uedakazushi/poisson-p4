#!/usr/bin/env python3
"""SOURCE ONLY: eight cached-kappa tail contractions and their global chart family."""
import argparse, collections, hashlib, json, math, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TAIL8_COMPLEMENT_SOURCE_20261105_v63'
TAILS=[1,2,3,4,6,7,8,9]
ROWS=[2,4,5,8]
KROOTS=[179,198,217,237]
B=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]

def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def accepted(st,flag=None):
  r=load(st,'report.json');assert r['all_passed']
  if flag:assert r[flag]
  return r
 def dump(name,v):(out/name).write_text(json.dumps(v,indent=2)+'\n')
 accepted('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW')
 accepted('S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW')
 accepted('S61_NONSQUARE_PREPOWER_M12_AXIS_ONE_COEFFICIENT_20261105_v63_REVIEW')
 accepted('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW')
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63')
 cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json')
 pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63')
 oldsource=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63','source-binding.json')
 assert oldsource['source_current15_rows']==ROWS and oldsource['cofactor_roots']==KROOTS
 kappas=[load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63',f'root{i}.json')['polynomial'] for i in KROOTS]
 T=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63','root239.json')['polynomial']
 assert len(T)==7208 and min(m[0] for m,c in T)==4
 Tmeta=dict(terms=len(T),box=[[min(m[k]for m,c in T),max(m[k]for m,c in T)]for k in range(3)],
   degree=[min(sum(m)for m,c in T),max(sum(m)for m,c in T)],named_unit_a_power4_removable_by_offset_only=True,
   no_C_or_other_polynomial_factorization_performed=True)
 del T
 cpterms=sum(sum(len(p)for p in v)for row in cp['matrix25x32']for v in row)
 max_source_terms=cpterms+sum(len(p)for p in kappas)+7208
 assert max_source_terms<150000
 profiles=[];operators=[]
 def lattice(box,lo,hi):
  return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)
    for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
 for target in TAILS:
  N=[];cache={};atomkeys={};inventory=[]
  def put(op,args,bound,D,L,**kw):
   i=len(N);N.append(dict(id=i,op=op,args=args,bound=bound,QQ_D=str(D),integer_l1_upper=str(L),**kw));return i
  def atom(p,ref,unit=Q(1),shift=0):
   vals=[([m[0]+shift,m[1],m[2]],Q(c)/unit)for m,c in p if Q(c)]
   key=tuple((tuple(m),str(c))for m,c in vals)
   if key in atomkeys:return atomkeys[key]
   D=math.lcm(*(c.denominator for m,c in vals))if vals else 1
   L=sum(abs(c.numerator)*(D//c.denominator)for m,c in vals)
   bd=None if not vals else dict(support=len(vals),box=[[min(m[k]for m,c in vals),max(m[k]for m,c in vals)]for k in range(3)],degree=[min(sum(m)for m,c in vals),max(sum(m)for m,c in vals)])
   i=put('external_atom',[],bd,D,L,reference=ref,QQ_divisor=str(unit),a_shift=shift)
   atomkeys[key]=i;return i
  z=atom([],dict(QQ='0'));one=atom([[[0,0,0],'1']],dict(QQ='1'))
  bn=atom(B,dict(named_FULLU_unit='B'))
  def bd(i):return N[i]['bound']
  def dl(i):return int(N[i]['QQ_D']),int(N[i]['integer_l1_upper'])
  def mul(i,j):
   if i==z or j==z:return z
   if i==one:return j
   if j==one:return i
   A,C=bd(i),bd(j);d,n=dl(i);e,m=dl(j)
   box=[[A['box'][k][0]+C['box'][k][0],A['box'][k][1]+C['box'][k][1]]for k in range(3)]
   lo=A['degree'][0]+C['degree'][0];hi=A['degree'][1]+C['degree'][1]
   return put('mul_bound_only',[i,j],dict(support=min(A['support']*C['support'],lattice(box,lo,hi),n*m),box=box,degree=[lo,hi]),d*e,n*m)
  def add(i,j,sgn=1):
   if j==z:return i
   if sgn<0 and i==j:return z
   if i==z:return j if sgn==1 else mul(atom([[[0,0,0],'-1']],dict(QQ='-1')),j)
   A,C=bd(i),bd(j);d,n=dl(i);e,m=dl(j);D=math.lcm(d,e);L=n*(D//d)+m*(D//e)
   box=[[min(A['box'][k][0],C['box'][k][0]),max(A['box'][k][1],C['box'][k][1])]for k in range(3)]
   lo=min(A['degree'][0],C['degree'][0]);hi=max(A['degree'][1],C['degree'][1])
   return put('add_bound_only'if sgn==1 else'sub_bound_only',[i,j],dict(support=min(A['support']+C['support'],lattice(box,lo,hi),L),box=box,degree=[lo,hi]),D,L)
  def timesB(i,k):
   for _ in range(k):i=mul(i,bn)
   return i
  def imp(i):
   if i in cache:return cache[i]
   n=q['scalar_DAG']['nodes'][i];op=n[0]
   if op=='z':v=(z,0)
   elif op=='one':v=(one,0)
   elif op in ['a','m','s']:
    x,b=imp(n[1]);y,c=imp(n[2])
    if op=='m':v=(mul(x,y),b+c)
    else:
     bb=max(b,c);v=(add(timesB(x,bb-b),timesB(y,bb-c),-1 if op=='s'else 1),bb)
   else:
    num,den=cp['matrix25x32'][n[1]][n[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
    if len(den)==1:
     assert den[0][0][1:]==[0,0] and Q(den[0][1]);unit=Q(den[0][1]);shift=-den[0][0][0];bpower=0
    else:
     dd={tuple(m):Q(c)for m,c in den};bb={tuple(m):Q(c)for m,c in B}
     assert dd.keys()==bb.keys();unit=dd[next(iter(bb))]/bb[next(iter(bb))]
     assert unit and all(dd[m]==unit*bb[m]for m in bb);shift=0;bpower=1
    inventory.append(dict(old_node=i,definition=n,QQ_unit=str(unit),a_shift=shift,B_power=bpower))
    v=(atom(num,dict(kind=op,old_node=i,definition=n),unit,shift),bpower)
   cache[i]=v;return v
  j=q['new22_target_original32_columns'].index(target)
  source_roots=[q['remaining15x22_scalar_roots'][r][j]for r in ROWS]
  vv=[imp(i)for i in source_roots];bpower=max(b for i,b in vv)
  cols=[timesB(i,bpower-b)for i,b in vv]
  col_needed=set();todo=cols[:]
  while todo:
   i=todo.pop()
   if i in col_needed:continue
   col_needed.add(i);todo+=N[i]['args']
  uses={i:0 for i in col_needed}
  for i in col_needed:
   for k in N[i]['args']:uses[k]+=1
  for i in cols:uses[i]+=1
  live=peak=0
  for i in sorted(col_needed):
   ss=bd(i)['support']if bd(i)else 0;live+=ss
   peak=max(peak,live+(ss if N[i]['args']else 0))
   for k in N[i]['args']:
    uses[k]-=1
    if not uses[k]:live-=bd(k)['support']if bd(k)else 0
  ks=[atom(p,dict(cached_kappa_root=r,path=f'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root{r}.json'))for p,r in zip(kappas,KROOTS)]
  termnodes=[mul(k,v)for k,v in zip(ks,cols)];output=z
  for t in termnodes:output=add(output,t)
  D,L=dl(output);ob=bd(output)
  actual_literal=sum(bd(i)['support']for i in col_needed if not N[i]['args'] and bd(i))
  cached_k=sum(len(p)for p in kappas);col_support=sum(bd(i)['support']if bd(i)else 0 for i in cols)
  max_input_atom=max([len(p)for p in kappas]+[bd(i)['support']for i in col_needed if not N[i]['args'] and bd(i)])
  # A future changed direct accumulator has no full product buffer. Count BOTH
  # coefficient+packed-key cells and one transient JSON atom with3 exponents.
  column_phase=2*(peak+actual_literal+cached_k)+5*max_input_atom+15000
  contraction_phase=2*(cached_k+col_support+ob['support'])+15000+100
  readback_phase=2*ob['support']+15000+100
  owned=max(column_phase,contraction_phase,readback_phase)
  atom_bytes=sum((bd(i)['support']if bd(i)else 0)*(len(str(dl(i)[0]))+len(str(dl(i)[1]))+100)for i in col_needed if not N[i]['args'])
  atom_bytes+=sum(len(p)*(160+max(len(c)for m,c in p))for p in kappas)
  raw=2_000_000+atom_bytes+ob['support']*(len(str(D))+len(str(L))+80)
  pairs=sum((bd(k)['support']if bd(k)else 0)*(bd(v)['support']if bd(v)else 0)for k,v in zip(ks,cols))
  profile=dict(original32_target=target,already_evaluated_target2=(target==2),source_q2_roots=source_roots,
   typed_nodes=len(N),column_source_nodes=len(col_needed),literal_column_atom_terms=actual_literal,
   cached4_kappa_terms=cached_k,cleared_column_support_sum_upper=col_support,common_B_power=bpower,
   output_bound=ob,QQ_D_bits=D.bit_length(),output_integer_l1_bits=L.bit_length(),
   direct_accumulator_pair_upper=pairs,prospective_owned_numeric_cells=owned,
   phase_owned=dict(column=column_phase,dot=contraction_phase,readback=readback_phase),
   prospective_RAW_bytes_upper=raw,source_gate_fits=(owned<=150000 and raw<=32*1024*1024 and L.bit_length()<=2048),
   no_actual_contraction_or_column_evaluation=True)
  profiles.append(profile)
  operators.append(dict(target=target,nodes=N,cached_kappa_nodes=ks,column_roots=cols,
    source_current15_rows=ROWS,source_original32_columns=[27,31,22,target],
    common_B_power=bpower,term_roots=termnodes,output_root=output,atom_denominator_inventory=inventory,
    all_bound_nodes_UNEVALUATED=True))
 candidates=[p for p in profiles if p['original32_target']!=2]
 chosen=min(candidates,key=lambda p:(p['prospective_owned_numeric_cells'],p['direct_accumulator_pair_upper'],p['original32_target']))
 cov=load('S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63','coverage.json')
 accepted('S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63_REVIEW')
 cover_reports=[('S61_NONSQUARE_THETA_P80_REPLACEMENT_20261104_v63_REVIEW','ALL_P80_full_cubic57_proved'),
  ('S61_NONSQUARE_THETA_P100_REPLACEMENT_20261104_v63_REVIEW','ALL_P100_full_cubic57_proved'),
  ('S61_NONSQUARE_RAW_H_B_K13_CUBIC_REGULAR_REP_20261104_v63_REVIEW','whole_actual_Ka13_Bbad_full_cubic57')]
 for st,flag in cover_reports:accepted(st,flag)
 accepted('S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63_REVIEW','actual_B_bad_algebra_exact_Ka13_length13_nonempty')
 kop=load('S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63','operator.json')
 accepted('S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63_REVIEW')
 assert 'kappa3*detA8=cprod*Saux*T' in kop['identities']
 dump('all8-profiles.json',profiles)
 dump('typed-family-and-source-plans.json',dict(
   rows=ROWS,anchor=22,pair_H2=[27,31],pair_H3_literal_zero=[28,30,29],
   fixed9_raw_pure_rank_at_most1=[22,11,13,25,15,19,14,18,17],
   other8=TAILS,cached_kappa_roots=KROOTS,source_frame='EXACT same frame as accepted KAPPA_EXTRA2_ENTRY; no extra a on the anchor',
   cofactor_sign='kappa_i=(-1)^i det(A_without_row_i); det([A,v])=-kappa*v',
   Tj_definition='T_j=B^b_j * sum(kappa_i * M15[ROWS[i],original32_column_j])',
   determinant_family='D15_j uses the SAME normalized15 source rows as DIRECT_M12, columns[J11,27,31,j,17], with its j column scaledB^b_j',
   polynomial_identity='A_anchor*D15_j=-T_j*D12; no A/kappa3 inverse chart',
   identity_proof='same three row clearings/Laplace; detHcorr_j=A_anchor^2*detA4_j=-A_anchor^2*T_j; cancel A_anchor^2 in the source domain',
   source_plans=operators,
   sufficient_chart_family='FULLU D(D12) intersect union_j D(T_j)',
   T2_is_exact_existing_T=True,no_further_source_values_evaluated=True))
 dump('coverage-and-frontier.json',dict(existing_T_support=Tmeta,accepted_heavy_coverage=cov,
   covered_full57=['Theta P80 on its permitted source open','Theta P100 on its permitted source open','the actual finite reduced/unreduced Ka13 Bbad scheme (no etale premise needed)'],
   global_Q_L_role='global two Q/L unit completions only; G0=Q*L^4 exclusions are used in actual Bbad13 scheme proof, not an additional full-cubic cover of arbitrary Q/L zero loci',
   known_factor_identities=['kappa3*D8=cprod*Saux*T','A_anchor*D15_2=-T*D12','T=a^4*Tbar by saved support; a known FULLU unit'],
   no_known_factorization_or_gcd_between_T_D12_and_P80_P100_Theta=True,
   unresolved='actual cubic defect is contained in H7 and in V(D12) union (V(T_1,...,T_9) with only indices in the8-tail list); this is only a necessary containment derived from sufficient charts',
   candidate_chosen_from_source_cost_only=chosen['original32_target'],
   no_new_polynomial_gcd_factorization_minor_rank_point_prime_or_source572_replay=True))
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum(p.stat().st_size for p in out.iterdir());assert raw<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,typed_tail_family_global_identity=True,
   all8_profiled_no_new_values=True,selected_source_candidate=chosen['original32_target'],
   chosen_candidate_profile=chosen,largest_loaded_source_terms=max_source_terms,
   no_scalar_job_authorized_or_launched=True,independent_review_pending=True,
   output_RAW_bytes=raw,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump('report.json',report);signal.alarm(0)
 print(json.dumps(dict(id=ID,source_only=True,chosen=chosen,output_RAW_bytes=raw,seconds=report['seconds']),indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
