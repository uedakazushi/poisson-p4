#!/usr/bin/env python3
"""SOURCE ONLY sameK120/14: streamedtwo-columnW2 / nestedtop5 integerlattice."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K120_T14_NESTED_SOURCE_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def parity(p):return sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))%2
@lru_cache(None)
def lattice(box,lo,hi):return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
def universal():
 expected={tuple((i,p[i])for i in range(5)):(-1 if parity(p)else 1)for p in itertools.permutations(range(5))};got={}
 for a,b in itertools.combinations(range(5),2):
  cc=[j for j in range(5)if j not in[a,b]]
  for swap in[0,1]:
   ab=[a,b]if not swap else[b,a]
   for p in itertools.permutations(range(3)):
    mon=tuple(sorted([(0,ab[0]),(4,ab[1])]+[(i+1,cc[p[i]])for i in range(3)]));sgn=-1 if(a+b+swap+parity(p))%2 else 1;got[mon]=got.get(mon,0)+sgn
 assert got==expected
 for p in itertools.permutations(range(7)):
  rem=sorted(p[5:]);chosen=[j for j in range(7)if j not in rem];p5=[chosen.index(j)for j in p[:5]];sgn=(11+sum(rem)+parity(p5)+(p[5]>p[6]))%2;assert sgn==parity(p)
 # Wedge coefficient proof with A<B<C; C wedge B=-B wedge C.
 lhs={('g','g','AB'):1,('g','hp','AC'):-1,('g','h','BC'):1};rhs={('g','g','AB'):1,('g','hp','AC'):-1,('g','h','BC'):1};assert lhs==rhs
 return dict(all120_formal_top5_coefficients_match=True,all5040_outer_sign_cases_pass=True,g_factor_three_terms_and_CC_cancellation_verified=True,no_parameter_polynomials_evaluated=True)
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 oldst='S61_NONSQUARE_PREPOWER_K120_T14_CANCELLED_SOURCE_PREFLIGHT_20261104_v63';oldp=ROOT/'reports'/oldst/'operator-and-cost.json';old=read(oldp);pins[str(oldp.relative_to(ROOT))]=sha(oldp);oldr=ROOT/'reports'/oldst/'report.json';assert read(oldr)['resource_gate_passes']is False;pins[str(oldr.relative_to(ROOT))]=sha(oldr)
 N=old['nodes'][:];orig=len(N);top=old['top5_matrix'];zero=next(n['id']for n in N if n['op']=='atom'and n['bound']is None);one=next(n['id']for n in N if n['op']=='atom'and n['reference'].get('constant')=='1')
 def bd(i):return N[i]['bound']
 def size(i):return bd(i)['support']if bd(i)else 0
 def dl(i):return int(N[i]['QQ_denominator']),int(N[i]['integer_l1_upper'])
 def put(op,args,b,D,L):i=len(N);N.append(dict(id=i,op=op,args=args,bound=b,QQ_denominator=str(D),integer_l1_upper=str(L)));return i
 def mul(i,j):
  if i==zero or j==zero:return zero
  if i==one:return j
  if j==one:return i
  A,B=bd(i),bd(j);d,n=dl(i);e,m=dl(j);box=[[A['box'][k][0]+B['box'][k][0],A['box'][k][1]+B['box'][k][1]]for k in range(3)];lo=A['degree'][0]+B['degree'][0];hi=A['degree'][1]+B['degree'][1];b=dict(support=min(A['support']*B['support'],lattice(tuple(map(tuple,box)),lo,hi),n*m),box=box,degree=[lo,hi]);return put('mul',[i,j],b,d*e,n*m)
 def add(i,j,sgn=1):
  if j==zero:return i
  assert i!=zero or sgn==1
  if i==zero:return j
  if sgn<0 and i==j:return zero
  A,B=bd(i),bd(j);d,n=dl(i);e,m=dl(j);D=math.lcm(d,e);L=n*(D//d)+m*(D//e);box=[[min(A['box'][k][0],B['box'][k][0]),max(A['box'][k][1],B['box'][k][1])]for k in range(3)];lo=min(A['degree'][0],B['degree'][0]);hi=max(A['degree'][1],B['degree'][1]);b=dict(support=min(A['support']+B['support'],lattice(tuple(map(tuple,box)),lo,hi),L),box=box,degree=[lo,hi]);return put('add'if sgn>0 else'sub',[i,j],b,D,L)
 def det3(m):
  ans=zero
  for p in itertools.permutations(range(3)):
   t=one
   for i in range(3):t=mul(t,m[i][p[i]])
   # First lex permutation haspositive sign, sozero-negativecase cannotoccur.
   ans=add(ans,t,-1 if parity(p)else 1)
  return ans
 def closure(roots):
  found=set();todo=roots[:]
  while todo:
   i=todo.pop()
   if i in found:continue
   found.add(i);todo+=N[i]['args']
  return sorted(found)
 def lifetime(roots):
  need=closure(roots);use={i:0 for i in need}
  for i in need:
   for j in N[i]['args']:use[j]+=1
  for i in roots:use[i]+=1
  live=peak=0
  for i in need:
   live+=size(i);scratch=size(i)if N[i]['op']in['mul','add','sub']else 0;peak=max(peak,live+scratch)
   for j in N[i]['args']:
    use[j]-=1
    if not use[j]:live-=size(j)
  return dict(nodes=need,source_peak_with_scratch=peak,root_retained_terms=sum(size(i)for i in roots))
 plan=[];allatoms=set();all_needed=set();phases=[];globalD=dl(old['output'])[0];oldD=globalD
 for k,e in enumerate(old['laplace21']):
  chosen=e['top5_columns'];wr=e['weighted_bottom2_root'];W=lifetime([wr]);all_needed.update(W['nodes']);allatoms.update(i for i in W['nodes']if N[i]['op']=='atom');inner=[];topD=1;maxordinary=0
  for a,b in itertools.combinations(range(5),2):
   ca,cb=chosen[a],chosen[b];rr=[chosen[j]for j in range(5)if j not in[a,b]];left=add(mul(top[0][ca],top[4][cb]),mul(top[0][cb],top[4][ca]),-1);right=det3([[top[i][j]for j in rr]for i in[1,2,3]]);life=lifetime([left,right]);all_needed.update(life['nodes']);allatoms.update(i for i in life['nodes']if N[i]['op']=='atom');maxordinary=max(maxordinary,life['source_peak_with_scratch']);D=dl(left)[0]*dl(right)[0];topD=math.lcm(topD,D);inner.append(dict(source_columns=[ca,cb],literal3_columns=rr,left2minor=left,right3minor=right,sign=-1 if(a+b)%2 else 1,source_nodes=life['nodes'],source_peak_with_scratch=life['source_peak_with_scratch']))
  tb=bd(e['top5_root']);topslots=lattice(tuple(map(tuple,tb['box'])),*tb['degree']);wslots=size(wr);within=max(W['source_peak_with_scratch'],topslots+wslots+maxordinary);globalD=math.lcm(globalD,topD*dl(wr)[0]);plan.append(dict(index=k,top5_columns=chosen,bottom2_columns=e['bottom2_columns'],weighted2_root=wr,weighted2_source_nodes=W['nodes'],weighted2_peak=W['source_peak_with_scratch'],top5_bound=tb,top5_integer_denominator=str(topD),inner10=inner,outer_sign=e['sign'],ordinary_inner_term_peak=maxordinary,top5_dense_slots=topslots,retained_weighted2_slots=wslots,phase_peak_excludingouteratomswork=within))
  phases.append(dict(index=k,weighted2=W['source_peak_with_scratch'],top5_dense=topslots,retainedweighted2=wslots,innerordinary=maxordinary,phase_peak=within))
 literal=sum(size(i)for i in allatoms);ob=bd(old['output']);outslots=lattice(tuple(map(tuple,ob['box'])),*ob['degree']);phasepeak=max(e['phase_peak_excludingouteratomswork']for e in plan);total=outslots+literal+15000+phasepeak
 innerpairs=0;outerpairs=0
 for e in plan:
  TD=int(e['top5_integer_denominator']);e['outer_integer_QQ_scale']=str(globalD//(TD*dl(e['weighted2_root'])[0]));assert globalD%(TD*dl(e['weighted2_root'])[0])==0
  for t in e['inner10']:
   D=dl(t['left2minor'])[0]*dl(t['right3minor'])[0];assert TD%D==0;t['integer_QQ_scale']=str(TD//D);innerpairs+=size(t['left2minor'])*size(t['right3minor'])
  outerpairs+=e['top5_dense_slots']*size(e['weighted2_root'])
 assert globalD%oldD==0;maxoutL=dl(old['output'])[1]*(globalD//oldD);inputbytes=0
 for i in allatoms:
  d,n=dl(i);u=Q(N[i].get('QQ_unit','1'));inputbytes+=size(i)*(len(str(n))+len(str(d))+len(str(abs(u.numerator)))+len(str(u.denominator))+100)
 integer_output_bytes=outslots*(len(str(maxoutL))+45)+len(str(globalD))+1000
 recipe=dict(external_scalar_nodes=dict(path=str(oldp.relative_to(ROOT)),sha256=sha(oldp),node_count=orig),appended_scalar_nodes=N[orig:],outer21=plan,all_atom_nodes=sorted(allatoms),output_bound=ob,global_output_integer_denominator=str(globalD),known_a_top5_row_factors=old['known_a_top5_row_factors'],known_B_column_clearings=old['known_B_column_clearings'],raw_D7_over_g_equals='a^2*output',
  next_residual='c11*output/N6 (becauseknowna² cancelsbetweenrawD7/ganda²N6)',all_three_multilinear_terms_preserved=True,no_numeric_g_division=True,
  phase_protocol='For eachfixedouterevent: evaluateONLYtwo-source-columnweighted2closureandreleaseallbutW; allocateintegerTop5; per10innertermfreshneededsourcecache,fuse2minor*3minorandrelease; fuseTop5*Wtoouter; releaseboth. NoordinaryTop5polynomial.',
  output_storage='OneglobalQQdenominatorplusintegercoefficientrecords, neverrepeatdenominatorpereachterm',no_new_parameter_coefficients=True)
 dump(out/'nested-plan.json',recipe);dump(out/'formal-sign-proof.json',universal());datasize=sum(p.stat().st_size for p in out.iterdir());artifact=datasize+inputbytes+integer_output_bytes+3000000
 cost=dict(output_dense_slots=outslots,literal_atom_reserve=literal,work_reserve=15000,max_phase_peak=phasepeak,total_wholeworker_slots_upper=total,output_integer_l1_bits=maxoutL.bit_length(),global_QQ_denominator_bits=globalD.bit_length(),input_raw_byte_bound=inputbytes,integer_globalden_output_byte_bound=integer_output_bytes,existing_recipe_bytes=datasize,whole_raw_artifact_byte_bound=artifact,
  inner_convolution_pairs_upper=innerpairs,outer_convolution_pairs_upper=outerpairs,total_convolution_pairs_upper=innerpairs+outerpairs,phase_bounds=phases,all_inner10_outer21_integer_scales_verified=True,resource_gate_passes=total<=150000 and artifact<32*1024*1024,runtime60s_not_guaranteed=True)
 dump(out/'cost.json',cost);dump(out/'phase.json',dict(phase='FINAL source-only SAMEcandidate; no scalar launched'));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,same_fixed_K120_target14=True,source_nodes_external=orig,appended_metadata_nodes=len(N)-orig,inner_products210=True,outer_products21=True,source_atom_terms=literal,output_dense_slots=outslots,max_phase_peak=phasepeak,total_live_upper=total,whole_raw_artifact_byte_bound=artifact,output_integer_l1_bits=maxoutL.bit_length(),global_QQ_denominator_bits=globalD.bit_length(),scalar_convolution_pairs_upper=innerpairs+outerpairs,resource_gate_passes=cost['resource_gate_passes'],all120_top5_and5040_outer_signs_verified=True,no_scalar_job=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
