#!/usr/bin/env python3
"""SOURCE ONLY finite simple-path supports and original Newton intersection."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_BORDER5_SUPPORT_INTERSECTION_20261105_v63';W=(1,2244,68)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 old=src('reports/S61_NONSQUARE_PREPOWER_K92_POLYNOMIAL_BORDER5_SOURCE_20261105_v63/polynomial-border5-source.json');z=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_FIRST_PATH_SOURCE_20261105_v63/report.json');assert z['all_passed']and not z['source_support_possible']
 closed=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_SOURCE_20261105_v63/closed-inverse-source-proof.json');typ=src('reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63/typed-next-residual-Schur.json');alpha=typ['alpha'];beta=typ['beta'];Ast=closed['Ainverse_minplus'];R=old['R2_normalized'];QQ=old['Q2'];DA={0,1,69};cycles=min(e['min_cost']for e in closed['all_nonloop1_cycles_cost_gt69']);assert cycles==2243;horizon=cycles-1
 source={};inputterms=0
 for rel in['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for line in f:
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column'];source[b,i,j]={int(g)for g,v in e['sparse_integer_polynomial']};inputterms+=len(source[b,i,j])
 assert inputterms==23932
 assert all(min((g for g in source['A',i,i]if g>0),default=10**12)>=cycles for i in[0,2,3])
 assert {g for g in source['A',1,1]if g<cycles}==DA
 paths={};pathproof=[];peak=0
 def guard(extra=0):
  nonlocal peak
  owned=inputterms+sum(map(len,paths.values()))+extra+15000;peak=max(peak,owned)
  if owned>150000:raise RuntimeError('SOURCE_owned_scalar_guard')
  if time.monotonic()-start>45:raise RuntimeError('SOURCE_45s_guard')
 for k in range(4):
  for l in range(4):
   S=set();found=[]
   vslist=[[k]]if k==l else([k]+list(mid)+[l]for n in range(3)for mid in itertools.permutations([x for x in range(4)if x not in[k,l]],n))
   for vs in vslist:
    edges=list(zip(vs,vs[1:]));lower=sum(min(source['A',a,b],default=10**12)for a,b in edges)
    if lower>Ast[k][l]+horizon:continue
    pp={0}
    for a,b in edges:pp={g+h for g in pp for h in source['A',a,b]if g+h<=Ast[k][l]+horizon}
    if 1 not in vs:pp={g+h for g in pp for h in DA if g+h<=Ast[k][l]+horizon}
    S|=pp;found.append(dict(vertices=vs,minimum=lower,sign=(-1)**len(edges),through1=1 in vs,support=sorted(pp)))
   paths[k,l]=S;pathproof.append(dict(start=k,end=l,paths=found));guard()
 HB={};maxscratch=0
 for i in range(7):
  for j in range(7):
   sh=alpha[i]+beta[j];S={d+z-sh for d in source['D',i,j]for z in DA if 0<=d+z-sh<=horizon}
   for k in range(4):
    for l in range(4):
     A=paths[k,l];C=source['C',i,k];B=source['B',l,j]
     if not A or not C or not B or min(A)+min(C)+min(B)-sh>horizon:continue
     for c in C:
      if c+min(A)+min(B)-sh>horizon:continue
      for b in B:
       if c+b+min(A)-sh>horizon:continue
       S.update(c+a+b-sh for a in A if 0<=c+a+b-sh<=horizon)
   HB[i,j]=S;guard(sum(map(len,HB.values()))+len(S))
 TB={}
 for u in range(7):
  for v in range(7):
   S=set().union(*(HB[i,j]for i in range(7)for j in range(7)if R[u][i]!='0'and QQ[j][v]!='0'));S.discard(0)
   if u==v==0:S.add(0)
   TB[u,v]=S;guard(sum(map(len,HB.values()))+sum(map(len,TB.values()))+len(S))
 del HB
 maxshift=max(e['removed_t_power']for e in old['entries']);uniform=horizon-maxshift;detcap=uniform+4*2;assert uniform==2172 and detcap==2180
 def conv(p,q):
  S={a+b for a in p for b in q if a+b<=horizon};return S
 cof0=conv(TB[0,0],TB[6,6])|conv(TB[0,6],TB[6,0]);cofA={};cofB={}
 for q in old['freeQ']:
  j=1+q;cofA[q]=conv(TB[0,6],TB[6,j])|conv(TB[6,6],TB[0,j]);cofB[q]=conv(TB[6,0],TB[0,j])|conv(TB[0,0],TB[6,j])
 NB={};records=[]
 for u,p in enumerate(old['freeP']):
  for v,q in enumerate(old['freeQ']):
   i=1+p;j=1+q;sh=next(e['removed_t_power']for e in old['entries']if(e['row'],e['column'])==(u,v));full=conv(cof0,TB[i,j])|conv(cofA[q],TB[i,0])|conv(cofB[q],TB[i,6]);S={g-sh for g in full if 2<=g-sh<=uniform};NB[u,v]=S
   records.append(dict(row=u,column=v,support=sorted(S),first_above6=min((g for g in S if g>6),default=None),shift=sh));guard(sum(map(len,TB.values()))+sum(map(len,NB.values()))+len(full)*2+len(cof0)+sum(map(len,cofA.values()))+sum(map(len,cofB.values())))
 dp={0:{0}};stages=[]
 for i in range(5):
  nxt={}
  for mask,p in dp.items():
   for j in range(5):
    if mask>>j&1:continue
    k=mask|1<<j;qq={a+b for a in p for b in NB[i,j]if a+b<=detcap-2*(4-i)};nxt.setdefault(k,set()).update(qq);guard(sum(map(len,TB.values()))+sum(map(len,NB.values()))+sum(map(len,dp.values()))+sum(map(len,nxt.values()))+len(qq)*2)
  dp=nxt;stages.append(sum(map(len,dp.values())))
 allowed=dp.get(31,set());wb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');box=wb['output_metadata']['box'];weights=wb['weights'];linear=wb['output_metadata']['linear'];matches=[];outer=26653
 for deg in sorted(allowed):
  if deg<73:continue
  order=outer+deg;v=232227-order;n=v-5;a=5+n%68;n//=68;d=n%33;e=n//33;m=(a,e,d)
  if all(lo<=x<=hi for x,(lo,hi)in zip(m,box))and all(lo<=sum(a*b for a,b in zip(m,w))<=hi for w,(lo,hi)in zip(weights,linear)):matches.append(dict(border_det_degree=deg,whole_order=order,original_monomial=m))
 first=matches[0]if matches else None
 dump(out/'all-simple-path-support-proof.json',dict(horizon=horizon,first_omitted_cycle=cycles,path_records=pathproof,
  proof='allrepeatwalks have anaddednonloopcycle>=2243 or omitteddiagonal term>=2243; allsimplepaths retainvertex1 scalarloop exactly viaDA',no_signed_new_values=True))
 dump(out/'border5-supports-and-Newton-intersection.json',dict(Tbar=[dict(row=i,column=j,support=sorted(p))for(i,j),p in TB.items()],border5=records,uniform_known_entry_horizon=uniform,determinant_support_horizon=detcap,
  determinant_possible_support=sorted(allowed),DP_stages=stages,all_allowed_original_Newton_leading_matches=matches,first_possible_match=first,
  leading_unit_argument='at firstnonzero coefficient all formalunit factors contributeonly constants; any lowerdet support notmatching originalNewton mustvanish inductively',known_low_det_identicallyZERO=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins)
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_signed_target_values=True,first_omitted_cycle=cycles,uniform_border_entry_horizon=uniform,determinant_support_horizon=detcap,
  first_possible_leading_Newton_match=first,all_match_count=len(matches),whole_owned_scalar_upper_measured=peak,no_det11_zero_nonzero_or_global57_claim=True,
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 assert new+external<33554432;dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
