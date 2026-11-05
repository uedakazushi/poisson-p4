#!/usr/bin/env python3
"""SOURCE-only rank-one shift and support intersection from saved actual prefixes."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_BORDER5_ACTUAL_SUPPORT_RANK_20261105_v63';W=(1,2244,68)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 actual='reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_FOUR_MATRIX_ENTRY_20261105_v63';r=src(actual+'/report.json');assert r['all_passed']and not r['target_NONZERO'];a=src(actual+'/literal-source-and-all-finite-path-identities.json');m=src(actual+'/all100-four-matrix-and50-knownzero-identities.json');bound=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_SUPPORT_INTERSECTION_20261105_v63/border5-supports-and-Newton-intersection.json');plan=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_FOUR_MATRIX_SOURCE_20261105_v63/polynomial-border5-source.json')
 prefix={(e['row'],e['column']):{int(g)for g,v in e['integer_polynomial']if int(v)}for e in a['all_Tbar']};cofprefix={k:{int(g)for g,v in p if int(v)}for k,p in m['all11_cofactor_integer_polynomials'].items()};knownup=138;horizon=2242;uniform=2172;detcap=2180
 def conv(p,q):return{x+y for x in p for y in q if x+y<=horizon}
 lowproof=[]
 for e in plan['entries']:
  i=e['Tbar_row'];j=e['Tbar_column'];q=plan['freeQ'][e['column']];sh=e['removed_t_power'];s=conv(cofprefix['pivot2'],prefix[i,j])|conv(cofprefix['A'+str(q)],prefix[i,0])|conv(cofprefix['B'+str(q)],prefix[i,6]);bad={g-sh for g in s if 0<=g-sh<68};assert not bad
  # All requested absolute degrees <138, so unknown prefixes cannot enter.
  assert sh+67<=knownup;lowproof.append(dict(row=e['row'],column=e['column'],shift=sh,all_possible_prefix_border_degrees=sorted(g-sh for g in s if 0<=g-sh<=68),all_coefficients_below68_ZERO=True))
 n68=m['matrices']['68'];nonzero=[(i,j,int(v))for i,row in enumerate(n68)for j,v in enumerate(row)if int(v)];assert len(nonzero)==1;pi,pj,pv=nonzero[0];TD=int(plan['Tbar_denominator']);p=Q(pv,TD**3);assert p
 roworder=[pi]+[i for i in range(5)if i!=pi];colorder=[pj]+[j for j in range(5)if j!=pj];RR=[[Q(0)for j in range(5)]for i in range(5)];QQ=[[Q(0)for j in range(5)]for i in range(5)]
 for i,k in enumerate(roworder):RR[i][k]=1/p if i==0 else Q(1)
 for j,k in enumerate(colorder):QQ[k][j]=Q(1)
 assert all(sum(RR[i][k]*Q(n68[k][l],TD**3)*QQ[l][j]for k in range(5)for l in range(5))==Q(i==j==0)for i in range(5)for j in range(5))
 knownorder=26653+5*68+4;assert knownorder==26997
 TB={(e['row'],e['column']):set(e['support'])for e in bound['Tbar']}
 for k in TB:TB[k]={g for g in TB[k]if g>knownup}|prefix[k]
 cof0=conv(TB[0,0],TB[6,6])|conv(TB[0,6],TB[6,0]);cofA={};cofB={}
 def tighten(s,k):return{g for g in s if g>knownup}|cofprefix[k]
 cof0=tighten(cof0,'pivot2')
 for q in plan['freeQ']:
  j=1+q;cofA[q]=tighten(conv(TB[0,6],TB[6,j])|conv(TB[6,6],TB[0,j]),'A'+str(q));cofB[q]=tighten(conv(TB[6,0],TB[0,j])|conv(TB[0,0],TB[6,j]),'B'+str(q))
 NB={};entries=[]
 for e in plan['entries']:
  i=e['Tbar_row'];j=e['Tbar_column'];u=e['row'];v=e['column'];q=plan['freeQ'][v];sh=e['removed_t_power'];S={g-sh for g in(conv(cof0,TB[i,j])|conv(cofA[q],TB[i,0])|conv(cofB[q],TB[i,6]))if 69<=g-sh<=uniform}
  if n68[u][v]:S.add(68)
  NB[u,v]=S;entries.append(dict(row=u,column=v,support=sorted(S)))
 dp={0:{0}};stages=[];peak=0
 for i in range(5):
  nxt={}
  for mask,pp in dp.items():
   for j in range(5):
    if mask>>j&1:continue
    k=mask|1<<j;S={a+b for a in pp for b in NB[i,j]if a+b<=detcap-68*(4-i)};nxt.setdefault(k,set()).update(S);peak=max(peak,sum(map(len,dp.values()))+sum(map(len,nxt.values()))+2*len(S)+15000)
    assert peak<150000
  dp=nxt;stages.append(sum(map(len,dp.values())))
 possible=dp.get(31,set());wb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');box=wb['output_metadata']['box'];weights=wb['weights'];linear=wb['output_metadata']['linear'];matches=[]
 for deg in sorted(possible):
  if 26653+deg<knownorder:continue
  value=232227-26653-deg;n=value-5;mon=(5+n%68,(n//68)//33,(n//68)%33)
  if all(lo<=x<=hi for x,(lo,hi)in zip(mon,box))and all(lo<=sum(x*y for x,y in zip(mon,w))<=hi for w,(lo,hi)in zip(weights,linear)):matches.append(dict(border_det_degree=deg,whole_order=26653+deg,original_monomial=mon))
 first=matches[0]if matches else None
 dump(out/'actual-prefix-support-and-rankone-proof.json',dict(all25_low_support_checks=lowproof,all_border_coefficients_below68_ZERO=True,known_matrix68_integer=n68,integer_scale=str(TD**3),exact_nonzero_entry=dict(row=pi,column=pj,integer=str(pv),QQ=str(p)),rank68=1,
  whole_order_lower_bound=knownorder,proof='B5=t68*(B68+t...); rankB68=1, dimension5 =>orddetB5>=340+4=344; outer26653'))
 dump(out/'typed-residual4-and-next-source.json',dict(source='B5=t68*L5; R4*L5*Q4=[[U,W],[X,J]],U0=1,W0=X0=J0=0; residual4=J-X*U^-1*W',R4=RR,Q4=QQ,row_order=roworder,column_order=colorder,
  fixed_unit_pivot_QQ=str(p),only_constant_QQ_inverse=True,new_residual_size=4,outer_order_before_residual4=26993,known_whole_order=knownorder,
  first_remaining_source_Newton_match=first,target_residual4_degree=first['whole_order']-26993 if first else None,all_Schur_corrections_retained=True,no_new_coefficient_or_inverse_series=True))
 dump(out/'tightened-source-supports.json',dict(Tbar_prefix_bound=knownup,actual15entry21term_prefix=True,all11_cofactor_prefixes_reused=True,border5=entries,all_determinant_support=sorted(possible),DP_stages=stages,all_remaining_Newton_matches=matches,
  next_coefficient_NOT_evaluated=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_signed_coefficients=True,all_border_coefficients_below68_ZERO=True,actual_matrix68_rank=1,new_residual_size=4,whole_order_lower_bound=knownorder,
  next_source_Newton_match=first,no_full_det11_zero_nonzero_or_global57_claim=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
