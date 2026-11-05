#!/usr/bin/env python3
"""SOURCE ONLY: P02 completed group0 versus unfinished-group Newton oversets."""
import argparse,gc,hashlib,importlib.util,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_P02_NEWTON_SCREEN_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def mod(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
 prof=mod('code/s61_nonsquare_prePOWER_K92_70_laplace_source_preflight_20261104_v63.py','profile');parser=mod('code/s61_nonsquare_prePOWER_K92_restore_p02_sparse_state_20261105_v63.py','p02')
 reuse=src('reports/S61_NONSQUARE_PREPOWER_K92_70_LAPLACE_SOURCE_PREFLIGHT_20261104_v63/state-reuse.json');meta=reuse['retained_metadata']
 for p in meta.values():
  if p:p['c']=Q(p['c']);p['box']=tuple(map(tuple,p['box']));p['degree']=tuple(p['degree']);p['linear']=tuple(map(tuple,p['linear']))
 M=[[meta[str(i)]for i in r]for r in reuse['matrix_state_ids']];weights=[meta[str(i)]for i in reuse['weight_state_ids']];schedule=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json')['normative_sign_schedule']
 # Single unfinished group1 suffices to show interval isolation fails if it
 # already covers every nonzero coefficient of the exact group0 sum.
 cc=(0,1,2,4);assert list(itertools.combinations(range(8),4))[1]==cc;other=[i for i in range(8)if i not in cc];top,_=prof.det_profile([[M[i][j]for j in cc]for i in range(4)]);tail=None
 for ti,t in enumerate(schedule):
  rr=[4]+[5+i for i in t['pure_remaining_lower6']];bottom,_=prof.det_profile([[M[i][j]for j in other]for i in rr]);tail=prof.add(tail,prof.mul(weights[ti],bottom))
 bound=prof.mul(top,tail);dump(out/'unfinished-group1-bound.json',dict(group=1,columns=cc,upper=bound,weights=prof.W,no_polynomial_coefficients_evaluated=True))
 del reuse,meta,M,weights,schedule,top,tail,bottom;gc.collect()
 st='reports/S61_NONSQUARE_PREPOWER_K92_70_TAIL_GUARDED_CONTINUATION_20261105_v63';r=src(st+'/report.json');assert r['safety_and_roundtrip_checks_passed']and r['groups_completed']==1 and not r['full_scalar_complete'];p=ROOT/st/'state.bin';pins[str(p.relative_to(ROOT))]=sha(p);assert pins[str(p.relative_to(ROOT))]=='c2775cdfb3225397e517b1af19e7feadb3bb04e18aa68e0665aa9d3acb67e6ce'
 s=parser.read_state(p);h=s['header'];assert h['stage']=='bottom4'and h['group']==1 and h['groups_completed']==1 and h['terms_completed']==10;assert h['tail_id']>=0 and s['pending']['dst']!=h['acc_id'];acc=s['polynomials'][h['acc_id']];outside=[];checked=0
 for mon,c in acc['coefficients']:
  assert c;raw=(mon[0]+h['known_output_a_shift'],mon[1],mon[2]);inside=all(lo<=v<=hi for v,(lo,hi)in zip(raw,bound['box']))and all(lo<=sum(v*w for v,w in zip(raw,weight))<=hi for weight,(lo,hi)in zip(prof.W,bound['linear']));checked+=1
  if not inside:outside.append(mon)
 assert checked==48946
 dump(out/'screen.json',dict(exact_global_accumulator_is_completed_group0_ONLY=True,partial_group1_tail_and_DP_retained_unchanged=True,coefficients_checked=checked,all_nonzero_coefficients_inside_single_unfinished_group1_upper=not outside,outside_count=len(outside),first_outside=outside[:10],
  conclusion='NO coefficient can be isolated by these fixed group-level Newton oversets'if not outside else'INCONCLUSIVE: candidates require checking every other future group',nonzero_det11_proved=False,NO_new_polynomial_products=True))
 note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));total=sum(x.stat().st_size for x in out.iterdir())
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,NO_new_polynomial_coefficients=True,completed_group0_coefficients=checked,unfinished_group1_upper_covers_all=len(outside)==0,Newton_isolation_witness=False,no_det_or_rank_claim=True,
  source_state_unchanged=sha(p)==pins[str(p.relative_to(ROOT))],total_RAW_before_report=total,input_sha256=pins,evidence_sha256={x.name:sha(x)for x in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
