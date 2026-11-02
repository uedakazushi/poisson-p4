"""Bounded QQ first-base-order probe using exact scalar normal forms."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement as cwr
from s61_flag_kernel_scalar_export_20261102 import export_scalar
from s61_flag_kernel_first_linear_raw_20261102 import raw_first_linear
ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/S61_FLAG_ZERO12_FIRST_LINEAR_T_ACTION_20261102.md'

def verify_in(directory):
 start=time.monotonic();root=ROOT;outdir=directory
 scalar_export=outdir/'scalar-full.txt';raw_export=outdir/'raw7.txt'
 scalar_report=export_scalar(scalar_export);raw_report=raw_first_linear(raw_export)
 def rowread(line):
  z=line.split();assert len(z)==1+2*int(z[0]);return{int(i):Q(c)for i,c in zip(z[1::2],z[2::2])}
 def line(p):return str(len(p))+' '+' '.join(f'{i} {c}'for i,c in sorted(p.items()))+'\n'
 def add(ans,p,c=Q(1)):
  for i,d in p.items():
   v=ans.get(i,Q(0))+c*d
   if v:ans[i]=v
   elif i in ans:del ans[i]
 scalar_m=[m for d in range(7)for m in cwr(range(2,9),d)];scalar_ix={m:i for i,m in enumerate(scalar_m)}
 z=scalar_export.read_text().splitlines();L,count,n=map(int,z[0].split());assert(L,count,n)==(142,1716,7)
 sd=list(map(int,z[1].split()));nf=[rowread(q)for q in z[2:2+count]];actions=[rowread(q)for q in z[2+count:]];assert len(actions)==7*142
 reps={}
 for m,r in zip(scalar_m,nf):
  if len(r)==1:
   k,c=next(iter(r.items()))
   if c==1 and len(m)==sd[k]:reps.setdefault(k,m)
 assert len(reps)==142 and all(nf[scalar_ix[m]]=={k:Q(1)}for k,m in reps.items())
 raw_m=[m for d in range(8)for m in cwr(range(9),d)if sum(j<2 for j in m)<=1];rix={m:i for i,m in enumerate(raw_m)};assert len(raw_m)==6864
 top_m=[m for m in raw_m if len(m)==7]
 cache=outdir/'top7-projection.json'
 if not cache.exists():
  initial=json.loads((root/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json').read_text());polys=[]
  for f in initial['initial_relations']:
   polys.append('+'.join('('+c+')*'+'*'.join('x_'+str(j)for j in m)for m,c in f))
  program='R=QQ[x_0..x_8];J=ideal('+','.join(polys)+',x_0^2,x_0*x_1,x_1^2);G=gb J;assert(numColumns basis(7,R/J)==35);out=openOut "'+str(outdir/'top7-projection.txt')+'";\n'
  for m in top_m:
   expr='*'.join('x_'+str(j)for j in m)
   program+='a=listForm(('+expr+')%G);out<<toString(#a)<<" ";scan(a,z->{out<<toString(last z)<<" ";scan(first z,v->out<<toString(v)<<" ");});out<<endl;\n'
  program+='close out;';src=outdir/'top7-projection.m2';src.write_text(program)
  env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(root/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
  r=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(src)],env=env,capture_output=True,text=True,check=True,timeout=35)
  data=[]
  for s in (outdir/'top7-projection.txt').read_text().splitlines():
   v=s.split();assert len(v)==1+10*int(v[0]);data.append([[list(map(int,v[2+10*i:11+10*i])),v[1+10*i]]for i in range(int(v[0]))])
  cache.write_text(json.dumps(data))
 data=json.loads(cache.read_text());assert len(data)==len(top_m)
 topkeys=sorted({tuple(e)for row in data for e,c in row});assert len(topkeys)==35
 topreps=[tuple(j for j,c in enumerate(e)for _ in range(c))for e in topkeys]
 lowrecords=[(len(m),-1,m)for m in scalar_m]+[(sd[k]+1,b,reps[k])for b in range(2)for k in range(142)if sd[k]<6]
 lowrecords.sort();lowN=len(lowrecords);assert lowN==1970
 lowrep=[m if b<0 else tuple(sorted((b,)+m))for _,b,m in lowrecords];assert len(set(lowrep))==lowN
 lowix={m:i for i,m in enumerate(lowrep)};N=lowN+35;ki={k:lowN+i for i,k in enumerate(topkeys)}
 topproj={m:{ki[tuple(e)]:Q(c)for e,c in row}for m,row in zip(top_m,data)}
 def scalar_column(b,k):
  m=tuple(sorted((b,)+reps[k]));return{lowix[m]:Q(1)}if sd[k]<6 else topproj[m]
 def projectmon(m):
  if len(m)>7 or sum(j<2 for j in m)>1:return{}
  if len(m)==7:return topproj[m]
  if not m or m[0]>1:return{lowix[m]:Q(1)}
  b=m[0];ans={}
  for k,c in nf[scalar_ix[m[1:]]].items():add(ans,scalar_column(b,k),c)
  return ans
 proj=[projectmon(m)for m in raw_m];rep=lowrep+topreps
 assert all(projectmon(m)=={i:Q(1)}for i,m in enumerate(rep))
 degrees=[len(m)for m in rep];assert degrees==sorted(degrees)
 shift=[[projectmon(tuple(sorted(m+(j,))))for m in rep]for j in range(9)]
 assert all(all(degrees[k]>degrees[i]for k in r)for a in shift for i,r in enumerate(a))
 def project(p):
  ans={}
  for j,c in p.items():add(ans,proj[j],c)
  return ans
 raw=[]
 for z in raw_export.read_text().splitlines():
  a=z.split();assert len(a)==6854;raw.append({j:Q(c)for j,c in enumerate(a,10)if c!='0'})
 assert len(raw)==136
 pr=[project(r)for r in raw]
 print('COMPRESSED',N,'low',lowN,'projected row terms',sum(map(len,pr)),time.monotonic()-start,flush=True)
 maps=outdir/'maps.txt'
 with maps.open('w')as f:
  f.write(f'9 {N} {lowN} 135\n'+' '.join(map(str,degrees))+'\n')
  for a in shift:
   for r in a:f.write(line(r))
 qq=outdir/'qq.txt';qq.write_text('\n'.join(' '.join(str(r.get(j,0))for j in range(10,N))for r in pr[:135])+'\n')
 p=32003
 def modrow(r):return{k:int(c.numerator*pow(c.denominator,-1,p)%p)for k,c in r.items()if k<lowN and c.numerator%p}
 mod=outdir/'mod.txt'
 with mod.open('w')as f:
  f.write(f'{p} 9 {lowN} 135\n'+' '.join(map(str,degrees[:lowN]))+'\n')
  for a in shift:
   for r in a[:lowN]:f.write(line(modrow(r)))
  for r in pr[:135]:f.write(line(modrow(r)))
 ev=outdir/'evaluation.txt'
 with ev.open('w')as f:
  f.write(f'{len(raw_m)}\n')
  for m,r in zip(raw_m,proj):f.write((f'{rix[m[:-1]]} {m[-1]} 'if m else '-1 -1 ')+line(r))
 lifts=json.loads((root/'reports/S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102.json').read_text())['actual_generator_lifts']
 def lift(ff):
  return{rix[m]:Q(c)for e,c in ff if sum(e[:2])<=1 for m in [tuple(j for j,k in enumerate(e)for _ in range(k))]}
 generators=[lift(ff)for ff in lifts['kernel_generators_sixjet']]
 def product(f,g):
  ans={}
  for i,c in f.items():
   for j,d in g.items():
    m=tuple(sorted(raw_m[i]+raw_m[j]))
    if m in rix:ans[rix[m]]=ans.get(rix[m],Q(0))+c*d
  return{j:c for j,c in ans.items()if c}
 queries=[project(product(raw[-1],g))for g in generators]+[pr[-1]]
 query=outdir/'queries.txt';query.write_text(str(len(queries))+'\n'+''.join(line(r)for r in queries))
 (outdir/'metadata.json').write_text(json.dumps(dict(N=N,lowN=lowN,representatives=rep,scalar_basis_representatives=[reps[k]for k in range(142)],query_count=len(queries))))
 with tempfile.TemporaryDirectory(prefix='zero12-compressed-linear-')as td:
  td=Path(td);binaries=[]
  for src in [root/'code/s61_flag_kernel_sparse_trace_20261102.cpp',root/'code/s61_flag_kernel_first_linear_replay_20261102.cpp']:
   binary=td/src.stem;subprocess.run(['c++','-O3','-w','-std=c++17','-I/opt/homebrew/include',str(src),'-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(binary)],check=True);binaries.append(binary)
  trace=td/'trace.txt';order=td/'order.txt'
  r=subprocess.run([str(binaries[0]),str(mod),str(trace),str(order)],capture_output=True,text=True,check=True,timeout=30);print('GF TRACE',r.stdout,time.monotonic()-start,flush=True)
  assert r.stdout.strip()=='1574';assert trace.exists()
  (outdir/'trace.txt').write_bytes(trace.read_bytes());(outdir/'order.txt').write_bytes(order.read_bytes())
  r=subprocess.run([str(binaries[1]),str(maps),str(qq),str(trace),str(ev),str(query),str(td/'answers.txt')],capture_output=True,text=True,timeout=90)
  print('QQ RESULT',r.returncode,r.stdout,r.stderr,time.monotonic()-start,flush=True)
  assert r.returncode==0
  (outdir/'answers.txt').write_bytes((td/'answers.txt').read_bytes())
  answers=(outdir/'answers.txt').read_text().splitlines();assert answers[0]=='431 13';assert all(not rowread(z)for z in answers[1:13])
  report=dict(all_passed=True,field='QQ',length=431,raw_monomials=6864,scalar_algebra_length=142,scalar_monomials=1716,scalar_lower_trace_rank=665,compressed_columns=N,lower_columns=lowN,top_quotient_columns=35,exact_trace_rank=1574,t_times_actual12_generators_zero=True,actual_t_M_contained_in_b_squared_M=True,first_linear_T_matrices_zero=True,kernel_first_base_length=333,full_kernel_free_claimed=False,full_ring_nilpotence_claimed=False,relative_ambient_removal_claimed=False,all_original135_checked=True,all_monomial_evaluations_checked=True,all_coordinate_commutators_checked=True,base_square_operators_zero=True,query_normal_form_lengths=[len(rowread(z))for z in answers[1:]],raw_report=raw_report,scalar_report=scalar_report,replay_log=r.stderr,seconds=time.monotonic()-start)
  report['id']='S61_FLAG_ZERO12_FIRST_LINEAR_T_ACTION_20261102';report['proof_note']=NOTE
  paths=[Path(__file__),ROOT/NOTE]+[ROOT/'code'/f for f in ['s61_flag_kernel_scalar_export_20261102.py','s61_flag_kernel_residue_export_engine_20261102.cpp','s61_flag_kernel_first_linear_raw_20261102.py','s61_flag_kernel_sparse_trace_20261102.cpp','s61_flag_kernel_first_linear_replay_20261102.cpp','s61_flag_kernel_residue_raw_engine_20261102.py','monomial_quotient_implicit_jets_20261102.cpp','truncated_local_ideal_filtered_trace_20261101.cpp','s61_model.py','s61_minimal_infinity_local_data_20261101.py','s61_flag_transverse_obstruction_20261101.py']]
  paths += [ROOT/'reports'/f for f in ['S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json','S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101.json','S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102.json','S61_FLAG_ZERO12_KERNEL_FIRST_ORDER_20261102.json','S61_FLAG_ZERO12_RESIDUE_ACTIONS_20261102.json']]
  for path in paths:
   if path.suffix=='.json':
    primary=json.loads(path.read_text());assert primary['all_passed']
    for f,h in primary.get('input_sha256',{}).items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h,f
  report['input_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths}
  return report

def verify():
 with tempfile.TemporaryDirectory(prefix='zero12-first-linear-T-')as directory:
  return verify_in(Path(directory))

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();result=verify()
 with args.output.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
