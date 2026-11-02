#!/usr/bin/env python3
"""Two full cubic minors and exact finite-field codimension-one support exclusion."""
import argparse,json,os,subprocess,tempfile,time
from pathlib import Path
import numpy as np
from s61_binary_hyperplane_cubic_matrix_20261102 import ROOT,build_matrix,compile_native,digest
from s61_binary_hyperplane_cubic_horizontal_20261102 import PRIME,read_polys,mod
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod

def parse_poly(z,offset=0):
 nt=int(z[offset]);assert len(z)==offset+1+4*nt
 return [[list(map(int,z[offset+2+4*i:offset+5+4*i])),int(z[offset+1+4*i])] for i in range(nt)]
def bag(path,pointdet):
 lines=iter(Path(path).read_text().splitlines());h=next(lines).split();assert h[:2]==['MINOR','0'];factors=[];value=1;lead=[0,0,0];lc=1;degrees=[0,0,0];total=0
 for _ in range(int(h[3])):
  z=next(lines).split();assert z[0]=='EXP';exponent=int(z[1]);assert exponent>0;p=parse_poly(next(lines).split());leading=max(p,key=lambda v:(v[0][2],v[0][1],v[0][0]));val=sum(c*pow(2,m[1],PRIME)*pow(3,m[2],PRIME) for m,c in p)%PRIME;assert val
  value=value*pow(val,exponent,PRIME)%PRIME;lc=lc*pow(leading[1],exponent,PRIME)%PRIME;lead=[a+exponent*b for a,b in zip(lead,leading[0])];degrees=[a+exponent*max(m[j] for m,c in p) for j,a in enumerate(degrees)];total+=exponent*max(sum(m) for m,c in p)
  factors.append(dict(exponent=exponent,terms=len(p),leading_a_e_d=leading[0],leading_coefficient=leading[1]))
 scalar=pointdet*pow(value,-1,PRIME)%PRIME
 return dict(factor_count=len(factors),factors=factors,constant_scalar_modp=scalar,lex_order='d,e,a',leading_exponents_a_e_d=lead,leading_coefficient_modp=lc*scalar%PRIME,individual_degrees_a_e_d=degrees,total_degree=total)

def run_modular(wd,matrix=None):
 start=time.time();wd=Path(wd);wd.mkdir(exist_ok=True,parents=True);source=build_matrix(wd/'source') if matrix is None else None;matrix=wd/'source/schur' if matrix is None else Path(matrix)
 certpath=ROOT/'code/certificates/s61_binary_hyperplane_cubic_first_minor_20261102.json';cert=json.loads(certpath.read_text());assert digest(matrix)==cert['source_schur_sha256'];(nr,nc,k),q=read_polys(matrix,3);assert (nr,nc)==(42,70)
 q=[[(tuple(m),mod(c)) for m,c in poly if mod(c)] for poly in q]
 def matrix_value(vals):return np.array([sum(c*pow(vals[0],m[0],PRIME)*pow(vals[1],m[1],PRIME)*pow(vals[2],m[2],PRIME) for m,c in poly)%PRIME for poly in q],dtype=np.int64).reshape(nr,nc)
 M=matrix_value([1,2,3]);first=dict(rows=cert['rows'],columns=cert['columns']);adjacent=dict(rows=first['rows'].copy(),columns=first['columns'].copy());adjacent['rows'][1]=19;adjacent['columns'][22]=35
 choices=[first,adjacent]
 for i,choice in enumerate(choices):
  value=int(determinant_mod(M[np.ix_(choice['rows'],choice['columns'])],PRIME));assert value==[19715,5757][i];choice['determinant_at_a1_e2_d3']=value
  with (wd/f'input{i}').open('w') as f:
   f.write(f'{nr} {nc} 36 {PRIME}\n')
   for poly in q:f.write(str(len(poly))+''.join(' '+str(c)+' '+' '.join(map(str,m)) for m,c in poly)+'\n')
   f.write('1\n'+' '.join(map(str,choice['rows']+choice['columns']))+'\n')
 native={name:ROOT/f'code/s61_binary_hyperplane_modular_{name}_20261102.cpp' for name in ['determinants','bag_gcd','small_factors','known_units']}
 for name,p in native.items():compile_native(p,wd/name)
 env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'};logs={};timings={}
 def run(name,args,label,timeout=150):
  t=time.time();cp=subprocess.run([str(wd/name),*map(str,args)],env=env,capture_output=True,text=True,check=True,timeout=timeout);timings[label]=time.time()-t;logs[label]=cp.stderr;print(label,round(timings[label],3),cp.stderr[-120:],flush=True)
 for i in range(2):run('determinants',[wd/f'input{i}',wd/f'bag{i}'],f'determinant{i}')
 bags=[bag(wd/f'bag{i}',choice['determinant_at_a1_e2_d3']) for i,choice in enumerate(choices)]
 assert bags[0]['leading_exponents_a_e_d']==[101,129,79] and bags[0]['leading_coefficient_modp']==8492
 run('known_units',[wd/'bag0',wd/'known_orders'],'known_orders');totals=list(map(int,(wd/'known_orders').read_text().splitlines()[-1].split()[1:]));assert totals==[37,10,12,26]
 run('bag_gcd',[wd/'bag0',wd/'bag1',wd/'gcd'],'bag_gcd');lines=(wd/'gcd').read_text().splitlines();assert lines[-1]=='HITS 3';gcd=[]
 for line in lines[:-1]:
  z=line.split();assert z[0]=='HIT';gcd.append(dict(factor_pair=list(map(int,z[1:3])),polynomial_sparse=parse_poly(z,3)))
 run('small_factors',[wd/'gcd',wd/'irreducible'],'factorization');factorlines=iter((wd/'irreducible').read_text().splitlines())
 for i in range(3):
  z=next(factorlines).split();assert z[0]=='HIT' and z[-2:]==['FACTORS','1'];z=next(factorlines).split();assert z[:2]==['EXP','1'];assert parse_poly(z,2)==gcd[i]['polynomial_sparse']
 witnesses=[]
 for i,vals in enumerate([[15781,20631,1240],[26448,2,3],[5990,2,3]]):
  assert sum(c*pow(vals[0],m[0],PRIME)*pow(vals[1],m[1],PRIME)*pow(vals[2],m[2],PRIME) for m,c in gcd[i]['polynomial_sparse'])%PRIME==0
  a,e,d=vals;units=[a,(500*a+90*e*e+3*e-24)%PRIME,(50*a*d+18*d*e*e-9*d*e-2*d+60*e-40)%PRIME,(100*a+18*e*e-9*e-2)%PRIME];assert all(units)
  A=matrix_value(vals);rank,cols=elimination_rank(A,PRIME);assert rank==36;_,rows=elimination_rank(A[:,cols].T,PRIME);det=int(determinant_mod(A[np.ix_(rows,cols)],PRIME));assert det
  witnesses.append(dict(factor_index=i,point=vals,rows=rows,columns=cols,determinant_modp=det,unit_values=units))
 paths=[Path(__file__),certpath,*native.values(),ROOT/'code/s61_binary_hyperplane_cubic_horizontal_20261102.py',ROOT/'code/s61_binary_hyperplane_cubic_matrix_20261102.py',ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py'];hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
 if source:hashes.update(source['input_sha256'])
 return dict(all_passed=True,scope='Finite field F32003 only: maximal-minor gcd is a unit after inverting aFBC. Codimension-two rank drops remain possible; no QQ or S61 exhaustion conclusion.',prime=PRIME,source_matrix=source,source_matrix_sha256=digest(matrix),matrix_shape=[42,70],choices=choices,determinant_bags=bags,first_known_orders_a_F_B_C=totals,unknown_common_factors=gcd,each_unknown_factor_irreducible_modp=True,full_rank_witnesses=witnesses,maximal_minor_gcd_unit_on_known_open=True,generated_sha256={p.name:digest(p) for p in wd.iterdir() if p.is_file() and p.name not in native},native_logs=logs,stage_seconds=timings,input_sha256=hashes,seconds=time.time()-start)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists()
 with tempfile.TemporaryDirectory(prefix='poisson-hyperplane-modular-') as temp:result=run_modular(temp)
 with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:result[k] for k in ['all_passed','maximal_minor_gcd_unit_on_known_open','seconds']}))
