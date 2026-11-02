#!/usr/bin/env python3
"""Direct point witness and geometric generic weighted rank bounds on K=0."""
import sys,json,time,subprocess,os,argparse,hashlib,tempfile
from pathlib import Path
from functools import lru_cache
import sympy as s,numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'code'))
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod

def verify():
 start=time.time();data=json.loads((root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json').read_text());old=json.loads((root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
 a,e,d=s.symbols('a e d');point={a:1,e:2,d:3};value=lambda x:s.Rational(s.sympify(x).subs(point))
 T=s.Matrix([[value(row[j]) for j in data['normal34']] for row in data['T34']]);U=s.Matrix([[value(v) for v in row] for row in data['U34']])
 m=S61();beta={int(k[1:])-1:value(v) for k,v in data['beta'].items()};two,F=specialize(m,beta,True)
 fixed=old['fixed'];cols=old['cols'];J=s.zeros(135,45);H=[]
 for i,f in enumerate(F):
  row=[]
  for mon,c in f.items():
   cc=s.Rational(c)
   if len(mon)==2:
    x,y=mon
    if x!=fixed and y==fixed:J[i,cols.index(x)]+=cc
    if y!=fixed and x==fixed:J[i,cols.index(y)]+=cc
   elif len(mon)==1 and mon[0]==fixed:J[i,44]+=cc
   if fixed not in mon:
    ij=[cols.index(j) for j in mon]+[44]*(2-len(mon))
    if cc:row.append((ij,cc))
  H.append(row)
 assert J*T==s.zeros(135,8);assert J.extract(data['selected34_rows'],range(45))*U==s.eye(34)
 temporary=tempfile.TemporaryDirectory(prefix='poisson-hyperplane-ranks-');wd=Path(temporary.name)
 def write(f,c):
  if c:f.write('1\n'+str(c)+' '+' '.join(['0']*8)+'\n')
  else:f.write('0\n')
 with (wd/'input').open('w') as f:
  f.write('135 0 8 4 34\n'+' '.join(map(str,data['selected34_rows']))+'\n');write(f,1)
  for mat in [J,T,U]:
   for row in mat.tolist():
    for c in row:write(f,c)
  for row in H:
   f.write(str(len(row))+'\n')
   for (i,j),c in row:f.write(f'{i} {j}\n');write(f,c)
 env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(root/'work/m2-runtime/lib')+':/opt/homebrew/lib'
 native=root/'code/parametric_implicit_jets_20261101.cpp'
 command=['clang++','-O3','-std=c++17','-I'+str(root/'work/flint-headers/flint/3.5.0/include'),'-I'+str(root/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(native),str(root/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(wd/'native')]
 subprocess.run(command,check=True,capture_output=True,text=True)
 r=subprocess.run([str(wd/'native'),str(wd/'input'),str(wd/'raw')],env=env,text=True,capture_output=True,check=True)
 print(r.stderr,flush=True)
 p=32003;nw=data['normal_weights34'];rows=[{} for _ in range(135)];tpoly={}
 @lru_cache(maxsize=10000)
 def mod(c):
  q=s.Rational(c);assert int(q.q)%p;return int(q.p)%p*pow(int(q.q),-1,p)%p
 with (wd/'raw').open() as f:
  assert list(map(int,f.readline().split()))==[135,0,8,4]
  for degree in range(2,5):
   for i in range(45+135):
    ct=int(f.readline())
    for _ in range(ct):
     z=f.readline().split();mon=tuple(map(int,z[1:]));c=mod(z[0]);assert sum(mon)==degree
     weight=sum(x*w for x,w in zip(mon,nw))
     expected=data['local_coordinate_weights'][i] if i<45 else data['equation_weights'][i-45]
     assert weight==expected
     if i>=45:rows[i-45][mon]=(rows[i-45].get(mon,0)+c)%p
     elif i==44:tpoly[mon]=(tpoly.get(mon,0)+c)%p

 def monomials(weight):
  def rec(i,left,mon):
   if i==8:
    if left==0:yield tuple(mon)
    return
   for n in range(left//nw[i]+1):yield from rec(i+1,left-n*nw[i],mon+[n])
  return list(rec(0,weight,[]))
 results={}
 for weight in [2,3,4]:
  ms=monomials(weight);lookup={v:i for i,v in enumerate(ms)};mat=[];src=[]
  for i,row in enumerate(rows):
   if not row or data['equation_weights'][i]>weight:continue
   for shift in monomials(weight-data['equation_weights'][i]):
    v=np.zeros(len(ms),dtype=np.int64)
    for mon,c in row.items():v[lookup[tuple(x+y for x,y in zip(mon,shift))]]=c
    mat.append(v);src.append((i,shift))
  mat=np.array(mat,dtype=np.int64);rank,mr=elimination_rank(mat.T,p);_,mc=elimination_rank(mat[mr],p);det=determinant_mod(mat[mr][:,mc],p)
  results[weight]=dict(shape=list(mat.shape),rank=rank,minor_rows=mr,minor_columns=mc,determinant=det)
  print('RANK',weight,rank,'DET',det,flush=True)
 assert [results[i]['rank'] for i in [2,3,4]]==[6,57,239]

 contraction_path=root/'reports/S61_BINARY_NUISANCE_CONTRACTION_20261102.json'
 seven_path=root/'reports/S61_BINARY_DEEPER_SEVEN_SLICE_20261101.json'
 lower_path=root/'reports/S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101.json'
 contraction=json.loads(contraction_path.read_text());seven=json.loads(seven_path.read_text());lower=json.loads(lower_path.read_text())
 assert lower['all_passed'] and lower['source_beta_rank']==16
 assert set(contraction['source_beta'])==set(seven['normalized_beta'])
 for key,expr in contraction['source_beta'].items():
  assert s.expand(s.sympify(expr).subs(a,(-90*e**2-3*e+24)/500)-s.sympify(seven['normalized_beta'][key]))==0
 note=root/'research/notes/S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.md'
 paths=[Path(__file__),note,native,root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json',root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',root/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json',contraction_path,seven_path,lower_path,root/'reports/S61_BINARY_DEEPER_FULL_CORE_REMOVAL_20261102.json',root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'code/s61_binary_deeper_uniform_quadratic_20261101.py']
 for path in paths:
  if path.suffix=='.json' and path.parent.name=='reports':
   dep=json.loads(path.read_text());assert dep['all_passed']
   for name,expected in dep.get('input_sha256',{}).items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected
 out=dict(id='S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102',all_passed=True,point=[1,2,3],prime=p,weighted_rank_witnesses=results,generic_weighted_ranks=[6,57,239],generic_weighted_Hilbert_function=[1,7,23,34],generic_transverse_length=65,all_core_cubic_rank_upper=57,deep_source_inclusion_checked=True,source_dominance_rank=16,native_log=r.stderr,proof_note=str(note.relative_to(root)),input_sha256={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},limitations=['Nuisance zero on a*F*B*C nonzero.','The universal cubic upper bound uses the geometric same-component lower65 proof in the note.','No all-point cubic lower bound, quartic surjectivity, or ambient removal.'],seconds=time.time()-start)
 temporary.cleanup()
 return out

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();assert not args.output.exists()
 result=verify()
 with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
 print(json.dumps({k:result[k] for k in ['all_passed','generic_weighted_ranks','all_core_cubic_rank_upper','seconds']}))
