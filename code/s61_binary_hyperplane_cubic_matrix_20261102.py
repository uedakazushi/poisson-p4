#!/usr/bin/env python3
"""Reconstruct the complete weighted cubic matrix directly from the 135 equations."""
import argparse,hashlib,json,os,subprocess,sys,tempfile,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def compile_native(source,target):
 subprocess.run(['clang++','-w','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(source),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(target)],check=True,capture_output=True,text=True)
def build_matrix(wd):
 start=time.time();wd=Path(wd);wd.mkdir(exist_ok=True,parents=True)
 chart_path=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';old_path=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json';quad_path=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
 data=json.loads(chart_path.read_text());old=json.loads(old_path.read_text())['linear'];quad=json.loads(quad_path.read_text());assert data['all_passed'] and quad['all_passed']
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();pol=lambda x:R.from_expr(s.sympify(x));frac=lambda x:K.from_expr(s.sympify(x));known=[a,500*a+90*e**2+3*e-24,50*a*d+18*d*e**2-9*d*e-2*d+60*e-40,100*a+18*e**2-9*e-2]
 def known_unit(p):
  if not p:return False
  for f in known:
   while not p.rem(f):p=p.exquo(f)
  return all(not any(m) for m in p)
 T=[[frac(row[j]) for j in data['normal34']] for row in data['T34']];U=[[frac(v) for v in row] for row in data['U34']];td=[]
 for j in range(8):
  den=R.one
  for row in T:den=den.lcm(row[j].denom)
  assert known_unit(den);td.append(den)
 ud=R.one
 for row in U:
  for v in row:ud=ud.lcm(v.denom)
 assert known_unit(ud)
 Tn=[[v.numer*td[j].exquo(v.denom) for j,v in enumerate(row)] for row in T];Un=[[v.numer*ud.exquo(v.denom) for v in row] for row in U]
 m=S61();beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()};two,F=specialize(m,beta,True);fixed=old['fixed'];cols=old['cols'];J=[[R.zero for _ in range(45)] for _ in range(135)];H=[]
 for i,f in enumerate(F):
  row=[]
  for mon,c in f.items():
   cc=pol(c)
   if len(mon)==2:
    x,y=mon
    if x!=fixed and y==fixed:J[i][cols.index(x)]+=cc
    if y!=fixed and x==fixed:J[i][cols.index(y)]+=cc
   elif len(mon)==1 and mon[0]==fixed:J[i][44]+=cc
   if fixed not in mon:
    ij=[cols.index(j) for j in mon]+[44]*(2-len(mon))
    if cc:row.append((ij,cc))
  H.append(row)
 assert all(not sum((J[i][k]*Tn[k][j] for k in range(45)),R.zero) for i in range(135) for j in range(8))
 assert all(sum((J[i][k]*Un[k][j] for k in range(45)),R.zero)==(ud if r==j else R.zero) for r,i in enumerate(data['selected34_rows']) for j in range(34))
 def write(f,p,normal=0):
  f.write(str(len(p))+'\n')
  for mon,c in p.items():f.write(str(c)+' '+' '.join(map(str,mon+(0,)*normal))+'\n')
 with (wd/'input').open('w') as f:
  f.write('135 3 8 3 34\n'+' '.join(map(str,data['selected34_rows']))+'\n');write(f,ud,8)
  for mat in [J,Tn,Un]:
   for row in mat:
    for p in row:write(f,p,8)
  for row in H:
   f.write(str(len(row))+'\n')
   for (i,j),p in row:f.write(f'{i} {j}\n');write(f,p,8)
 (wd/'weights').write_text(' '.join(map(str,data['local_coordinate_weights']+data['equation_weights']+data['normal_weights34']))+'\n')
 with (wd/'coeff').open('w') as f:
  write(f,ud);f.write('3\n')
  for p in known[:3]:write(f,p)
 env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
 logs={}
 native=[ROOT/'code/parametric_implicit_jets_20261101.cpp',ROOT/'code/s61_binary_hyperplane_cubic_extract_20261102.cpp',ROOT/'code/s61_binary_hyperplane_cubic_unit_schur_20261102.cpp']
 for i,source in enumerate(native):compile_native(source,wd/f'native{i}')
 for name,cmd in [('implicit',[wd/'native0',wd/'input',wd/'raw']),('extract',[wd/'native1',wd/'raw',wd/'coeff',wd/'rows',wd/'row_meta',wd/'weights'])]:
  cp=subprocess.run(list(map(str,cmd)),env=env,text=True,capture_output=True,check=True);logs[name]=cp.stderr;print(cp.stderr,flush=True)
 qs=[{tuple(mon):R.from_dict({tuple(mm):s.QQ(c) for mm,c in pp}) for mon,pp in row} for row in quad['quadratic_rows_sparse']];rows=[];labels=[]
 for i,q in enumerate(qs):
  for j in [0,1,2,3,4,6,7]:
   rows.append({tuple(x+(1 if h==j else 0) for h,x in enumerate(mon)):p for mon,p in q.items()});labels.append(['product',i,j])
 with (wd/'rows').open() as f:
  count=int(f.readline());assert count==27
  for _ in range(count):
   idx,w,ct=map(int,f.readline().split());row={}
   for j in range(ct):
    mon=tuple(map(int,f.readline().split()));z=f.readline().split();p=R.from_dict({tuple(map(int,z[2+4*k:5+4*k])):s.QQ(z[1+4*k]) for k in range(int(z[0]))});row[mon]=p
   if w==2:
    ref=qs[quad['quadratic_rows'].index(idx)];assert set(row)==set(ref)
    mon=next(iter(row));num=row[mon];den=ref[mon];g=num.gcd(den);assert known_unit(num.exquo(g)) and known_unit(den.exquo(g))
    assert all(row[m]*den==ref[m]*num for m in row)
    continue
   assert w==3
   g=R.zero
   for p in row.values():g=g.gcd(p)
   div=R.one
   for p in known:
    while not g.rem(p):g=g.exquo(p);div*=p
   rows.append({mon:p.exquo(div) for mon,p in row.items()});labels.append(['cubic',idx])
 mons=sorted(set().union(*(row for row in rows)));assert len(rows)==63 and len(mons)==91
 with (wd/'full').open('w') as f:
  f.write('63 91 57\n')
  for row in rows:
   for mon in mons:
    p=row.get(mon,R.zero);f.write(str(len(p)))
    for mm,c in p.items():f.write(' '+str(c)+' '+' '.join(map(str,mm)))
    f.write('\n')
 cp=subprocess.run([str(wd/'native2'),str(wd/'full'),str(wd/'schur')],env=env,capture_output=True,text=True,check=True);logs['unit_schur']=cp.stderr;print(cp.stderr,flush=True)
 head=(wd/'schur').open().readline().split();assert head==['42','70','21']
 paths=[Path(__file__),chart_path,old_path,quad_path,*native,ROOT/'code/s61_binary_hyperplane_cubic_qq3_20261102.hpp',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py']
 return dict(all_passed=True,ordinary_jet_order=3,full_matrix_shape=[63,91],known_unit_pivots=21,residual_shape=[42,70],full_matrix_sha256=digest(wd/'full'),schur_matrix_sha256=digest(wd/'schur'),labels=labels,normal_monomials=mons,Tden=list(map(str,td)),Uden=str(ud),native_logs=logs,input_sha256={str(p.relative_to(ROOT)):digest(p) for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists()
 with tempfile.TemporaryDirectory(prefix='poisson-hyperplane-cubic-matrix-') as temp:result=build_matrix(temp)
 with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:result[k] for k in ['all_passed','full_matrix_sha256','schur_matrix_sha256','seconds']}))
