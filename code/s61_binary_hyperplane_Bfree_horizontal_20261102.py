#!/usr/bin/env python3
"""Reconstruct the exact three Schur kernels of the B-free quartic module.

All intermediate files are regenerated in a private directory. The small step
sources are immutable inputs; none of the earlier work/ caches is read.
"""
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_hyperplane_cubic_matrix_20261102 import digest
STEPS=ROOT/'code/s61_binary_Bfree_horizontal_steps_20261102'
DATE='_20261102'

def compile_native(src,dst):
 subprocess.run(['clang++','-w','-O3','-std=c++17','-I'+str(ROOT/'code'),'-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(src),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(dst)],check=True,capture_output=True,text=True)

def build_source(wd,env,run):
 """Original 134 equations plus G6, with the certified B-free frame."""
 wd.mkdir(parents=True,exist_ok=True)
 files=[ROOT/'reports/S61_BINARY_HYPERPLANE_BFREE_FRAME_20261102.json',ROOT/'reports/S61_BINARY_HYPERPLANE_LINEAR_LIGHT_FRAME_20261102.json',ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json',ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json']
 bf,lf,ch,old=[json.loads(f.read_text())for f in files];old=old['linear']
 R,a,e,d=ring('a,e,d',sp.QQ);K=R.to_field();pol=lambda x:R.from_expr(sp.sympify(x));frac=lambda x:K.from_expr(sp.sympify(x))
 U=[[frac(v)for v in row]for row in bf['U_new']];ud=pol(bf['U_common_denominator']);T=[[pol(v)for v in row]for row in lf['T_new']];Un=[[v.numer*ud.exquo(v.denom)for v in row]for row in U]
 beta={int(k[1:])-1:sp.sympify(v)for k,v in ch['beta'].items()};m=S61();two,F=specialize(m,beta,True);fixed=old['fixed'];cols=old['cols'];J=[[R.zero for _ in range(45)]for _ in range(135)];H=[]
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
    ij=[cols.index(j)for j in mon]+[44]*(2-len(mon))
    if cc:row.append((ij,cc))
  H.append(row)
 assert all(not sum((J[i][k]*T[k][j]for k in range(45)),R.zero)for i in range(135)for j in range(8))
 assert all(sum((J[i][k]*Un[k][j]for k in range(45)),R.zero)==(ud if r==j else R.zero)for r,i in enumerate(ch['selected34_rows'])for j in range(34))
 def write(f,p,extra=0):
  f.write(str(len(p))+'\n')
  for mon,c in p.items():f.write(str(c)+' '+' '.join(map(str,mon+(0,)*extra))+'\n')
 with (wd/'input').open('w')as f:
  f.write('135 3 8 3 34\n'+' '.join(map(str,ch['selected34_rows']))+'\n');write(f,ud,8)
  for mat in [J,T,Un]:
   for row in mat:
    for p in row:write(f,p,8)
  for row in H:
   f.write(str(len(row))+'\n')
   for (i,j),p in row:f.write(f'{i} {j}\n');write(f,p,8)
 (wd/'weights').write_text(' '.join(map(str,ch['local_coordinate_weights']+ch['equation_weights']+ch['normal_weights34']))+'\n')
 with (wd/'coeff').open('w')as f:
  write(f,ud);f.write('4\n')
  for p in [a,500*a+90*e**2+3*e-24,50*a+9*e**2+3*e-6,100*a+18*e**2-9*e-2]:write(f,p)
 for label,source in [('raw_native',ROOT/'code/parametric_implicit_jets_20261101.cpp'),('extract_native',ROOT/'code/s61_binary_hyperplane_cubic_extract_20261102.cpp')]:compile_native(source,wd/label);files.append(source)
 run('source_implicit',[wd/'raw_native',wd/'input',wd/'raw'])
 run('source_extract',[wd/'extract_native',wd/'raw',wd/'coeff',wd/'rows',wd/'row_meta',wd/'weights'])
 files += [ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py']
 return dict(all_passed=True,original_equations=135,all_JT_and_selected_JU_checked=True,denominator=str(ud),source_rows_sha256=digest(wd/'rows'),input_sha256={str(f.relative_to(ROOT)):digest(f)for f in files})

def check_Bzero_low_ranks(rows_path):
 """A good-reduction minor on the remaining prime B in QQ(a,e)[d]."""
 import functools
 from fractions import Fraction
 import numpy as np
 from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
 p=32003;point=[1,2,-40*pow(51,-1,p)%p];weights=[1,1,1,1,1,2,1,1];rows=[]
 @functools.lru_cache(None)
 def co(c):
  c=Fraction(c);assert c.denominator%p;return c.numerator*pow(c.denominator,-1,p)%p
 with Path(rows_path).open()as stream:
  for _ in range(int(stream.readline())):
   source,w,n=map(int,stream.readline().split());row={}
   for _ in range(n):
    m=tuple(map(int,stream.readline().split()));z=stream.readline().split();row[m]=sum(co(z[1+4*k])*pow(point[0],int(z[2+4*k]),p)*pow(point[1],int(z[3+4*k]),p)*pow(point[2],int(z[4+4*k]),p)for k in range(int(z[0])))%p
   rows.append((source,w,row))
 def monomials(w,i=0,t=()):
  if i==8:
   if w==0:yield t
   return
  for k in range(w//weights[i]+1):yield from monomials(w-k*weights[i],i+1,t+(k,))
 result=[]
 for target in [2,3]:
  mons=list(monomials(target));mat=[];labels=[]
  for src,w,row in rows:
   if w>target:continue
   for shift in monomials(target-w):
    shifted={tuple(x+y for x,y in zip(m,shift)):c for m,c in row.items()};mat.append([shifted.get(m,0)for m in mons]);labels.append([src,shift])
  A=np.array(mat,dtype=np.int64);rank,cc=elimination_rank(A,p);_,rr=elimination_rank(A[:,cc].T,p);det=int(determinant_mod(A[np.ix_(rr,cc)],p));result.append(dict(weight=target,shape=list(A.shape),rank=rank,pivot_rows=rr,pivot_columns=cc,minor=det,labels=labels))
 assert [(x['rank'],x['minor'])for x in result]==[(6,9925),(57,24692)]
 aa,ee,dd=point;F=500*aa+90*ee**2+3*ee-24;D0=50*aa+9*ee**2+3*ee-6;C=100*aa+18*ee**2-9*ee-2;B=50*aa*dd+18*dd*ee**2-9*dd*ee-2*dd+60*ee-40
 assert B%p==0 and aa*F*D0*C%p
 return dict(all_passed=True,prime=p,point_rational=['1','2','-40/51'],point_residues=point,known_open_residue=aa*F*D0*C%p,matrices=result)

def run_certificate(wd, resume=False):
 st=time.time();wd=Path(wd).resolve();wd.mkdir(parents=True,exist_ok=True);work=wd/'work';tmp=wd/'tmp';work.mkdir(exist_ok=True);tmp.mkdir(exist_ok=True)
 for name,target in [('code',ROOT/'code'),('reports',ROOT/'reports'),('research',ROOT/'research')]:
  if not (wd/name).exists():(wd/name).symlink_to(target,target_is_directory=True)
 if not (work/'m2-runtime').exists():(work/'m2-runtime').symlink_to(ROOT/'work/m2-runtime',target_is_directory=True)
 env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'};logs={};timings={}
 def run(label,args,timeout=240):
  if resume and (work/(label+'.log')).exists():
   logs[label]=(work/(label+'.log')).read_text();timings[label]=0;print('RESUME',label,flush=True);return
  begin=time.time();cp=subprocess.run(list(map(str,args)),cwd=wd,env=env,capture_output=True,text=True,check=False,timeout=timeout);logs[label]=cp.stdout+cp.stderr;timings[label]=time.time()-begin;(work/(label+'.log')).write_text(logs[label]);cp.check_returncode();print(label,round(timings[label],3),logs[label][-160:].strip(),flush=True)
 scripts={}
 for src in sorted(STEPS.glob('*.py')):
  dst=work/src.name;dst.write_text(src.read_text().replace('/tmp/',str(tmp)+'/'));scripts[src.stem]=dst
 def step(name,timeout=240):run(name,[sys.executable,scripts[name+DATE]],timeout)
 def data(s,ext='out'):return tmp/f'hyperplane-Bfree-infinity-{s}-20261102.{ext}'
 def native(name,args,timeout=240):run(name,args,timeout)
 low_dependencies={}
 for stem in ['QUADRATIC_CODIMENSION','CUBIC_CODIMENSION','GENERIC_RANKS']:
  q=json.loads((ROOT/f'reports/S61_BINARY_HYPERPLANE_{stem}_20261102.json').read_text());assert q['all_passed']
  for key in ['input_sha256','primary_input_sha256']:
   for path,h in q.get(key,{}).items():assert digest(ROOT/path)==h
  low_dependencies[stem]=q
 assert low_dependencies['QUADRATIC_CODIMENSION']['bad_locus_dimension_upper_on_three_core_open']==1
 assert low_dependencies['CUBIC_CODIMENSION']['cubic_rank_drop_dimension_upper']==1
 assert low_dependencies['CUBIC_CODIMENSION']['all_core_cubic_rank_upper']==57
 source=build_source(tmp/'poisson-hyperplane-Bfree-cubic-20261102',env,run)
 step('binary_hyperplane_Bfree_t_targets_univariate')
 step('binary_hyperplane_Bfree_targets_exact_core_input')
 from s61_binary_hyperplane_Bfree_modular_certificate_20261102 import run_minor_module
 modular=run_minor_module(tmp/'modular',tmp/'hyperplane-Bfree-t-targets-exact-core-20261102.in',tmp/'hyperplane-Bfree-t-targets-univariate-20261102.in')
 print('modular minor/module',modular.get('seconds'),flush=True)
 from s61_binary_hyperplane_Bfree_heavy_square_20261102 import run_heavy_square
 heavy=run_heavy_square(tmp/'heavy_square',tmp/'poisson-hyperplane-Bfree-cubic-20261102/raw',tmp/'poisson-hyperplane-Bfree-cubic-20261102/coeff')
 low_ranks=check_Bzero_low_ranks(tmp/'poisson-hyperplane-Bfree-cubic-20261102/rows')
 print('heavy-square and B=0 low ranks PASS',flush=True)
 # Frozen coordinate report is read as data; all six input quadrics are matched to fresh source below.
 shutil.copyfile(ROOT/'reports/S61_BINARY_HYPERPLANE_INFINITY_QUADRATIC_NORMAL_FORM_20261102.json',work/'generic_Bfree_infinity_constant_form_20261102.json')
 natives={}
 for src in sorted(STEPS.glob('*.cpp')):
  dst=tmp/src.stem;compile_native(src,dst);natives[src.stem]=dst
 kernel=tmp/'kernel';compile_native(ROOT/'code/s61_binary_hyperplane_leading_kernel_20261102.cpp',kernel)
 def kr(inp,out,label):run(label,[kernel,inp,out])
 def transpose(inp,out):
  z=inp.read_text().splitlines();nr,nc,r,k=map(int,z[0].split());rr=list(map(int,z[1].split()));cc=list(map(int,z[2].split()));free=[i for i in range(nr)if i not in rr]
  with out.open('w')as f:
   f.write(f'{nc} {nr} {r} {len(free)}\n'+' '.join(map(str,cc))+'\n'+' '.join(map(str,rr))+'\n'+' '.join(map(str,free))+'\n')
   for j in range(nc):
    for i in range(nr):f.write(z[4+i*nc+j]+'\n')
 step('binary_hyperplane_Bfree_infinity_initial_forms')
 for name in ['constant_quadratic_projection','constant_cubic_projection','traced_projection']:step('binary_hyperplane_infinity_'+name)
 step('binary_hyperplane_Bfree_infinity_cubic_reduction');kr(data('reduced-kernel','in'),data('reduced-kernels'),'right75')
 step('binary_hyperplane_Bfree_infinity_kernel_audit');transpose(data('reduced-kernel','in'),data('reduced-left-kernel','in'));kr(data('reduced-left-kernel','in'),data('reduced-left-kernels'),'left116')
 for name in ['quadratic_delta1','quadratic_syzygy_delta1','quadratic_delta2_blocks']:step('generic_Bfree_infinity_'+name)
 step('binary_hyperplane_Bfree_infinity_cubic_second_jet')
 # The second-jet file also contains the identical first-jet fields.
 shutil.copyfile(work/'binary_hyperplane_Bfree_infinity_cubic_second_jet_20261102.json',work/'binary_hyperplane_Bfree_infinity_cubic_first_jet_20261102.json')
 for name in ['corrected_cubic_delta1','corrected_cubic_delta2','cubic_product_delta1','cubic_product_delta2','cubic_product_full_delta1','first_obstruction']:step('binary_hyperplane_Bfree_infinity_'+name)
 obs=json.loads((work/'binary_hyperplane_Bfree_infinity_first_obstruction_20261102.json').read_text())
 def enc(z):return str(len(z))+''.join(' '+c+' '+' '.join(map(str,m))for m,c in z)+'\n'
 with data('first-obstruction','in').open('w')as f:
  f.write('158 106\n')
  for row in obs['rows']:
   row=dict(row)
   for j in range(106):f.write(enc(row.get(j,[])))
 with data('first-obstruction-denominators','in').open('w')as f:
  for q in obs['denominators']:f.write(enc(q))
 native('first_product',[natives['binary_hyperplane_obstruction_product'+DATE],data('first-obstruction','in'),data('reduced-kernels'),data('first-obstruction-product')])
 step('binary_hyperplane_first_obstruction_kernel_input');kr(data('second-kernel','in'),data('second-kernels'),'right34');step('binary_hyperplane_Bfree_infinity_second_kernel_audit')
 transpose(data('second-kernel','in'),data('second-left-kernel','in'));kr(data('second-left-kernel','in'),data('second-left-kernels'),'left117')
 native('second_left_operator',[natives['binary_hyperplane_second_left_operator'+DATE],data('first-obstruction-denominators','in'),data('first-obstruction-product'),data('reduced-left-kernels'),data('second-left-kernels'),data('second-left-operator')])
 step('global_blocks_input')
 native('global_blocks',[natives['global_blocks'+DATE],data('second-obstruction-uncleared','in'),data('second-obstruction-blocks','in')])
 native('second_dual_blocks',[natives['binary_hyperplane_second_dual_blocks'+DATE],data('second-obstruction-blocks','in'),data('second-left-operator'),data('second-dual-blocks'),data('second-dual-cfree')])
 z=data('reduced-kernel','in').read_text().splitlines();rr=list(map(int,z[1].split()));cc=list(map(int,z[2].split()));A=[z[4+i*106:4+(i+1)*106]for i in range(147)];c=data('second-dual-cfree').read_text().splitlines();assert c[0]=='117 106'
 with data('second-dual-lift-solve','in').open('w')as f:
  f.write('106 148 31 117\n'+' '.join(map(str,cc))+'\n'+' '.join(map(str,range(31)))+'\n'+' '.join(map(str,range(31,148)))+'\n')
  for j in range(106):
   for i in rr:f.write(A[i][j]+'\n')
   for i in range(117):f.write(c[1+i*106+j]+'\n')
 kr(data('second-dual-lift-solve','in'),data('second-dual-lift-solutions'),'dual_lift117')
 native('second_dual_effective',[natives['binary_hyperplane_second_dual_effective'+DATE],data('second-obstruction-blocks','in'),data('second-dual-blocks'),data('second-dual-lift-solutions'),data('reduced-kernel','in'),data('second-dual-effective')])
 def strip(src,dst):
  with src.open()as f,dst.open('w')as out:
   n,m=map(int,f.readline().split());out.write(f'{n} {m}\n')
   for i in range(n):
    f.readline()
    for j in range(m):out.write(f.readline())
   assert not f.read().strip()
 current=data('second-dual-effective')
 for name,ker in [('first','reduced-kernels'),('second','second-kernels')]:
  inp=data('second-dual-'+name+'-product','in');out=data('second-dual-'+name+'-product');strip(current,inp)
  native('dual_product_'+name,[natives['binary_hyperplane_second_dual_'+name+'_kernel_product'+DATE],inp,data(ker),out]);current=out
 step('binary_hyperplane_second_obstruction_kernel_input');kr(data('third-kernel','in'),data('third-kernels'),'right16');step('binary_hyperplane_Bfree_infinity_third_kernel_audit')
 audits=[json.loads((work/('binary_hyperplane_Bfree_infinity_'+n+'_20261102.json')).read_text())for n in ['kernel_audit','second_kernel_audit','third_kernel_audit']]
 assert [x['kernel_columns']for x in audits]==[75,34,16] and all(x['all_passed']for x in audits)
 paths=[Path(__file__),*STEPS.glob('*'),ROOT/'code/s61_binary_hyperplane_leading_kernel_20261102.cpp',ROOT/'code/s61_binary_deeper_cubic_field_20261101.hpp',ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py',ROOT/'reports/S61_BINARY_HYPERPLANE_INFINITY_QUADRATIC_NORMAL_FORM_20261102.json',ROOT/'code/certificates/s61_binary_hyperplane_infinity_quadrics_20261102.json',ROOT/'research/notes/S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102.md',ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json',ROOT/'reports/S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102.json',ROOT/'reports/S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.json']
 hashes={str(q.relative_to(ROOT)):digest(q)for q in paths if q.is_file()};hashes.update(source['input_sha256']);hashes.update(modular['input_sha256']);hashes.update(heavy['input_sha256'])
 return dict(all_passed=True,resumed_run=resume,scope='Over K[d], K=QQ(a,e), the B-free weighted transverse algebra is free65 with weighted Hilbert function(1,7,23,34), allweight>=4 zero and t^2=0. The321-row prolongation spans238 columns; raw H^2 Bezout supplies the239th. Vertical(a,e) exceptions, nuisance and ambient saturation are separate.',horizontal_row_module_surjective=True,modular_certificate=modular,heavy_square_certificate=heavy,Bzero_low_ranks=low_ranks,weighted_Hilbert_function=[1,7,23,34],weighted_algebra_free_rank_over_Kd=65,t_square_zero_over_Kd=True,finite_local_base='Z_(32003)[a,e]_(32003,a-1,e-2)',minor_degree_upper=259,source=source,leading_constant_quadratic_rank=132,leading_cubic_rank=31,first_obstruction_rank=41,second_obstruction_rank=18,kernel_dimensions=[75,34,16],universal_maximal_minor_valuation_lower=125,all_full_polynomial_kernel_identities_checked=True,kernel_audits=audits,coordinate_convention='old_light = D*Ti(new_light), old_heavy = H; equivalently D^w*f(Ti(new_light),H/D^2)',stage_seconds=timings,native_logs=logs,input_sha256=hashes,generated_sha256={q.name:digest(q)for q in tmp.iterdir()if q.is_file()},seconds=time.time()-st)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workdir',type=Path);ap.add_argument('--resume',action='store_true',help='Development only: reuse already logged stages, marked in report');args=ap.parse_args();assert not args.output.exists()
 if args.workdir:result=run_certificate(args.workdir,args.resume)
 else:
  with tempfile.TemporaryDirectory(prefix='poisson-Bfree-horizontal-')as wd:result=run_certificate(wd)
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:result[k]for k in ['all_passed','kernel_dimensions','universal_maximal_minor_valuation_lower','seconds']}))
