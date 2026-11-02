"""Exact finite positive-weight models at three ordinary-flag frame artifacts."""
from pathlib import Path
import argparse,json,subprocess,tempfile,os,time,hashlib
import sympy as s,numpy as np
from functools import lru_cache
from s61_flag_linear_artifact_weighted_engine_20261102 import build
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FLAG_LINEAR_ARTIFACT_WEIGHTED_REMOVAL_20261102';NOTE='research/notes/'+STEM+'.md'
W=(1,1,2,2,1,2,2);P=32003
@lru_cache(None)
def mons(d,n=7):
 if n==0:return((),)if d==0 else()
 return tuple(t+(e,)for e in range(d//W[n-1]+1)for t in mons(d-e*W[n-1],n-1))
def mod(q):
 x=s.Rational(q);assert x.q%P;return int(x.p)%P*pow(int(x.q)%P,-1,P)%P

def verify():
 start=time.monotonic();result=[];base_path='reports/S61_FLAG_ORDINARY_WEIGHTED_FINITE_KERNEL_20261102.json';base=json.loads((ROOT/base_path).read_text())
 for p,h in base['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 dpath='reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json';dr=json.loads((ROOT/dpath).read_text());us=s.symbols('u1:62')
 with tempfile.TemporaryDirectory(prefix='flag-artifact-points-')as td:
  td=Path(td);env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'};exe=td/'rows'
  cpp='code/s61_flag_ordinary_weighted_rows_20261102.cpp'
  subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(ROOT/cpp),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(exe)],check=True)
  for num,a in enumerate([-6,-s.Rational(12,5),-s.Rational(119,20)]):
   a=s.Rational(a);meta,m,eq,J,T,U=build(td/f'raw{num}',a);assert a!=2
   inp=td/f'den{num}';out=td/f'out{num}';inp.write_text('\n'.join('1 '+meta[k]for k in['ud','td'])+'\n')
   subprocess.run([str(exe),str(td/f'raw{num}/raw.txt'),str(inp),str(out)],env=env,check=True,capture_output=True,text=True,timeout=60)
   def readp(f):
    ts=f.readline().split();assert len(ts)==int(ts[0])+1 and len(ts)<=2
    return mod(ts[1])if len(ts)>1 else 0
   rows=[]
   with out.open()as f:
    assert int(f.readline())==136
    for i in range(136):
     ii,w,n=map(int,f.readline().split());assert ii==i;dd=readp(f);cc=readp(f);assert dd
     ts=[]
     for j in range(n):
      aa=tuple(map(int,f.readline().split()));c=readp(f);assert sum(x*y for x,y in zip(aa,W))==w;ts.append((aa,c))
     if ts:assert cc
     if i in meta['selected']:assert not ts
     rows.append((w,ts))
   ranks=[]
   for degree in[7,8]:
    mm=mons(degree);idx={v:i for i,v in enumerate(mm)};basis=[];piv=[];det=1;des=[];nr=0
    for i,(w,ts)in enumerate(rows[:135]):
     if not ts or not 0<=w<=degree:continue
     for sh in mons(degree-w):
      nr+=1;v=np.zeros(len(mm),dtype=np.int64)
      for aa,c in ts:v[idx[tuple(x+y for x,y in zip(aa,sh))]]=c
      for p,b in zip(piv,basis):
       if v[p]:v=(v-int(v[p])*b)%P
      nz=np.flatnonzero(v)
      if len(nz):
       p=int(nz[0]);c=int(v[p]);det=det*c%P;v=v*pow(c,-1,P)%P;basis.append(v);piv.append(p);des.append([i,sh])
    assert len(piv)==len(mm)-5 and det
    ranks.append(dict(degree=degree,rows=nr,rank=len(piv),columns=len(mm),minor_mod=det,selected_rows=des,selected_columns=piv,monomials=mm))
    print('CORE',a,'WEIGHT',degree,'RANK',len(piv),'SECONDS',time.monotonic()-start,flush=True)
   qw=base['raw_coordinate_weights'];assert all(not T[i,j]or qw[i]==W[j]for i in range(45)for j in range(7))
   rw=base['row_weights'];assert all(not U[i,j]or qw[i]==rw[meta['selected'][j]]for i in range(45)for j in range(36))
   point=s.zeros(61,1);point[4]=point[48]=1;point[12]=a;point[37]=-1;sub=dict(zip(us,point))
   K=s.Matrix(dr['kernel']).applyfunc(s.sympify).subs(sub);C=s.Matrix(dr['actual_C_tangent']).applyfunc(s.sympify).subs(sub);N=K[:,dr['squarezero_normal_columns']];B=C.row_join(N);cut=[19,26]
   assert C[cut,:]==s.Matrix([[0,0,1],[-2,0,0]])
   H=B*s.Matrix.hstack(*B[cut,:].nullspace());Cc=C*s.Matrix.hstack(*C[cut,:].nullspace());qweights=[0 if i==25 else qw[i-(i>25)]for i in range(45)]
   def chars(M):return {w:M[[i for i,v in enumerate(qweights)if v==w],:].rank()for w in sorted(set(qweights))if M[[i for i,v in enumerate(qweights)if v==w],:].rank()}
   assert chars(H)=={1:3,2:2}and chars(Cc)=={1:1}
   comp=s.Matrix([[-1,2],[-1,a]]).det();assert comp
   result.append(dict(core=str(a),frame=meta,rank_witnesses=ranks,coordinate_weights=qw,normal_weights=W,D_cut_characters={'1':3,'2':2},C_cut_characters={'1':1},relative_compensator_minor=str(comp),ordinary_derivative=str(a-2),actual_kernel_weights=[2,6],relative_t_squared_zero=True,all_nuisance_by_contraction=True))
 paths=[str(Path(__file__).relative_to(ROOT)),NOTE,'code/s61_flag_linear_artifact_weighted_engine_20261102.py',base_path,*base['input_sha256']]
 hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in dict.fromkeys(paths)}
 return dict(id=STEM,all_passed=True,prime=P,models=result,full_Jacobian_rank=36,whole_relative_not_only_transverse=True,all_general_beta_points_on_three_cores=True,other_cores_claimed=False,all_flags_or_affine_exhaustion_claimed=False,proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();result=verify()
 with a.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=result['seconds'],hashes=len(result['input_sha256']))))
