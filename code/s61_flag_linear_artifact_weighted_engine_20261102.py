"""Regenerate the three rational frame-artifact models from all 135 equations."""
import json,time,subprocess,os,tempfile
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]

def build(wd,core):
 start=time.monotonic()
 m=S61();r=s.Symbol('r');one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1]
 beta=dict.fromkeys(one,s.S.Zero);beta[4]=beta[48]=1;beta[12]=core
 _,eq=specialize(m,beta,True);fixed=25;chart=[i for i in range(45)if i!=fixed]
 J=s.Matrix([[-f.get(tuple(sorted((fixed,j))),0)for j in chart]+[-f.get((fixed,),0)]for f in eq])
 cut=[chart.index(19),chart.index(26)];columns=[i for i in range(45)if i not in cut]
 Jc=J[:,columns];RR,piv=Jc.to_DM().to_field().rref();RR=RR.to_Matrix();assert len(piv)==36
 free=[j for j in range(43)if j not in piv];assert len(free)==7
 T=s.zeros(45,7)
 for k,j in enumerate(free):T[columns[j],k]=1
 for i,j in enumerate(piv):
  for k,z in enumerate(free):T[columns[j],k]=-RR[i,z]
 sel=list(J[:,[columns[j]for j in piv]].T.to_DM().to_field().rref()[1]);assert len(sel)==36
 inv=J[sel,[columns[j]for j in piv]].to_DM().to_field().inv().to_Matrix();U=s.zeros(45,36)
 for i,j in enumerate(piv):U[columns[j],:]=inv[i,:]
 T=T.applyfunc(s.cancel);U=U.applyfunc(s.cancel)
 assert (J*T).applyfunc(s.cancel)==s.zeros(135,7)
 assert (J[sel,:]*U).applyfunc(s.cancel)==s.eye(36)
 r0=json.load(open(ROOT/'reports/S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_20261102.json'));qw=r0['raw_coordinate_weights']
 nw=[]
 for j in range(7):
  w={qw[i]for i in range(45)if T[i,j]};assert len(w)==1;nw.append(next(iter(w)))
 assert nw==[1,1,2,2,1,2,2]
 def denominator(M):
  out=s.S.One
  for v in M:out=s.lcm(out,s.denom(v))
  return s.factor(out)
 td=denominator(T);ud=denominator(U);print('FRAME',time.monotonic()-start,'td',td,'ud',ud,flush=True)
 meta=dict(core=str(core),chart=chart,fixed=fixed,selected=sel,td=str(td),ud=str(ud),normal_weights=nw,T=[[str(z)for z in row]for row in T.tolist()],U=[[str(z)for z in row]for row in U.tolist()])
 wd.mkdir(exist_ok=True)
 (wd/'frame.json').write_text(json.dumps(meta,indent=2)+'\n')
 # Native engine uses Wd_num/(ud^(d-1)*td^d).
 cpp=ROOT/'code/parametric_implicit_jets_20261101.cpp';native_temp=tempfile.TemporaryDirectory(prefix='flag-weighted-core-build-');exe=Path(native_temp.name)/'native'
 subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(cpp),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(exe)],check=True)
 def writepoly(out,v):
  p=s.Poly(s.cancel(v),r,domain=s.QQ);ts=[(a,c)for a,c in p.terms()if c]
  out.write(str(len(ts))+'\n')
  for a,c in ts:out.write(str(c)+' '+str(a[0])+' '+' '.join(['0']*7)+'\n')
 with tempfile.TemporaryDirectory(prefix='flag-weighted-core-native-')as temp:
  inp=Path(temp)/'input.txt';out=Path(temp)/'output.txt'
  with inp.open('w')as f:
   f.write('135 1 7 10 36\n'+' '.join(map(str,sel))+'\n');writepoly(f,ud)
   for M in [J,td*T,ud*U]:
    for row in M.tolist():
     for v in row:writepoly(f,v)
   for pol in eq:
    terms=[(mon,c)for mon,c in pol.items()if fixed not in mon]
    f.write(str(len(terms))+'\n')
    for mon,c in terms:
     a,b=[chart.index(j)for j in mon]+[44]*(2-len(mon));f.write(f'{a} {b}\n');writepoly(f,c)
  env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
  print('NATIVE_INPUT',inp.stat().st_size,time.monotonic()-start,flush=True)
  pr=subprocess.run([str(exe),str(inp),str(out)],env=env,check=True,timeout=100)
  (wd/'raw.txt').write_bytes(out.read_bytes())
  print('RAW_COMPLETE',out.stat().st_size,time.monotonic()-start,flush=True)
 return meta,m,eq,J,T,U
