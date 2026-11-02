"""Finite two-weight certificate for the actual ordinary flag kernel over QQ(r)."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
from collections import Counter
from functools import lru_cache
import sympy as s
import numpy as np
from s61_flag_ordinary_weighted_core_engine_20261102 import build
from s61_model import BASIS,X,J as J5
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FLAG_ORDINARY_WEIGHTED_FINITE_KERNEL_20261102'
NOTE='research/notes/'+STEM+'.md'
P=32003;weights=(1,1,2,2,1,2,2)

def verify():
 start=time.monotonic();r=s.Symbol('r')
 with tempfile.TemporaryDirectory(prefix='flag-ordinary-weighted-')as tmp:
  temp=Path(tmp);meta,m,eq,J,T,U=build(temp/'raw');sel=meta['selected']
  env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'}
  exe=temp/'rows';inp=temp/'den';out=temp/'rows.txt'
  inp.write_text('\n'.join(str(s.Poly(s.sympify(meta[k]),r).degree()+1)+' '+' '.join(str(s.Poly(s.sympify(meta[k]),r).nth(i))for i in range(s.Poly(s.sympify(meta[k]),r).degree()+1))for k in ['ud','td'])+'\n')
  cpp=ROOT/'code/s61_flag_ordinary_weighted_rows_20261102.cpp'
  subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(cpp),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(exe)],check=True)
  result=subprocess.run([str(exe),str(temp/'raw/raw.txt'),str(inp),str(out)],env=env,capture_output=True,text=True,check=True,timeout=25)
  def readp(f):
   ts=f.readline().split();assert len(ts)-1==int(ts[0]);return ts[1:]
  def rat(t):
   a,*b=t.split('/');assert not b or int(b[0])%P
   return int(a)%P*(pow(int(b[0])%P,-1,P)if b else 1)%P
  def atone(ts):return sum(map(rat,ts))%P
  D=s.Poly((r+6)*(5*r+12)*(20*r+119)*(25*r+4),r,domain=s.QQ);rows=[];rowinfo=[]
  with out.open()as f:
   assert int(f.readline())==136
   for i in range(136):
    ii,w,n=map(int,f.readline().split());assert ii==i
    dd=readp(f);cc=readp(f);assert len(cc)<=1 and atone(dd)
    den=s.Poly.from_list([s.Rational(z)for z in dd[::-1]],r,domain=s.QQ)
    rest=den
    while rest.degree()>0:
     g=s.gcd(rest,D);assert g.degree()>0;rest=rest.exquo(g)
    ts=[]
    for j in range(n):
     a=tuple(map(int,f.readline().split()));cs=readp(f);assert sum(x*y for x,y in zip(a,weights))==w
     value=atone(cs)*atone(cc)%P
     if value:ts.append((a,value))
    if i in sel:assert not ts and w==-1
    rows.append((w,ts));rowinfo.append(dict(weight=w,terms=n,common_denominator=str(den.as_expr()),rational_content=cc))
  @lru_cache(None)
  def mons(w,n=7):
   if n==0:return ((),)if w==0 else()
   return tuple(t+(e,)for e in range(w//weights[n-1]+1)for t in mons(w-e*weights[n-1],n-1))
  certificates=[]
  for degree in [7,8]:
   mon=mons(degree);index={a:i for i,a in enumerate(mon)};basis=[];piv=[];desc=[];pv=[];count=0
   for i,(w,ts) in enumerate(rows[:135]):
    if not ts or not 0<=w<=degree:continue
    for shift in mons(degree-w):
     count+=1;v=np.zeros(len(mon),dtype=np.int64)
     for a,c in ts:v[index[tuple(x+y for x,y in zip(a,shift))]]=c
     for p,b in zip(piv,basis):
      if v[p]:v=(v-int(v[p])*b)%P
     nz=np.flatnonzero(v)
     if len(nz):
      p=int(nz[0]);c=int(v[p]);v=v*pow(c,-1,P)%P;basis.append(v);piv.append(p);desc.append([i,list(shift)]);pv.append(c)
   det=1
   for z in pv:det=det*z%P
   assert len(mon)in[280,462]and len(basis)==len(mon)-5 and det
   certificates.append(dict(weight=degree,rows=count,columns=len(mon),rank=len(basis),minor_mod=det,selected_rows=desc,selected_columns=piv,monomials=mon))
  assert [c['minor_mod']for c in certificates]==[31741,12056]
  # The positive grading is computed from the original homogeneous monomials.
  one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2]
  torus=[0,1,-1,0,0];ww=[]
  for j,index in enumerate(m.certificate['raw_free_indices']):
   pair,mon=BASIS[index];ep=s.Poly(mon,*X).monoms()[0];raw=sum(torus[i]for i in pair)-sum(a*b for a,b in zip(torus,ep));ww.append(m.reduced['parameter_weights'][j]-raw)
  qw=[5-ww[two[j]]for j in meta['chart']]+[5];assert max(qw)==6
  rw=[]
  for f in eq:
   ws={sum(5-ww[two[j]]for j in mon)+5*(2-len(mon))for mon,c in f.items()if c};assert len(ws)<=1;rw.append(next(iter(ws))if ws else None)
  assert max(w for w in rw if w is not None)==10
  assert all(not T[i,j]or qw[i]==weights[j]for i in range(45)for j in range(7))
  assert all(not U[i,j]or qw[i]==rw[sel[j]]for i in range(45)for j in range(36))
  # Actual D and C, rather than an initial-ideal surrogate.
  pd='reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json';dr=json.loads((ROOT/pd).read_text());assert dr['all_passed']
  point=s.zeros(61,1);point[4]=point[48]=1;point[12]=r;point[37]=-1;us=s.symbols('u1:62');sub=dict(zip(us,point))
  K=s.Matrix(dr['kernel']).applyfunc(s.sympify).subs(sub);C=s.Matrix(dr['actual_C_tangent']).applyfunc(s.sympify).subs(sub);N=K[:,dr['squarezero_normal_columns']];B=C.row_join(N);cut=[19,26]
  assert C[cut,:]==s.Matrix([[0,0,1],[-2,0,0]])
  H=B*s.Matrix.hstack(*B[cut,:].nullspace());Cc=C*s.Matrix.hstack(*C[cut,:].nullspace());qweights=[0 if i==25 else qw[i-(i>25)]for i in range(45)]
  chars=lambda M:{w:M[[i for i,a in enumerate(qweights)if a==w],:].rank()for w in sorted(set(qweights))if M[[i for i,a in enumerate(qweights)if a==w],:].rank()}
  assert chars(H)=={1:3,2:2}and chars(Cc)=={1:1}
  z,b,beta=primitive_beta(m);v=s.Matrix([beta[i,j].diff(z[k])for(i,j),k in COORDS]);base=v.jacobian(b);br=list(base.T.to_DM().rref()[1]);left=base[br,:].inv();sp=sp_basis();act=s.Matrix.hstack(*[left*infinitesimal_action(v,a)[br,:]for a in sp]);ar=act.subs(dict(zip(b,[point[i]for i in one])))
  coeff=s.Matrix.hstack(s.eye(10)[:,4],-r*s.eye(10)[:,7]+s.eye(10)[:,8]/2+s.eye(10)[:,9]);assert ar*coeff==s.zeros(16,2)
  normal=ar[[one.index(3),one.index(11)],:];assert normal[:,[5,6]]==s.Matrix([[-1,2],[-1,r]])
  gs=[sp[4],-r*sp[7]+sp[8]/2+sp[9]];assert gs[0]*gs[1]==gs[1]*gs[0]
  p=s.zeros(61,1);p[37]=-1;Q0=m.reconstruct(list(p))-X[0]**2*J5;orbit=[]
  for a in gs+[sp[6]+2*sp[5]]:
   G=s.zeros(5);G[1:,1:]=a;Ax=G*s.Matrix(X);dQ=G*Q0+Q0*G.T-sum((Ax[i]*Q0.diff(X[i])for i in range(5)),s.zeros(5));nn=m.normalize(X[0]**2*J5+dQ);q=s.Matrix([nn['parameters'][i]for i in two]);q+=q[25]*s.Matrix([-int(j==25)for j in range(45)]);orbit.append(q)
  orbit=s.Matrix.hstack(*orbit);assert orbit[cut,:2].det()==2
  assert (orbit[:,2]+C[:,1]).applyfunc(s.cancel)==s.zeros(45,1)
  assert (normal*s.Matrix([0,0,0,0,0,2,1,0,0,0]))==s.Matrix([0,r-2])
  assert max(-ww[i]for i in one if i not in[3,11])==0
  nu=[24,33,34,35,36,47,49];assert[ww[i]for i in nu]==[3,3,1,2,2,2,1]
  sample=point.subs(r,1)
  for k,i in enumerate(nu):sample[i]=k+1
  source=s.eye(16)[:,[one.index(i)for i in [12]+nu]].row_join(act.subs(dict(zip(b,[sample[i]for i in one]))))
  sourcecols=list(source.to_DM().rref()[1]);assert len(sourcecols)==16
  sourcedet=source[:,sourcecols].det(method='domain-ge');assert sourcedet
  paths=[str(Path(__file__).relative_to(ROOT)),NOTE,'code/s61_flag_ordinary_weighted_core_engine_20261102.py','code/s61_flag_ordinary_weighted_rows_20261102.cpp','code/parametric_implicit_jets_20261101.cpp','code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py',pd,'research/notes/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.md','research/notes/S61_FLAG_UNIFORM_OPEN_20261101.md','research/notes/S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102.md','research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md','research/notes/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.md','reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json','reports/S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_20261102.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
  paths+=list(dr['input_sha256']);paths+=[str((m.source/f).relative_to(ROOT))for f in['symplectic_slice_certificate.json','symplectic_slice_equations.json','symplectic_slice_reduced_equations.json']]
  hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in dict.fromkeys(paths)}
  for p,h in dr['input_sha256'].items():assert hashes[p]==h,p
  return dict(id=STEM,all_passed=True,field='QQ(r)',core_beta='u5=u49=1,u13=r,other weight-one coordinates zero',frame=meta,raw_coordinate_weights=qw,row_weights=rw,normal_weights=weights,all_implicit_coordinates_above_degree6_zero=True,full_positive_weight_model=True,residual_degree_bound=10,primitive_row_data=rowinfo,rank_witnesses=certificates,prime=P,witness_core=1,D_cut_characters={'1':3,'2':2},C_cut_characters={'1':1},D_ideal_generator_weights=[2,2,2,2,2,3,3,3,3,4,4,4],actual_kernel_generated_weight_at_most=4,actual_kernel_zero_weight7_and8_over_function_field=True,actual_kernel_weight_bounds=[2,6],ordinary_C_derivative='r-2 up to inverse-action sign',two_orbit_cut_minor='2',relative_compensation_minor='2-r',source_sample_beta=[str(sample[i])for i in one],source_minor_columns=sourcecols,source_minor=str(sourcedet),relative_t_kernel_zero_on_explicit_minor_open=True,relative_full_t_squared_zero_on_minor_open=True,all_nuisance_by_positive_contraction=True,all_core_claimed=False,minor_open='(r-2)*(r+6)*(5r+12)*(20r+119)*(25r+4)*Delta7(r)*Delta8(r)!=0; Delta_d is the selected actual weighted-row minor recorded here',minors_expanded_over_QQ=False,finite_field_support_exploration_used_as_proof=False,proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();result=verify()
 with a.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=result['seconds'],ranks=[c['rank']for c in result['rank_witnesses']])))
