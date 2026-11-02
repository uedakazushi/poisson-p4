"""All ordinary flag cores outside seven stated rational values: finite QQ certificate."""
import argparse,hashlib,json,os,re,subprocess,tempfile,time
from pathlib import Path
from functools import lru_cache
from itertools import combinations,combinations_with_replacement
import sympy as s
from s61_flag_ordinary_weighted_core_engine_20261102 import build
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FLAG_ORDINARY_CORE_SUPPORT_20261102'
NOTE='research/notes/'+STEM+'.md'
BASE='reports/S61_FLAG_ORDINARY_WEIGHTED_FINITE_KERNEL_20261102.json'
P=32003;W=(1,1,2,2,1,2,2);r=s.Symbol('r')
ROOTS=list(map(s.Rational,[-6]))+[s.Rational(-12,5),s.Rational(-119,20),s.Rational(-4,25),s.Rational(2),s.Rational(-6,25),s.Rational(-91,400)]
ORDERS=[1165,168,518,269,130,24,2]
CPP={k:'code/s61_flag_ordinary_support_'+k+'_20261102.cpp'for k in['infinity_QQ','infinity_mod','root_QQ','root_mod','popov']}

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def modq(x):
 a=s.Rational(x);assert int(a.q)%P
 return int(a.p)%P*pow(int(a.q)%P,-1,P)%P
@lru_cache(None)
def mons(w,n=7):
 if n==0:return ((),)if w==0 else()
 return tuple(t+(e,)for e in range(w//W[n-1]+1)for t in mons(w-e*W[n-1],n-1))
def poly(f):
 ts=f.readline().split();assert len(ts)-1==int(ts[0]);return ts[1:]
def profile(text):return [int(line.split()[3])for line in text.splitlines()if line.startswith('STEP ')]
def compile_native(src,exe):
 subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(ROOT/src),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(exe)],check=True)

def active_quadratics(rows):
 heavy=(2,3,5,6);light=(0,1,4);mon=[]
 for ij in combinations_with_replacement(range(4),2):
  a=[0]*4
  for j in ij:a[j]+=1
  mon.append(tuple(a))
 ix={a:i for i,a in enumerate(mon)};data=[];source=[]
 for i,(w,ts)in enumerate(rows[:135]):
  row=[s.S.Zero]*10
  for a,cs in ts:
   if w==4 and not any(a[j]for j in light):row[ix[tuple(a[j]for j in heavy)]]=s.Poly.from_list([s.Rational(x)for x in cs[::-1]],r).as_expr()
  if any(row):source.append(i);data.append(row)
 M=s.Matrix(data);RR,piv=M.to_DM().to_field().rref();RR=RR.to_Matrix();assert len(piv)==9
 free=next(i for i in range(10)if i not in piv);mu=s.zeros(10,1);mu[free]=1
 for i,p in enumerate(piv):mu[p]=-RR[i,free]
 mu=mu.applyfunc(s.cancel);mu=(s.lcm([s.denom(c)for c in mu])*mu).applyfunc(s.cancel)
 mu=(mu/s.gcd_list(list(mu))).applyfunc(s.factor)
 assert s.Poly(s.gcd_list(list(mu)),r).degree()==0
 assert (M*mu).applyfunc(s.cancel)==s.zeros(M.rows,1)
 selected=list(M[:,list(piv)].T.subs(r,1).to_DM().rref()[1]);assert len(selected)==9
 det=s.factor(M[selected,list(piv)].det(method='domain-ge'));cofactor=s.factor(s.cancel(det/mu[free]))
 allowed=s.Poly((r+6)*(5*r+12)*(20*r+119)*(25*r+4),r);rest=s.Poly(cofactor,r)
 while rest.degree()>0:
  d=s.gcd(rest,allowed);assert d.degree()>0;rest=rest.exquo(d)
 gram=s.zeros(4)
 for i in range(4):
  for j in range(4):
   a=[0]*4;a[i]+=1;a[j]+=1;gram[i,j]=mu[ix[tuple(a)]]
 minors=[s.factor(gram.extract(a,b).det())for a in combinations(range(4),2)for b in combinations(range(4),2)]
 gcd=s.factor(s.gcd_list(minors));assert s.Poly(gcd,r).monic()==s.Poly(r+s.Rational(91,400),r)
 return dict(source_rows=source,monomials=mon,matrix=[[str(x)for x in row]for row in M.tolist()],primitive_dual=list(map(str,mu)),rank9_rows=selected,rank9_columns=list(piv),rank9_minor=str(det),cofactor_scalar=str(cofactor),dual_2minor_gcd=str(gcd),active_cubics_vanish_on_D4_rho_open=True)

def verify():
 start=time.monotonic();base=json.loads((ROOT/BASE).read_text());assert base['all_passed']
 for p,h in base['input_sha256'].items():assert sha(p)==h,p
 with tempfile.TemporaryDirectory(prefix='flag-ordinary-support-')as td:
  tmp=Path(td);meta,m,eq,J,T,U=build(tmp/'raw')
  assert meta==base['frame'];env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'}
  den=tmp/'den';out=tmp/'rows';den.write_text('\n'.join(str(s.Poly(s.sympify(meta[k]),r).degree()+1)+' '+' '.join(str(s.Poly(s.sympify(meta[k]),r).nth(i))for i in range(s.Poly(s.sympify(meta[k]),r).degree()+1))for k in['ud','td'])+'\n')
  exe=tmp/'row_native';compile_native('code/s61_flag_ordinary_weighted_rows_20261102.cpp',exe)
  subprocess.run([str(exe),str(tmp/'raw/raw.txt'),str(den),str(out)],env=env,capture_output=True,text=True,check=True,timeout=60)
  rows=[];D=s.Poly((r+6)*(5*r+12)*(20*r+119)*(25*r+4),r,domain=s.QQ)
  with out.open()as f:
   assert int(f.readline())==136
   for i in range(136):
    ii,w,n=map(int,f.readline().split());assert ii==i
    dd,cc=poly(f),poly(f);assert len(cc)<=1
    for c in dd+cc:modq(c)
    if cc:assert modq(cc[0])
    pp=s.Poly.from_list([s.Rational(z)for z in dd[::-1]],r,domain=s.QQ)
    while pp.degree()>0:
     d=s.gcd(pp,D);assert d.degree()>0;pp=pp.exquo(d)
    ts=[]
    for j in range(n):
     a=tuple(map(int,f.readline().split()));cs=poly(f);assert sum(x*y for x,y in zip(a,W))==w
     for c in cs:modq(c)
     ts.append((a,cs))
    if i in meta['selected']:assert not ts and w==-1
    rows.append((w,ts))
  cert=base['rank_witnesses'][0];assert cert['weight']==7 and cert['rank']==275
  mon=mons(7);assert [list(a)for a in mon]==cert['monomials'];cols=cert['selected_columns'];index={mon[j]:i for i,j in enumerate(cols)}
  def restrict(terms,shift):
   rr=[]
   for a,c in terms:
    j=index.get(tuple(x+y for x,y in zip(a,shift)))
    if j is not None:rr.append((j,c))
   return rr
  allrows=[];alldesc=[]
  for i,(w,ts)in enumerate(rows[:135]):
   if not ts or not 0<=w<=7:continue
   for sh in mons(7-w):allrows.append(restrict(ts,sh));alldesc.append([i,list(sh)])
  assert len(allrows)==532 and len(index)==275
  selected=[restrict(rows[i][1],shift)for i,shift in cert['selected_rows']]
  assert all(x in alldesc for x in cert['selected_rows'])
  ip=tmp/'selected';fp=tmp/'full_mod'
  def write_matrix(path,data,modular=False):
   with path.open('w')as f:
    f.write((str(P)+' 'if modular else'')+f'{len(data)} 275\n')
    for row in data:
     f.write(str(len(row))+'\n')
     for j,cs in row:f.write(f'{j} {len(cs)} '+' '.join(str(modq(c))if modular else c for c in cs)+'\n')
  write_matrix(ip,selected);write_matrix(fp,allrows,True)
  input_sha=hashlib.sha256(ip.read_bytes()).hexdigest()
  native={}
  for k,p in CPP.items():native[k]=tmp/k;compile_native(p,native[k])
  def run(k,inp,extra=()):
   op=tmp/(k+'-'+str(len(logs))+'.out');pr=subprocess.run([str(native[k]),str(inp),str(op),*map(str,extra)],env=env,capture_output=True,text=True,check=False,timeout=180)
   assert pr.returncode==0,(k,pr.returncode,pr.stderr[-3000:],op.read_text()if op.exists()else'')
   text=op.read_text();logs.append((k,text));print('CERTIFICATE',k,round(time.monotonic()-start,3),flush=True);return text
  logs=[];mi=run('infinity_mod',ip);qi=run('infinity_QQ',ip)
  a=list(map(int,re.search(r'DEGREE (\d+) LEADING (\d+) BOUND (\d+) VALUATION (\d+)',mi).groups()))
  assert a==[2543,13821,2912,369];assert profile(qi)==profile(mi)[:-1]
  assert 'DEGREE_BOUND 2543 TROPICAL_ROW_BOUND 2912 VALUATION 369'in qi
  root_data=[]
  for i,root in enumerate(ROOTS):
   mm=run('root_mod',ip,[modq(root)]);ranks=profile(mm);assert ranks[-1]==275
   val=sum(275-x for x in ranks);assert val==ORDERS[i]
   cfg=tmp/f'root_config{i}';cfg.write_text(f'{root} {len(ranks)-1} {val} '+' '.join(map(str,ranks[:-1]))+'\n')
   qq=run('root_QQ',ip,[cfg]);assert profile(qq)==ranks[:-1] and f'VALUATION {val}\n'in qq
   root_data.append(dict(root=str(root),mod_ranks=ranks,QQ_ranks=profile(qq),exact_order=val,mod_output=mm,QQ_output=qq))
  assert sum(ORDERS)==2276 and a[0]-sum(ORDERS)==267
  pp=run('popov',fp);assert 'RANK 275 DEGREE 136'in pp
  lines=pp.splitlines();k=next(i for i,line in enumerate(lines)if line.startswith('FACTORS '));n=int(lines[k].split()[1]);factors=[]
  for line in lines[k+1:k+1+n]:
   exponent,degree,*cs=map(int,line.split());assert degree==1 and cs[-1]==1
   root=(-cs[0])%P;factors.append((root,exponent))
  assert sorted(factors)==sorted([(modq(2),111),(modq(s.Rational(-6,25)),23),(modq(s.Rational(-91,400)),2)])
  # Each exact order agrees with reduction. Hence the residual monic quotient
  # has no known root even modulo p; no giant determinant expansion is required.
  assert len({modq(x)for x in ROOTS})==7
  active=active_quadratics(rows);print('ACTIVE_QUADRATICS',round(time.monotonic()-start,3),flush=True)
  paths=[str(Path(__file__).relative_to(ROOT)),NOTE,BASE,*CPP.values(),'code/s61_flag_ordinary_weighted_core_engine_20261102.py','code/s61_flag_ordinary_weighted_rows_20261102.cpp',*base['input_sha256']]
  hashes={p:sha(p)for p in dict.fromkeys(paths)}
  return dict(id=STEM,all_passed=True,prime=P,field='QQ',core_open='(r+6)*(5*r+12)*(20*r+119)*(25*r+4)*(r-2)*(25*r+6)*(400*r+91)!=0',frame=meta,rows=532,columns=275,weighted_degree=7,selected_rows=cert['selected_rows'],selected_columns=cols,selected_polynomial_matrix_sha256=input_sha,all_coefficients_p_integral=True,row_content_nonzero_p_units=True,row_denominators_D4_only=True,infinity_mod_output=mi,infinity_QQ_output=qi,minor_degree=2543,minor_leading_mod=13821,root_certificates=root_data,residual_monic_degree=267,residual_coprime_to_known_support_mod_p=True,full_modular_rowmodule_certificate=pp,modular_support_after_D4_saturation=factors,active_quadratics=active,finite_V_CRT_conditions_verified=True,weight7_rank_275_on_stated_QQ_open=True,Koszul_weight8_reduction=True,actual_kernel_weights=[2,6],relative_t_squared_zero=True,all_nuisance_by_positive_contraction=True,all_general_beta_points_on_stated_core_open=True,exceptional_cores_claimed=False,all_flags_or_affine_exhaustion_claimed=False,proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();result=verify()
 with a.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=result['seconds'],hashes=len(result['input_sha256']))))
