#!/usr/bin/env python3
"""ONE exact existing-cache worker: m0, radical X0/X1, and constant W84.
No new source/jet solve and no evaluation of original regular r/g quotients.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_native128_first_20261104_v63 import Algebra
from ambient_B18_T_native_S2_B11_20261104_v63 import Finite
from ambient_B18_honest_nil_T_first_values_preflight_20261104_v63 import PATHS
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def nqq(x):
 if isinstance(x,Q):return 1
 if isinstance(x,dict):return sum(nqq(v)for v in x.values())
 if isinstance(x,(list,tuple)):return sum(nqq(v)for v in x)
 return 0
def invq(a):
 n=len(a);m=[r[:]+[Q(i==j)for j in range(n)]for i,r in enumerate(a)];det=Q(1)
 for j in range(n):
  p=next(i for i in range(j,n)if m[i][j])
  if p!=j:m[j],m[p]=m[p],m[j];det=-det
  c=m[j][j];det*=c;m[j]=[x/c for x in m[j]]
  for i in range(n):
   if i!=j and m[i][j]:
    c=m[i][j];m[i]=[x-c*y for x,y in zip(m[i],m[j])]
 return [r[n:]for r in m],det
def pivots(cols):
 basis={};chosen=[]
 for j,c in enumerate(cols):
  v={i:x for i,x in enumerate(c)if x}
  for p,b in sorted(basis.items()):
   if p in v:
    a=v[p]
    for i,x in b.items():
     y=v.get(i,Q(0))-a*x
     if y:v[i]=y
     elif i in v:del v[i]
  if v:
   p=min(v);a=v[p];basis[p]={i:x/a for i,x in v.items()};chosen.append(j)
 return chosen,list(sorted(basis))
def mv(a,v):return [sum((x*y for x,y in zip(r,v)),Q(0))for r in a]
def polynomial(a,p,vs):
 vals=[]
 for m,c in p:
  v=a.c(Q(c))
  for x,e in zip(vs,m):
   if e:v=a.mul(v,a.power(x,e))
  vals.append(v)
 return a.add(*vals)
def matrix_inverse(a,M):
 n=len(M);B=[r[:]+[a.c(i==j)for j in range(n)]for i,r in enumerate(M)]
 for j in range(n):
  p=next(i for i in range(j,n)if B[i][j][0])
  B[j],B[p]=B[p],B[j];t=a.inv(B[j][j]);B[j]=[a.mul(x,t)for x in B[j]]
  for i in range(n):
   if i!=j and B[i][j]!=a.Z:
    t=B[i][j];B[i]=[a.add(x,a.neg(a.mul(t,y)))for x,y in zip(B[i],B[j])]
 return [r[n:]for r in B]
def mm(a,A,B):return [[a.add(*(a.mul(x,y)for x,y in zip(r,c)))for c in zip(*B)]for r in A]
def comp(v,i):return tuple(v[4*k+i]for k in range(len(v)//4))
def sparse(v):return [[i,x]for i,x in enumerate(v)if x]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);saved=0;bytecount=0
 def save(name,data):
  nonlocal saved,bytecount
  count=nqq(data);assert saved+count<150000,('saved_QQ_guard',saved,count)
  txt=json.dumps(enc(data),separators=(',',':'))+'\n';assert bytecount+len(txt.encode())<32<<20
  (out/name).write_text(txt);saved+=count;bytecount+=len(txt.encode())
  (out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ_slots':saved,'seconds':time.monotonic()-start})+'\n')
  print(json.dumps({'phase':name,'seconds':time.monotonic()-start}),flush=True)
 ds={k:read(ROOT/p)for k,p in PATHS.items()}
 pfpath='reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_PREFLIGHT_20261104_v63/report.json';pf=read(ROOT/pfpath)
 assert pf['all_passed']and pf['source_visible_QQ_bound']==128812
 for p,h in pf['input_sha256'].items():assert sha(ROOT/p)==h
 mono1path='reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase1-source-firsts.json';N1=read(ROOT/mono1path)
 R=ds['radical'];V=ds['nil'];G=ds['native1'];K=ds['Koszul']
 a=Algebra(8,2,[list(map(Q,v))for v in N1['nu_first9']])
 b=Finite(2,[list(map(Q,v))for v in ds['monicB11']['nu_A_B_B11']])
 def actual_variables(alg,kappa,G8):
  core=[tuple(map(Q,kappa[i]))for i in [0,1,2,6]]
  GL={s:polynomial(alg,p,core)for s,p in V['forward_GL4'].items()}
  return [tuple(map(Q,v))for v in G8[:5]]+[GL[s]for s in ['h','c','r','l']]
 va=actual_variables(a,ds['proper1']['kappa9'],G['G8'])
 vb=actual_variables(b,ds['properB11']['kappa9_C_B44'],ds['properB11']['actualG8_C_B44'])
 def qcols(alg,vs):
  d=polynomial(alg,R['d_polynomial'],vs)
  native=[alg.mul(d,polynomial(alg,p,vs))for p in R['R9_generators_polynomials9']]
  return [alg.mul(v,basis(alg,j))for v in native for j in range(4)]
 def basis(alg,j):v=[Q(0)]*alg.n;v[j]=Q(1);return tuple(v)
 qa=qcols(a,va);qb=qcols(b,vb)
 # Actual complete source U8=beta8*tauA; gA order3 vanishes in both finite coefficient rings.
 for i in range(8):
  za=[Q(0)]*a.n;zb=[Q(0)]*b.n
  if i<7:
   za[4*(i+1)+1]=-Q(1)
   mon=[0]*8;mon[i]=1
   mi=ds['B11_basis']['B11_basis_monomials8']
   if mon in mi:zb[4*mi.index(mon)+1]=-Q(1)
  qa.append(tuple(za));qb.append(tuple(zb))
 gamma=lambda v:[comp(v,1),comp(v,3),tuple(x-Q(11,12)*y for x,y in zip(comp(v,2),comp(v,1)))]
 gv=[gamma(v)for v in qa];gb=[gamma(v)for v in qb]
 assert all(not q[0]for col in gv for q in col)
 # Direct source m0 values, without evaluating any original regular division.
 m0=[[Q(0)]*44 for _ in range(24)]
 for j,(v,w)in enumerate(zip(gv,gb)):
  for sector in range(3):
   for k in range(7):m0[8*sector+k][j]=v[sector][k+1]
   if sector<2:m0[8*sector+7][j]=v[sector][8]
   else:
    assert not v[sector][8];m0[23][j]=w[sector][2]
 mcols=list(map(list,zip(*m0)));E,mrows=pivots(mcols);assert len(E)==24
 Mi,mdet=invq([[m0[i][j]for j in E]for i in range(24)])
 assert all(sum((m0[i][E[k]]*Mi[k][j]for k in range(24)),Q(0))==Q(i==j)for i in range(24)for j in range(24))
 save('phase1-actual-m0-and-section.json',{'actual_frame':K['actual_S_variable_order'],'m0_24_by44':m0,'rank':24,
  'selected_source_columns24':E,'selected_minor_determinant':mdet,'selected_minor_inverse':Mi,
  'Vgamma_firsts_all44':[[sparse(v)for v in col]for col in gv],
  'third_gamma_r_squared_all44':[q[2][2]for q in gb],'n0_is_E_times_inverse':True,
  'all_gA_order3_and_no_quotient_derivative_evaluation':True})
 # Exact actual 12x9 source adapter in C36, from SAVED N12/G8/coordinates.
 N=[[tuple(map(Q,v))for v in row]for row in G['N']];Ni=matrix_inverse(a,N)
 ident=mm(a,N,Ni);assert all(ident[i][j]==a.c(i==j)for i in range(12)for j in range(12))
 D=[[polynomial(a,p,va)for p in row]for row in R['G8_polynomial_factor_D8_by9']]
 A=mm(a,[r[:8]for r in Ni],D)
 # Prove all old/new12 row identities with no coefficient solve.
 oldF=[tuple(map(Q,v))for v in ds['proper1']['physicalF12']]
 newQ=qa[:36:4]
 assert all(x==y for x,y in zip([r[0]for r in mm(a,A,[[x]for x in newQ])],oldF))
 del Ni,N,ident,D
 zero=(Q(0),)*9
 def add9(x,y):return tuple(u+v for u,v in zip(x,y))
 def mul9(x,y):return (x[0]*y[0],)+tuple(x[0]*y[k]+x[k]*y[0]for k in range(1,9))
 X=[{}for _ in range(640)]
 index=collections.defaultdict(list)
 for j,col in enumerate(ds['Bfirst']['source_columns640_rows56_sparse_dual9']):
  for row,sp in col:
   v=[Q(0)]*9
   for k,c in sp:v[k]=Q(c)
   if row<48:index[row].append((j,tuple(v)))
   else:X[j][row-12]=tuple(v)
 for row,entries in sorted(index.items()):
  i,j=divmod(row,4);res=[]
  for k in range(9):
   v=a.mul(A[i][k],basis(a,j));res.extend(comp(v,t)for t in range(4))
  for c,v in entries:
   for k,w in enumerate(res):
    value=mul9(w,v)
    if value!=zero:
     q=add9(X[c].get(k,zero),value)
     if q!=zero:X[c][k]=q
     elif k in X[c]:del X[c][k]
 # Complete 140 radical syzygies, not the old128 native syzygies.
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  nc={i:polynomial(a,p,va)for i,p in col}
  for j in range(4):
   dest={}
   for i,v in nc.items():
    w=a.mul(v,basis(a,j))
    for k in range(4):
     p=comp(w,k)
     if p!=zero:dest[4*i+k]=p
   X.append(dest)
 assert len(X)==780
 xsave=[[[i,sparse(v)]for i,v in sorted(c.items())]for c in X]
 assert sum(len(v)for c in xsave for i,v in c)<=20356
 X0=[[c.get(i,zero)[0]for i in range(44)]for c in X]
 # Complete source cycle forces this stronger equality, including the r² row.
 assert all(not x for col in X0 for x in mv(m0,col))
 xE,xrows=pivots(X0);r=len(xE);assert r<=20
 Xi,xdet=invq([[X0[j][i]for j in xE]for i in xrows])if r else([],Q(1))
 fixed=[X0[j]for j in xE]
 coeffs=[]
 for col in X0:
  c=mv(Xi,[col[i]for i in xrows])if r else[]
  assert all(sum((fixed[k][i]*c[k]for k in range(r)),Q(0))==col[i]for i in range(44))
  coeffs.append(c)
 save('phase2-actual-radical-X-firsts.json',{'all780_source_columns44_sparse_dual9':xsave,
  'new_A12_by9_C36':A,'actual_X0_rank':r,'X0_pivot_columns':xE,'X0_pivot_rows':xrows,
  'X0_fixed_minor_determinant':xdet,'X0_fixed_minor_inverse':Xi,'all780_X0_factorization_coefficients':coeffs,
  'm0_X0_exact_zero_including_third_r_squared_row':True,'all12_F12_equals_A_dR9_exact_C36':True})
 # After X0 cancellation, m1*X0 vanishes. Read constant wedge from m0*X1 only.
 pair_basis=K['basis_sector_wedge'][2];assert len(pair_basis)==84
 W=[];original=[];verified=0
 for j in range(780):
  if j in xE:continue
  c=coeffs[j];res={}
  for i in range(44):
   v=[X[j].get(i,zero)[k+1]-sum((c[t]*X[xE[t]].get(i,zero)[k+1]for t in range(r)),Q(0))for k in range(8)]
   if any(v):res[i]=v
  v1=[[sum((m0[i][k]*v[t]for k,v in res.items()),Q(0))for t in range(8)]for i in range(24)]
  w=[]
  for spec in pair_basis:
   sector,I=spec if isinstance(spec,list)else(spec['sector'],spec['wedge'])
   i,k=I;w.append(v1[8*sector+k][i])
  replay=[[Q(0)]*8 for _ in range(24)]
  for z,col in zip(w,K['all8_differentials_sparse_columns'][1]):
   if z:
    for row,mon,sign in col:
     if sum(mon)==1:replay[row][mon.index(1)]+=Q(sign)*z
  assert replay==v1,('CONSTANT_KOSZUL_LIFT_LINEAR_IDENTITY',j)
  # Explicitly includes component7*r in sector2, which must be zero.
  assert not v1[23][7]
  W.append(w);original.append(j);verified+=24*8
 wp,wr=pivots(W);rankW=len(wp)
 Wi,wdet=invq([[W[j][i]for j in wp]for i in wr])if rankW else([],Q(1))
 factorhash=hashlib.sha256()
 for col in W:
  c=mv(Wi,[col[i]for i in wr])if rankW else[]
  assert all(sum((W[wp[k]][i]*c[k]for k in range(rankW)),Q(0))==col[i]for i in range(84))
  factorhash.update((json.dumps(enc(sparse(c)),separators=(',',':'))+'\n').encode())
 mu=104-r-rankW
 save('phase3-actual-W84-rank-and-mu.json',{'residual_original_X_columns':original,'all_W84_columns_sparse':[sparse(v)for v in W],
  'all_Koszul_linear_coefficients_verified':verified,'all_third_sector_component7_r_coefficients_zero':True,
  'W_rank':rankW,'W_pivot_residual_columns':wp,'W_pivot_rows':wr,'W_minor_determinant':wdet,
  'W_fixed_minor_inverse':Wi,'all_W_column_factorizations_verified':True,
  'W_factorization_stream_sha256':factorhash.hexdigest(),'mu_S_Tnil':mu,
  'formula':'104-rankX0-rankW0','m1_and_original_regular_quotient_firsts_not_evaluated':True})
 paths=list(PATHS.values())+[pfpath,mono1path];note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_honest_Tnil_minimal_generator_count_from_existing_cache',
  'm0_rank':24,'actual_X0_rank':r,'actual_W0_rank':rankW,'mu_R0_Tnil':mu,
  'all780_source_X0_X1_built_and_complete':True,'all24x8_Koszul_linear_equations_checked_each_residual':True,
  'r_squared_sector_r_coefficient_checked_not_assumed':True,
  'new_source_jet_IFT_GB_or_old_r_g_division':False,'exact_Tnil_length_Nil_generators_or_cube_computed':False,
  'visible_saved_QQ_slots':saved,'output_bytes':bytecount,'input_sha256':{p:sha(ROOT/p)for p in paths},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{p.name:sha(p)for p in out.glob('phase[123]-*.json')},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report),flush=True)
if __name__=='__main__':main()
