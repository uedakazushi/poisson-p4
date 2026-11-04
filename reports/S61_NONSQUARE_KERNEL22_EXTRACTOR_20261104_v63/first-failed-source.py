#!/usr/bin/env python3
"""One120s constructive extractor job; no coefficient/minor/GB search."""
import argparse,gzip,hashlib,itertools,json,resource,shlex,signal,subprocess,tarfile,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for z in iter(lambda:f.read(1024*1024),b''):h.update(z)
 return h.hexdigest()
def mons(n,d):return[(d,)]if n==1 else[(i,)+v for i in range(d+1)for v in mons(n-1,d-i)]

class Circuit:
 def __init__(self,R,K,known):
  self.R,self.K,self.known=R,K,known;self.atoms=[];self.ai={};self.s=[['z']];self.si={('z',):0};self.l=[['z']];self.li={('z',):0};self.one=self.c(K.one)
 def c(self,c):
  c=self.K(c)
  if not c:return 0
  p=c.denom
  for f in self.known:
   while not p.rem(f):p=p.exquo(f)
  assert p.is_ground and p
  key=str(c)
  if key not in self.ai:self.ai[key]=len(self.atoms);self.atoms.append(c)
  return self.sc(['c',self.ai[key]])
 def sc(self,node):
  key=tuple(node)
  if key not in self.si:self.si[key]=len(self.s);self.s.append(node)
  return self.si[key]
 def add(self,a,b):
  if not a:return b
  if not b:return a
  return self.sc(['a',a,b])
 def sub(self,a,b):
  if not b:return a
  if a==b:return 0
  return self.sc(['s',a,b])
 def mul(self,a,b):
  if not a or not b:return 0
  if a==self.one:return b
  if b==self.one:return a
  return self.sc(['m',a,b])
 def sum(self,ts):
  z=0
  for a,b in ts:z=self.add(z,self.mul(a,b))
  return z
 def var(self,i):
  key=('v',i)
  if key not in self.li:self.li[key]=len(self.l);self.l.append(['v',i])
  return self.li[key]
 def lin(self,terms):
  z={}
  for c,n in terms:
   if c and n:z[n]=self.add(z.get(n,0),c)
  if not z:return 0
  if len(z)==1:
   n,c=next(iter(z.items()))
   if c==self.one:return n
  ts=tuple((c,n)for n,c in sorted(z.items()));key=('l',ts)
  if key not in self.li:self.li[key]=len(self.l);self.l.append(['l',[list(t)for t in ts]])
  return self.li[key]
 def encoded(self):
  enc=lambda p:[[list(ex),str(c)]for ex,c in p.items()]
  return dict(atoms_num_den=[[enc(c.numer),enc(c.denom)]for c in self.atoms],scalar_nodes=self.s,linear_nodes=self.l)

def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(120);out.mkdir(parents=True,exist_ok=False)
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4
 o=Circuit(R,K,[a,d,F,B,C,J]);frac=lambda x:K.from_expr(s.sympify(x));dec=lambda p:K(R.from_dict({tuple(ex):s.QQ(c)for ex,c in p}))
 paths=dict(K22=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json',
  K22review=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63_REVIEW/report.json',
  q=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json',
  anchors=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json',
  pairarchive=ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63/original_ZW_rank19_generic146_exact.tar.gz',
  pairreview=ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW/report.json',
  ellchart=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json',
  shadow=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/right_inverse.txt.gz',
  shadowreview=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v63_REVIEW/report.json')
 assert read(paths['K22review'])['all_passed'];kd=read(paths['K22']);ad=read(paths['anchors']);qd=read(paths['q'])
 T,*x=ring('x0,x1,x2,x3,x4,Z,W',K);x0,x1,x2,x3,x4,Z,W=x
 load=lambda p:T.from_dict({tuple(m):dec(c)for m,c in p})
 L,ell,TZ,TW=[load(ad[k])for k in ['L','ell','TZ','TW']]
 xm=mons(7,2);mi={m:i for i,m in enumerate(xm)};unitmon=lambda i:tuple(int(i==j)for j in range(7))
 coef=lambda p,i:p.get(unitmon(i),K.zero)
 L1,ell1=coef(L,1),coef(ell,1);assert L1==K(-40*C)and ell1==K(-16*C*J/243)
 g=K(-80*a*d**2/81);Delta=ell-T(K(2*J)/1215)*L
 D=[K.zero,-(3*coef(L,2)-K(7)*coef(L,4)/2)/L1,K(3),K.zero,-K(7)/2]
 assert sum((D[i]*coef(ell,i)for i in range(5)),K.zero)==g
 imgs=[x0,x1,x2,W,x3,T.zero,Z+x1-3*x0,x4];qs=[]
 for row in qd['quadratic_rows_sparse']:
  q=T.zero
  for m,c in row:
   z=T(dec(c))
   for im,k in zip(imgs,m):
    if k:z*=im**k
   q+=z
  qs.append(q)
 ch=[[K.zero]*6 for _ in range(6)];ch[0][0]=ch[1][1]=ch[2][3]=ch[3][5]=K.one
 ch[4][1]=K(10*a+3*e**2-3*e)/K(a);ch[4][3]=-K(3*(e-1))/K(a);ch[4][4]=K(-10)
 ch[5][1]=K(-(d*e+8));ch[5][2]=K(8);ch[5][3]=K(d)
 qp=[qs[0],qs[1],qs[3],qs[5],ell*Z,ell*W]
 assert all(sum((T(c)*q for c,q in zip(row,qs)),T.zero)==p for row,p in zip(ch,qp))
 inv=[[K.zero]*6 for _ in range(6)];inv[0][0]=inv[1][1]=inv[3][2]=inv[5][3]=K.one
 inv[2][1]=K(d*e+8)/8;inv[2][2]=-K(d)/8;inv[2][5]=K.one/8
 inv[4][1]=K(10*a+3*e**2-3*e)/K(10*a);inv[4][2]=-K(3*(e-1))/K(10*a);inv[4][4]=-K.one/10
 assert all(sum((ch[i][k]*inv[k][j]for k in range(6)),K.zero)==(K.one if i==j else K.zero)for i in range(6)for j in range(6))
 # Frozen2pair, with the literal116-nonzero-term unit certificate.
 with tarfile.open(paths['pairarchive'])as tf:
  def member(suffix):
   matches=[m for m in tf.getmembers()if m.name.endswith('/'+suffix)];assert len(matches)==1
   data=tf.extractfile(matches[0]).read();(out/('input-'+suffix)).write_bytes(data);return data
  pairbytes=member('pair-native.txt');bez=json.loads(member('literal-Bezout.json'))
 def natpoly(line):
  t=line.split();n=int(t[0]);assert len(t)==1+4*n
  return K(R.from_dict({tuple(map(int,t[2+4*i:5+4*i])):s.QQ(t[1+4*i])for i in range(n)}))
 lines=pairbytes.decode().splitlines();assert lines[0]=='2 6 0';M=[[natpoly(lines[1+i*6+j])for j in range(6)]for i in range(2)]
 fs=[dec(p)for p in bez['source_scalar5']];bs=[dec(p)for p in bez['scalar_unit_cofactors']]
 assert sum((f*b for f,b in zip(fs,bs)),K.zero)==K.one
 rem=[0,1,3,4,5];pivot=M[0][2];alpha=-K(a*C);assert pivot==K(5*a)
 assert all(M[1][j]-M[1][2]*M[0][j]/pivot==alpha*fs[k]for k,j in enumerate(rem))
 U3,*z=ring('u2,u3,u4',K);zimg=[U3.zero,-sum((U3(coef(L,j))*z[j-2]for j in range(2,5)),U3.zero)/L1,*z]
 def substitute(p,images,target):
  r=target.zero
  for m,c in p.items():
   v=target(c)
   for im,k in zip(images,m):
    if k:v*=im**k
   r+=v
  return r
 pairmons=mons(3,2);bar=[T.from_dict({m:c for m,c in p.items()if not m[5]and not m[6]})for p in qp[:4]]
 pairactual=[substitute(p,[*zimg,U3.zero,U3.zero],U3)for p in [bar[0],bar[3]]]
 assert all(pairactual[i].get(m,K.zero)*(K(C)if i else K.one)==M[i][j]for i in range(2)for j,m in enumerate(pairmons))
 u=[0]*6
 for k,j in enumerate(rem):u[j]=o.c(bs[k]/alpha)
 u[2]=o.sub(0,o.sum([(o.c(M[0][j]/pivot),u[j])for j in rem]))
 pairdual=[[o.sub(o.c(K.one/pivot)if j==2 else 0,o.mul(o.c(M[1][2]/pivot),u[j])),o.mul(o.c(C),u[j])]for j in range(6)]
 checks=[]
 for i in range(2):
  for j in range(2):checks.append([o.sum([(o.c(pairactual[i].get(m,K.zero)),pairdual[k][j])for k,m in enumerate(pairmons)]),int(i==j)])
 # lambda on5-variable quadrics: derivative after projecting offq0,q5.
 m5=mons(5,2);lam=[]
 beta=[]
 for mon in m5:
  b=K.zero
  if mon[1]:
   for j in range(5):
    if mon[j]-int(j==1):b+=mon[1]*(mon[j]-int(j==1))*D[j]/(g*L1)
  beta.append(b)
 bq=[sum((beta[k]*bar[i].get(tuple([*m,0,0]),K.zero)for k,m in enumerate(m5)),K.zero)for i in[0,3]]
 for k,m in enumerate(m5):
  p=T.from_dict({tuple([*m,0,0]):K.one});cut=substitute(p,[*zimg,U3.zero,U3.zero],U3)
  val=o.c(beta[k])
  for t in range(2):val=o.sub(val,o.mul(o.c(bq[t]),o.sum([(o.c(cut.get(mm,K.zero)),pairdual[h][t])for h,mm in enumerate(pairmons)])))
  lam.append(val)
 for p,want in [(p,0)for p in bar]+[(T.from_dict({m:c for m,c in (ell*L).items()if not m[5]and not m[6]}),1)]:
  checks.append([o.sum([(o.c(p.get(tuple([*m,0,0]),K.zero)),lam[k])for k,m in enumerate(m5)]),want])
 # Actual original cofactor174 variables -> prime6 polynomial coefficients.
 X=[{}for _ in range(6)]
 def plin(ts):
  keys=set().union(*(p.keys()for c,p in ts))if ts else set()
  return{m:o.lin([(c,p.get(m,0))for c,p in ts])for m in keys}
 for row,definition in enumerate(kd['actual_source174_definitions']):
  _,qi,m=definition
  if m[5]:continue
  p=T.one
  for im,k in zip(imgs,m):
   if k:p*=im**k
  for mon,c in p.items():X[qi][mon]=o.lin([(o.one,X[qi].get(mon,0)),(o.c(c),o.var(row))])
 Xp=[plin([(o.c(inv[i][j]),X[i])for i in range(6)if inv[i][j]])for j in range(6)]
 n=o.lin([(o.sub(0,lam[k]),Xp[2].get(tuple([*m,0,0]),0))for k,m in enumerate(m5)])
 source=ROOT/'code/s61_nonsquare_kernel22_extractor_20261104_v63.cpp';exe=out/'native';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(source),*flags,'-lgmp','-o',str(exe)],check=True,capture_output=True,text=True)
 km=[[None]*22 for _ in range(174)]
 for j,col in enumerate(kd['sparse_columns']):
  for i,num,den in col:km[i][j]=(num,den)
 def putpoly(f,p):
  f.write(str(len(p)))
  for ex,c in p:f.write(' '+c+' '+' '.join(map(str,ex)))
  f.write('\n')
 enc=lambda p:[[list(ex),str(c)]for ex,c in p.items()]
 def runstage(name,outputs):
  data=o.encoded();(out/(name+'.json')).write_text(json.dumps(dict(circuit=data,scalar_checks=checks,linear_output_expected_basis_columns=outputs),separators=(',',':'))+'\n')
  path=out/(name+'.native')
  with path.open('w')as f:
   f.write(f'{len(o.atoms)} {len(o.s)} {len(o.l)} {len(checks)} {len(outputs)}\n')
   for c in o.atoms:putpoly(f,enc(c.numer));putpoly(f,enc(c.denom))
   for nd in o.s:f.write(' '.join(map(str,nd))+'\n')
   for nd in o.l:
    if nd[0]=='l':f.write('l '+str(len(nd[1]))+''.join(' '+str(c)+' '+str(k)for c,k in nd[1])+'\n')
    else:f.write(' '.join(map(str,nd))+'\n')
   for k,v in checks:f.write(f'{k} {v}\n')
   for k,v in outputs:f.write(f'{k} {v}\n')
   f.write('174 22\n')
   for row in km:
    for c in row:
     if c:putpoly(f,c[0]);putpoly(f,c[1])
     else:putpoly(f,[]);putpoly(f,[[[0,0,0],'1']])
  left=117-(time.monotonic()-began);assert left>0
  cp=subprocess.run([str(exe),str(path),str(out/(name+'.log'))],capture_output=True,text=True,timeout=left)
  (out/(name+'-process.log')).write_text(cp.stdout+cp.stderr);assert cp.returncode==0,(name,cp.returncode)
  print(json.dumps(dict(stage=name,passed=True,seconds=time.monotonic()-began,scalar_nodes=len(o.s),linear_nodes=len(o.l))),flush=True)
 runstage('new22-extractor',[[n,21]])
 (out/'new22-checkpoint.json').write_text(json.dumps(dict(new22_coefficient_extractor_passed_all22_frozen_columns=True,new_linear_output=n,lambda5_coefficients=lam,pairdual6x2=pairdual,core_anchors=dict(L_x1=str(L1),ell_x1=str(ell1),g=str(g)),seconds=time.monotonic()-began),indent=2)+'\n')
 # Construct the global6-quad dual from the frozen4-quad shadow rightinverse.
 V6,*w=ring('w0,w1,w2,w3,Z,W',K);ell0=-sum((V6(coef(ell,i))*w[k]for k,i in enumerate([0,2,3,4,5,6])),V6.zero)/ell1
 ellimages=[w[0],ell0,w[1],w[2],w[3],w[4],w[5]];wm=mons(6,2)
 Qbar=[substitute(q,ellimages,V6)for q in qp[:4]];ec=read(paths['ellchart']);mult=[frac(v)for v in ec['exact_clearing_unit_multipliers']]
 for i,row in enumerate(ec['cleared_qfour_sparse']):
  saved={tuple(mon[:6]):dec(c)for mon,c in row};assert all(Qbar[i].get(m,K.zero)*mult[i]==saved.get(m,K.zero)for m in wm)
 with gzip.open(paths['shadow'],'rt')as f:
  assert f.readline().split()==['21','4','0'];raw=[]
  for _ in range(21):
   row=[]
   for __ in range(4):
    t=f.readline().split();ct=int(t[0]);assert len(t)==1+4*ct;terms={}
    for k in range(ct):
     ex=list(map(int,t[2+4*k:5+4*k]));terms[(ex[0]+ex[2],ex[1],ex[2])]=s.QQ(t[1+4*k])
    row.append(K(R.from_dict(terms)))
   raw.append(row)
  assert not f.read().strip()
 gshadow=-K(a**4*d**2*J);shift=[0,2,2,1]
 v4=[[o.c(raw[m][i]*K(a)**shift[i]*mult[i]/gshadow)for i in range(4)]for m in range(21)]
 restriction=[]
 for mon in xm:restriction.append(substitute(T.from_dict({mon:K.one}),ellimages,V6))
 dual=[]
 for i in range(4):dual.append([o.sum([(o.c(restriction[k].get(m,K.zero)),v4[h][i])for h,m in enumerate(wm)])for k in range(28)])
 for zindex in[5,6]:
  psi=[]
  for mon in xm:
   z=K.zero
   if mon==tuple(2*int(j==1)for j in range(7)):z=-2*coef(ell,zindex)/ell1**2
   if mon==tuple(int(j==1)+int(j==zindex)for j in range(7)):z=K.one/ell1
   psi.append(z)
  pq=[sum((psi[k]*q.get(mon,K.zero)for k,mon in enumerate(xm)),K.zero)for q in qp[:4]]
  dual.append([o.sub(o.c(psi[k]),o.sum([(o.c(pq[i]),dual[i][k])for i in range(4)]))for k in range(28)])
 for j,q in enumerate(qp):
  for i in range(6):checks.append([o.sum([(o.c(q.get(mon,K.zero)),dual[i][k])for k,mon in enumerate(xm)]),int(i==j)])
 N=[T.zero,ell*(T(K(e))*L-T(K(400*a))*x0),-ell*L,T.zero,-TZ,-TW]
 current=[]
 for j in range(6):current.append(plin([(o.one,Xp[j]),(o.sub(0,o.one),{m:o.lin([(o.c(c),n)])for m,c in N[j].items()})]))
 apr={}
 for i in range(4):
  for j in range(i+1,6):apr[i,j]=o.lin([(dual[j][k],current[i].get(mon,0))for k,mon in enumerate(xm)])
 residual=[p.copy()for p in current]
 for(i,j),v in apr.items():
  residual[i]=plin([(o.one,residual[i]),(o.sub(0,o.one),{m:o.lin([(o.c(c),v)])for m,c in qp[j].items()})])
  residual[j]=plin([(o.one,residual[j]),(o.one,{m:o.lin([(o.c(c),v)])for m,c in qp[i].items()})])
 monWX1=tuple(int(i in[1,6])for i in range(7));apr[4,5]=o.lin([(o.c(K.one/ell1),residual[4].get(monWX1,0))])
 lp=[]
 for j in[0,2,3,4,5,6]:
  mon=tuple(int(i==6)+int(i==j)for i in range(7));lp.append(o.lin([(o.one,residual[4].get(mon,0)),(o.c(-coef(ell,j)),apr[4,5])]))
 orig=[];pairs=list(itertools.combinations(range(6),2))
 for i,j in pairs:orig.append(o.lin([(o.c(ch[p][i]*ch[q][j]-ch[p][j]*ch[q][i]),apr[p,q])for p,q in pairs]))
 outputs=orig+lp+[n];assert len(outputs)==22
 runstage('full-P22',[[v,i]for i,v in enumerate(outputs)])
 art=dict(circuit=o.encoded(),P22_linear_output_roots_in_frozen_original_kernel_order=outputs,source_variables=kd['actual_source174_definitions'],
  global_quad_dual6x28_scalar_roots=dual,quadric_monomials7=xm,pair_dual6x2_scalar_roots=pairdual,
  coefficient_extractor_new22=n,prime_quadric_source_change=[[str(c)for c in row]for row in ch],
  inverse_prime_cofactor_transform=[[str(c)for c in row]for row in inv],
  global282_source_formula='Use acceptedM304=[A174 B22 C108]. Actualglobal282 map=[A+B*P22,C], in the sameTn/Qhatframe. No quarticcoefficients expanded.',
  fixed282_submatrix_not_claimed=True)
 (out/'global-P22.json').write_text(json.dumps(art,separators=(',',':'))+'\n')
 inputs=[Path(__file__),source,ROOT/'code/s61_binary_hyperplane_cubic_qq3_20261102.hpp',*paths.values()]
 result=dict(id=ID,all_exact_extractor_checks_passed=True,new22_extractor_on_all22_basis_columns=True,
  actual_pair_dual2_and_global_quad_dual6_verified=True,literal_global_P22_times_frozen_K_equals_I22=True,
  global282_basis_formula_constructed=True,fixed282_submatrix_claimed=False,pointwise_quad_rank_without_dual_not_assumed=True,
  open='FULLU=D(a*d*F*B*C*J)',unknown_factor_inverted=False,no_new_search_minor_GB_field_point_or_fullquartic_expansion=True,
  no_full239_C122rank64_or_additional_kernel_directions_claim=True,scalar_nodes=len(o.s),linear_nodes=len(o.l),core_atoms=len(o.atoms),
  seconds=time.monotonic()-began,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name not in['native','report.json']},frozen=True,independent_review_pending=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 for p in out.iterdir():assert p.stat().st_size<32*1024*1024
 print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
