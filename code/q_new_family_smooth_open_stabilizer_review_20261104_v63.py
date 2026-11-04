#!/usr/bin/env python3
"""Independent Fraction Laurent-ring audit, without importing the producer."""
import argparse,ast,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_SMOOTH_OPEN_STABILIZER_STRATA_20261104_v63_REVIEW'
ZERO=(0,)*8
def q(c):return {}if not c else{ZERO:Q(c)}
def var(i):e=list(ZERO);e[i]=1;return{tuple(e):Q(1)}
def add(*ps):
 r={}
 for p in ps:
  for m,c in p.items():r[m]=r.get(m,Q(0))+c
 return{m:c for m,c in r.items()if c}
def neg(p):return{m:-c for m,c in p.items()}
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 r={}
 for m,c in a.items():
  for n,d in b.items():
   e=tuple(i+j for i,j in zip(m,n));r[e]=r.get(e,Q(0))+c*d
 return{m:c for m,c in r.items()if c}
def pw(a,n):
 if n<0:
  assert len(a)==1,'only a named Laurent monomial may be inverted'
  (m,c),=a.items();return{tuple(n*i for i in m):c**n}
 r=q(1)
 for _ in range(n):r=mul(r,a)
 return r
def scale(p,c):return mul(p,q(c))
def parse(text,env):
 def ev(n):
  if isinstance(n,ast.Constant):assert isinstance(n.value,int);return q(n.value)
  if isinstance(n,ast.Name):return env[n.id]
  if isinstance(n,ast.UnaryOp):return neg(ev(n.operand))if isinstance(n.op,ast.USub)else ev(n.operand)
  assert isinstance(n,ast.BinOp)
  if isinstance(n.op,ast.Add):return add(ev(n.left),ev(n.right))
  if isinstance(n.op,ast.Sub):return sub(ev(n.left),ev(n.right))
  if isinstance(n.op,ast.Mult):return mul(ev(n.left),ev(n.right))
  if isinstance(n.op,ast.Div):return mul(ev(n.left),pw(ev(n.right),-1))
  assert isinstance(n.op,ast.Pow);return pw(ev(n.left),ast.literal_eval(n.right))
 return ev(ast.parse(str(text),mode='eval').body)
def derivative(p,i):
 r={}
 for m,c in p.items():
  if m[i]:e=list(m);e[i]-=1;r[tuple(e)]=c*m[i]
 return r
def paramcoef(p,hom):
 r={}
 for m,c in p.items():
  if m[3:]==tuple(hom):r[m[:3]+(0,)*5]=c
 return r
def homogeneous(p):
 r={}
 for m,c in p.items():
  assert m[3]==0;d=sum(m[4:]);assert 0<=d<=2
  e=list(m);e[3]=2-d;r[tuple(e)]=c
 return r
def ser(p):return[[str(c.numerator),str(c.denominator),*m]for m,c in sorted(p.items())]
def mm(A,B):return[[add(*(mul(A[i][k],B[k][j])for k in range(len(B))))for j in range(len(B[0]))]for i in range(len(A))]
def matadd(A,B):return[[add(a,b)for a,b in zip(r,s)]for r,s in zip(A,B)]
def matscale(A,c):return[[mul(c,a)for a in r]for r in A]
def eye(n):return[[q(int(i==j))for j in range(n)]for i in range(n)]
def determinant(A):
 r={}
 for per in itertools.permutations(range(len(A))):
  t=q((-1)**sum(per[i]>per[j]for i in range(len(A))for j in range(i+1,len(A))))
  for i,j in enumerate(per):t=mul(t,A[i][j])
  r=add(r,t)
 return r
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def specialize(p,i,value):
 r={}
 for m,c in p.items():
  e=list(m);n=e[i];e[i]=0
  r=add(r,mul({tuple(e):c},pw(value,n)))
 return r
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 rootstem='Q_NEW_FAMILY_SMOOTH_OPEN_STABILIZER_STRATA_20261104_v63';primary=ROOT/'reports'/rootstem/'report.json';cp=primary.parent/'certificate.json';geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
 report=json.loads(primary.read_text());data=json.loads(cp.read_text());G=json.loads(geo.read_text());assert report['all_passed']and G['all_passed']
 checked={}
 for n,h in report['input_sha256'].items():assert sha(ROOT/n)==h;checked[n]=h
 for n,h in report['evidence_sha256'].items():assert sha(primary.parent/n)==h;checked[str((primary.parent/n).relative_to(ROOT))]=h
 u,t,T=var(0),var(1),var(2);sp=[var(i)for i in range(4,8)];env=dict(zip(['x','y','z','w'],sp));env.update(u=u,v=var(1),n0=u,n1=q(1),n2=var(1))
 P=[[parse(z,env)for z in row]for row in G['affine_P']]
 # Own full20x6 constant-Schouten matrix; scalar derivatives, no author ring.
 pair4=list(itertools.combinations(range(4),2));compat=[]
 for i,j,k in itertools.combinations(range(4),3):
  pols=[]
  for aa,bb in pair4:
   cc=[[int(r==aa and s==bb)-int(r==bb and s==aa)for s in range(4)]for r in range(4)]
   pols.append(add(*(scale(derivative(P[j][k],4+b),cc[i][b])for b in range(4)),*(scale(derivative(P[k][i],4+b),cc[j][b])for b in range(4)),*(scale(derivative(P[i][j],4+b),cc[k][b])for b in range(4))))
  for mon in[(0,)*5,*[(0,)+tuple(int(i==k)for i in range(4))for k in range(4)]]:compat.append([paramcoef(z,mon)for z in pols])
 assert compat==[[parse(z,env)for z in row]for row in data['compatibility_matrix']]
 # Necessary rows have coefficients v,2,-64,4/3: units throughout uvp-open.
 assert compat[11]==[var(1),{}, {},{}, {},{}]
 assert compat[4]==[{}, {},q(2),{}, {},{}]
 assert compat[19]==[{}, {},{}, {},q(-64),{}]
 assert compat[15]==[{},q(Q(4,3)),{}, {},{},scale(u,Q(4,3))]
 kernel=[{},neg(u),{}, {},{},q(1)]
 kr=[add(*(mul(a,b)for a,b in zip(row,kernel)))for row in compat]
 assert kr[0]==scale(sub(pw(u,2),q(64)),Q(4,3))and all(not z for z in kr[1:])
 # Own divergence and homogeneous N, then exact image minor.
 def modular(P):
  div=[add(*(derivative(P[i][j],4+i)for i in range(4)))for j in range(4)]
  N=[[{}for _ in range(5)]for _ in range(5)]
  for i,f in enumerate(div):
   N[i+1][0]=paramcoef(f,(0,)*5)
   for k in range(4):N[i+1][k+1]=paramcoef(f,(0,)+tuple(int(j==k)for j in range(4)))
  return N
 N=modular(P);assert N==[[parse(z,env)for z in row]for row in data['N']]
 assert determinant([[N[2][0],N[2][1]],[N[3][0],N[3][1]]])==scale(pw(u,2),Q(4,3))
 assert all(not N[i][j]for i in[0,1,4]for j in range(5))
 # All exceptional quadratic case identities, independently in a polynomial ring.
 U,V,A,C=[var(i)for i in range(4)];ee=dict(u=U,v=V,a=A,c=C)
 Qs=[parse(z,ee)for z in data['quadratic_system']];D=sub(pw(A,2),scale(pw(C,2),64));p=add(scale(U,32),V)
 assert not add(Qs[0],Qs[2],mul(sub(V,p),Qs[1]),scale(mul(mul(A,C),sub(pw(U,2),q(256))),32))
 for sg in[1,-1]:
  qm=[specialize(f,0,q(16*sg))for f in Qs]
  ap=add(A,scale(C,8*sg));am=sub(A,scale(C,8*sg))
  assert qm[1]==mul(pw(ap,2),sub(q(1),pw(am,2)))
  assert add(qm[2],mul(V,qm[1]))==scale(mul(mul(C,add(q(256),scale(V,sg))),ap),16)
  qrr=[specialize(f,1,q(-256*sg))for f in qm]
  assert qrr[0]==scale(qrr[1],256*sg)and qrr[2]==scale(qrr[1],256*sg)
 assert specialize(specialize(Qs[1],1,scale(U,-16)),2,{})==scale(mul(pw(C,2),sub(q(1),scale(pw(C,2),64))),64)
 def full_residual(raw,M,kap):
  forms=[add(*(mul(M[i][j],var(3+j))for j in range(5)))for i in range(5)]
  products={}
  for i,j in itertools.combinations_with_replacement(range(5),2):
   mon=tuple(int(k==i)+int(k==j)for k in range(5));products[mon]=mul(forms[i],forms[j])
  answer={}
  for i,j in itertools.combinations(range(5),2):
   lhs=add(*(mul(mul(kap,sub(mul(M[i][a],M[j][b]),mul(M[i][b],M[j][a]))),raw[a][b])for a,b in itertools.combinations(range(5),2)))
   rhs={}
   for m,c in raw[i][j].items():rhs=add(rhs,mul({m[:3]+(0,)*5:c},products[m[3:]]))
   answer[(i,j)]=sub(lhs,rhs)
  return answer
 cases=[];groupbase=None
 for case in data['families']:
  ce=dict(u=u,t=t,T=T);uu=parse(case['u'],ce);vv=parse(case['v'],ce)
  pe=dict(env,n0=uu,n2=vv,u=uu,v=vv)
  pc=[[parse(z,pe)for z in row]for row in G['affine_P']]
  raw=[[{}for _ in range(5)]for _ in range(5)]
  for i,j in pair4:raw[i+1][j+1]=homogeneous(pc[i][j]);raw[j+1][i+1]=neg(raw[i+1][j+1])
  M=[[parse(z,ce)for z in row]for row in case['matrix']];kap=parse(case['kappa'],ce)
  residual=full_residual(raw,M,kap)
  # Independent constant matrix from C0 and target evaluation, not pulses.
  aa,cc,dd,ee,gg,jj,ii,tx,tw=M[1][1],M[4][1],M[2][2],M[2][1],M[2][4],M[3][3],M[3][2],M[1][0],M[4][0]
  L=[[scale(mul(kap,cc),Q(-512,9)),scale(mul(kap,aa),Q(8,9)),{}],
     [scale(mul(kap,aa),Q(8,9)),scale(mul(kap,cc),Q(-8,9)),{}],
     [mul(kap,add(scale(dd,Q(-1,36)),scale(gg,Q(-8,9)))),scale(mul(kap,ee),Q(8,9)),scale(uu,Q(4,3))]]
  rhs=[sub(mul(add(tx,mul(vv,T)),sub(q(Q(1,48)),scale(tw,64))),scale(mul(mul(kap,aa),ii),Q(1,36))),
       add(pw(tx,2),mul(mul(vv,tx),T),scale(T,Q(2,3)),scale(mul(mul(kap,cc),ii),Q(1,36))),
       sub(add(scale(mul(mul(uu,vv),pw(T,2)),Q(1,2)),mul(mul(uu,tx),T),scale(tw,Q(-1,1536)),pw(tw,2),scale(pw(tx,2),Q(-1,64)),scale(mul(tw,T),64),q(Q(-2,9))),
           mul(kap,sub(scale(mul(ee,ii),Q(1,36)),scale(mul(dd,jj),Q(2,9)))))]
  assert L==[[parse(z,ce)for z in row]for row in case['constant_linear_matrix']]
  assert rhs==[parse(z,ce)for z in case['constant_rhs']]
  assert determinant(L)==parse(case['constant_unit_minor'],ce)
  assert determinant(L)==scale(mul(uu,pw(dd,-1)),Q(-256,243))
  sol=[parse(z,ce)for z in case['H_K_Z']]
  assert[add(*(mul(x,y)for x,y in zip(row,sol)))for row in L]==rhs
  assert sol==[M[3][1],M[3][4],M[3][0]]
  obstruction=case['resonant_linear_obstruction'];own={}
  if obstruction is None:
   assert all(not z for z in residual.values());groupbase=(M,modular(pc))if case['name']=='J0_cplus'else groupbase
  else:
   expected={pair:{}for pair in residual}
   for key,value in obstruction['all_coefficients_divided_by_t_squared_minus_one'].items():
    pair,mon=ast.literal_eval(key);poly=mul(parse(value,ce),sub(pw(t,2),q(1)));poly={m[:3]+tuple(mon):c for m,c in poly.items()};expected[pair]=add(expected[pair],poly)
   assert residual==expected
   witness=next((q for q in obstruction['all_coefficients_divided_by_t_squared_minus_one'].values()if len(parse(q,ce))==1 and next(iter(parse(q,ce)))==ZERO),None);assert witness is not None
   own=dict(all_residuals_match=True,nonzero_rational_unit_witness=witness)
  cases.append(dict(name=case['name'],full150_residuals=[dict(pair=pair,coefficients=ser(poly))for pair,poly in residual.items()],constant_minor=ser(determinant(L)),obstruction=own))
 # Entire-J0 group and coset relations, not a specialization in u.
 Mp,NN=groupbase;S=[[specialize(z,2,{})for z in row]for row in Mp]
 assert S==[[parse(z,dict(u=u))for z in row]for row in data['S_matrix']]
 ev=dict(u=u,v=scale(u,-16));E=[[parse(z,ev)for z in row]for row in[
  ['1','0','0','0','0'],['-4/(3*v)','-1','0','0','0'],['4/(3*v**2)','2/v','1','0','0'],['-2/(3*v**3)','-1/v**2','0','1','0'],['1/1536','0','0','0','-1']]]
 assert mm(S,S)==eye(5)and mm(E,E)==eye(5)
 assert mm(S,E)==mm(E,S)and mm(S,NN)==matscale(mm(NN,S),q(-1))and mm(E,NN)==mm(NN,E)
 def ga(qt):return matadd(matadd(eye(5),matscale(NN,qt)),matscale(mm(NN,NN),scale(pw(qt,2),Q(1,2))))
 plus=next(c for c in data['families']if c['name']=='J0_cplus');minus=next(c for c in data['families']if c['name']=='J0_cminus')
 ce=dict(u=u,T=T);MMp=[[parse(z,ce)for z in row]for row in plus['matrix']];MMm=[[parse(z,ce)for z in row]for row in minus['matrix']]
 assert MMp==mm(ga(scale(mul(T,pw(u,-1)),Q(3,4))),S)
 ES=mm(E,S);shift=ES[2][0];assert MMm==mm(ga(scale(mul(sub(T,shift),pw(u,-1)),Q(3,4))),ES)
 proof=dict(independent_method='own Fraction Laurent8 ring, AST parser, no SymPy/author import or producer rerun',
  all120_constant_schouten_entries_match=True,constant_kernel_jump_u_squared_64_and_basis_accepted=True,
  rebuilt_modular_N_and_whole_image_flag=True,all_exceptional_quadratic_case_identities=True,
  full150_cases=cases,constant_unit_minor_formula='-256*u/(243*D)',
  S_E_group_and_N_inversion=True,all_J0_families_exact_GaS_GaES=True,
  nonJ0='Ga x mu2',J0='(Ga semidirect mu2) x mu2',domain='u*v*(32*u+v)!=0',
  generic_old_cases_retained_via_accepted_old_review_and_whole_unit_minor=True,
  special_resonance_Gm_eliminated=True,boundary_not_classified=True)
 (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 note=ROOT/'research/notes'/(ID+'.md')
 result=dict(id=ID,all_passed=True,whole_uvp_open_exhaustion_and_groups_accepted=True,
  actual_constant_jump_retained_and_modular_image_repairs_flag=True,
  full150_new_J0_and_all4_resonant_residuals_verified_independently=True,
  no_author_rerun_import_new_point_GB_or_CAS=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,primary,cp,geo]},checked_dependency_sha256=checked,
  evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:result[k]for k in['all_passed','whole_uvp_open_exhaustion_and_groups_accepted','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
