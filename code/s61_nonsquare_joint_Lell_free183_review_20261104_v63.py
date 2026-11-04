#!/usr/bin/env python3
"""Independent typed circuit audit of unit98 Schur and free183 transport."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_free281_circuit_review_20261104_v63 import ProofExpressions,polynomial,sha,read
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63_REVIEW'
def matrix(p):
 ts=iter(p.read_text().split());r,c=int(next(ts)),int(next(ts));out=[]
 for _ in range(r):
  row=[]
  for _ in range(c):
   terms=[]
   for _ in range(int(next(ts))):
    co=next(ts);terms.append([[int(next(ts))for j in range(3)],co])
   row.append(terms)
  out.append(row)
 assert next(ts,None)is None;return out
class LinearProof:
 def __init__(self,sc):self.sc=sc;self.ids={('zero',):0};self.data=[('zero',)]
 def intern(self,key):
  if key not in self.ids:self.ids[key]=len(self.data);self.data.append(key)
  return self.ids[key]
 def variable(self,i):return self.intern(('variable',i))
 def combine(self,terms):
  values={}
  for c,j in terms:
   if c and j:values[j]=self.sc.plus(values.get(j,0),c)
  if not values:return 0
  if len(values)==1:
   j,c=next(iter(values.items()))
   if c==self.sc.one:return j
  return self.intern(('linear',tuple((c,j)for j,c in sorted(values.items()))))
 def substitute(self,roots,values):
  cache={0:0}
  def rec(i):
   if i in cache:return cache[i]
   nd=self.data[i]
   if nd[0]=='variable':v=values[nd[1]]
   else:v=self.sc.total((c,rec(j))for c,j in nd[1])
   cache[i]=v;return v
  return[rec(i)for i in roots]
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63';r=read(b/'report.json');c=read(b/'certificate.json');assert r['frozen']and r['actual_free_source_target_shape']==[183,141]
 assert sha(b/'certificate.json')==r['certificate_sha256']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
 pre=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';pf=read(pre/'preflight.json');stage=pf['stages'][1];oldb=ROOT/'reports/S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63';fd=read(oldb/'certificate.json')
 assert c['scalar_circuit']['nodes'][:len(fd['circuit']['nodes'])]==fd['circuit']['nodes']and c['scalar_circuit']['atoms_num_den'][:len(fd['circuit']['atoms_num_den'])]==fd['circuit']['atoms_num_den']
 assert c['old_source319_definitions']==fd['source319_definitions']and c['only_quadratic174_source_basis_changed']and c['cubic112_keep_old_light_multipliers']
 sc=ProofExpressions();sg=c['scalar_circuit'];sm=[]
 for i,nd in enumerate(sg['nodes']):
  if nd[0]=='z':assert i==0;v=0
  elif nd[0]=='c':v=sc.atom(*sg['atoms_num_den'][nd[1]])
  else:assert nd[0]in['a','s','m']and all(0<=j<i for j in nd[1:]);v={'a':sc.plus,'s':sc.minus,'m':sc.times}[nd[0]](sm[nd[1]],sm[nd[2]])
  sm.append(v)
 lin=LinearProof(sc);lm=[];masks=[]
 for i,nd in enumerate(c['linear_nodes']):
  if nd[0]=='z':assert i==0;v=0;mask=0
  elif nd[0]=='v':assert 0<=nd[1]<685;v=lin.variable(nd[1]);mask=1<<next(j for j,(lo,hi)in enumerate([(0,174),(174,348),(348,424),(424,446),(446,685)])if lo<=nd[1]<hi)
  else:
   assert nd[0]=='l'and all(0<=k<len(sm)and 0<=j<i for k,j in nd[1]);v=lin.combine((sm[k],lm[j])for k,j in nd[1]);mask=0
   for k,j in nd[1]:mask|=masks[j]
  lm.append(v);masks.append(mask)
 def checkS(key,expected):assert[sm[i]for i in c[key]]==expected,key
 def checkL(record,key,expected,group):
  assert[lm[i]for i in c[record][key]]==expected,record
  assert all(masks[i]in[0,1<<group]for i in c[record][key]),record
 R,a,e,d=ring('a,e,d',s.QQ);field=R.to_field();enc=lambda p:[[list(m),str(v)]for m,v in p.items()]
 def atom(v):
  v=field(v);return sc.atom(enc(v.numer),enc(v.denom))
 C=100*a+18*e**2-9*e-2;delta=field(3200)/81*field(a*C*d**2)
 Nraw=matrix(pre/'coordinate-inverse-numerators7x7.txt');Fraw=matrix(pre/'coordinate-forward7x7.txt');Qraw=matrix(pre/'six-prime-quadrics-numerators.txt');D=matrix(pre/'six-prime-quadrics-denominators.txt')
 unit=[[[0,0,0],'1']];N=[[sc.times(sc.atom(v,unit),atom(1/delta))for v in row]for row in Nraw];F=[[sc.atom(v,unit)for v in row]for row in Fraw];quad=[[sc.atom(Qraw[i][j],D[i][0])for j in range(29)]for i in range(6)]
 assert[[sm[j]for j in row]for row in c['quadric_generator_coefficient_roots6x29']]==quad
 # Coordinate matrices in physicalv8 and newY8; H changes its index5->7 only.
 B=[[0]*8 for _ in range(8)];A=[[0]*8 for _ in range(8)]
 for v,x in[(0,0),(1,1),(2,2),(3,6),(4,3),(7,4)]:B[v][:7]=N[x]
 for j in range(7):B[6][j]=sc.minus(sc.plus(N[5][j],N[1][j]),sc.times(atom(3),N[0][j]))
 B[5][7]=sc.one
 for y in range(7):
  A[y][0]=sc.plus(F[y][0],sc.times(atom(3),F[y][5]));A[y][1]=sc.minus(F[y][1],F[y][5]);A[y][2]=F[y][2];A[y][3]=F[y][6];A[y][4]=F[y][3];A[y][6]=F[y][5];A[y][7]=F[y][4]
 A[7][5]=sc.one
 oldm2=[d[2]for d in fd['source319_definitions'][:29]];newm2=[m for _,m in stage['four_quadric_product_rows'][:29]];assert c['new_prime_quadric_source174_definitions']==[[i,m]for i in range(6)for m in newm2]
 def sym2(m,images,index):
  # Degree2 cofactor is either two light linear forms or a single H.
  variables=[j for j,n in enumerate(m)for _ in range(n)];poly={(0,)*8:sc.one}
  for v in variables:
   tmp={}
   for mon,co in poly.items():
    for j,im in enumerate(images[v]):
     if im:
      dest=list(mon);dest[j]+=1;dest=tuple(dest);tmp[dest]=sc.plus(tmp.get(dest,0),sc.times(co,im))
   poly=tmp
  return{index[m]:v for m,v in poly.items()}
 TB=[sym2(m,B,{tuple(v):j for j,v in enumerate(newm2)})for m in oldm2];TA=[sym2(m,A,{tuple(v):j for j,v in enumerate(oldm2)})for m in newm2]
 for key,expected in[('monomial_old_to_new_scalar_roots',TB),('monomial_new_to_old_scalar_roots',TA)]:
  assert[{j:sm[v]for j,v in row}for row in c[key]]==expected,key
 change=[[field.from_expr(s.sympify(v))for v in row]for row in c['six_quadric_change_K6']];inverse=[[field.from_expr(s.sympify(v))for v in row]for row in c['inverse_six_quadric_change']]
 expected=[[0,8100,0,0,0,0],[0,field(81*e)/field(4*a),0,-field(81)/field(4*a),0,0],[0,field(10*a+3*e**2-3*e)/field(a),0,-field(3*(e-1))/field(a),-10,0],[0,-d*e-8,8,d,0,0],[1,0,0,0,0,0],[0,0,0,0,0,1]]
 assert change==[[field(v)for v in row]for row in expected]
 for left,right in[(change,inverse),(inverse,change)]:assert all(sum((left[i][k]*right[k][j]for k in range(6)),field.zero)==int(i==j)for i in range(6)for j in range(6))
 cm=[[atom(v)for v in row]for row in change];im=[[atom(v)for v in row]for row in inverse]
 def transform(values,forward):
  monmap,coeff=(TB,im)if forward else(TA,cm);tmp=[]
  for qi in range(6):
   lists=[[]for _ in range(29)]
   for j,row in enumerate(monmap):
    for k,co in row.items():lists[k].append((co,values[qi*29+j]))
   tmp.append([lin.combine(v)for v in lists])
  return[lin.combine((coeff[i][j],tmp[i][k])for i in range(6)if coeff[i][j])for j in range(6)for k in range(29)]
 oldrange=list(range(174));primerange=list(range(174,348));erange=list(range(348,424));krange=list(range(424,446));trange=list(range(446,685))
 for key,var in[('old_to_prime_quad_cofactor_operator',oldrange),('prime_to_old_quad_cofactor_operator',primerange),('quad_restoration_operator',erange),('descended_P22_operator',erange),('descended_K22_product_operator',krange),('target_projection_operator',trange)]:assert c[key]['input_variables']==var
 old2prime=transform([lin.variable(i)for i in oldrange],True);prime2old=transform([lin.variable(i)for i in primerange],False)
 checkL('old_to_prime_quad_cofactor_operator','output174_linear_roots',old2prime,0);checkL('prime_to_old_quad_cofactor_operator','output174_linear_roots',prime2old,1)
 mons=stage['target_monomials'];mi={tuple(m):j for j,m in enumerate(mons)};qrows=[]
 for qi in range(6):
  for cof in newm2:qrows.append({mi[tuple(x+y for x,y in zip(cof,m))]:co for m,co in zip(newm2,quad[qi])if co})
 assert[{j:sm[v]for j,v in row}for row in c['quadric_product_coefficient_sparse_rows174']]==qrows
 selected=stage['monic_pivot_source_rows'];targets=stage['monic_pivot_target_columns'];targetD=stage['retained_target_columns'];rest=[i for i in range(174)if i not in selected]
 assert c['unit98_source_rows']==selected and c['unit98_target_columns']==targets and c['retained76_quadratic_source_rows']==rest and c['retained141_target_columns']==targetD
 assert len(rest)==76 and len(targetD)==141 and c['all18_overlap_source_rows_retained']==stage['all_unselected_overlap_rows_retained']and set(stage['all_unselected_overlap_rows_retained'])<=set(rest)
 U=[[qrows[i].get(j,0)for j in targets]for i in selected];RD=[[qrows[i].get(j,0)for j in targetD]for i in selected];NC=[[qrows[i].get(j,0)for j in targets]for i in rest]
 for key,expected in[('U98_scalar_roots',U),('RD98x141_scalar_roots',RD),('NC76x98_scalar_roots',NC)]:assert[[sm[j]for j in row]for row in c[key]]==expected,key
 assert all(U[i][i]==sc.one for i in range(98))and all(U[i][j]==0 for i in range(98)for j in range(i))
 xlin=[lin.variable(i)for i in erange];beta=[];bcol=[lin.combine((NC[e][j],xlin[e])for e in range(76)if NC[e][j])for j in range(98)]
 for j in range(98):beta.append(lin.combine([(atom(-1),bcol[j])]+[(sc.minus(0,U[i][j]),beta[i])for i in range(j)if U[i][j]]))
 restored=[0]*174
 for i,row in enumerate(selected):restored[row]=beta[i]
 for i,row in enumerate(rest):restored[row]=xlin[i]
 restoredold=transform(restored,False);checkL('quad_restoration_operator','prime174_linear_roots',restored,2);checkL('quad_restoration_operator','old174_linear_roots',restoredold,2)
 pp=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63/global-P22.json';pd=read(pp);pc=pd['circuit'];psc={0:0};plc={0:0}
 def ps(i):
  if i in psc:return psc[i]
  nd=pc['scalar_nodes'][i]
  if nd[0]=='c':v=sc.atom(*pc['atoms_num_den'][nd[1]])
  else:v={'a':sc.plus,'s':sc.minus,'m':sc.times}[nd[0]](ps(nd[1]),ps(nd[2]))
  psc[i]=v;return v
 def pl(i):
  if i in plc:return plc[i]
  nd=pc['linear_nodes'][i]
  if nd[0]=='v':v=restoredold[nd[1]]
  else:v=lin.combine((ps(k),pl(j))for k,j in nd[1])
  plc[i]=v;return v
 Prest=[pl(i)for i in pd['P22_linear_output_roots_in_frozen_original_kernel_order']];checkL('descended_P22_operator','output22_linear_roots',Prest,2)
 kp=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json';kd=read(kp);rows=[[]for _ in range(174)]
 for j,col in enumerate(kd['sparse_columns']):
  for i,num,den in col:rows[i].append((sc.atom(num,den),lin.variable(krange[j])))
 Kprime=transform([lin.combine(row)for row in rows],True);Krest=[Kprime[i]for i in rest];checkL('descended_K22_product_operator','output76_linear_roots',Krest,3)
 def sameScalar(key,values):assert[sm[j]for j in c[key]]==values,key
 oldx=[sm[j]for j in fd['new_quadratic174_roots']];xprime=lin.substitute(old2prime,dict(zip(oldrange,oldx)));x76=[xprime[i]for i in rest];sameScalar('new_quad76_kernel_after_old15_roots',x76)
 u=[sm[j]for j in fd['new_B22_roots']];h=[sm[j]for j in fd['new_C108_roots']];sameScalar('new_B22_kernel_after_old15_roots',u);sameScalar('new_C108_kernel_after_old15_roots',h)
 Px=lin.substitute(Prest,dict(zip(erange,x76)));sameScalar('Pnew_on_actualx76_roots',Px);k=[sc.minus(Px[j],u[j])for j in range(22)];sameScalar('kernel22_coordinate_roots',k)
 Kk=lin.substitute(Krest,dict(zip(krange,k)));sameScalar('Knew_times_actualk22_roots',Kk);kernel=[sc.minus(x76[i],Kk[i])for i in range(76)]+h;sameScalar('actual184_kernel_roots',kernel)
 assert c['Bnew22_original319_labels']==fd['source304_new_block_B22']and c['Cnew108_original319_labels']==fd['source304_new_block_C108']
 p,q=[76+fd['source304_new_block_C108'].index(i)for i in[182,186]];assert[p,q]==c['separator184_indices']==[76,77]
 theta=sm[c['separator_known_unit_root']];assert c['separator_known_unit_root']==fd['known_unit_separator_root']and c['separator_zero_check_root']==fd['separator_identity_zero_root']
 assert sc.minus(sc.plus(kernel[p],kernel[q]),theta)==sm[c['separator_zero_check_root']]
 w=kernel[:];w[p]=theta;sameScalar('changed184_kernel_roots',w);inverseTheta=atom(-field(243)/(field(1024000)*field(a*C**3*(5*a*d-6*e-4))));normalized=[sc.times(v,inverseTheta)for v in w];normalized[p]=sc.one;sameScalar('normalized_unit_kernel_roots',normalized)
 assert c['actual_free183_shape']==[183,141]and c['source183_changed184_indices']==[i for i in range(184)if i!=p]
 mult=[]
 for j in range(98):mult.append(lin.combine([(sc.one,lin.variable(trange[targets[j]]))]+[(sc.minus(0,U[i][j]),mult[i])for i in range(j)if U[i][j]]))
 NF=[lin.combine([(sc.one,lin.variable(trange[d]))]+[(sc.minus(0,RD[i][j]),mult[i])for i in range(98)if RD[i][j]])for j,d in enumerate(targetD)];checkL('target_projection_operator','output141_linear_roots',NF,4)
 # Denominator and scalar-alias checks; all large roots stay factored.
 with(out/'scalar-proof-input.txt').open('w')as f:
  f.write(f'{len(sg["atoms_num_den"])} {len(sg["nodes"])} {c["separator_zero_check_root"]} {c["separator_known_unit_root"]}\n')
  for num,den in sg['atoms_num_den']:
   for terms in[num,den]:
    poly=polynomial(terms);f.write(str(len(poly))+''.join(' '+str(co)+' '+' '.join(map(str,m))for m,co in poly)+'\n')
  for nd in sg['nodes']:
   vals=[0,0,0]if nd[0]=='z'else[1,nd[1],0]if nd[0]=='c'else[{'a':2,'s':3,'m':4}[nd[0]],nd[1],nd[2]];f.write(' '.join(map(str,vals))+'\n')
 cpp=ROOT/'code/s61_nonsquare_joint_Lell_free183_review_20261104_v63.cpp';exe=out/'review-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();rem=lambda:max(.1,60-(time.monotonic()-began))
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=rem())
 run=subprocess.run([str(exe),str(out/'scalar-proof-input.txt')],capture_output=True,text=True,timeout=rem());exe.unlink();(out/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0,run.stderr;native=json.loads(run.stdout);assert native['all_passed']
 deps=[ROOT/'reports'/j/'report.json'for j in['S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63_REVIEW','S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63_REVIEW','S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW']];assert all(read(p)['all_passed']for p in deps)
 note=ROOT/'research/notes'/(ID+'.md');utility=ROOT/'code/s61_nonsquare_free281_circuit_review_20261104_v63.py';cpputility=ROOT/'code/s61_nonsquare_ell_Lframe_review_20261104_v63.cpp'
 result=dict(id=ID,all_passed=True,actual174_cofactor_basis_contragredients_both_directions_bound=True,
  all174_quad_products_U98_RD_NC_and_18overlap_indices_checked=True,
  all2705_linear_nodes_typed_topological_and_operator_variable_groups_checked=True,
  exact_beta_minus_U_inverse_transpose_NC_transpose_recursion_checked=True,
  Pnew_Knew_identity_accepted_by_pure_quad_kernel_restoration_and_pinned_old_PK=True,
  all_actual184_kernel_roots_and_same_reviewed_unit_separator_checked=True,
  actual_full319_cokernel_equals_free183_to141_on_FULLU=True,
  genuine112_OLD_light_multipliers_and_all33_quartics_preserved=True,
  no_rank141_full239_expanded183matrix_newpoints_GB_solver_or_author_rerun=True,
  target141_full_cokernel_analysis_open=True,native=native,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),cpp,utility,cpputility,note,b/'report.json',b/'certificate.json',pre/'report.json',pre/'preflight.json',oldb/'certificate.json',pp,kp,*deps]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-began,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
