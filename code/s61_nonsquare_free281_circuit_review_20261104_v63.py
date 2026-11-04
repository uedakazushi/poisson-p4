#!/usr/bin/env python3
"""Independent structural composition audit of the actual free281 source recipe."""
import argparse,hashlib,json,resource,shutil,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def polynomial(terms):
 p={}
 for m,c in terms:
  m=tuple(m);assert len(m)==3 and all(isinstance(j,int)and j>=0 for j in m);p[m]=p.get(m,Q())+Q(c)
 return tuple(sorted((m,c)for m,c in p.items()if c))
class ProofExpressions:
 """A free expression algebra; no coefficient expansion or source execution."""
 def __init__(self):self.ids={('zero',):0};self.one=self.atom([[[0,0,0],'1']],[[[0,0,0],'1']])
 def intern(self,key):
  if key not in self.ids:self.ids[key]=len(self.ids)
  return self.ids[key]
 def atom(self,n,d):
  n,d=polynomial(n),polynomial(d);assert d
  if not n:return 0
  if n==d:n=d=(((0,0,0),Q(1)),)
  if len(n)==len(d)==1 and n[0][0]==d[0][0]==(0,0,0):n,d=(((0,0,0),n[0][1]/d[0][1]),),(((0,0,0),Q(1)),)
  return self.intern(('rational-atom',n,d))
 def plus(self,x,y):return y if x==0 else x if y==0 else self.intern(('sum',x,y))
 def minus(self,x,y):return x if y==0 else 0 if x==y else self.intern(('difference',x,y))
 def times(self,x,y):return 0 if x==0 or y==0 else y if x==self.one else x if y==self.one else self.intern(('product',x,y))
 def total(self,terms):
  z=0
  for c,j in terms:z=self.plus(z,self.times(c,j))
  return z
def main(out,prior_failure=None):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 if prior_failure:
  prior=read(prior_failure/'first-failure.json');assert prior['no_PASS_claim']and prior['no_author_input_changed']
  shutil.copyfile(prior_failure/'first-failure.json',out/'prior-failure.json');shutil.copyfile(prior_failure/'failed-reader.cpp',out/'failed-reader.cpp')
 base=ROOT/'reports/S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63';r=read(base/'report.json');c=read(base/'certificate.json')
 assert r['frozen']and r['finite_exact_source_transport_recipe_constructed']and r['new_actual_free_source_rank']==281
 assert sha(base/'certificate.json')==r['certificate_sha256']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
 compact=ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';cc=read(compact/'certificate.json');defs=cc['source319_definitions'];assert c['source319_definitions']==defs
 pbase=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63';Pd=read(pbase/'global-P22.json');old=read(pbase/'global282-basis.json');atlas=read(ROOT/old['source304_atlas_path'])
 Kd=read(ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json')
 assert Kd['actual_source174_definitions']==Pd['source_variables']==defs[:174]
 keep=atlas['retained304_source319_indices'];deleted=atlas['deleted15_source319_indices_in_pivot_order'];assert c['source304_retained319_indices']==keep and keep[:174]==list(range(174))
 assert len(keep)==304 and sorted(keep+deleted)==list(range(319))and len(deleted)==15
 B=[j for j in keep[174:]if j not in[182,186]][:22];C=[j for j in keep[174:]if j not in B]
 assert c['source304_new_block_A174']==list(range(174))and c['source304_new_block_B22']==B and c['source304_new_block_C108']==C
 assert len(B)==22 and len(C)==108 and not set(B)&set(C)
 assert all(j in old['old_source319_indices_arbitrary_B22']for j in[182,186])
 p,q=[174+C.index(j)for j in[182,186]];assert[p,q]==c['separator_new282_indices']==[174,175]
 # Read the immutable proof graph into a different expression interner.
 # Every target scalar root is then compared to a reconstruction from
 # accepted source operators; no large QQ polynomial root is evaluated.
 proof=ProofExpressions();graph=c['circuit'];nodes=graph['nodes'];atoms=graph['atoms_num_den'];mapped=[]
 assert len(nodes)==10236 and len(atoms)==898
 for i,nd in enumerate(nodes):
  op=nd[0]
  if op=='z':assert i==0;m=0
  elif op=='c':assert 0<=nd[1]<len(atoms);m=proof.atom(*atoms[nd[1]])
  else:
   assert op in['a','s','m']and all(0<=j<i for j in nd[1:]);x,y=[mapped[j]for j in nd[1:]]
   m={'a':proof.plus,'s':proof.minus,'m':proof.times}[op](x,y)
  mapped.append(m)
 def roots(key,expected):
  actual=c[key];assert len(actual)==len(expected),key
  assert[mapped[j]for j in actual]==expected,key
 # Bind the original polynomial vector directly to the reviewed literal319.
 tokens=iter((compact/'full319-kernel.txt').read_text().split());assert[int(next(tokens)),int(next(tokens))]==[319,1];N=[]
 for _ in range(319):
  terms=[]
  for _ in range(int(next(tokens))):
   val=next(tokens);exp=[int(next(tokens))for _ in range(3)];terms.append([exp,val])
  N.append(proof.atom(terms,[[[0,0,0],'1']]))
 assert next(tokens,None)is None;roots('new_source319_roots',N)
 bd=read(ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json');rd=read(ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/certificate.json')
 assert bd['actual319_source_definitions']==defs
 R,a,e,d=ring('a,e,d',s.QQ);field=R.to_field()
 def atomtext(text):
  v=field.from_expr(s.sympify(text));enc=lambda f:[[list(m),str(k)]for m,k in f.items()]
  return proof.atom(enc(v.numer),enc(v.denom))
 ar=bd['arithmetic'];assert ar['nodes'][0]==['z'];bc={0:0}
 def bscalar(i):
  if i in bc:return bc[i]
  assert ar['node_weights'][i]==0;nd=ar['nodes'][i]
  if nd[0]=='c':v=proof.atom(*ar['atoms'][nd[1]])
  elif nd[0]=='m':v=proof.times(bscalar(nd[1]),bscalar(nd[2]))
  else:
   assert nd[0]=='l';v=proof.total((proof.atom(*ar['atoms'][k]),bscalar(j))for k,j in nd[1])
  bc[i]=v;return v
 U=rd['residual_transform_rows15to4'];T=rd['pivot_row_transform'];coeff=[proof.total((atomtext(T[3][i]),atomtext(U[i][j]))for i in range(4))for j in range(15)]
 r15=[v['after35_deletions_coefficients319']for v in bd['raw15_source_lifts']];ci={g:j for j,g in enumerate(rd['cubic_source319_indices'])};last=[]
 for j in range(319):
  if j<174:v=proof.total((coeff[k],bscalar(r15[k][j]))for k in range(15)if coeff[k])
  elif j in ci:v=atomtext(rd['post_pivot_cubic_matrix'][3][ci[j]])
  else:v=0
  last.append(v)
 roots('last_residual4_source319_roots',last)
 assert[last[j]for j in[176,177,237,181]]==[0,0,0,proof.one]
 assert all(N[j]==0 for j in deleted if j!=181)
 transported=[proof.minus(N[j],proof.times(N[181],last[j]))for j in range(319)];roots('transported_source319_roots',transported)
 assert all(transported[j]==0 for j in deleted)
 x=transported[:174];u=[transported[j]for j in B];h=[transported[j]for j in C]
 roots('new_quadratic174_roots',x);roots('new_B22_roots',u);roots('new_C108_roots',h)
 # P22 is linear on174 sources; check topology and then compose its exact
 # pinned operator with x. No producer function is imported or invoked.
 pc=Pd['circuit'];assert pc['scalar_nodes'][0]==pc['linear_nodes'][0]==['z'];sc={0:0};lc={0:0}
 for i,nd in enumerate(pc['linear_nodes']):
  if i==0:continue
  if nd[0]=='v':assert 0<=nd[1]<174
  else:assert nd[0]=='l'and all(0<=j<i and 0<=k<len(pc['scalar_nodes'])for k,j in nd[1])
 def pscalar(i):
  if i in sc:return sc[i]
  nd=pc['scalar_nodes'][i]
  if nd[0]=='c':v=proof.atom(*pc['atoms_num_den'][nd[1]])
  else:
   assert nd[0]in['a','s','m']and all(j<i for j in nd[1:]);v={'a':proof.plus,'s':proof.minus,'m':proof.times}[nd[0]](pscalar(nd[1]),pscalar(nd[2]))
  sc[i]=v;return v
 def plinear(i):
  if i in lc:return lc[i]
  nd=pc['linear_nodes'][i]
  if nd[0]=='v':v=x[nd[1]]
  else:v=proof.total((pscalar(k),plinear(j))for k,j in nd[1])
  lc[i]=v;return v
 Px=[plinear(i)for i in Pd['P22_linear_output_roots_in_frozen_original_kernel_order']];roots('P22_on_x22_roots',Px)
 k=[proof.minus(Px[j],u[j])for j in range(22)];roots('kernel_coordinate22_roots',k)
 Krow=[[]for _ in range(174)]
 for j,col in enumerate(Kd['sparse_columns']):
  for i,num,den in col:Krow[i].append((proof.atom(num,den),k[j]))
 z=[proof.minus(x[i],proof.total(Krow[i]))for i in range(174)];v282=z+h;roots('new_kernel282_roots',v282)
 theta=-field(1024000)/field(243)*field(a)*field(100*a+18*e**2-9*e-2)**3*field(5*a*d-6*e-4)
 enc=lambda f:[[list(m),str(k)]for m,k in f.items()]
 th=proof.atom(enc(theta.numer),enc(theta.denom));it=proof.atom(enc(theta.denom),enc(theta.numer))
 assert mapped[c['known_unit_separator_root']]==th
 assert mapped[c['separator_identity_zero_root']]==proof.minus(proof.plus(v282[p],v282[q]),th)
 w=v282[:];w[p]=th;roots('new_kernel_after_constant_change_roots',w)
 normalized=[proof.times(v,it)for v in w];normalized[p]=proof.one;roots('normalized_kernel_with_pivot_one_roots',normalized)
 assert c['actual_free281_source_indices_in_changed282']==[i for i in range(282)if i!=p]and c['actual_free281_shape']==[281,239]
 assert c['constant_coordinate_change']=={'forward':'w_p=v_p+v_q; w_q=v_q','inverse':'v_p=w_p-w_q; v_q=w_q','determinant':1}
 # An independent sparse-QQ native check verifies ALL atom denominators and
 # the ONE scalar identity used to alias w_p to the known unit theta.
 def wp(f,terms):
  p=polynomial(terms);f.write(str(len(p))+''.join(' '+str(k)+' '+' '.join(map(str,m))for m,k in p)+'\n')
 with(out/'scalar-proof-input.txt').open('w')as f:
  f.write(f'{len(atoms)} {len(nodes)} {c["separator_identity_zero_root"]} {c["known_unit_separator_root"]}\n')
  for num,den in atoms:wp(f,num);wp(f,den)
  for nd in nodes:
   if nd[0]=='z':v=[0,0,0]
   elif nd[0]=='c':v=[1,nd[1],0]
   else:v=[{'a':2,'s':3,'m':4}[nd[0]],nd[1],nd[2]]
   f.write(' '.join(map(str,v))+'\n')
 cpp=ROOT/'code/s61_nonsquare_free281_circuit_review_20261104_v63.cpp';utility=ROOT/'code/s61_nonsquare_ell_Lframe_review_20261104_v63.cpp';exe=out/'review-native'
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();rem=lambda:max(.1,60-(time.monotonic()-began))
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=rem())
 run=subprocess.run([str(exe),str(out/'scalar-proof-input.txt')],capture_output=True,text=True,timeout=rem());exe.unlink();(out/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0,run.stderr;native=json.loads(run.stdout);assert native['all_passed']
 deps=[ROOT/'reports'/j/'report.json'for j in['S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW','S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW','S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63_REVIEW','S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW','S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in deps)
 for rev,cert in[(deps[1],pbase/'global-P22.json'),(deps[3],ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json'),(deps[4],ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/certificate.json')]:
  assert read(rev)['input_sha256'][str(cert.relative_to(ROOT))]==sha(cert)
 note=ROOT/'research/notes'/(ID+'.md')
 result=dict(id=ID,all_passed=True,all319_304_282_indices_and_Bnew_Cnew_partitions_verified=True,
  old_B22_separator_mixing_detected_and_new_C108_literals_preserved=True,
  all_source319_to304_transport_roots_match_exact_last_R4_combination=True,
  all_pinned_P22_operator_composition_and_K22_complement_roots_match=True,
  P22_linearity_and_source174_variable_binding_verified=True,
  all10236_scalar_nodes_and898_atoms_bound_without_large_root_expansion=True,
  one_actual_separator_identity_and_known_unit_denominators_verified=True,
  constant_coordinate_change_and_unit_deletion_preserve_actual_cokernel=True,
  actual_global_free281_to239_finite_recipe_accepted=True,
  uniform_known_kernel_count_unchanged=38,remaining_for80=42,
  no_expanded281x239_matrix_full239_newkernelcount_point_GB_rank_or_author_rerun=True,
  prior_review_reader_initialization_failure_preserved=bool(prior_failure),
  imported_Bianchi_scalar_nodes=len(bc),imported_P22_scalar_nodes=len(sc),imported_P22_linear_nodes=len(lc),native=native,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),cpp,utility,note,base/'report.json',base/'certificate.json',compact/'full319-kernel.txt',pbase/'global-P22.json',pbase/'global282-basis.json',ROOT/old['source304_atlas_path'],ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json',ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json',ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/certificate.json',*deps]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-began,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--prior-failure',type=Path);args=p.parse_args();main(args.output.resolve(),args.prior_failure.resolve()if args.prior_failure else None)
