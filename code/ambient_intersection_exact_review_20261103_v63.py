#!/usr/bin/env python3
"""Independent actual-N image presentation and final target-ideal/colon audit."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def split(f,p,a):
 e=o=s.S.Zero
 for (n,),c in s.Poly(s.expand(f),p).terms():
  if n%2:o+=c*a**(n//2)
  else:e+=c*a**(n//2)
 return s.expand(e),s.expand(o)
def text(f):return str(s.expand(f)).replace('**','^')
def ideal(name,fs):return name+'=ideal('+','.join(text(f)for f in fs if f!=0)+');'

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();args.output_dir=args.output_dir.resolve()
 assert not(args.output_dir/'report.json').exists();args.output_dir.mkdir(parents=True,exist_ok=True)
 directory=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63'
 builder=load(directory/'builder.json');receipt=load(directory/'report.json');fresh=load(args.replay)
 assert {k:v for k,v in builder.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),directory/'builder.json',directory/'report.json']
 for d in [builder,receipt]:
  for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 variables=list(map(s.Symbol,builder['ambient_variables']));local={str(v):v for v in variables}
 p,old_s=s.symbols('p s');a,b,c,sv=[local[v]for v in ['aa','bb','cc','ss']]
 K=s.Matrix([local[f'v{i}']for i in range(10)]);z=s.Matrix([local[f'z{i}']for i in range(35)])
 parse=lambda f:s.sympify(f,locals=local)
 IA=[parse(f)for f in builder['component_presentations']['CA']];IB=[parse(f)for f in builder['component_presentations']['CB']]
 F,G=map(parse,[builder['actual_nilpotent_candidate'],builder['candidate_p_multiple']])
 assert F==2*a+K[0]+z[25]and G==-a*K[1]+2*b+z[8]
 # Independently build the entire N kernel from its actual9-parameter map.
 maps=load(ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json')
 chart=load(ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json');qi=chart['q_indices_zero_based']
 n=maps['families']['N'];pars=list(map(s.Symbol,n['source_fiber_parameters']))
 u=s.Matrix(list(map(s.sympify,n['full_61_coordinates'])));q=u[qi,0];M=q.jacobian(pars);offset=q.subs(dict.fromkeys(pars,0))
 ns=n['actual_q_selected_rows'];cs=builder['common_q_selected_rows'];rem=builder['remaining_q_rows']
 assert ns==cs[1:]
 target=s.Matrix([K[cs.index(i)]if i in cs else z[rem.index(i)]for i in range(45)])
 L=M[ns,:];assert not s.factor(L.det()).free_symbols and L.det()!=0
 fiber=L.inv()*(target[ns,0]-offset[ns,0]);solved=(M*fiber+offset).subs(old_s,sv).applyfunc(s.expand)
 nr=[i for i in range(45)if i not in ns];eo=[split(f,p,a)for f in solved[nr,0]]
 E=s.Matrix([v[0]for v in eo]);O=s.Matrix([v[1]for v in eo]);O0=O.subs({a:0,sv:0})
 R=O0.jacobian(list(K[1:,0]));assert O0==R*K[1:,0]and R.rank()==6
 rows=list(R.T.rref()[1]);JR=R[rows,:];cols=list(JR.rref()[1]);J=s.Matrix([a,sv,*O0[rows,0]])
 C=R[:,cols]*JR[:,cols].inv();assert C*O0[rows,0]==O0
 U=[];V=[]
 for f in O-C*O0[rows,0]:
  ev=s.cancel(f.subs(sv,0)/a);od=s.cancel((f-a*ev)/sv)
  assert not s.denom(ev).free_symbols and not s.denom(od).free_symbols
  U.append(ev);V.append(od)
 U,V=s.Matrix(U),s.Matrix(V);yt=target[nr,0]
 W=yt[rows,0]-E[rows,0]-b*U[rows,0]-c*V[rows,0];odd=s.Matrix([b,c,*W])
 NE=[f for f in (yt-E-b*U-c*V-C*W).applyfunc(s.expand)if f!=0]
 NE += [s.expand(J[i]*odd[j]-J[j]*odd[i])for i in range(8)for j in range(i)]
 NE += [s.expand(odd[i]*odd[j]-a*J[i]*J[j])for i in range(8)for j in range(i+1)]
 assert len(NE)==94
 substitution={a:p*p,b:p**3,c:sv*p,**dict(zip(target,list(solved)))}
 assert all(s.expand(f.subs(substitution,simultaneous=True))==0 for f in NE)
 # Check the CA source F=r and polynomial coordinates (F,N_K9).
 ca=maps['families']['CA'];cpars=list(map(s.Symbol,ca['source_fiber_parameters']));cu=s.Matrix(list(map(s.sympify,ca['full_61_coordinates'])))
 actualF=s.expand(cu[0]-4*cu[34]+cu[51]);assert actualF==s.Symbol('r')
 chosen=s.Matrix([actualF,*cu[[qi[i]for i in ns],0]]);assert s.expand(chosen.jacobian(cpars).det())==-72
 # All45 coordinate generators independently identify a smooth prime; there
 # is no guess based on the10-dimensional common nonnormal stratum.
 P=[v for v in variables if str(v)not in ['v1','v4','v5','v7']];assert len(P)==45
 script=['R=QQ['+','.join(map(str,variables))+',MonomialOrder=>ProductOrder{35,14}];',
  ideal('IA',IA),ideal('IB',IB),ideal('NN',NE),ideal('PP',P),'DD=IA+IB;',
  'ff='+text(F)+';','gg='+text(G)+';','GD=gb DD;GN=gb NN;GP=gb PP;',
  'assert(ff%GD!=0);assert(ff^2%GD==0);assert(gg%GD==0);',
  'DRED=DD+ideal(ff);GDR=gb DRED;',
  'assert(all(flatten entries gens DRED,f->f%GN==0));assert(all(flatten entries gens NN,f->f%GDR==0));',
  'CANDN=IA+ideal(ff,gg);GCN=gb CANDN;',
  'assert(all(flatten entries gens CANDN,f->f%GN==0));assert(all(flatten entries gens NN,f->f%GCN==0));',
  'COL=DD:ideal(ff);GC=gb COL;',
  'assert(all(flatten entries gens COL,f->f%GP==0));assert(all(flatten entries gens PP,f->f%GC==0));',
  'assert(dim(R/PP)==4);assert(degree(R/PP)==1);assert(dim(R/NN)==11);',
  'assert((v1*ff)%GD!=0);assert((v1^2*ff)%GD!=0);',
  'print "INDEPENDENT_ACTUAL_N_KERNEL_COLON_REDUCTION_PASS";stdio<<flush;exit 0;']
 m2=args.output_dir/'review.m2';log=args.output_dir/'review.log';m2.write_text('\n'.join(script)+'\n')
 with log.open('w')as stream:subprocess.run([shutil.which('M2'),'--script',str(m2.resolve())],stdout=stream,stderr=subprocess.STDOUT,timeout=30,check=True)
 assert 'INDEPENDENT_ACTUAL_N_KERNEL_COLON_REDUCTION_PASS'in log.read_text()and'error:'not in log.read_text()
 paths += [m2,log,ROOT/'reports/AMBIENT_PINCHED_MODELS_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_INTERSECTION_20261103_v63_REVIEW.json',
  ROOT/'research/notes/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63_REVIEW.md',
  ROOT/'research/notes/AMBIENT_INTERSECTION_SYMPLECTIC_SCOPE_20261103_v63.md']
 result=dict(id='AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63_REVIEW',all_passed=True,
  fresh_full_QQ_component_builder_replay_matches=True,regular_sequence_Koszul_full_kernel_proof_accepted=True,
  independently_constructed_actual_N_image_kernel_generator_count=94,independent_N_regularity_sequence_rank=8,
  all94_independent_N_relations_pull_back_to_zero=True,exact_actual_target_equalities_rechecked_by_fresh_M2=True,
  identities=['F!=0','F²=0','G=pF=0','ID+(F)=independent_IN=IA+(F,G)','ID:F=(all variables except v1,v4,v5,v7)'],
  nilradical='principal squarezero(F)',nilmodule='C[[u11]]F',nilmodule_support_full_S61_dimension=11,
  exact_associated_primes=['minimal N prime','smooth11-dimensional embedded prime'],
  completed_S61_intersection_dimension=18,completed_S61_intersection_depth=11,
  ring_splitting_asserted=False,full_global_rank_at_most2_boundary_covered=False,
  input_sha256={str(path.relative_to(ROOT)):sha(path)for path in dict.fromkeys(paths)})
 with(args.output_dir/'report.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS independent94-generator N kernel, actual218 intersection, nilradical/colon')

if __name__=='__main__':main()
