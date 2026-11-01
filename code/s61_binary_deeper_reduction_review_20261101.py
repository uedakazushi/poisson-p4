#!/usr/bin/env python3
"""Review the deeper B-rank locus and independently rebuild primitive group action."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
import sympy as s
from itertools import combinations
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def review():
 start=time.time();path=ROOT/'reports/S61_BINARY_DEEPER_REDUCTION_20261101.json';r=json.loads(path.read_text());assert r['all_passed']
 # Replay the explicit ideal certificate; its matrix is a saved input from the reviewed block decomposition.
 env=os.environ.copy();env['DYLD_LIBRARY_PATH']=':'.join([str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib','/opt/homebrew/opt/givaro/lib'])
 with tempfile.TemporaryDirectory(dir=ROOT/'work') as temp:
  f=Path(temp)/'review.m2';f.write_text(r['m2_certificate']);run=subprocess.run([r['m2_binary'],'--script',str(f)],env=env,capture_output=True,text=True,timeout=90)
 assert run.returncode==0 and 'PASS' in run.stdout,(run.stdout,run.stderr)
 m=S61();z,bu,beta=primitive_beta(m);v=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS]);B=v.jacobian(bu);rows=list(B.T.to_DM().rref()[1]);left=B[rows,:].inv()
 action=s.Matrix.hstack(*[left*infinitesimal_action(v,A)[rows,:] for A in sp_basis()]).applyfunc(s.expand)
 u=s.symbols('u1:62');out=[]
 rank_source=json.loads((ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json').read_text())
 basis=s.Matrix.hstack(*[A.reshape(16,1) for A in sp_basis()]);br=list(basis.T.to_DM().rref()[1]);bleft=basis[br,:].inv()
 stab=s.Matrix.hstack(*[bleft*s.Matrix(A).applyfunc(s.sympify).reshape(16,1)[br,:] for A in rank_source['stabilizer_generators']])
 comp=[2,4,5,6,7];assert stab.row_join(s.eye(10)[:,comp]).det()!=0
 for c in r['excluded_components']:
  ideal=[s.sympify(f) for f in r['reduced_prime_ideals'][c['ideal_index']]]
  cons=s.Matrix([u[11],u[12]]+ideal);sub={s.Symbol(k):s.sympify(v) for k,v in c['dense_parameterization'].items()}
  assert cons.subs(sub,simultaneous=True).applyfunc(s.cancel)==s.zeros(len(cons),1)
  N=(cons.jacobian(bu)*action).subs(sub,simultaneous=True).applyfunc(s.cancel)
  assert (N*stab).applyfunc(s.cancel)==s.zeros(N.rows,5)
  V=N[:,comp];den=s.lcm([s.denom(v) for v in V]);V=V.applyfunc(lambda v:s.cancel(v*den))
  size=c['group_normal_rank_upper_bound']+1;zero_minors=0
  for rr in combinations(range(V.rows),size):
   for cc in combinations(range(5),size):
    M=V[list(rr),list(cc)].to_DM();assert M.det()==M.domain.zero;zero_minors+=1
  vals={x:i+2 for i,x in enumerate(bu)};sample=V.subs(vals,simultaneous=True);rank=sample.rank()
  assert rank==size-1
  image=c['source_dimension']+rank
  assert image<=15
  out.append(dict(ideal_index=c['ideal_index'],all_ten_group_normal_rank=rank,beta_image_dimension=image));print('component',c['ideal_index'],'image',image,flush=True)
 paths=[Path(__file__),path,ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json']
 return dict(id='S61_BINARY_DEEPER_REDUCTION_REVIEW_20261101',all_passed=True,m2_certificate_replayed=True,independent_primitive_action_results=out,proof_audit=['Four stated primes and square containment give the radical equality.','The dense charts cover dense opens of the prime sources; their boundary does not create a separate source component.','All ten group directions were recomputed from primitive beta.','Only B-rank<=7 on K=0 is classified; C-block exceptions remain separate.'],limitations=['The B block and Macaulay2 membership script are supplied certificates, not independently reconstructed from all135 equations.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=review()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','independent_primitive_action_results','seconds']}))
