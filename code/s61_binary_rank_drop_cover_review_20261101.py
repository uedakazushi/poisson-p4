#!/usr/bin/env python3
"""Independent finite ideal and geometric audit of the binary rank-drop cover."""
import argparse,json,hashlib,time,subprocess,os,tempfile
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def verify():
 start=time.time();rp=ROOT/'reports/S61_BINARY_RANK_DROP_COVER_20261101.json';r=json.loads(rp.read_text());hp=ROOT/'reports/S61_BINARY_HYPERSURFACE_REVIEW_20261101.json';assert r['all_passed'] and json.loads(hp.read_text())['all_passed']
 m=S61();u=s.symbols('u1:62');one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two,F=specialize(m,{i:(s.S.Zero if i in [11,12] else u[i]) for i in one},True);fixed=two.index(29);D=s.zeros(135,45);tc=s.zeros(135,1)
 for row,f in enumerate(F):
  for mon,c in f.items():
   if len(mon)==2:
    for k in [0,1]:
     if mon[1-k]==fixed:D[row,mon[k]]+=c
   elif len(mon)==1 and mon[0]==fixed:tc[row]=c
 K=s.Matrix.hstack(*D[:70,:].nullspace());E=(D[70:,:]*K).applyfunc(s.expand);block=json.loads((ROOT/'reports/S61_BINARY_HYPERSURFACE_RANK_DROP_20261101.json').read_text());rr=block['residual_block_rows'];cc=block['residual_block_columns'];A,B,C=[E[rr[i],cc[i]] for i in range(3)];C=C[[0,1,2,3,4,5,6,7,8,10,12,14],:]
 assert [A,B,C]==[s.Matrix(M).applyfunc(s.sympify) for M in r['block_matrices']]
 core=list(map(s.Symbol,r['core_variables']));a,b,c,d,e,f=core;kk=5*c+f;ll=5*b-e;P=600*c*c+25*b*d-155*d*e+230*c*f-18*f*f;Q=600*b*b-25*a*c-255*b*e-63*e*e+95*a*f
 assert s.expand(P.subs(b,31*e/5)-2*(15*c-f)*(20*c+9*f))==0
 assert s.expand(Q.subs(c,19*f/5)-3*(5*b-3*e)*(40*b+7*e))==0
 assert s.gcd(P,Q)==1
 prime=[[kk],[ll],[d,15*c-f],[5*b-3*e,a],[P,Q]];assert prime==[list(map(s.sympify,g)) for g in r['prime_ideals']]
 expr=lambda v:str(v).replace('**','^')
 mat=lambda M:'matrix{'+','.join('{'+','.join(expr(v) for v in row)+'}' for row in M.tolist())+'}'
 script='R=QQ['+','.join(map(str,core))+'];\nA='+mat(A)+';B='+mat(B)+';C='+mat(C)+';\n'
 script+='IB=minors(9,B);IA=minors(4,A);IC=minors(5,C);\n'
 for i,g in enumerate(prime):script+='J'+str(i)+'=ideal('+','.join(map(expr,g))+');\n'
 script+='J=intersect{J4,J3,J2,J1,J0}; G=groebnerBasis IB;\n'
 script+='assert all(flatten entries gens IB,g->((g%J)==0));\nassert all(flatten entries gens J,g->(((g^2)%G)==0));\n'
 script+='kl='+expr(kk*ll)+';assert (((kl^2)%IA)==0);assert (((kl^3)%IC)==0);\nprint(\"INDEPENDENT_MEMBERSHIPS_PASS\");\n'
 env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
 with tempfile.TemporaryDirectory(dir=ROOT/'work') as tmp:
  path=Path(tmp)/'audit.m2';path.write_text(script);run=subprocess.run([r['m2_binary'],'--script',str(path)],env=env,capture_output=True,text=True,timeout=60)
 assert run.returncode==0 and 'INDEPENDENT_MEMBERSHIPS_PASS' in run.stdout,(run.stdout,run.stderr)
 left=s.Matrix.hstack(*D[:70,:].T.nullspace()).T;F0=(left*tc[:70,:]).applyfunc(s.expand);assert s.expand(10*F0[8]-20*F0[28]-ll**2)==0 and s.expand(10*F0[13]+50*F0[29]-kk**2)==0
 assert all(s.expand(v.subs({e:5*b,f:-5*c}))==0 for v in F0)
 for case in r['orbit_strata']:
  M=s.Matrix(case['beta_normal_action']).applyfunc(s.Rational);assert M[:,case['beta_minor_columns']].det()==s.Rational(case['beta_minor'])!=0
 paths=[Path(__file__),rp,hp,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'research/notes/S61_BINARY_RANK_DROP_COVER_20261101.md']
 return dict(id='S61_BINARY_RANK_DROP_COVER_REVIEW_20261101',all_passed=True,raw_blocks_reconstructed=True,independent_m2_script=script,independent_m2_output=run.stdout,geometric_CI_primality_audited=True,first70_t_radical_verified=True,accepted_conclusion='All nonsquare binary full rank<=34 points lie in exactly the displayed five prime loci, with three symmetry orbit families.',limitations=['No saturation exhaustion follows merely from the rank-drop cover.','Beta normal action entries use the primary report; their nonzero minors are recomputed.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','seconds']}))
