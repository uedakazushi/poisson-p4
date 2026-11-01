#!/usr/bin/env python3
"""Independent raw double-line check and replay/proof audit of the C3 formal ring."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
ROOT=Path(__file__).resolve().parents[1]
def review():
 start=time.time();path=ROOT/'reports/S61_BINARY_MIXED_FORMAL_20261101.json';r=json.loads(path.read_text());assert r['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];p=s.zeros(61,1)
 for k,v in r['point']['sample'].items():p[int(k[1:])-1]=s.Rational(v)
 p[29]=1;local=calculate(m,True,p);assert local['homogenized_jacobian_rank']==33
 cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols);T=s.Matrix(r['point']['tangent_basis']).applyfunc(s.Rational);J=s.Matrix(local['jacobian']).applyfunc(s.Rational);assert J*T==s.zeros(135,12)
 delta=s.Matrix(r['point']['actual_double_line_tangent']).applyfunc(s.Rational);nil=s.Matrix(r['point']['actual_double_line_nil_direction']).applyfunc(s.Rational);ss,ee=s.symbols('s e');q=s.zeros(45,1);q[fixed]=1
 for i,c in enumerate(cols):q[c]=ss*delta[i]+ee*nil[i]
 two,F=specialize(m,{i:p[i] for i in one},True)
 for f in F:
  expr=s.Poly(s.expand(sum(co*s.prod(q[i] for i in mon)*(ss*delta[-1]+ee*nil[-1])**(2-len(mon)) for mon,co in f.items())),ss,ee)
  assert all(c==0 or mon[1]>=2 for mon,c in expr.terms())
 assert s.Matrix.hstack(delta,nil).rank()==2 and delta[-1]==nil[-1]==0
 source=r['macaulay2_program'];assert hashlib.sha256(source.encode()).hexdigest()==r['source_program_sha256']
 env=os.environ.copy();env['DYLD_LIBRARY_PATH']=':'.join([str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib','/opt/homebrew/opt/givaro/lib'])
 with tempfile.TemporaryDirectory(dir=ROOT/'work') as tmp:
  f=Path(tmp)/'replay.m2';f.write_text(source);run=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(f)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=180)
 assert run.returncode==0 and 'CERTIFICATE PASS' in run.stdout,(run.stdout+run.stderr)[-3000:]
 paths=[Path(__file__),path,ROOT/'research/notes/S61_BINARY_MIXED_FORMAL_20261101.md',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py']
 return dict(id='S61_BINARY_MIXED_FORMAL_REVIEW_20261101',all_passed=True,independent_raw_local_rank=33,independent_actual_double_line=True,qq_macaulay_certificate_replayed=True,proof_audit=['The homogeneous uniform tangent scaling preserves the initial ideals and linear double-line ideal.','v1 saturation and coefficient saturation with a monic normal Groebner basis prove base torsion-freeness, hence C-primaryness.','The kernel of A/J -> grS lies in D/J, a C-primary module; equality of top multiplicity65 forces that kernel to vanish.','Strict induced filtration gives grN=D/J; its finite torsion-free base module structure lifts to N.','Good filtration gives rankN=65; the same actual binary component exhausts the generic part, where I^4=0 and t in I^2.','t^2 lies in N and vanishes generically, so N torsion-freeness gives full t^2=0.','QQ low rank226 and modular augmented227 prove t nonzero modulo m^4.'],limitations=['The large QQ Schur/primary program is replayed from the supplied certificate, not independently regenerated.','Only the specified fixed-beta local ring; no relative or whole-stratum removal.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=review()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','seconds']}))
