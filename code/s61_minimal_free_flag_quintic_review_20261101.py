#!/usr/bin/env python3
"""Independent unreduced numeric fifth-jet check and filtered-proof review."""
import sys,json,time,os,subprocess,tempfile
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import numpy as np
import sympy as s
from sympy.polys.rings import PolyElement
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_free_flag_filtered_data_20261101 import build_four
from s61_minimal_free_flag_universal_quintic_20261101 import fifth_quotient,raw_fifth

import argparse,hashlib
ROOT=Path(__file__).resolve().parents[1]

def verify():
    start=time.time();p=32003;a=s.symbols('a0:5');sample=dict(zip(a,range(1,6)))
    q=json.load(open(ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101.json'));nf=json.load(open(ROOT/'reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json'));pt=json.load(open(ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json'))['strata']['free_M0']['local_data'];cr=json.load(open(ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101.json'));qr=json.load(open(ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101.json'))
    def mod(x):
     x=s.Rational(x);return int(x.p)*pow(int(x.q),-1,p)%p
    def array(M):return np.array([[mod(x)for x in row]for row in M.tolist()],dtype=np.int64)
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1];beta=s.Matrix(nf['normalized_beta16']).applyfunc(s.sympify).subs(sample);two,F=specialize(model,dict(zip(one,beta)),True);fixed=two.index(pt['chart_fixed_u']-1);cols=pt['chart_columns'];q0=s.Matrix(pt['point61']).applyfunc(s.Rational)[two,:];J=s.zeros(135,45);HH=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if len(mon)==2:
       i,j=mon
       if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
       if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
      elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
      if fixed not in mon:
       i,j=[cols.index(k)for k in mon]+[44]*(2-len(mon));cc=mod(c)
       if i==j:HH[row,i,j]=(HH[row,i,j]+cc)%p
       else:HH[row,i,j]=(HH[row,i,j]+cc*pow(2,-1,p))%p;HH[row,j,i]=HH[row,i,j]
    T=s.Matrix(q['tangent_basis_after_shift']).applyfunc(s.sympify);T[:,1]-=a[1]*T[:,0]/6;T=T.subs(sample);U=s.Matrix(q['right_inverse']).applyfunc(s.sympify).subs(sample);assert J*T==s.zeros(135,11);rows=q['pivot_rows'];assert J[rows,:]*U==s.eye(34);j=array(J);uu=array(U);JU=j@uu%p;W={1:array(T)};A={};mons={d:list(cwr(range(11),d))for d in range(1,6)};ix={d:{m:i for i,m in enumerate(mons[d])}for d in mons}
    for degree in range(2,6):
     Q=np.zeros((135,len(mons[degree])),dtype=np.int64)
     for da in range(1,degree//2+1):
      db=degree-da;fac=1 if da==db else 2;indices=np.array([ix[degree][tuple(sorted(x+y))]for x in mons[da]for y in mons[db]])
      for row in range(135):
       cross=W[da].T@(HH[row]@W[db]%p)%p;acc=np.bincount(indices,weights=cross.ravel(),minlength=len(mons[degree])).astype(np.int64);Q[row]=(Q[row]+fac*acc)%p
     W[degree]=-uu@Q[rows,:]%p;A[degree]=(Q-JU@Q[rows,:])%p
     print('RAW',degree,time.time()-start,flush=True)
    state=build_four();R=state['R'];vals=list(range(1,6));memo={}
    def cv(x):
     if isinstance(x,PolyElement):
      k=tuple(x.items())
      if k not in memo:memo[k]=R.ground_new(sum(c*s.prod(t**e for t,e in zip(vals,mon))for mon,c in x.items()))
      return memo[k]
     if isinstance(x,s.MatrixBase):return x.xreplace(sample)
     if isinstance(x,dict):return{k:cv(v)for k,v in x.items()}
     if isinstance(x,list):return[cv(v)for v in x]
     if isinstance(x,tuple):return tuple(cv(v)for v in x)
     if isinstance(x,s.Basic):return x.xreplace(sample)
     return x
    state=cv(state);print('STATE',time.time()-start,flush=True)
    # Verify numeric unreduced W2 and W3 against the source state before testing highest-order reduction.
    for d in[2,3]:
     got=np.array([[mod(v.as_expr())*pow(pow(2,d-1,p),-1,p)%p for v in row]for row in state['W'+str(d)]],dtype=np.int64)
     assert np.array_equal(got,W[d]),d
    project,free5,_=fifth_quotient(state,mons);got=raw_fifth(state,mons,project);difference=A[5].copy();den=pow(pow(2,4,p),-1,p)
    for i,row in enumerate(got):
     for k,v in enumerate(row):difference[i,free5[k]]=(difference[i,free5[k]]-mod(v.as_expr())*den)%p
    print('REDUCED_FIFTH',time.time()-start,flush=True)
    def poly(vals,ms):return '+'.join(str(int(c))+'*'+'*'.join('v_'+str(k)for k in mon)for c,mon in zip(vals,ms)if c)or'0'
    lower=[poly([mod(s.Rational(c))for c in row],list(map(tuple,q['quadratic_monomials'])))for row in q['constant_quadratic_rref']]
    v=s.symbols('v0:11')
    for expression in cr['constant_cubic_polynomials']:
     pp=s.Poly(s.sympify(expression),*v);lower.append('+'.join(str(mod(c))+'*'+'*'.join('v_'+str(j)for j,e in enumerate(ex)for _ in range(e))for ex,c in pp.terms())or'0')
    for row in qr['new_quartic_rref']:
     lower.append(poly([mod(s.sympify(c).subs(sample))for c in row],[tuple(qr['monomials'][j])for j in qr['quotient_free_columns']]))
    program='R=ZZ/32003[v_0..v_10];I=ideal('+','.join(lower)+');G=gb I;assert(hilbertFunction(5,R/I)==763);\n'
    program+='D=matrix{{'+','.join(poly(row,mons[5])for row in difference)+'}};assert(D%G==0);print("RAW UNREDUCED FIFTH VS HIGHEST-ONLY REDUCTION ALL135 PASS");\n'
    with tempfile.TemporaryDirectory()as td:
     f=Path(td)/'verify.m2';f.write_text(program);env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib';run=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(f)],env=env,capture_output=True,text=True,timeout=60);print(run.stdout,run.stderr,time.time()-start,flush=True);assert run.returncode==0
    paths=[Path(__file__)]+[ROOT/name for name in [
     'code/s61_minimal_free_flag_filtered_data_20261101.py',
     'code/s61_minimal_free_flag_universal_quintic_20261101.py',
     'code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py',
     'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101.json',
     'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101.json',
     'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101.json',
     'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101.json',
     'reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json',
     'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json',
     'research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101.md']]
    primary=json.loads((ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101.json').read_text())
    for name,h in primary['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    return dict(id='S61_MINIMAL_FREE_FLAG_QUINTIC_REVIEW_20261101',all_passed=True,prime=p,
     sample={str(k):str(v)for k,v in sample.items()},raw5_residual_entries=135*3003,
     modulus_degree5_quotient_dimension=763,all135_differences_zero=True,
     independent_unreduced_raw135_numeric_jets=True,primary_reduced_formula_specialized=True,
     lower_W2_W3_agree=True,accepted_true_hilbert_function=[1,11,54,168,391,740],
     proof_note='research/notes/S61_MINIMAL_FREE_FLAG_QUINTIC_REVIEW_20261101.md',
     limitations=['The new coefficient comparison is at one sample over GF(32003), not an independent universal QQ contraction.',
     'The full1180-row QQ identity and rank23 remain the primary producer certificate.',
     'No full initial-ideal equality, flatness, t-nilpotence, or endpoint removal.'],
     m2_log=run.stdout+run.stderr,input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths},seconds=time.time()-start)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
