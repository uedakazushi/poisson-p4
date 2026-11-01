#!/usr/bin/env python3
"""Modular truncated standard-basis exploration at the free-M0 flag endpoint."""
import argparse, base64, gzip, hashlib, json, subprocess, tempfile, time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
import numpy as np
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from truncated_local_ideal_io_20261101 import write_input,read_initial
ROOT=Path(__file__).resolve().parents[1]
def build_jets(maximum_degree):
    r=json.load(open(ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json'))['strata']['free_M0']['local_data'];p=32003;m=S61();pt=s.Matrix(r['point61']).applyfunc(s.Rational);b={i:pt[i]for i,wg in enumerate(m.reduced['parameter_weights'])if wg==1};two,F=specialize(m,b,True);fixed=two.index(r['chart_fixed_u']-1);cols=r['chart_columns'];J=s.Matrix(r['jacobian']);T=s.Matrix(r['tangent_basis']);T=s.Matrix.hstack(T[:,3],T[:,4],T[:,5]-T[:,1],T[:,6:11],T[:,0],T[:,2],T[:,1]);n=T.cols

    def mod(v):v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,p)%p

    def arr(M):return np.array([[mod(v)for v in row]for row in M.tolist()],dtype=np.int64)

    def rref(A):
     A=A.copy()%p;row=0;piv=[]
     for col in range(A.shape[1]):
      choices=np.flatnonzero(A[row:,col])
      if len(choices)==0:continue
      j=row+int(choices[0]);A[[row,j]]=A[[j,row]];A[row]=A[row]*pow(int(A[row,col]),-1,p)%p;mask=np.flatnonzero(A[:,col]);mask=mask[mask!=row]
      if len(mask):A[mask]=(A[mask]-A[mask,col,None]*A[row,None,:])%p
      piv.append(col);row+=1
      if row==A.shape[0]:break
     return A,piv
    j=arr(J);t=arr(T);_,C=rref(j);_,R=rref(j[:,C].T);rank=len(C);rr,_=rref(np.concatenate([j[np.ix_(R,C)],np.eye(rank,dtype=np.int64)],axis=1));U=np.zeros((45,rank),dtype=np.int64);U[C,:]=rr[:,rank:];JU=j@U%p;H=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      i,l=[cols.index(v)for v in mon]+[44]*(2-len(mon));v=mod(c)
      if i==l:H[row,i,l]+=v
      else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p
    H%=p;W={1:t};A={};ms={d:list(combinations_with_replacement(range(n),d))for d in range(1,maximum_degree+1)};idx={d:{m:i for i,m in enumerate(ms[d])}for d in range(1,maximum_degree+1)}
    for degree in range(2,maximum_degree+1):
     Qd=np.zeros((135,len(ms[degree])),dtype=np.int64)
     for da in range(1,degree//2+1):
      db=degree-da;fac=1 if da==db else 2;indices=np.array([idx[degree][tuple(sorted(ma+mb))]for ma in ms[da]for mb in ms[db]])
      for row in range(135):
       cross=(W[da].T@(H[row]@W[db]%p))%p;acc=np.bincount(indices,weights=cross.ravel(),minlength=len(ms[degree])).astype(np.int64);Qd[row]=(Qd[row]+fac*acc)%p
     W[degree]=-U@Qd[R,:]%p;A[degree]=(Qd-JU@Qd[R,:])%p
     print('degree',degree,flush=True)
    jet={d:v.tolist()for d,v in A.items()};ts={d:W[d][-1,:].tolist()for d in A}
    return n,jet,ts

def parse_lengths(text):
    return [int(line.split()[-1]) for line in text.splitlines() if line.startswith('modm ')]
def run(maximum_degree,m2_binary=None):
    start=time.time();n,jets,tjets=build_jets(maximum_degree);assert n==11
    with tempfile.TemporaryDirectory(prefix='poisson-local-standard-') as temp:
        td=Path(temp);exe=td/'standard';ref=td/'closure'
        subprocess.run(['clang++','-O3','-std=c++17',str(ROOT/'code/truncated_local_standard_basis_20261101.cpp'),'-o',str(exe)],check=True,capture_output=True)
        subprocess.run(['clang++','-O3','-std=c++17',str(ROOT/'code/truncated_local_ideal_20261101.cpp'),'-o',str(ref)],check=True,capture_output=True)
        six=td/'six.txt';write_input(six,32003,n,jets,tjets,maximum_degree=6)
        old=subprocess.run([str(ref),str(six),str(td/'old.txt')],check=True,capture_output=True,text=True,timeout=120)
        new=subprocess.run([str(exe),str(six),str(td/'new.txt')],check=True,capture_output=True,text=True,timeout=60)
        lengths=parse_lengths(old.stdout);assert lengths==parse_lengths(new.stdout)==[1,12,66,234,625,1365,2568]
        # The reference returns a homogeneous echelon basis in every degree.
        oldrows=read_initial(td/'old.txt',n,6);newrows=read_initial(td/'new.txt',n,6)
        piv={min(row):row for row in oldrows};assert len(piv)==len(oldrows)
        for row in newrows:
            row=dict(row)
            while row:
                lead=min(row);assert lead in piv;a=row[lead]*pow(piv[lead][lead],-1,32003)%32003
                for mon,c in piv[lead].items():
                    v=(row.get(mon,0)-a*c)%32003
                    if v:row[mon]=v
                    else:row.pop(mon,None)
        inp=td/'input.txt';write_input(inp,32003,n,jets,tjets,maximum_degree=maximum_degree)
        done=subprocess.run([str(exe),str(inp),str(td/'initial.txt')],check=True,capture_output=True,text=True,timeout=300)
        initial=read_initial(td/'initial.txt',n,maximum_degree);hf=parse_lengths(done.stdout)
        assert hf[:7]==lengths
        raw=(td/'initial.txt').read_bytes();encoded=base64.b64encode(gzip.compress(raw,mtime=0)).decode();cas=None
        if m2_binary is not None:
            def expression(row):return '+'.join(str(c)+'*'+'*'.join('v_'+str(j)for j in mon) for mon,c in row.items())
            program='R=ZZ/32003[v_0..v_10];J=ideal('+','.join(expression(row)for row in initial)+');\n'
            program+='print("DIMDEG",dim(R/J),degree(R/J));print("HILBERT",hilbertSeries(R/J));stdio<<flush;\n'
            if maximum_degree==10:
                f=json.loads((ROOT/'reports/S61_MINIMAL_FREE_FLAG_FLAT4_PROBE_20261101.json').read_text());V=s.Matrix(f['endpoint_tangent']).applyfunc(s.Rational);C=s.zeros(11)
                for j,i in enumerate([3,4,5,6,7,8,9,10,0,2,1]):C[i,j]=1
                C[1,2]=-1;V=C.inv()*V;v=s.Matrix(s.symbols('v_0:11'));ann=s.Matrix.hstack(*V.T.nullspace()).T;rows=list(V.T.to_DM().rref()[1]);e=V[rows,:].inv()*v[rows,:]
                forms=list(ann*v)+[e[i]*e[j]for i in range(1,4)for j in range(i,4)]
                program+='Dactual=ideal('+','.join(str(s.expand(f)).replace('**','^')for f in forms)+');\n'
                program+='BF=saturate(J,ideal(v_9*(8*v_8+97*v_4)));BB=saturate(J,ideal(v_9*v_2));DD=saturate(J,ideal(v_0+2*v_4));\n'
                program+='assert(dim(R/J)==3 and degree(R/J)==70);assert(degree(R/BF)==5 and degree(R/BB)==65 and degree(R/DD)==4);assert(DD==Dactual and isSubset(J,Dactual));\n'
                program+='E=intersect(BF,BB,Dactual)/J;assert(dim E==2 and degree E==31);print("PARTS",5,65,4,"EXTRA",dim E,degree E);print("EXTRA_HILBERT",hilbertSeries E);stdio<<flush;\n'
                program+='PF=ideal(v_7,v_6,v_5,v_4+4*v_10,v_3,v_2,v_1,v_0-8*v_10);PB=ideal(v_7,v_6,v_5,97*v_4+8*v_8,v_3,97*v_2+2*v_8-97*v_10,v_1,97*v_0-16*v_8);PL=ideal(v_10,v_9,v_7,v_6,v_5,v_4,v_3,v_2,v_1,9*v_0+v_8);PS=minimalPrimes J;assert(#PS==3 and any(PS,P->P==PF) and any(PS,P->P==PB) and any(PS,P->P==PL));print("ALL_CHECKS_PASS");\n'
            script=td/'probe.m2';script.write_text(program);rr=subprocess.run([str(m2_binary),'--script',str(script)],check=True,text=True,capture_output=True,timeout=180);cas=dict(output=rr.stdout,script=program)
            if maximum_degree==10:assert 'ALL_CHECKS_PASS' in rr.stdout
    deps=['code/truncated_local_standard_basis_20261101.cpp','code/truncated_local_ideal_20261101.cpp','code/truncated_local_ideal_io_20261101.py','code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json','reports/S61_MINIMAL_FREE_FLAG_FLAT4_PROBE_20261101.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    paths=[Path(__file__)]+[ROOT/p for p in deps]
    return dict(id='S61_MINIMAL_FREE_FLAG_INITIAL_PROBE_20261101',all_passed=True,prime=32003,variables=n,maximum_degree=maximum_degree,cumulative_lengths=hf,initial_rows_encoding='base64 of gzip of native sparse initial-row format',initial_rows_compressed=encoded,initial_rows_sha256=hashlib.sha256(raw).hexdigest(),initial_rows_count=len(initial),cas=cas,reference_six_jet=dict(algorithm='independent full sparse Macaulay row closure',cumulative_lengths=lengths,all_new_initial_rows_in_old_span=True,old_echelon_rows=len(oldrows),new_standard_basis_rows=len(newrows)),native_stdout=done.stdout.splitlines(),scope=['Exact finite-field truncated ideals only.','The exported homogeneous ideal is contained in the full initial ideal; equality with the full initial ideal is not asserted.','No characteristic-zero t-nilpotence or saturation-removal conclusion.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.time()-start)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--maximum-degree',type=int,default=10,choices=range(6,11));parser.add_argument('--m2-binary',type=Path);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();result=run(args.maximum_degree,args.m2_binary)
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'],cumulative_lengths=result['cumulative_lengths'])))
