"""Fixed-beta exact weighted rings at r=2, -6/25 and -91/400; no relative assertion."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,tempfile,time
from fractions import Fraction
import sympy as s
from sympy.polys.rings import ring
from sympy import QQ
from s61_flag_linear_artifact_weighted_engine_20261102 import build
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FLAG_REMAINING_CORE_WEIGHTED_POINTS_20261102';NOTE='research/notes/'+STEM+'.md'
BASE='reports/S61_FLAG_ORDINARY_WEIGHTED_FINITE_KERNEL_20261102.json';WEIGHTS=(1,1,2,2,1,2,2)
CORES=[s.Rational(2),s.Rational(-6,25),s.Rational(-91,400)]
NUMS=[[1,2,7,12,18,9],[1,2,7,12,18,9,-34,-3],[1,2,7,12,18,9,-34,-9]]

def verify():
 start=time.monotonic();base=json.loads((ROOT/BASE).read_text());results=[]
 for p,h in base['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 K,*xs=ring(','.join('x'+str(i)for i in range(7)),QQ)
 def q(s0):
  f=Fraction(s0);return QQ(f.numerator,f.denominator)
 with tempfile.TemporaryDirectory(prefix='flag-three-special-points-')as td:
  td=Path(td);env=os.environ|{'DYLD_LIBRARY_PATH':str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'};exe=td/'rows';cpp='code/s61_flag_ordinary_weighted_rows_20261102.cpp'
  subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(ROOT/cpp),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(exe)],check=True)
  for k,core in enumerate(CORES):
   meta,m,eq,J,T,U=build(td/f'raw{k}',core);assert meta['normal_weights']==list(WEIGHTS)
   inp=td/f'den{k}';out=td/f'rows{k}';inp.write_text('\n'.join('1 '+meta[x]for x in['ud','td'])+'\n');subprocess.run([str(exe),str(td/f'raw{k}/raw.txt'),str(inp),str(out)],env=env,capture_output=True,text=True,check=True,timeout=60)
   def readp(f):
    ts=f.readline().split();assert len(ts)-1==int(ts[0])and len(ts)<=2
    return ts[1]if len(ts)>1 else'0'
   rows=[];expr=[];scales=[]
   with out.open()as f:
    assert int(f.readline())==136
    for i in range(136):
     ii,w,n=map(int,f.readline().split());assert ii==i;dd=readp(f);cc=readp(f);assert Fraction(dd)
     terms={};parts=[]
     for j in range(n):
      a=tuple(map(int,f.readline().split()));c=readp(f);assert sum(v*z for v,z in zip(a,WEIGHTS))==w;terms[a]=q(c)
      parts.append('('+c+')'+''.join('*x'+str(z)+('^'+str(e)if e>1 else'')for z,e in enumerate(a)if e))
     if terms:assert Fraction(cc)
     if i in meta['selected']:assert not terms
     rows.append(K.from_dict(terms));expr.append('+'.join(parts)or'0');scales.append(str(Fraction(cc)/Fraction(dd)))
   member=td/f'cof{k}';metadata=td/f'meta{k}';script=td/f'check{k}.m2'
   # Export the exact QQ standard basis and its leading monomials.  An
   # independent sparse division below checks the two observables without
   # constructing unused change-of-basis rows in all higher degrees.
   program='R=QQ[x0,x1,x2,x3,x4,x5,x6,Degrees=>{1,1,2,2,1,2,2}];F=matrix{{'+','.join(expr[:135])+'}};tt='+expr[-1]+';G=gb F;gg=gens G;assert(tt^2%G==0);assert(tt%G!=0);out=openOut "'+str(member)+'";out<<numColumns gg<<endl;scan(numColumns gg,i->{f=gg_(0,i);lm=first exponents leadMonomial f;scan(lm,e->out<<e<<" ");out<<endl;a=listForm f;out<<#a;scan(a,z->{out<<" "<<toString(last z);scan(first z,e->out<<" "<<e);});out<<endl;});close out;out=openOut "'+str(metadata)+'";out<<numColumns gg<<endl;out<<dim coker F<<endl;out<<toString(tt%G)<<endl;out<<toString hilbertSeries(coker F,Reduce=>true)<<endl;out<<toString apply(0..12,i->hilbertFunction(i,coker F))<<endl;L=ideal(F)+ideal(x0+2*x1+3*x4);out<<dim(R/L)<<endl;out<<toString hilbertSeries(R/L,Reduce=>true)<<endl;close out;'
   script.write_text(program);p=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(script)],env=env,capture_output=True,text=True,timeout=80);assert p.returncode==0,(p.stdout,p.stderr)
   lines=member.read_text().splitlines();ng=int(lines[0]);assert len(lines)==1+2*ng;basis=[]
   for i in range(ng):
    lm=tuple(map(int,lines[1+2*i].split()));ss=lines[2+2*i].split();n=int(ss[0]);assert len(ss)==1+8*n
    bp=K.from_dict({tuple(map(int,ss[2+8*j:9+8*j])):q(ss[1+8*j])for j in range(n)});assert lm in bp;basis.append((lm,bp))
   # Use the acyclic rewrite graph itself rather than silently assuming
   # Macaulay2 and this Python ring use the same term order.  In a fixed
   # weight there are only finitely many monomials.
   def reduce_observable(f):
    order=[];done=set();active=set();rule={}
    def visit(a):
     if a in done:return
     assert a not in active;active.add(a)
     for i,(lm,bp)in enumerate(basis):
      if all(l<=x for l,x in zip(lm,a)):
       shift=tuple(x-l for x,l in zip(a,lm));rule[a]=(i,shift)
       for b in bp:
        if b!=lm:visit(tuple(x+y for x,y in zip(b,shift)))
       break
     active.remove(a);done.add(a);order.append(a)
    for a in f:visit(a)
    rem=dict(f);steps=0
    for a in reversed(order):
     if a not in rule or not rem.get(a):continue
     i,shift=rule[a];lm,bp=basis[i];c=rem[a]/bp[lm]
     for b,bv in bp.items():
      dest=tuple(x+y for x,y in zip(b,shift));rem[dest]=rem.get(dest,QQ.zero)-c*bv
     assert not rem[a];steps+=1
    return K.from_dict({a:c for a,c in rem.items()if c}),steps
   nt,steps1=reduce_observable(rows[-1]);nt2,steps2=reduce_observable(rows[-1]**2);assert nt and not nt2
   stats=metadata.read_text().splitlines();assert stats[1]=='1'and stats[2]!='0'
   # Check the reduced exact Hilbert rational function, not just its initial coefficients.
   z=s.Symbol('T');hil=s.sympify(stats[3].replace('^','**'),locals={'T':z});expected=sum(v*z**i for i,v in enumerate(NUMS[k]))/(1-z)
   assert s.cancel(hil-expected)==0
   hff=[sum(NUMS[k][:n+1])for n in range(13)];assert list(map(int,stats[4].strip('()').split(',')))==hff
   assert sum(NUMS[k])in[49,12,6]
   light_cut_hil=s.sympify(stats[6].replace('^','**'),locals={'T':z})
   if k<2:assert stats[5]=='0'and s.cancel(light_cut_hil-sum(v*z**i for i,v in enumerate(NUMS[0])))==0
   else:assert stats[5]=='1'and s.cancel(light_cut_hil-(1+2*z+6*z**2+10*z**3+11*z**4-3*z**5-17*z**6-9*z**7)/(1-z**2))==0
   dr=json.loads((ROOT/'reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json').read_text());us=s.symbols('u1:62');pt=s.zeros(61,1);pt[4]=pt[48]=1;pt[12]=core;pt[37]=-1
   C=s.Matrix(dr['actual_C_tangent']).applyfunc(s.sympify).subs(dict(zip(us,pt)));assert C[[19,26],:]==s.Matrix([[0,0,1],[-2,0,0]])
   cc=C*C[[19,26],:].nullspace()[0];cr=s.Matrix([cc[i]for i in meta['chart']]+[0]);pivs=list(T.T.to_DM().rref()[1]);cpar=T[pivs,:].inv()*cr[pivs,:];assert T*cpar==cr
   light=next(i for i,v in enumerate(cpar)if v and WEIGHTS[i]==1);assert all(not v or WEIGHTS[i]==1 for i,v in enumerate(cpar));assert cpar[0]+2*cpar[1]+3*cpar[4]==-5
   results.append(dict(core=str(core),frame=meta,standard_basis_output_sha256=hashlib.sha256(member.read_bytes()).hexdigest(),independent_sparse_QQ_reduction_steps=[steps1,steps2],independent_t_nonzero_t_square_zero=True,homogenizer_scale=scales[-1],t_nonzero_normal_form=stats[2],groebner_count=int(stats[0]),slice_dimension=1,weighted_Hilbert_numerator=NUMS[k],weighted_Hilbert_denominator='1-T',eventual_Hilbert_function=sum(NUMS[k]),eventual_actual_D_kernel_Hilbert_function=(sum(NUMS[k])-5 if k else None),h_parameter_test="h=x0+2*x1+3*x4",h_derivative_on_actual_C="-5",h_cut_dimension=int(stats[5]),h_cut_Hilbert_series=stats[6],S_over_h_rank=(sum(NUMS[k])if k<2 else None),actual_D_kernel_over_h_rank=(sum(NUMS[k])-5 if k<2 else None),actual_D_kernel_over_h_torsion_minimal_generators=([0,37,None][k]),actual_D_cut_regular_light_index=light,actual_D_cut_regular_light_derivative=str(cpar[light]),t_nonzero=True,t_squared_zero=True,full_fixed_beta_dimension=3,two_fixed_beta_orbit_minor='2'))
   print('POINT',core,'REDUCTIONS',[steps1,steps2],'HF',sum(NUMS[k]),'SECONDS',time.monotonic()-start,flush=True)
 paths=[str(Path(__file__).relative_to(ROOT)),NOTE,BASE,'code/s61_flag_linear_artifact_weighted_engine_20261102.py',cpp,*base['input_sha256']]
 hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in dict.fromkeys(paths)}
 return dict(id=STEM,all_passed=True,field='QQ',points=results,all_fixed_beta_zero_nuisance_only=True,actual_kernel_finite_length=False,relative_t_nilpotence_claimed=False,general_beta_removal_claimed=False,full_S61_exhaustion_claimed=False,proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();result=verify()
 with a.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=result['seconds'],hashes=len(result['input_sha256']))))
