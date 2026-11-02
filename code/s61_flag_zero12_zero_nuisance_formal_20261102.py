"""Bounded exact weighted finite model at the zero-nuisance flag point."""
import argparse,hashlib,json,subprocess,tempfile,time,os
from pathlib import Path
from fractions import Fraction
from itertools import combinations_with_replacement as cwr
import sympy as s
from s61_model import S61,BASIS,X
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import sp_basis,COORDS,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_20261102.md'

def verify():
 root=ROOT;started=time.monotonic();m=S61()
 one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1]
 initial=json.loads((root/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json').read_text());point=s.Matrix(initial['point61']).applyfunc(s.Rational)
 for i in one:point[i]=0
 point[4]=point[12]=1
 local=calculate(m,True,point);J=s.Matrix(local['jacobian']).applyfunc(s.Rational);Tfull=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);T=Tfull[:,list(range(2,9))]
 chart=local['chart_columns'];fixed=next(i for i in range(45)if i not in chart)
 torus=[0,1,-1,0,0];raw=[]
 for index in m.certificate['raw_free_indices']:
  pair,mon=BASIS[index];ex=s.Poly(mon,*X).monoms()[0];raw.append(sum(torus[i]for i in pair)-sum(a*b for a,b in zip(torus,ex)))
 weights=[a-b for a,b in zip(m.reduced['parameter_weights'],raw)];two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2]
 qw=[5-weights[two[j]]for j in chart]+[5];nw=[1,1,2,2,1,2,2]
 assert max(qw)==6
 assert all(T[i,j]==0 or qw[i]==nw[j]for i in range(45)for j in range(7))
 _,eq=specialize(m,dict(zip(one,[point[i]for i in one])),True);terms=[];rowweight=[]
 for f in eq:
  ws={sum(5-weights[two[i]]for i in mon)+5*(2-len(mon))for mon,c in f.items()if c};assert len(ws)<=1;rowweight.append(next(iter(ws))if ws else None)
  row=[]
  for mon,c in f.items():
   if fixed not in mon:
    a,b=[chart.index(i)for i in mon]+[44]*(2-len(mon));row.append((a,b,c))
  terms.append(row)
 assert max(z for z in rowweight if z is not None)==10
 piv=list(J.to_DM().rref()[1]);selected=list(J[:,piv].T.to_DM().rref()[1]);U=s.zeros(45,36);inverse=J[selected,piv].inv(method='DM')
 for i,col in enumerate(piv):U[col,:]=inverse[i,:]
 assert J*T==s.zeros(135,7)and J[selected,:]*U==s.eye(36)
 assert all(U[i,j]==0 or qw[i]==rowweight[selected[j]]for i in range(45)for j in range(36))
 scales=[s.ilcm(*[s.denom(c)for c in list(J[i,:])+[c for a,b,c in terms[i]]])for i in range(135)]
 Jn=s.diag(*scales)*J;Un=U*s.diag(*[1/s.Integer(scales[i])for i in selected]);den=s.ilcm(*[s.denom(c)for c in Un]);Un=(den*Un).applyfunc(s.cancel);td=s.ilcm(*[s.denom(c)for c in T]);Tn=td*T
 ms={d:[a for a in cwr(range(7),d)if sum(nw[j]for j in a)<=10]for d in range(11)};ix={d:{a:i for i,a in enumerate(v)}for d,v in ms.items()}
 mons=sum((ms[d]for d in range(11)),[]);assert len(mons)==2892
 with tempfile.TemporaryDirectory(prefix='zero12-weighted-finite-')as tmp:
  temp=Path(tmp);cpp=ROOT/'code/s61_flag_zero12_zero_nuisance_raw_engine_20261102.cpp';exe=temp/'raw';subprocess.run(['c++','-w','-O3','-std=c++17','-I/opt/homebrew/include',str(cpp),'-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(exe)],check=True)
  inp=temp/'raw.txt'
  with inp.open('w')as out:
   out.write(f'135 7 10 36\n{den} {td}\n'+' '.join(str(len(ms[d]))for d in range(11))+'\n'+' '.join(map(str,selected))+'\n')
   for M in [Jn,Tn,Un]:
    for row in M.tolist():out.write(' '.join(map(str,row))+'\n')
   for i,row in enumerate(terms):out.write(str(len(row))+' '+' '.join(f'{a} {b} {c*scales[i]}'for a,b,c in row)+'\n')
   for a in range(1,10):
    for b in range(1,11-a):out.write(' '.join(str(ix[a+b].get(tuple(sorted(x+y)),-1))for x in ms[a]for y in ms[b])+'\n')
   out.write(' '.join(map(str,qw))+'\n')
   for d in range(1,11):out.write(' '.join(str(sum(nw[j]for j in a))for a in ms[d])+'\n')
  rawout=temp/'rawout.txt';coordinate_file=temp/'coords.txt';r=subprocess.run([str(exe),str(inp),str(rawout),str(coordinate_file)],capture_output=True,text=True,check=True,timeout=45);print(r.stderr,flush=True)
  rows=[list(map(Fraction,line.split()))for line in rawout.read_text().splitlines()];assert len(rows)==136
  for i,row in enumerate(rows):
   for mon,c in zip(mons[8:],row):
    if c:assert sum(nw[j]for j in mon)==(rowweight[i]if i<135 else 5)
  expr=[]
  for row in rows:
   expr.append('+'.join('('+str(c)+')*'+'*'.join('x_'+str(j)for j in mon)for mon,c in zip(mons[8:],row)if c)or'0')

  # Reconstruct the whole polynomial substitution independently of the native residual recursion.
  K,*V=ring(','.join('x_'+str(i)for i in range(7)),QQ)
  def monomial(a):
   z=K.one
   for j in a:z*=V[j]
   return z
  monpolys=[monomial(a)for a in mons]
  coordlines=coordinate_file.read_text().splitlines();assert coordlines[0]=='45';coordinates=[]
  for row in coordlines[1:]:
   f=row.split();assert len(f)==1+2*int(f[0]);coordinates.append(sum((K.ground_new(QQ.from_sympy(s.Rational(c)))*monpolys[int(i)]for i,c in zip(f[1::2],f[2::2])),K.zero))
  residuals=[sum((K.ground_new(QQ.from_sympy(s.Rational(c.numerator,c.denominator)))*q for q,c in zip(monpolys[8:],row)if c),K.zero)for row in rows]
  for i in range(135):
   original=sum((K.ground_new(QQ.from_sympy(J[i,j]))*coordinates[j]for j in range(45)if J[i,j]),K.zero)
   original+=sum((K.ground_new(QQ.from_sympy(c))*coordinates[a]*coordinates[b]for a,b,c in terms[i]),K.zero)
   assert scales[i]*original==residuals[i]
  assert coordinates[44]==residuals[-1]
  # A finite polynomial combination certifies t^2, in addition to the Groebner normal-form test.
  membership=temp/'membership.txt';metadata=temp/'metadata.txt'
  program='R=QQ[x_0..x_6,Degrees=>{1,1,2,2,1,2,2}];F=matrix{{'+','.join(expr[:135])+'}};t='+expr[-1]+';G=gb(F,ChangeMatrix=>true);C=t^2//G;assert((F*C)_(0,0)==t^2);assert(numRows C==135);assert(t%G!=0);out=openOut "'+str(membership)+'";out<<numRows C<<endl;scan(numRows C,i->{a=listForm C_(i,0);out<<#a;scan(a,z->{out<<" "<<toString(last z);scan(first z,e->out<<" "<<e);});out<<endl;});close out;out=openOut "'+str(metadata)+'";out<<numColumns gens G<<endl;out<<dim coker F<<endl;out<<toString(t%G)<<endl;close out;'
  script=temp/'membership.m2';script.write_text(program);env=os.environ|{'DYLD_LIBRARY_PATH':str(root/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'}
  rr=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(script)],capture_output=True,text=True,timeout=60,env=env)
  assert rr.returncode==0,(rr.stdout,rr.stderr)
  text=membership.read_text().splitlines();assert text[0]=='135';cof=[]
  for row in text[1:]:
   f=row.split();assert len(f)==1+8*int(f[0]);cof.append(K.from_dict({tuple(map(int,f[2+8*i:9+8*i])):QQ.from_sympy(s.Rational(f[1+8*i]))for i in range(int(f[0]))}))
  assert sum((a*b for a,b in zip(cof,residuals[:135])),K.zero)==residuals[-1]**2
  meta=metadata.read_text().splitlines();assert meta[:2]==['122','1'] and meta[2]!='0'
  print('EXACT WEIGHTED MEMBERSHIP PASS',time.monotonic()-started,flush=True)
  # Independent finite Lie action for the two genuine fixed-beta group directions.
  z,bu,beta=primitive_beta(m);vec=s.Matrix([beta[i,j].diff(z[k])for(i,j),k in COORDS]);basis=vec.jacobian(bu);br=list(basis.T.to_DM().rref()[1]);left=basis[br,:].inv();sp=sp_basis()
  bzero=dict(zip(bu,[point[i]for i in one]));action=s.Matrix.hstack(*[left*infinitesimal_action(vec,a)[br,:]for a in sp]).subs(bzero)
  assert action.rank()==8 and action[:,4]==s.zeros(16,1)and action[:,9]==action[:,7]
  generators=[sp[4],sp[9]-sp[7]];assert generators[0]*generators[1]==generators[1]*generators[0]
  qpoint=s.zeros(61,1);qpoint[37]=-1
  from s61_model import J as J5
  Q0=m.reconstruct(list(qpoint))-X[0]**2*J5;orbit=[]
  for a in generators:
   G=s.zeros(5);G[1:,1:]=a;Ax=G*s.Matrix(X);dQ=G*Q0+Q0*G.T-sum((Ax[i]*Q0.diff(X[i])for i in range(5)),s.zeros(5));normalized=m.normalize(X[0]**2*J5+dQ);q=s.Matrix([normalized['parameters'][i]for i in two]);q+=q[25]*s.Matrix([(-1 if j==25 else 0)for j in range(45)]);orbit.append(q)
  free=[19,26,27,28,33,34,35,37,38];orbit=s.Matrix.hstack(*orbit)[free,:]
  assert orbit[:2,:].det()==2 and orbit[2:,:]==s.zeros(7,2)
  # The full seven-nuisance source has a positive contracting action; this is not a relative nilpotence proof.
  indices=[i for i in one if i!=3];nuisance_b=[4,7,8,9,10,11,14];nuisance_u=[indices[j]for j in nuisance_b]
  assert [weights[i]for i in nuisance_u]==[3,3,1,2,2,2,1]and weights[4]==weights[12]==0
  extra=json.loads((root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text());all_equations=m.reduced['equations']+[extra['universal_new_equation']]
  assert len(all_equations)==135
  assert all(len({sum(weights[i]for i in mon)for mon,_,_ in polynomial})==1 for polynomial in all_equations)
  enc=lambda p:[[list(a),str(c)]for a,c in sorted(p.items())]
  report=dict(id='S61_FLAG_ZERO12_ZERO_NUISANCE_FORMAL_20261102',all_passed=True,field='QQ',point61=list(map(str,point)),primitive_stabilizer_dimension=2,stabilizer_coefficients=[[int(i==4)for i in range(10)],[int(i==9)-int(i==7)for i in range(10)]],orbit_tangents=[[str(z)for z in row]for row in orbit.tolist()],orbit_minor='2',raw_coordinate_weights=qw,normal_weights=nw,row_weights=rowweight,implicit_coordinate_degree_bound=6,residual_weight_bound=10,weighted_monomials=2892,all45_coordinate_weights_checked=True,all135_full_polynomial_substitutions_checked=True,all_implicit_coordinates_degree7_and_higher_zero=True,polynomial_residuals=[enc(p)for p in residuals[:135]],homogenizer_polynomial=enc(residuals[-1]),t_squared_cofactors=[enc(p)for p in cof],nonzero_cofactors=sum(bool(p)for p in cof),cofactor_terms=sum(map(len,cof)),t_squared_polynomial_identity_independently_multiplied=True,groebner_generators=122,slice_dimension=1,full_fixed_beta_dimension=3,t_nonzero_normal_form=meta[2],full_fixed_beta_t_nonzero=True,full_fixed_beta_t_squared_zero=True,nuisance_weights=[3,3,1,2,2,2,1],all135_global_weight_homogeneous=True,relative_nuisance_nilpotence_claimed=False,general_beta_removal_claimed=False,proof_note=NOTE,raw_log=r.stderr,seconds=time.monotonic()-started)
  paths=[Path(__file__),root/NOTE,root/'code/s61_flag_zero12_zero_nuisance_raw_engine_20261102.cpp',root/'code/s61_model.py',root/'code/s61_minimal_infinity_local_data_20261101.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'code/weight1_structure.py',root/'code/primitive_jet_elliptic_curve_20261101.py',root/'reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json',root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',root/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
  paths += [m.source/f for f in ['symplectic_slice_certificate.json','symplectic_slice_equations.json','symplectic_slice_reduced_equations.json']]
  report['input_sha256']={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths}
  return report

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();result=verify()
 with args.output.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
