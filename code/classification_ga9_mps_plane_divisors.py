#!/usr/bin/env python3
"""Exact divisor restrictions on finite explicit Lagrangian-plane families.
All sourceperiod/smoothing=1 tensors represent denseopenPGLorbits ofT36-38.
No exhaustive Grassmann assertion is made for omitted nonrational curves.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xs=s.symbols('x0:5');base=s.symbols('A B C');t,u,v=s.symbols('t u v');order=[0,3,4,1,2];records=[];inputs=[Path(__file__)]
 cases={
 'T36':[("012",'main',[0,t,-t*u/5,0,u,-25/t]),("013",'prime1',[-3*t*t/100,-10/t,3*t**4/2000,t,0]),("013",'prime2',[0,t,0,u,(t*u-5)/(5*t**3)]),("014",'prime1',[0,t,0,-5*t*t]),("014",'prime2',[3*t**3,t,15*t**4,10*t*t]),("014",'prime3',[0,0,t,u]),("023",'main',[0,t,u,0]),("034",'axis1',[t,0]),("034",'axis2',[0,t]),("134",'main',[t]),("234",'point',[])],
 'T38':[("012",'main',[0,t,u,0,-2*u/t,-4/t]),("013",'prime1',[-1/(2*t),1/t**2,u,v,t]),("013",'prime2',[0,t,0,u,(t*u-2)/(2*t**3)]),("014",'prime1',[0,t,0,-2*t*t]),("014",'prime2',[0,0,t,u]),("023",'main',[0,t,u,0]),("034",'main',[t,u]),("134",'main',[t]),("234",'point',[])],
 'T37':[("034",'main',[t,-1/(4*t)]),("134",'main',[t]),("234",'point',[]),("013",'prime1',[-(t**5+1)/(2*t**2),t,(t**5+1)**2/(12*t**4),-(t**5+1)/(3*t),-t*t])]
 }
 # The second rationalcurve in cell013, derived directly from the smalloriginalGB.
 a4=t;a2=(3-t**5)/(12*t*t);a1=s.cancel((t**3+12*a2)/3);a0=s.cancel((a1*t*t+1)/(2*t));a3=s.cancel((-6*a0**3+a1*a1*t+6*a1*a2*t+18*a0*a2+t*t)/(3*a0*t));cases['T37'].append(('013','prime2',[a0,a1,a2,a3,t]))
 for name in ['T36','T37','T38']:
  path=ROOT/f'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/{name}.json';inputs.append(path);source=json.loads(path.read_text());env=dict(zip(map(str,xs),xs));env.update({'p0':1,**{'t'+str(i):1 for i in range(4)}});P=s.zeros(5)
  for ij,f in source['homogeneous_bivector'].items():
   i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
  pf=[]
  for missing in range(5):
   i,j,k,l=[a for a in range(5)if a!=missing];pf.append((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
  D=s.expand(sum(xs[i]*pf[i]for i in range(5)))
  for label,branch,params in cases[name]:
   pivot=list(map(int,label));M=s.zeros(5,3);count=0
   for c,r in enumerate(pivot):M[order[r],c]=1
   for r in range(5):
    if r in pivot:continue
    for c,pivotrow in enumerate(pivot):
     if r>pivotrow:M[order[r],c]=params[count];count+=1
   assert count==len(params);normals=[]
   for r in range(5):
    if r in pivot:continue
    n=s.zeros(5,1);n[order[r]]=1
    for c,pivotrow in enumerate(pivot):n[order[pivotrow]]=-M[order[r],c]
    normals.append(n)
   sub=dict(zip(xs,M*s.Matrix(base)));lag=s.Poly((normals[0].T*P*normals[1])[0].subs(sub,simultaneous=True),*base);assert all(s.cancel(a)==0 for a in lag.coeffs()),(name,label,branch,lag)
   restriction=s.cancel(D.subs(sub,simultaneous=True));co,factors=s.factor_list(restriction,*base);data=[dict(factor=str(f),multiplicity=int(m),degree=int(s.Poly(f,*base).total_degree()))for f,m in factors];degrees=sorted([q['degree']for q in data for _ in range(q['multiplicity'])]);assert s.cancel(restriction-co*s.prod(f**m for f,m in factors))==0
   row=dict(component=name,cell=label,branch=branch,parameters=sorted(str(a)for a in set().union(*(a.free_symbols if hasattr(a,'free_symbols')else set()for a in params))),plane_matrix=[[str(a)for a in rr]for rr in M.tolist()],exact_Lagrangian_identity=True,restricted_actual_quintic=str(restriction),factor_coefficient=str(co),factors=data,degrees=degrees,possible_Ga9_irreducible_quartic=degrees==[1,4]);records.append(row);print(name,label,branch,degrees,round(time.monotonic()-started,3),flush=True)
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),explicit_families=records,scope='Exact generic divisor restrictions for the listed explicitLagrangianplane families. A positiveH+irreducquartic candidate is necessaryforuniform2normal scaling withGa9endpoint. No exhaustiveGrassmannor MPSnonmembership claim; T37cell012 nonrationalcurve not included.',possible_Ga9_families=[dict(component=r['component'],cell=r['cell'],branch=r['branch'])for r in records if r['possible_Ga9_irreducible_quartic']],seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_MPS_PLANE_DIVISORS_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS',result['seconds'])

def t37_open_lagrangian_locus(output):
 """Whole previously omitted cell012, fixed actual accepted T37 source."""
 started=time.monotonic();xs=s.symbols('x0:5');ABC=s.symbols('A B C');av=s.symbols('a0:6');sourcepath=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T37.json';source=json.loads(sourcepath.read_text());env=dict(zip(map(str,xs),xs));env.update({'p0':1,**{'t'+str(i):1 for i in range(4)}});P=s.zeros(5)
 for ij,text in source['homogeneous_bivector'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(text,locals=env);P[j,i]=-P[i,j]
 M=s.Matrix([[1,0,0],[*av[:3]],[*av[3:]],[0,1,0],[0,0,1]]);n1=s.Matrix([-av[0],1,0,-av[1],-av[2]]);n2=s.Matrix([-av[3],0,1,-av[4],-av[5]]);sub=dict(zip(xs,M*s.Matrix(ABC)));lag=s.Poly(s.expand((n1.T*P*n2)[0].subs(sub,simultaneous=True)),*ABC);equations=[s.factor(v)for v in lag.coeffs()];assert len(equations)==6
 pf=[]
 for missing in range(5):
  i,j,k,l=[a for a in range(5)if a!=missing];pf.append((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
 D=s.expand(sum(xs[i]*pf[i]for i in range(5)));restriction=s.expand(D.subs(sub,simultaneous=True));paths=[Path(__file__),sourcepath,ROOT/'reports/CLASSIFICATION_GA9_MPS_PLANE_DIVISORS_20261107.json']
 report=dict(task_id='E_reduced_nonnormal',phase='EXACT_WHOLE_CELL_INPUT_SAVED',all_passed=False,source='T37 allactualperiods/smoothing=1 accepted source point',source_period_inversions=['p0=1 nonzero'],whole_cell_plane_matrix=[[str(v)for v in M.row(i)]for i in range(5)],whole_cell_6_Lagrangian_equations=list(map(str,equations)),restricted_actual_quintic=str(restriction),actual_full_quintic=str(D),coefficient_ring='QQ[a0,...,a5]',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_bytes=Path(__file__).read_text(),scope='One entire Schubertcell012 for one actualknownT37source, not allperiods/PGLparents; no point sampling or higherjets. Only divisor/marked recognition entrance, not Poisson-source containment until actualregularfamily/allmag matched.')
 assert not output.exists();output.write_text(json.dumps(report,indent=2)+'\n');print('exact cell input saved',flush=True)
 G=s.groebner(equations,*av,order='lex');basis=[s.factor(v.as_expr())for v in G.polys];assert all(G.reduce(v)[1]==0 for v in equations)
 coeff=s.Poly(restriction,*ABC);reduced=s.expand(sum(s.factor(G.reduce(c)[1])*s.prod(q**k for q,k in zip(ABC,ex))for ex,c in coeff.terms()));report.update(phase='COMPLETE_WHOLE_CELL_GROEBNER',all_passed=True,lex_Groebner_basis=list(map(str,basis)),Groebner_basis_length=len(basis),all6_original_equations_reduce_zero=True,restricted_quintic_mod_whole_cell_ideal=str(reduced),seconds=time.monotonic()-started);output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,basis_length=len(basis),basis=list(map(str,basis)),seconds=report['seconds'])),flush=True)


def t37_nonrational_factor(output):
 """Exact actualD factorization on the entire previously omitted curve."""
 import subprocess,resource
 started=time.monotonic();path=ROOT/'reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_OPEN_LOCUS_20261107.json';data=json.loads(path.read_text());av=s.symbols('a0:6');A,B,C=s.symbols('A B C');env=dict(zip(map(str,av),av));env.update(A=A,B=B,C=C);gb=[s.sympify(v,locals=env)for v in data['lex_Groebner_basis']]
 factors=s.factor_list(gb[-1],*av)[1];curve=next(f for f,power in factors if s.Poly(f,*av).total_degree()>2);assert s.factor_list(curve,*av)[1]==[(curve,1)]
 a3=s.cancel(s.solve(gb[4],av[3])[0]);bindings={av[3]:a3}
 for index in [0,1,2]:bindings[av[index]]=s.cancel(s.solve(gb[index],av[index])[0].subs(bindings))
 den=s.factor(s.lcm([s.denom(v)for v in bindings.values()]));assert s.rem(den,curve,av[5])!=0
 # All original equations are checked over the exact curve function field.
 identities=[]
 for text in data['whole_cell_6_Lagrangian_equations']:
  num=s.fraction(s.cancel(s.sympify(text,locals=env).subs(bindings,simultaneous=True)))[0];remainder=s.rem(num,curve,av[5]);assert s.expand(remainder)==0;identities.append(True)
 def native(expr):return str(expr).replace('**','^')
 script='ring R=(0,a4,a5),(A,B,C),dp;\nminpoly='+native(curve)+';\n'
 for index in [3,0,1,2]:script+='number a'+str(index)+'='+native(bindings[av[index]])+';\n'
 for i,text in enumerate(data['whole_cell_6_Lagrangian_equations']):script+='number test'+str(i)+'='+native(s.sympify(text,locals=env))+';\nprint("LAG'+str(i)+'="+string(test'+str(i)+'));\n'
 script+='poly f='+native(s.sympify(data['restricted_actual_quintic'],locals=env))+';\nprint("ACTUAL_D=");print(f);\nlist F=factorize(f);\nprint("FACTOR_LIST=");F;\nfor(int j=1;j<=size(F[1]);j++){print("FACTOR_"+string(j)+"_DEG="+string(deg(F[1][j])));print(F[1][j]);print("MULT="+string(F[2][j]));}\nquit;\n'
 scriptpath=ROOT/'research/classification/inputs/t37_nonrational_actual_D_factor_20261107.sing';assert not scriptpath.exists();scriptpath.write_text(script)
 run=subprocess.run(['bash',str(ROOT/'work/nonsquare-v62/singular-runtime/run-singular-optimized-v6.sh'),'-q',str(scriptpath)],cwd=ROOT,text=True,capture_output=True,timeout=90)
 result=dict(task_id='E_reduced_nonnormal',phase='NATIVE_RESULT_SAVED',all_passed=False,input_report=str(path.relative_to(ROOT)),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),nonrational_curve=str(curve),generic_parameter_bindings={str(k):str(v)for k,v in bindings.items()},coefficient_denominator=str(den),all6_Lagrangian_identities_mod_curve_zero=identities,native_script=str(scriptpath.relative_to(ROOT)),native_script_sha256=hashlib.sha256(scriptpath.read_bytes()).hexdigest(),native_returncode=run.returncode,native_stdout=run.stdout,native_stderr=run.stderr,producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='Exact generic functionfield factorization on whole omittedcurve, without newpoint. If a properfactorpartition identity is obtained, its closedmultiplicationcondition extends to every curvepoint. Any actualpositive membership still needs fullPoisson/PGL/source bridge.')
 assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n');print(run.stdout);print(run.stderr);print('saved',result['seconds'])


def t37_curve_double_conic(output):
 """One whole-curve identity via exact function-field gcd and linear lift."""
 from sympy.polys.agca.extensions import MonogenicFiniteExtension
 import resource
 started=time.monotonic();path=ROOT/'reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_OPEN_LOCUS_20261107.json';data=json.loads(path.read_text());av=s.symbols('a0:6');A,B,C,T=s.symbols('A B C T');env=dict(zip(map(str,av),av));env.update(A=A,B=B,C=C);gb=[s.sympify(v,locals=env)for v in data['lex_Groebner_basis']];curve=next(f for f,power in s.factor_list(gb[-1],*av)[1]if s.Poly(f,*av).total_degree()>2)
 K=MonogenicFiniteExtension(s.Poly(curve,av[5],domain=s.QQ.frac_field(av[4])));zero,one=K.zero,K.one;vals={av[4]:K.convert(av[4]),av[5]:K.generator}
 def evaluate(expr):
  poly=s.Poly(expr,*av);total=zero
  for ex,c in poly.terms():
   term=K.convert(c)
   for symbol,n in zip(av,ex):
    if n:term*=vals[symbol]**n
   total+=term
  return total
 for index,rowindex in [(3,4),(0,0),(1,1),(2,2)]:
  expression=s.Poly(gb[rowindex],av[index]);assert expression.degree()==1;vals[av[index]]=-evaluate(expression.nth(0))/evaluate(expression.nth(1))
 assert all(evaluate(s.sympify(text,locals=env))==zero for text in data['whole_cell_6_Lagrangian_equations'])
 D={ex:evaluate(c)for ex,c in s.Poly(s.sympify(data['restricted_actual_quintic'],locals=env),A,B,C).terms()};D={ex:c for ex,c in D.items()if c}
 def trim(f):
  while f and not f[-1]:f.pop()
  return f
 def divrem(f,g):
  f=f.copy();q=[zero]*max(0,len(f)-len(g)+1)
  while f and len(f)>=len(g):
   j=len(f)-len(g);c=f[-1]/g[-1];q[j]=c
   for k,d in enumerate(g):f[j+k]-=c*d
   trim(f)
  return trim(q),f
 def pmul(f,g):
  out=[zero]*(len(f)+len(g)-1)
  for i,c in enumerate(f):
   for j,d in enumerate(g):out[i+j]+=c*d
  return trim(out)
 f=[D.get((0,b,5-b),zero)for b in range(6)];fp=[K.convert(i)*f[i]for i in range(1,len(f))];aa,bb=f,fp
 while bb:aa,bb=bb,divrem(aa,bb)[1]
 g=[c/aa[-1]for c in aa];degree=len(g)-1;report=dict(task_id='E_reduced_nonnormal',phase='EXACT_BINARY_GCD',all_passed=False,curve=str(curve),function_field=str(K),all6_Lagrangian_identities_zero=True,binary_restriction_gcd_degree=degree,input_report_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),scope='One entirefixedT37 omittedcurve exactfactoridentity. No newpoint/GB/jet. Properclosedfactorlocus conclusion onlyafterfullternarymultiplicationidentity.')
 assert not output.exists()
 if degree!=2:
  report.update(status='NO_QUADRATIC_BINARY_GCD_ROUTE_STOP',seconds=time.monotonic()-started);output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['status','binary_restriction_gcd_degree','seconds']}));return
 l,remainder=divrem(f,pmul(g,g));assert not remainder and len(l)==2
 def clean(P):return{ex:c for ex,c in P.items()if c}
 def add(P,Q):
  R=P.copy()
  for ex,c in Q.items():R[ex]=R.get(ex,zero)+c
  return clean(R)
 def mul(P,Q):
  R={}
  for e,c in P.items():
   for f,d in Q.items():
    ex=tuple(a+b for a,b in zip(e,f));R[ex]=R.get(ex,zero)+c*d
  return clean(R)
 def mon(ex,c=one):return {ex:c}
 def dictpoly(poly):return {ex[0]:c for ex,c in poly.rep.to_dict().items()}
 gd={j:c for j,c in enumerate(g)if c};ld={j:c for j,c in enumerate(l)if c};Q0={(0,b,2-b):c for b,c in gd.items()};L0={(0,b,1-b):c for b,c in ld.items()};q02=mul(Q0,Q0);lq=mul(L0,Q0)
 columns=[mul(mon((1,0,0)),q02),mul(mon((1,1,0),K.convert(2)),lq),mul(mon((1,0,1),K.convert(2)),lq)];rows=[(1,b,4-b)for b in range(4,-1,-1)];matrix=[[col.get(ex,zero)for col in columns]+[D.get(ex,zero)]for ex in rows];pivotrows=[];at=0
 for col in range(3):
  k=next(j for j in range(at,len(matrix))if matrix[j][col]);matrix[k],matrix[at]=matrix[at],matrix[k];pivotrows.append(k);pivot=matrix[at][col];matrix[at]=[v/pivot for v in matrix[at]]
  for j in range(len(matrix)):
   if j!=at:
    scalar=matrix[j][col];matrix[j]=[u-scalar*v for u,v in zip(matrix[j],matrix[at])]
  at+=1
 assert all(not any(row)for row in matrix[3:]);lA,qAB,qAC=[matrix[j][3]for j in range(3)];L=add(L0,mon((1,0,0),lA));Qp=add(Q0,{(1,1,0):qAB,(1,0,1):qAC});known=mul(L,mul(Qp,Qp));column=mul(mon((2,0,0),K.convert(2)),lq);rows2=[(2,b,3-b)for b in range(3,-1,-1)];ex=next(ex for ex in rows2 if column.get(ex,zero));qAA=(D.get(ex,zero)-known.get(ex,zero))/column[ex];Q=add(Qp,mon((2,0,0),qAA));product=mul(L,mul(Q,Q));assert add(D,{ex:-c for ex,c in product.items()})=={}
 def enc(P):return {','.join(map(str,ex)):str(K.to_sympy(c))for ex,c in P.items()}
 report.update(phase='COMPLETE_WHOLE_CURVE_IDENTITY',all_passed=True,status='WHOLE_CURVE_D_EQUALS_LINEAR_TIMES_CONIC_SQUARED',all21_ternary_degree5_coefficients_exactly_equal=True,linear_factor_coefficients=enc(L),quadratic_factor_coefficients=enc(Q),generic_plane_bindings={str(k):str(K.to_sympy(v))for k,v in vals.items()},binary_gcd_coefficients={str(k):str(K.to_sympy(v))for k,v in gd.items()},linear_lift_pivot_rows=pivotrows,proper_closed_factor_partition=[1,2,2],no_E6_A5_reduced_H_irreducible_quartic_normalcone_on_whole_curve=True,closure_logic='Allfactorcoefficients are functionfieldelements overQQ(a4)[a5]/B. Clearingfinitelymanydenominators givesa denseopenon everygeometriccomponent (QQirreducibility/characteristiczero); imagebelongsproperprojectivemultiplicationlocusLQ², whichisclosed, soentirecurveimagebelongsit.',seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','status','seconds','max_rss_kib']}))

def a5_smooth_hyperflex_graph(output):
 from s61_model import S61,X,pushforward
 import ast,resource,itertools as it
 h,x,y,z,w=coords=s.symbols('h x y z w')
 started=time.monotonic();root=Path(__file__).resolve().parents[1]
 sourcepath=Path('reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_SOURCES_20261107.json');normalpath=Path('reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_NORMALFORMS_20261107.json')
 source=next(v for v in json.loads(sourcepath.read_text())['cases']if v['name']=='A5_smooth_hyperflex');normal=next(v for v in json.loads(normalpath.read_text())['cases']if v['name']==source['name'])
 values={s.sympify(k):s.Rational(v)for k,v in normal['canonical_parameter_values'].items()};P=s.zeros(5)
 for ij,text in source['all150_source_brackets'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(text).subs(values);P[j,i]=-P[i,j]
 P=P.subs(dict(zip(coords,X)),simultaneous=True);F=s.sympify(source['transformed_irreducible_quartic']);assert F.subs({h:1,y:0,w:1})==1
 shift=s.eye(5);shift[4,0]=-1;translated=pushforward(P,shift);K=translated.subs(dict(zip(X,[1,0,0,0,0])))[1:,1:]
 pf=s.factor(K[0,1]*K[2,3]-K[0,2]*K[1,3]+K[0,3]*K[1,2]);assert pf!=0 and not pf.free_symbols
 pivot=next((i,j)for i,j in it.combinations(range(4),2)if K[i,j]!=0);i,j=pivot;unit=[s.eye(4)[:,k]for k in range(4)];r0=unit[i];r1=unit[j]/K[i,j];other=[]
 for k in range(4):
  if k in pivot:continue
  rr=unit[k];other.append(rr-(rr.T*K*r1)[0]*r0+(rr.T*K*r0)[0]*r1)
 r2,r3=other;r3=r3/s.cancel((r2.T*K*r3)[0]);frame=s.eye(5);frame[1:,1:]=s.Matrix.vstack(*(rr.T for rr in [r0,r1,r2,r3]));J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1
 assert (frame[1:,1:]*K*frame[1:,1:].T-J).applyfunc(s.cancel)==s.zeros(4)
 G=frame*shift;assert G.det()!=0;GI=G.inv();assert G*GI==s.eye(5)
 model=S61();framed=pushforward(P,G);normalized=model.normalize(framed);verification=model.verify_point(normalized,check_equations=True);point=normalized['parameters'];us=s.symbols('u1:62');results={}
 for name,path in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
  packet=json.loads(Path(path).read_text());indices=packet['weight2_u_indices_zero_based'];residuals=[s.factor(s.cancel(point[k]-s.sympify(expression).subs(dict(zip(us,point)),simultaneous=True)))for k,expression in zip(indices,packet['quadratic_section_q'])]
  assert len(residuals)==45
  results[name]=dict(all45_graph_residuals=list(map(str,residuals)),canonical_point_is_in_graph=all(v==0 for v in residuals),nonzero_residuals={str(k):str(v)for k,v in zip(indices,residuals)if v!=0})
 paths=[Path(__file__),sourcepath,normalpath,Path('code/s61_model.py'),Path('reports/Q51_UNIVERSAL_SECTION_20261101.json'),Path('reports/Q52_UNIVERSAL_SECTION_20261101.json'),*[Path('legacy/v02/data')/name for name in ['symplectic_slice_certificate.json','symplectic_slice_equations.json','symplectic_slice_reduced_equations.json']]]
 old=json.loads(Path('research/classification/jobs/plane_hyperplane_quartic_a5_hyperflex_graph.json').read_text())['old_function_AST_sha256'];current={n.name:hashlib.sha256(ast.dump(n,include_attributes=False).encode()).hexdigest()for n in ast.parse(Path(__file__).read_text()).body if isinstance(n,ast.FunctionDef)};assert all(current[k]==v for k,v in old.items())
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,source_case=source['name'],canonical_magnetic_values=normal['canonical_parameter_values'],structural_rank4_point=[1,0,0,0,1],rank4_point_Pf=str(pf),constant_value_matrix=[[str(v)for v in K.row(i)]for i in range(4)],actual_shift=[[str(v)for v in shift.row(i)]for i in range(5)],Darboux_frame=[[str(v)for v in frame.row(i)]for i in range(5)],full_PGL_matrix=[[str(v)for v in G.row(i)]for i in range(5)],full_inverse_PGL_matrix=[[str(v)for v in GI.row(i)]for i in range(5)],PGL_determinant=str(G.det()),all61_parameters=list(map(str,point)),perspective=list(map(str,normalized['perspective'])),perspective_matrix=[[str(v)for v in normalized['matrix'].row(i)]for i in range(5)],faithful150_full134_verification=verification,exact_graph_membership=results,whole_original_six_source_membership_if_positive=any(v['canonical_point_is_in_graph']for v in results.values()),old_all_FunctionDef_AST_preserved=True,input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='Exactly one structurally selected rank4 chart and its faithful actual PGL/S61 bridge. Positive graph membership propagates through the accepted full6mag canonicalorbit closure. The new global evaluation-fiber lemma promotes negative residuals to fullT39/T40parent exclusion once independently accepted; unrelated components remain unassigned. No pointfarm, no higherjets.')
 assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,membership={k:v['canonical_point_is_in_graph']for k,v in results.items()},seconds=result['seconds'],max_rss_kib=result['max_rss_kib'])))



def t37_curve_patch_cover(output):
 """Whole lifted curve coverage by the two existing a3-linear GB rows."""
 begun=time.monotonic();path=ROOT/'reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_OPEN_LOCUS_20261107.json';data=json.loads(path.read_text());av=s.symbols('a0:6');env=dict(zip(map(str,av),av));gb=[s.sympify(v,locals=env)for v in data['lex_Groebner_basis']];curve=next(f for f,power in s.factor_list(gb[-1],*av)[1]if s.Poly(f,*av).total_degree()>2);u,v=av[4:];rows=[s.Poly(gb[k],av[3])for k in [4,5]];c4,c5=[row.nth(1)for row in rows];assert not c5.has(u);phi=s.solve(c4,u)[0];E=s.expand(curve.subs(u,phi));A,B,g=s.gcdex(E,c5,v);g=s.factor(g);passed=g==1;result=dict(task_id='E_reduced_nonnormal',phase='COMPLETE_UNIVARIATE_PATCH_CHECK',all_passed=True,whole_B_covered_by_two_existing_linear_rows=passed,curve=str(curve),a3_linear_row_indices=[4,5],row4_coefficient=str(c4),row5_coefficient=str(c5),row4_zero_substitution=str(phi),curve_on_row4_zero=str(E),univariate_gcd=str(g),producer_bytes=Path(__file__).read_text(),input_report_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.monotonic()-begun,scope='No newGB/point. Checks whether thetwoexactexistinglinearrows covereveryactualplane aboveB; ifnot, retainsfinitecommonzero locus anddoesnotclaimfullcell.' )
 if passed:
  Q,remainder=s.div(s.expand(curve-E),c4,u);assert s.expand(remainder)==0;U=s.factor(A);V=s.factor(-A*Q);W=s.factor(B);assert s.expand(U*curve+V*c4+W*c5-1)==0
  result.update(unit_identity_coefficients_for_B_c4_c5=list(map(str,[U,V,W])),exact_unit_identity_verified=True,all_actual_lifted_planes_on_B_are_covered=True,no_extra_isolated_or_vertical_lifts=True,closure_conclusion='On eachU4/U5open a3 isuniquely determined; a0,a1,a2arepolynomiallydeterminedbyfirst3GBrows. Theseopens coverB byunitidentity. Exactfieldidentitiesyieldtheactualcurvebridge; allparametercomplements areincluded. ThealreadyprovedproperclosedLQ² condition extendsovereveryliftedBpoint.')
 assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['all_passed','whole_B_covered_by_two_existing_linear_rows','univariate_gcd','seconds']}))


def t37_exceptional_linear_symmetry(output):
 """ONE whole finite-algebra support proof, never a list of P4 points."""
 from sympy.polys.matrices import DomainMatrix
 import resource
 begun=time.monotonic();path=ROOT/'reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_OPEN_LOCUS_20261107.json';data=json.loads(path.read_text());patch=ROOT/'reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_CURVE_PATCH_COVER_20261107.json';pd=json.loads(patch.read_text());av=s.symbols('a0:6');A,B,C=s.symbols('A B C');aa,v=av[3],av[5];env=dict(zip(map(str,av),av));env.update(A=A,B=B,C=C);gb=[s.sympify(t,locals=env)for t in data['lex_Groebner_basis']];rr=s.Poly(s.sympify(pd['univariate_gcd'],locals=env),v).monic().as_expr();phi=s.sympify(pd['row4_zero_substitution'],locals=env);bindings={av[4]:phi}
 for i in range(3):bindings[av[i]]=s.expand(s.solve(gb[i],av[i])[0].subs(bindings,simultaneous=True))
 rawq=s.Poly(s.expand(gb[3].subs(bindings,simultaneous=True)),aa);lc=rawq.nth(2);inv=s.invert(lc,rr,v);qq=aa**2+sum(s.rem(s.expand(inv*rawq.nth(j)),rr,v)*aa**j for j in range(2));basis=[aa**j*v**i for j in range(2)for i in range(10)]
 def normal(expr):
  z=s.Poly(s.rem(s.expand(expr),qq,aa),aa);out=s.zeros(20,1)
  for j in range(2):
   q=s.Poly(s.rem(z.nth(j),rr,v),v)
   for i in range(10):out[10*j+i]=q.nth(i)
  return out
 Mv=s.Matrix.hstack(*[normal(b*v)for b in basis]);Ma=s.Matrix.hstack(*[normal(b*aa)for b in basis]);assert Mv*Ma==Ma*Mv;mult=[Mv**i for i in range(10)]+[Ma*Mv**i for i in range(10)];e=s.eye(20)[:,0]
 def multiply_matrix(c):return sum([q*M for q,M in zip(c,mult)],s.zeros(20))
 remaining=[normal(s.sympify(t,locals=env).subs(bindings,simultaneous=True))for t in data['whole_cell_6_Lagrangian_equations']]+[normal(g.subs(bindings,simultaneous=True))for g in gb]
 ideal=s.Matrix.hstack(*[multiply_matrix(c)for c in remaining]);dom=lambda X:DomainMatrix.from_Matrix(X).convert_to(s.QQ);ic=list(dom(ideal).rref()[1]);I=ideal[:,ic];joint=I;comp=[]
 for k in range(20):
  trial=joint.row_join(s.eye(20)[:,k])
  if len(dom(trial).rref()[1])>joint.cols:joint=trial;comp.append(k)
 assert joint.shape==(20,20);E=s.eye(20)[:,comp];dim=len(comp);projection=dom(joint).inv().to_Matrix()[len(ic):,:];report=dict(task_id='E_reduced_nonnormal',phase='FINITE_SOURCE_ALGEBRA_CONSTRUCTED',all_passed=False,exceptional_polynomial=str(rr),a4_binding=str(phi),a0_a1_a2_bindings={str(k):str(z)for k,z in bindings.items()if k!=av[4]},monic_a3_quadratic=str(qq),all6_original_Lagrangian_equations_imposed=True,full20_algebra_basis=list(map(str,basis)),full_original_ideal_rank=len(ic),actual_finite_quotient_dimension=dim,quotient_standard_basis_indices=comp,source_plane_parameter_map_preserved=True,producer_bytes=Path(__file__).read_text(),input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [path,patch]},scope='Wholefiniteexceptional algebra, not20P4anchors. Weakesttarget: nonzero strictlinearinfinitesimal symmetry at everygeometricpoint excludesfiniteAut markedE6/A5; no blanketproperfactorclaim.')
 assert not output.exists();output.write_text(json.dumps(report,indent=2)+'\n')
 if dim==0:
  report.update(phase='COMPLETE',all_passed=True,all_exceptional_sources_empty=True,seconds=time.monotonic()-begun);output.write_text(json.dumps(report,indent=2)+'\n');print('EMPTY actualfinitealgebra');return
 unit=projection*e;qm=[projection*mult[k]*E for k in comp]
 def qmul(c):return sum([z*M for z,M in zip(c,qm)],s.zeros(dim))
 restricted=s.sympify(data['restricted_actual_quintic'],locals=env).subs(bindings,simultaneous=True);terms={ex:projection*normal(c)for ex,c in s.Poly(s.expand(restricted),A,B,C).terms()};assert terms[(0,0,5)]==unit/3
 exponents=[(i,j,5-i-j)for i in range(5,-1,-1)for j in range(5-i,-1,-1)];polys=[]
 for i in range(3):
  for j in range(3):
   entry={}
   for ex,c in terms.items():
    if not ex[i]:continue
    out=list(ex);out[i]-=1;out[j]+=1;out=tuple(out);entry[out]=entry.get(out,s.zeros(dim,1))+ex[i]*c
   polys.append(entry)
 matrix=s.Matrix.vstack(*[s.Matrix.hstack(*[qmul(poly.get(ex,s.zeros(dim,1)))for poly in polys])for ex in exponents]);kernel=dom(matrix).nullspace().to_Matrix().T;assert matrix*kernel==s.zeros(matrix.rows,kernel.cols)
 cols=[];labels=[]
 for k in range(kernel.cols):
  for j in range(9):
   c=kernel[j*dim:(j+1)*dim,k];M=qmul(c)
   for b in range(dim):cols.append(M[:,b]);labels.append([k,j,b])
 idealK=s.Matrix.hstack(*cols);cc=list(dom(idealK).rref()[1]);covered=len(cc)==dim;report.update(strict_linear_action_QQ_matrix_shape=list(matrix.shape),all9_linear_fields_retained=True,whole_QQ_kernel_dimension=kernel.cols,kernel_basis=[[str(z)for z in kernel[:,k]]for k in range(kernel.cols)],all_strict_linear_derivative_products_zero=True,kernel_coefficient_ideal_rank=len(cc),all_exceptional_geometric_sources_have_nonzero_strict_linear_symmetry=covered,no_E6_A5_target_in_exceptional_source=covered,finite_Aut_target_argument='Nonzero strictlinearfield cannot beEuler sinceD_C5=1/3. Itsconnectedprojectiveflow preservesuniqueH+irreducF4markedpair; E6/A5smoothhyperflexpairs havefiniteAut, soexcluded.',quotient_projection=[[str(z)for z in projection.row(i)]for i in range(dim)])
 if covered:
  selected=idealK[:,cc];weights=dom(selected).inv().to_Matrix()*unit;assert selected*weights==unit
  report.update(unit_ideal_selected_coefficients=[labels[k]for k in cc],unit_identity_weights=list(map(str,weights)),unit_identity_exactly_verified=True)
 report.update(phase='COMPLETE',all_passed=True,seconds=time.monotonic()-begun,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','actual_finite_quotient_dimension','whole_QQ_kernel_dimension','kernel_coefficient_ideal_rank','no_E6_A5_target_in_exceptional_source','seconds','max_rss_kib']}))

def t37_exceptional_algebra_first(output):
 """ONE whole finite-algebra support proof, never a list of P4 points."""
 from sympy.polys.matrices import DomainMatrix
 import resource
 begun=time.monotonic();path=ROOT/'reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_OPEN_LOCUS_20261107.json';data=json.loads(path.read_text());patch=ROOT/'reports/CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_CURVE_PATCH_COVER_20261107.json';pd=json.loads(patch.read_text());av=s.symbols('a0:6');A,B,C=s.symbols('A B C');aa,v=av[3],av[5];env=dict(zip(map(str,av),av));env.update(A=A,B=B,C=C);gb=[s.sympify(t,locals=env)for t in data['lex_Groebner_basis']];rr=s.Poly(s.sympify(pd['univariate_gcd'],locals=env),v).monic().as_expr();phi=s.sympify(pd['row4_zero_substitution'],locals=env);bindings={av[4]:phi}
 for i in range(3):bindings[av[i]]=s.expand(s.solve(gb[i],av[i])[0].subs(bindings,simultaneous=True))
 rawq=s.Poly(s.expand(gb[3].subs(bindings,simultaneous=True)),aa);lc=rawq.nth(2);inv=s.invert(lc,rr,v);qq=aa**2+sum(s.rem(s.expand(inv*rawq.nth(j)),rr,v)*aa**j for j in range(2));basis=[aa**j*v**i for j in range(2)for i in range(10)]
 def normal(expr):
  z=s.Poly(s.rem(s.expand(expr),qq,aa),aa);out=s.zeros(20,1)
  for j in range(2):
   q=s.Poly(s.rem(z.nth(j),rr,v),v)
   for i in range(10):out[10*j+i]=q.nth(i)
  return out
 Mv=s.Matrix.hstack(*[normal(b*v)for b in basis]);Ma=s.Matrix.hstack(*[normal(b*aa)for b in basis]);assert Mv*Ma==Ma*Mv;mult=[Mv**i for i in range(10)]+[Ma*Mv**i for i in range(10)];e=s.eye(20)[:,0]
 def multiply_matrix(c):return sum([q*M for q,M in zip(c,mult)],s.zeros(20))
 remaining=[normal(s.sympify(t,locals=env).subs(bindings,simultaneous=True))for t in data['whole_cell_6_Lagrangian_equations']]+[normal(g.subs(bindings,simultaneous=True))for g in gb]
 ideal=s.Matrix.hstack(*[multiply_matrix(c)for c in remaining]);dom=lambda X:DomainMatrix.from_Matrix(X).convert_to(s.QQ);ic=list(dom(ideal).rref()[1]);I=ideal[:,ic];joint=I;comp=[]
 for k in range(20):
  trial=joint.row_join(s.eye(20)[:,k])
  if len(dom(trial).rref()[1])>joint.cols:joint=trial;comp.append(k)
 assert joint.shape==(20,20);E=s.eye(20)[:,comp];dim=len(comp);projection=dom(joint).inv().to_Matrix()[len(ic):,:];report=dict(task_id='E_reduced_nonnormal',phase='FINITE_SOURCE_ALGEBRA_CONSTRUCTED',all_passed=False,exceptional_polynomial=str(rr),a4_binding=str(phi),a0_a1_a2_bindings={str(k):str(z)for k,z in bindings.items()if k!=av[4]},monic_a3_quadratic=str(qq),all6_original_Lagrangian_equations_imposed=True,full20_algebra_basis=list(map(str,basis)),full_original_ideal_rank=len(ic),actual_finite_quotient_dimension=dim,quotient_standard_basis_indices=comp,source_plane_parameter_map_preserved=True,producer_bytes=Path(__file__).read_text(),input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [path,patch]},scope='Wholefiniteexceptional algebra, not20P4anchors. Weakesttarget: nonzero strictlinearinfinitesimal symmetry at everygeometricpoint excludesfiniteAut markedE6/A5; no blanketproperfactorclaim.')
 assert not output.exists();output.write_text(json.dumps(report,indent=2)+'\n')
 if dim==0:
  report.update(phase='COMPLETE',all_passed=True,all_exceptional_sources_empty=True,seconds=time.monotonic()-begun);output.write_text(json.dumps(report,indent=2)+'\n');print('EMPTY actualfinitealgebra');return
 unit=projection*e;qm=[projection*mult[k]*E for k in comp]
 def qmul(c):return sum([z*M for z,M in zip(c,qm)],s.zeros(dim))
 vals={symbol:projection*normal(bindings.get(symbol,symbol))for symbol in av};value_matrices={symbol:qmul(value)for symbol,value in vals.items()}
 def evaluate_coefficient(expr):
  result=s.zeros(dim,1)
  for ex,c in s.Poly(expr,*av).terms():
   term=c*unit
   for symbol,power in zip(av,ex):
    if power:term=value_matrices[symbol]**power*term
   result+=term
  return result
 restricted=s.Poly(s.sympify(data['restricted_actual_quintic'],locals=env),A,B,C);terms={ex:evaluate_coefficient(c)for ex,c in restricted.terms()};assert terms[(0,0,5)]==unit/3
 assert all(evaluate_coefficient(s.sympify(t,locals=env))==s.zeros(dim,1)for t in data['whole_cell_6_Lagrangian_equations'])
 report.update(phase='ALL_ACTUAL_D_COEFFICIENTS_EVALUATED_ALGEBRA_FIRST',all6_actual_Lag_evaluations_zero=True,algebra_first_original_D_coefficients={','.join(map(str,ex)):list(map(str,c))for ex,c in terms.items()});output.write_text(json.dumps(report,indent=2)+'\n')
 exponents=[(i,j,5-i-j)for i in range(5,-1,-1)for j in range(5-i,-1,-1)];polys=[]
 for i in range(3):
  for j in range(3):
   entry={}
   for ex,c in terms.items():
    if not ex[i]:continue
    out=list(ex);out[i]-=1;out[j]+=1;out=tuple(out);entry[out]=entry.get(out,s.zeros(dim,1))+ex[i]*c
   polys.append(entry)
 matrix=s.Matrix.vstack(*[s.Matrix.hstack(*[qmul(poly.get(ex,s.zeros(dim,1)))for poly in polys])for ex in exponents]);report.update(phase='EXACT_FINITE_LINEAR_MATRIX_READY',linear_matrix_shape=list(matrix.shape));output.write_text(json.dumps(report,indent=2)+'\n');kernel=dom(matrix).nullspace().to_Matrix().T;assert matrix*kernel==s.zeros(matrix.rows,kernel.cols)
 cols=[];labels=[]
 for k in range(kernel.cols):
  for j in range(9):
   c=kernel[j*dim:(j+1)*dim,k];M=qmul(c)
   for b in range(dim):cols.append(M[:,b]);labels.append([k,j,b])
 idealK=s.Matrix.hstack(*cols);cc=list(dom(idealK).rref()[1]);covered=len(cc)==dim;report.update(strict_linear_action_QQ_matrix_shape=list(matrix.shape),all9_linear_fields_retained=True,whole_QQ_kernel_dimension=kernel.cols,kernel_basis=[[str(z)for z in kernel[:,k]]for k in range(kernel.cols)],all_strict_linear_derivative_products_zero=True,kernel_coefficient_ideal_rank=len(cc),all_exceptional_geometric_sources_have_nonzero_strict_linear_symmetry=covered,no_E6_A5_target_in_exceptional_source=covered,finite_Aut_target_argument='Nonzero strictlinearfield cannot beEuler sinceD_C5=1/3. Itsconnectedprojectiveflow preservesuniqueH+irreducF4markedpair; E6/A5smoothhyperflexpairs havefiniteAut, soexcluded.',quotient_projection=[[str(z)for z in projection.row(i)]for i in range(dim)])
 if covered:
  selected=idealK[:,cc];weights=dom(selected).inv().to_Matrix()*unit;assert selected*weights==unit
  report.update(unit_ideal_selected_coefficients=[labels[k]for k in cc],unit_identity_weights=list(map(str,weights)),unit_identity_exactly_verified=True)
 report.update(phase='COMPLETE',all_passed=True,seconds=time.monotonic()-begun,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','actual_finite_quotient_dimension','whole_QQ_kernel_dimension','kernel_coefficient_ideal_rank','no_E6_A5_target_in_exceptional_source','seconds','max_rss_kib']}))


if __name__=='__main__':
 import sys,argparse
 if '--t37-exceptional-algebra-first' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--t37-exceptional-algebra-first',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();t37_exceptional_algebra_first(args.output)
 elif '--t37-exceptional-symmetry' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--t37-exceptional-symmetry',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();t37_exceptional_linear_symmetry(args.output)
 elif '--t37-curve-patch-cover' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--t37-curve-patch-cover',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();t37_curve_patch_cover(args.output)
 elif '--a5-hyperflex-graph' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--a5-hyperflex-graph',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();a5_smooth_hyperflex_graph(args.output)
 elif '--t37-double-conic' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--t37-double-conic',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();t37_curve_double_conic(args.output)
 elif '--t37-nonrational-factor' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--t37-nonrational-factor',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();t37_nonrational_factor(args.output)
 elif '--t37-open-locus' in sys.argv:
  parser=argparse.ArgumentParser();parser.add_argument('--t37-open-locus',action='store_true');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();t37_open_lagrangian_locus(args.output)
 else:main()

