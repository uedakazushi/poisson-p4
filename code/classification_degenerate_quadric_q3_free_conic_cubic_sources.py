#!/usr/bin/env python3
"""Seven exhaustively motivated conic+cubic configurations and full R6 sources."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();a,b,c,r,z=s.symbols('a b c r z');base=(a,b,c);xx=(a,b,c,r,z);T=s.Symbol('T');Q=s.symbols('q0:18');mm=s.symbols('m0:6');mons=[a*a,a*b,a*c,b*b,b*c,c*c]
cusp=b*b*c-a**3;node=b*b*c-a*a*(a+c);cases=[]
for m in [1,2,3]:
 f=s.Poly(s.expand((T-1)**m*(T+s.Rational(6-m,m))**(6-m)),T)
 q=f.nth(6)*b*b+f.nth(5)*a*b+f.nth(4)*a*a+f.nth(3)*b*c+f.nth(2)*a*c+f.nth(0)*c*c;assert f.nth(1)==0
 cases.append(('cusp_outside_'+str(m),q,cusp))
for m in [2,3]:
 f=s.Poly(s.expand(T**m*(T-1)**(6-m)),T);q=f.nth(6)*b*b+f.nth(5)*a*b+f.nth(4)*a*a+f.nth(3)*b*c+f.nth(2)*a*c;cases.append(('cusp_at_cusp_'+str(m),q,cusp))
cases.extend([('node_smooth_contact',2*a*a+b*b+3*a*c+c*c,node),('node_at_node',7*a*a+4*a*b+b*b+8*a*c+8*b*c,node)])
records=[]
for name,q,C in cases:
 D=s.expand(q*C);Y=s.Matrix([sum(Q[6*i+j]*m for j,m in enumerate(mons)) for i in range(3)]);poly=s.Poly(sum(s.diff(D,x)*Y[i] for i,x in enumerate(base)),*base);M,_=s.linear_eq_to_matrix(poly.coeffs(),Q);kernel=M.nullspace();rank=18-len(kernel);rec={'case':name,'quadric':str(q),'cubic':str(C),'quadric_hessian_determinant':str(s.hessian(q,base).det()),'quadratic_syzygy_matrix_shape':list(M.shape),'quadratic_syzygy_rank_QQ':rank,'quadratic_syzygy_kernel_dimension':len(kernel)}
 if len(kernel)==2:
  fields=[s.Matrix([s.expand(sum(v[6*i+j]*mon for j,mon in enumerate(mons))) for i in range(3)]) for v in kernel];U,V=fields;det=s.factor(s.Matrix.hstack(s.Matrix(base),U,V).det());ratio=s.cancel(det/D);assert ratio.is_number and ratio!=0
  du=sum(s.diff(U[i],base[i]) for i in range(3));dv=sum(s.diff(V[i],base[i]) for i in range(3));bracket=s.Matrix([s.expand(sum(U[j]*s.diff(V[i],base[j])-V[j]*s.diff(U[i],base[j]) for j in range(3))) for i in range(3)]);assert (bracket-(dv*U-du*V)).applyfunc(s.expand)==s.zeros(3,1)
  B=s.expand(-dv*r+du*z+sum(m*p for m,p in zip(mm,mons)));P=s.zeros(5)
  for i in range(3):P[i,3]=U[i];P[i,4]=V[i]
  P[3,4]=B;P-=P.T
  assert all(s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l]) for l in range(5)))==0 for i,j,k in itertools.combinations(range(5),3))
  assert all(s.expand(sum(s.diff(P[i,j],xx[i]) for i in range(5)))==0 for j in range(5))
  pf=[]
  for missing in range(5):
   i,j,k,l=[v for v in range(5) if v!=missing];pf.append((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
  actual=s.factor(sum(x*f for x,f in zip(xx,pf))/5);actualratio=s.cancel(actual/D);assert actualratio.is_number and actualratio!=0
  cofactors=[s.cancel(sum(s.diff(q,x)*field[i] for i,x in enumerate(base))/q) for field in fields];cofmat=s.Matrix([[s.diff(f,x) for x in base] for f in cofactors]);cofrank=cofmat.rank()
  rec.update(anchor_basis=[list(map(str,F)) for F in fields],full_six_magnetic_parameters=list(map(str,mm)),determinant_equals_D_times=str(ratio),full_unimodular_brackets={''.join(map(str,ij)):str(s.expand(P[ij])) for ij in itertools.combinations(range(5),2)},all_parameter_strict_Jacobi_zero=True,all_parameter_divergence_zero=True,actual_Pfaffian_equals_qc_times=str(actualratio),q_cofactor_components=['0','0','0']+list(map(str,cofactors)),q_cofactor_B_rank=cofrank)
 records.append(rec)
report={'task_id':'E_reduced_nonnormal','all_passed':True,'coordinate_order':list(map(str,xx)),'cases':records,'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start,'scope':'Finite exact anchor existence/freeness and full R6 source verification for the seven configurations. Geometric exhaustion of the configuration list is proved separately; component membership is not inferred.'}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_SOURCES_20261107.json')
with out.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps({'seconds':report['seconds'],'cases':[{k:r.get(k) for k in ['case','quadratic_syzygy_kernel_dimension','quadric_hessian_determinant','q_cofactor_B_rank']} for r in records]},indent=2))
