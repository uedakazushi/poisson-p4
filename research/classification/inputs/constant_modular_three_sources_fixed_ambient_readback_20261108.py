#!/usr/bin/env python3
import hashlib, itertools, json, resource, sys, time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
resource.setrlimit(resource.RLIMIT_CPU,(175,176))
ROOT=Path('/workspace/poisson-p4');sys.path.insert(0,str(ROOT/'code'))
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import COLUMNS,ROWS,M2,monomials,sparse,vector,jacobi
started=time.monotonic(); xx=s.symbols('h t u z w');h,t,u,z,w=xx
JOB=ROOT/'research/classification/jobs/constant_modular_three_sources_fixed_ambient.json'
job=json.loads(JOB.read_text())
for row in job['exact_inputs_and_hashes']:
 assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
packet=json.loads((ROOT/job['exact_inputs_and_hashes'][0]['path']).read_text())
sourcebytes=Path(__file__).read_text();sourcehash=hashlib.sha256(sourcebytes.encode()).hexdigest()

def build(case):
 values=dict(case['parameters']);values.update({'K_'+str(i):v for i,v in enumerate(case.get('K2',{}).values())})
 pp={k:s.Symbol('p_'+k) for k in values}
 K=sum(pp['K_'+str(i)]*s.sympify(m,locals={'u':u,'z':z,'w':w}) for i,m in enumerate(case.get('K2',{})))
 if case['id'].startswith('I_'):
  alpha,beta,gamma,delta=(pp[k] for k in ['alpha','beta','gamma','delta'])
  f=z*w;A=u*(-beta*z+gamma*w)+K
  B=z*(alpha*u+beta*z+gamma*w+delta);C=w*(1-alpha*u-beta*z-gamma*w-delta)
 elif case['id'].startswith('II_'):
  alpha,beta,gamma,lam,rho,sigma,b,c=(pp[k] for k in ['alpha','beta','gamma','lambda','rho','sigma','b','c'])
  f=z*w+u;L=alpha*u+gamma*z-beta*w
  A=u+z*w-beta*u*w-gamma*u*z-rho*w*w-sigma*z*z-b*w-c*z
  B=z*L+beta*u+lam*z+rho*w+b
  C=-w*L+gamma*u+sigma*z-lam*w+c
 else:
  beta,gamma,delta,e=(pp[k] for k in ['beta','gamma','delta','e'])
  f=z*z+u;B=beta*u+gamma*z+delta*w+e
  A=u+z*z-2*z*B;C=(2*beta*z-gamma)*w+K
 Y=[s.expand(q) for q in [A,B,C]]
 assert s.expand(sum(s.diff(q,x) for q,x in zip(Y,(u,z,w)))-1)==0
 assert s.expand(sum(q*s.diff(f,x) for q,x in zip(Y,(u,z,w)))-f)==0
 P=s.zeros(5)
 raw={(1,2):A,(1,3):B,(1,4):C,(2,3):s.diff(f,w),(2,4):-s.diff(f,z),(3,4):s.diff(f,u)}
 for ij,a in raw.items():
  a=s.expand(h*h*a.subs({u:u/h,z:z/h,w:w/h},simultaneous=True))
  P[ij]=a;P[ij[::-1]]=-a
 substitutions={pp[k]:s.Integer(v) for k,v in values.items()}
 Q=P.subs(substitutions);directions=[P.diff(pp[k]).subs(substitutions) for k in values]
 assert len(directions)==case['Y_source_dimension']
 # Recover every prescribed coefficient from an original150 raw coordinate.
 vpar=vector(sparse(Q,xx))
 baseline=P.subs({a:0 for a in pp.values()})
 bv=vector(sparse(baseline,xx));dv=s.Matrix.hstack(*(vector(sparse(d,xx)) for d in directions))
 recovery=[]
 for k in range(len(values)):
  found=None
  for row in range(150):
   if bv[row]==0 and dv[row,k] in [1,-1] and all(dv[row,j]==0 for j in range(len(values)) if j!=k):
    found=(row,int(dv[row,k]));break
  assert found is not None,k
  row,sign=found
  assert sign*vpar[row]==list(values.values())[k]
  recovery.append({'parameter':list(values)[k],'raw150_row':row,'sign':sign,'recovered_value':int(sign*vpar[row])})
 assert vpar==bv+dv*s.Matrix(list(values.values()))
 return P,Q,directions,values,recovery,bv,dv

def certify(M):
 dm=DomainMatrix.from_Matrix(M).convert_to(s.QQ)
 columns=list(dm.rref()[1]);rank=len(columns)
 rows=list(dm.extract(range(M.rows),columns).transpose().rref()[1])
 det=dm.extract(rows,columns).det();assert det
 return {'shape':list(M.shape),'rank_QQ':rank,'rows':rows,'columns':columns,'determinant':str(det)}

def matrix_json(M):return [[str(q) for q in M.row(i)] for i in range(M.rows)]

def signed_pf(P):
 a=0
 for omit in range(5):
  i,j,k,l=[v for v in range(5) if v!=omit]
  a+=(-1)**omit*xx[omit]*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])
 return s.expand(a)

summaries=[]
for case in packet['cases']:
 began=time.monotonic();P,Q,directions,values,recovery,bv,dv=build(case)
 ps=sparse(Q,xx);fullJ=jacobi(ps);assert fullJ==s.zeros(350,1)
 F=signed_pf(Q);expected=s.sympify(case['expected_actual_F'],locals=dict(zip(map(str,xx),xx)))
 assert s.expand(F-expected)==0
 full126=[str(s.Poly(F,*xx).coeff_monomial(e)) for e in monomials(5)]
 assert len(full126)==126
 D,V,G,T=ambient_matrices(Q,xx,directions)
 assert D*T==s.zeros(350,T.cols)
 certificates={name:certify(M) for name,M in [('D',D),('Euler',V),('GL_Euler_scale',G),('GL_Euler_scale_source',T)]}
 rd=certificates['D']['rank_QQ'];re=certificates['Euler']['rank_QQ'];rf=certificates['GL_Euler_scale_source']['rank_QQ']
 assert re==24 and rd<=101 and rf<=49
 status=('whole_source_component24_candidate' if rd==101 else 'ambient_tangent_larger_parent_unresolved')
 out=ROOT/('reports/CLASSIFICATION_CONSTANT_MODULAR_FIXED_AMBIENT_'+case['id'].split('_')[0]+'_20261108.json')
 result={'task_id':job['task_id'],'id':case['id'],'status':status,'original_parameters':values,'all_parameter_recovery':recovery,'universal_homogeneous_brackets':{str(i)+','+str(j):str(s.expand(P[i,j])) for i,j in itertools.combinations(range(5),2)},'raw150_point':[str(a) for a in vector(ps)],'source_baseline150':[str(a) for a in bv],'source_full150_parameter_matrix':matrix_json(dv),'full350_Jacobi':[str(a) for a in fullJ],'actual_signed_F':str(F),'full126_actual_F':full126,'full350x150_D':matrix_json(D),'full150x25_Euler':matrix_json(V),'full150_GL_Euler_scale':matrix_json(G),'full150_GL_Euler_scale_source':matrix_json(T),'all_D_source_products_zero':True,'rank_certificates':certificates,'projective_ambient_tangent_dimension':150-rd-re-1,'projective_source_tangent_rank':rf-re-1,'component_conclusion_scope':'Rank101 is positive only after independent certificate replay and prior distinction. Rank<=100 is not a component or parent recognition. No adaptive second point.','executed_source_sha256':sourcehash,'job_sha256':hashlib.sha256(JOB.read_bytes()).hexdigest(),'seconds':time.monotonic()-began,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 with out.open('x') as fp:json.dump(result,fp,ensure_ascii=False,indent=2);fp.write('\n')
 summary={'id':case['id'],'D':rd,'Euler':re,'group':certificates['GL_Euler_scale']['rank_QQ'],'source_tangent':rf,'ambient_projective_tangent':result['projective_ambient_tangent_dimension'],'status':status,'path':str(out.relative_to(ROOT)),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'seconds':result['seconds']}
 summaries.append(summary);print(json.dumps(summary),flush=True)
for row in job['exact_inputs_and_hashes']:
 assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256']
final={'task_id':job['task_id'],'all_three_prescribed_points_executed':True,'all_exact_original_source_Jacobi_Pf_checks_passed':True,'all_source_dimension_consistency_checks_passed':True,'new_component_count':0,'cases':summaries,'input_bindings':job['exact_inputs_and_hashes'],'coordinates':list(map(str,xx)),'full150_column_labels':[{'pair':list(p),'monomial':list(e)} for p,e in COLUMNS],'full350_row_labels':[{'quad':list(q),'monomial':list(e)} for q,e in ROWS],'full126_F_monomials':monomials(5),'executed_source_bytes':sourcebytes,'executed_source_sha256':sourcehash,'job_sha256':hashlib.sha256(JOB.read_bytes()).hexdigest(),'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'One exact fixed witness per entire irreducible24 source; no profile continuation and no new component count until independent review.'}
out=ROOT/'reports/CLASSIFICATION_CONSTANT_MODULAR_THREE_SOURCES_FIXED_AMBIENT_20261108.json'
with out.open('x') as fp:json.dump(final,fp,ensure_ascii=False,indent=2);fp.write('\n')
print(json.dumps({'all_passed':True,'seconds':final['seconds'],'maximum_RSS_KiB':final['maximum_RSS_KiB'],'report':str(out.relative_to(ROOT)),'sha256':hashlib.sha256(out.read_bytes()).hexdigest()}),flush=True)
