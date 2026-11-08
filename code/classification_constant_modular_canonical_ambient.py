#!/usr/bin/env python3
"""One predeclared whole-orbit canonical, full original ambient certificate.

Generic reader for literal constant-modular whole-source packets. No point loop,
profile, slice, adaptive minor sequence or higher jet is implemented.
"""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import COLUMNS,ROWS,monomials,sparse,vector,jacobi
ROOT=Path(__file__).resolve().parents[1]
PAIRS4=list(itertools.combinations(range(4),2))

def matrix_json(A):return [[str(q)for q in A.row(i)]for i in range(A.rows)]
def signed_F(P,xx):
 out=0
 for omit in range(5):
  i,j,k,l=[v for v in range(5)if v!=omit]
  out+=(-1)**omit*xx[omit]*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])
 return s.expand(out)
def rank_certificate(A):
 dm=DomainMatrix.from_Matrix(A).convert_to(s.QQ);rr,piv=dm.rref();rank=len(piv)
 rows=list(dm.extract(range(A.rows),piv).transpose().rref()[1]);det=dm.extract(rows,piv).det()
 assert det
 R=rr.to_Matrix();free=[j for j in range(A.cols)if j not in piv];K=s.zeros(A.cols,len(free))
 for at,j in enumerate(free):
  K[j,at]=1
  for i,k in enumerate(piv):K[k,at]=-R[i,j]
 assert A*K==s.zeros(A.rows,len(free))
 assert K.extract(free,range(len(free)))==s.eye(len(free))
 return {'shape':list(A.shape),'rank_QQ':rank,'rows':rows,'columns':list(piv),'determinant':str(det),'kernel':matrix_json(K),'free_columns':free}
def main(inp,outp,jobp):
 begun=time.monotonic();job=json.loads(jobp.read_text())
 resource.setrlimit(resource.RLIMIT_AS,(job['resource_limit']['memory_mb']*1024**2,)*2)
 cpu=job['resource_limit']['timeout_seconds']-2;resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu+1))
 for row in job['exact_inputs_and_hashes']:
  assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
 packet=json.loads(inp.read_text());assert packet['one_canonical_only']
 xx=s.symbols(' '.join(packet['homogeneous_coordinates']));h=xx[0];xa=xx[1:]
 pars={k:s.Symbol(k)for k in packet['parameters']};loc={str(q):q for q in xx};loc.update(pars)
 parse=lambda z:s.sympify(z,locals=loc)
 entries=[parse(q)for q in packet['universal_affine_brackets']]
 P4=s.zeros(4)
 for ij,a in zip(PAIRS4,entries):P4[ij]=a;P4[ij[::-1]]=-a
 for i,j,k in itertools.combinations(range(4),3):
  assert s.expand(sum(P4[a,l]*s.diff(P4[b,c],xa[l])for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(4)))==0
 div4=s.Matrix([s.expand(sum(s.diff(P4[i,j],xa[i])for i in range(4)))for j in range(4)])
 assert div4==s.Matrix([parse(q)for q in packet['universal_affine_divergence']])
 pf=s.expand(P4[0,1]*P4[2,3]-P4[0,2]*P4[1,3]+P4[0,3]*P4[1,2]);assert s.expand(pf-parse(packet['universal_affine_pf']))==0
 # Literal actual forward/inverse identities of each predeclared normalizer.
 normalizer_results=[]
 for ns in packet.get('normalizers',[]):
  nloc=dict(loc);nloc.update({q:s.Symbol(q)for q in ns['extra_parameters']})
  nparse=lambda z:s.sympify(z,locals=nloc)
  fw=s.Matrix([nparse(q)for q in ns['forward']]);iv=s.Matrix([nparse(q)for q in ns['inverse']])
  subf=dict(zip(xa,fw));subi=dict(zip(xa,iv))
  assert all(s.cancel(q)==0 for q in fw.subs(subi,simultaneous=True)-s.Matrix(xa))
  assert all(s.cancel(q)==0 for q in iv.subs(subf,simultaneous=True)-s.Matrix(xa))
  jac=fw.jacobian(xa);assert s.cancel(jac.det()-nparse(ns['determinant']))==0
  moved=(nparse(ns['tensor_scalar'])*jac*P4*jac.T).subs(subi,simultaneous=True)
  psub={pars[k]:nparse(v)for k,v in ns['parameter_image'].items()}
  expected=P4.subs(psub,simultaneous=True)
  assert all(s.cancel(a-b)==0 for a,b in zip(moved,expected))
  normalizer_results.append({'id':ns['id'],'all16_matrix_entries_and_both_inverse_products':True})
 P=s.zeros(5)
 for ij,a in zip(PAIRS4,entries):
  z=s.expand(h*h*a.subs({q:q/h for q in xa},simultaneous=True));assert s.Poly(z,*xx).total_degree()==2 or z==0
  P[ij[0]+1,ij[1]+1]=z;P[ij[1]+1,ij[0]+1]=-z
 sv={pars[k]:s.Rational(v)for k,v in packet['parameters'].items()};Q=P.subs(sv)
 assert not set().union(*(q.free_symbols for q in Q))-set(xx)
 directions=[P.diff(a).subs(sv)for a in pars.values()]
 # Each original whole-source coefficient is read back from a declared unit raw150 row.
 pv=s.Matrix([s.Poly(P[ij],*xx).coeff_monomial(e)for ij,e in COLUMNS]);base=pv.subs({q:0 for q in pars.values()})
 directions150=s.Matrix.hstack(*(pv.diff(q)for q in pars.values()))
 assert pv==base+directions150*s.Matrix(list(pars.values()))
 recovery=[]
 for name,row in packet['parameter_recovery'].items():
  index=COLUMNS.index((tuple(row['pair']),tuple(row['exponent'])))
  assert s.expand(pv[index]-pars[name])==0
  recovery.append({'parameter':name,'raw150_row':index,'coefficient':str(pv[index])})
 ps=sparse(Q,xx);assert jacobi(ps)==s.zeros(350,1)
 assert all(not q for q in jacobi(ps,raw=True).values())
 F=signed_F(Q,xx);assert s.expand(F-parse(packet['expected_actual_F']))==0
 div=s.Matrix([s.expand(sum(s.diff(Q[i,j],xx[i])for i in range(5)))for j in range(5)])
 assert div==s.Matrix([parse(q)for q in packet['expected_homogeneous_divergence']])
 E=s.Matrix(xx);canon=Q-(E*div.T-div*E.T)/5;Nh=s.simplify(5*canon[:,0]/h);assert Nh==div
 D,V,M,T=ambient_matrices(Q,xx,directions);assert D*M==s.zeros(350,M.cols)
 mats={'D350x150':D,'Euler150x25':V,'GL_Euler_scale150x51':M,'whole_source_tangent':T}
 certs={k:rank_certificate(A)for k,A in mats.items()}
 rd=certs['D350x150']['rank_QQ'];rm=certs['GL_Euler_scale150x51']['rank_QQ']
 assert certs['Euler150x25']['rank_QQ']==24 and rd+rm<=150
 assert certs['whole_source_tangent']['rank_QQ']==rm
 criterion=rd+rm==150
 result={'task_id':job['task_id'],'id':packet['id'],'status':'PASS','one_predeclared_canonical_only':True,'universal_affine_brackets':dict(zip(packet['pair_labels'],map(str,entries))),'all_universal_strict_Jacobi_and_div_Pf':True,'actual_normalizer_results':normalizer_results,'all_whole_parameter_recovery':recovery,'universal150':[str(q)for q in pv],'raw150_point':[str(q)for q in vector(ps)],'full350_Jacobi':['0']*350,'full126_actual_F':[str(s.Poly(F,*xx).coeff_monomial(e))for e in monomials(5)],'actual_F':str(F),'homogeneous_divergence':list(map(str,div)),'canonical_marked_Nh':list(map(str,Nh)),'column_order':[[list(ij),list(e)]for ij,e in COLUMNS],'row_order':[[list(ij),list(e)]for ij,e in ROWS],'matrices':{k:matrix_json(A)for k,A in mats.items()},'certificates':certs,'criterion_D_plus_M_150':criterion,'projective_orbit_dimension':rm-25,'projective_ambient_tangent_dimension':125-rd,'component_scope':'smooth orbit-closure component candidate only if criterion true; prior distinction and independent full replay still required'if criterion else 'canonical parent unresolved; no extra point/profile/jet authorized','input_bindings':job['exact_inputs_and_hashes'],'executed_source':Path(__file__).read_text(),'executed_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-begun,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256']
 with outp.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:result[k]for k in ['id','status','criterion_D_plus_M_150','projective_orbit_dimension','projective_ambient_tangent_dimension','seconds','maximum_RSS_KiB']}))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--job',type=Path,required=True);a=ap.parse_args();main(a.input,a.output,a.job)
