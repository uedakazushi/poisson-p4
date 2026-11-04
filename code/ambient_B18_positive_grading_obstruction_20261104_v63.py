#!/usr/bin/env python3
"""Physical stabilizer tori and a zero-weight actual B18 transverse class."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(c)for c in row]for row in M.tolist()]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 orbitp=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/report.json';op=read(orbitp);assert op['all_passed']and op['orbit_dimension']==19
 for rel,h in op['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in op['evidence_sha256'].items():assert sha(orbitp.parent/rel)==h,rel
 ocp=orbitp.parent/'certificate.json.gz';oc=read(ocp)
 tangentp=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/report.json';tp=read(tangentp);tcp=tangentp.parent/'certificate.json.gz';tc=read(tcp)
 assert tp['all_passed']and tp['actual_D_cone_embedding_dimension']==39
 assert sha(tcp)==tp['evidence_sha256']['certificate.json.gz']
 nilp=ROOT/'reports/AMBIENT_B18_NILCOTANGENT_BASIS_20261103_v63/report.json';nr=read(nilp);assert nr['all_passed']
 pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));keys=[(pair,mon)for pair in pairs for mon in mons];idx={k:i for i,k in enumerate(keys)}
 p=s.Matrix(list(map(s.Rational,oc['physical_point150'])))
 # The actual four-term B18 representative, order(h,x,y,z,w).
 expected={((0,2),(1,1)):-1,((2,3),(0,1)):1,((2,4),(0,0)):1,((2,4),(3,3)):1}
 assert {keys[i]:int(c)for i,c in enumerate(p)if c}==expected
 def wedge(a,b):
  return((a,b),1)if a<b else((b,a),-1)if a>b else(None,0)
 columns=[]
 for a,b in itertools.product(range(5),repeat=2):
  v=s.zeros(150,1)
  for k,c in enumerate(p):
   if not c:continue
   (i,j),mon=keys[k]
   if i==b:
    pair,sign=wedge(a,j)
    if sign:v[idx[pair,mon]]+=sign*c
   if j==b:
    pair,sign=wedge(i,a)
    if sign:v[idx[pair,mon]]+=sign*c
   for t,u in enumerate(mon):
    if u==a:
     shifted=list(mon);shifted[t]=b;v[idx[(i,j),tuple(sorted(shifted))]]-=c
  columns.append(v)
 O=s.Matrix.hstack(*columns);assert O.rank()==19
 O19=p.row_join(O[:,oc['GL18_generator_indices']]);assert O19==s.Matrix(oc['physical150_by19_scalar_then_GL'])
 glker=O.nullspace();assert len(glker)==6
 pure=[s.Matrix(5,5,v)for v in glker]
 I=s.eye(5);assert O*s.Matrix(I).reshape(25,1)==s.zeros(150,1)
 tf=[M-s.trace(M)*I/5 for M in pure];flat=s.Matrix.hstack(*(s.Matrix(M).reshape(25,1)for M in tf));chosen=list(flat.rref()[1]);assert len(chosen)==5
 N=[tf[j]for j in chosen];t=s.symbols('t0:5');M=sum((t[j]*N[j]for j in range(5)),s.zeros(5));powers=[];W=s.eye(5)
 for k in range(1,6):W=(W*M).applyfunc(s.expand);powers.append(W)
 assert all(s.expand(s.trace(powers[k-1]))==0 for k in range(1,6));assert powers[4]==s.zeros(5)
 for B in N:assert O*s.Matrix(B).reshape(25,1)==s.zeros(150,1)and s.trace(B)==0
 # Extended torus: GL weights(0,0,-1,0,0), tensor multiplier +1.
 T0=s.diag(0,0,-1,0,0);assert O*s.Matrix(T0).reshape(25,1)+p==s.zeros(150,1)
 extended=O.row_join(p);assert len(extended.nullspace())==7
 physicalweights=[int(1+(1 if 2 in pair else 0)*(-1)-sum(1 for v in mon if v==2)*(-1))for pair,mon in keys]
 assert all(physicalweights[i]==0 for i,c in enumerate(p)if c)
 T=s.Matrix(tc['common150_by39_normalized_tangent']);free=tc['common39_physical_coordinate_rows'];assert T.shape==(150,39)and T[free,:]==s.eye(39)
 assert T*O19[free,:]==O19
 def homogeneous_weight(v):
  ws={physicalweights[i]for i,c in enumerate(v)if c};assert len(ws)==1;return ws.pop()
 tw=[homogeneous_weight(T[:,j])for j in range(39)];ow=[homogeneous_weight(O19[:,j])for j in range(19)]
 from collections import Counter
 ct,co=Counter(tw),Counter(ow);quot={w:ct[w]-co[w]for w in sorted(set(ct)|set(co))};assert all(v>=0 for v in quot.values())and sum(quot.values())==20
 assert keys[63]==((1,2),(0,3))and physicalweights[63]==0
 assert any(T[63,j]for j in range(39))and all(O19[63,j]==0 for j in range(19))
 assert quot[0]>0
 # Diagonal fixing constraints in(H,X,Y,Z,W,tensor_sigma).
 fixing=s.Matrix([[1,-2,1,0,0,1],[-1,-1,1,1,0,1],[-2,0,1,0,1,1],[0,0,1,-2,1,1]])
 assert fixing.rank()==4 and fixing*s.Matrix([1,1,1,1,1,0])==s.zeros(4,1)and fixing*s.Matrix([0,0,-1,0,0,1])==s.zeros(4,1)
 evidence={'actual_representative':'e_y wedge (x^2 partial_h+h*x partial_z+(h^2+z^2)partial_w)','ambient_coordinate_order':['h','x','y','z','w'],'physical_nonzero_terms':[[list(pair),list(mon),str(c)]for(pair,mon),c in expected.items()],
 'GL25_kernel_basis':[[str(c)for c in v]for v in glker],'GL_stabilizer_lie_dimension':6,'tracefree_GL_stabilizer_basis':[enc(B)for B in N],'tracefree_dimension':5,
 'generic_tracefree_element':'sum_i t_i*N_i; all coefficient-polynomial traces1..5=0 and M^5=0 literally','all_five_generic_trace_polynomials_zero':True,'generic_matrix_fifth_power_all25_entries_zero':True,
 'central_GL_scalar_acts_identically_on_all_quadratic_bivectors':True,'extended_stabilizer_lie_dimension':7,'effective_extended_stabilizer_dimension':6,
 'extended_fixing_torus_GL_weights':[0,0,-1,0,0],'extended_tensor_scalar_weight':1,'physical150_torus_weights':physicalweights,'diagonal_fixing_constraint_matrix':enc(fixing),
 'actual_tangent39_column_weights':tw,'orbit19_column_weights':ow,'actual_transverse20_weight_multiplicities':{str(k):v for k,v in quot.items()},
 'zero_weight_actual_transverse_cotangent_class':{'raw_row':63,'pair':['x','y'],'quadratic_monomial':['h','z'],'function':'u63=P^{xy}_{hz}','nonzero_on_actual_T_D':True,'zero_on_full_orbit19':True,'torus_weight':0,'accepted_nil_cotangent_class_but_raw_function_not_claimed_nilpotent':True},
 'mathematical_gate':'tracefree pureGL stabilizer is nilpotent Lie; its connected SL5 subgroup is unipotent; modulo centralGL, any extended stabilizer torus injects into tensorGm, so rank<=1; displayed torus gives rank1; all maximal tori conjugate in connected stabilizer',
 'scope':'No cocharacter induced by GL5 times tensorGm fixing actualP0 can have all20 actual transverse cotangent weights strictly positive. Arbitrary abstract/formal local-ring torus actions are not classified.'}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),orbitp,ocp,tangentp,tcp,nilp]
 report={'id':ID,'all_exact_physical_tangent_and_stabilizer_weight_checks_passed':True,'actual_B18_representative_four_terms_not_square_constant_plane':True,'GL_stabilizer_lie_dimension':6,'tracefree_pure_GL_all_elements_nilpotent':True,'effective_extended_stabilizer_maximal_torus_rank':1,'actual_transverse_dimension':20,'transverse_weight_multiplicities':{str(k):v for k,v in quot.items()},'intrinsic_zero_weight_class':'u63 modulo orbit19','positive_stabilizer_induced_grading_on_all_transverse_coordinates_impossible':True,'full_nonreduced_D_used_not_Dred_or_normalizer':True,'no_full_ideal_nil_depth_or_arbitrary_formal_grading_claim':True,'no_GB_new_point_or_primary_solver':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
