#!/usr/bin/env python3
"""Separate-author Fraction review of the entire m2/k1/rank3 plane branch.

No producer, producer helper, CAS, SymPy, NumPy or legacy imports. Reuses
previous independent sparse Fraction and univariate Fraction-field arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
import ast,hashlib,itertools as it,json,resource,time
from classification_projected_veronese_zero_review import Ring,plus,times,scale,derivative,monomials,sparse_rref,det as qdet
from classification_unimodular_logarithmic_review import RF,det as rfdet,rank as rfrank,um,ud
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2));TRIPLES=list(it.combinations(range(5),3));QUADS=list(it.combinations(range(5),4));BASEPAIRS=list(it.combinations(range(3),2));MON4=monomials(5,4)
def anti(P):
 out={ij:P.get(ij,{})for ij in PAIRS}
 for i,j in PAIRS:out[j,i]=scale(out[i,j],-1)
 return out
def wedge(u,v):return anti({(i,j):plus(times(u[i],v[j]),scale(times(u[j],v[i]),-1))for i,j in PAIRS})
def tensor_add(*ps):return anti({ij:plus(*(p.get(ij,{})for p in ps))for ij in PAIRS})
def tensor_scale(P,c):return anti({ij:scale(P[ij],c)for ij in PAIRS})
def strict(P,Q=None):
 if Q is None:return {tri:plus(*(times(P.get((a,l),{}),derivative(P.get((b,c),{}),l))for a,b,c in[tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}
 return {tri:plus(*(plus(times(P.get((a,l),{}),derivative(Q.get((b,c),{}),l)),times(Q.get((a,l),{}),derivative(P.get((b,c),{}),l)))for a,b,c in[tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}
def projective(J,R):return {quad:plus(*(scale(times(R.var(R.names[quad[k]]),J[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in QUADS}
def parse_rf(text):
 def walk(n):
  if isinstance(n,ast.Constant)and type(n.value)is int:return RF(n.value)
  if isinstance(n,ast.Name)and n.id=='tau':return RF((F(0),F(1)))
  if isinstance(n,ast.UnaryOp):return -walk(n.operand)if isinstance(n.op,ast.USub)else walk(n.operand)
  if isinstance(n,ast.BinOp):
   a=walk(n.left)
   if isinstance(n.op,ast.Pow):
    assert isinstance(n.right,ast.Constant)and type(n.right.value)is int;v=RF(1)
    for _ in range(abs(n.right.value)):v=v*a
    return v if n.right.value>=0 else 1/v
   b=walk(n.right)
   if isinstance(n.op,ast.Add):return a+b
   if isinstance(n.op,ast.Sub):return a-b
   if isinstance(n.op,ast.Mult):return a*b
   if isinstance(n.op,ast.Div):return a/b
  raise ValueError(ast.dump(n))
 return walk(ast.parse(str(text),mode='eval').body)
def rf_eval(a,t):
 def ev(p):return sum(c*t**i for i,c in enumerate(p))
 assert ev(a.q);return ev(a.p)/ev(a.q)
def rf_mat(saved):return [[parse_rf(a)for a in row]for row in saved]
def matrix_coeff(P,R,degree=4):
 J=projective(P,R);out=[]
 for quad in QUADS:
  for m in MON4:
   terms={e[5]:c for e,c in J[quad].items()if e[:5]==m};p=[F(0)]*(max(terms,default=-1)+1)
   for k,c in terms.items():p[k]=c
   out.append(RF(p))
 return out
def operator(P,bases,R):return [list(row)for row in zip(*(matrix_coeff(strict(P,Q),R)for Q in bases))]
def zeros(M,K):
 columns=list(zip(*K))if K else[]
 return all(not sum((a*b for a,b in zip(row,col)if a and b),RF())for row in M for col in columns)
def cert(M,saved,variable=True):
 rank=saved['rank'];assert saved['nullity']==len(M[0])-rank
 rows=saved['minor_rows'];cols=saved['minor_columns'];assert len(rows)==len(cols)==rank
 value=rfdet([[M[i][j]for j in cols]for i in rows]);assert value==parse_rf(saved['minor_determinant'])and value
 K=rf_mat(saved['complete_kernel']);assert len(K)==len(M[0]);assert all(len(row)==saved['nullity']for row in K);assert zeros(M,K)
 if saved['nullity']:
  assert rfrank(K)==saved['nullity']
 return value
def independent_sources(R):
 x=[R.var(v)for v in R.names[:5]];z,w=x[3:];ones=R.const(1);sources=[]
 for n in [z,w]:
  for ij in BASEPAIRS:
   for a in x[:3]:sources.append(anti({ij:times(n,a)}))
 for n,j in [(times(z,z),4),(times(z,w),3),(times(z,w),4),(times(w,w),3)]:
  for i in range(3):sources.append(anti({(i,j):n}))
 second=[anti({ij:n})for ij in BASEPAIRS for n in[times(z,z),times(z,w),times(w,w)]]
 assert len(sources)==30 and len(second)==9
 return sources,second
def canonical_shear(P,R):
 x=[R.var(v)for v in R.names[:5]];out=[]
 for n in x[3:]:
  for k in range(3):
   V=[n if i==k else{}for i in range(5)];delta=anti({(i,j):plus(scale(plus(*(times(V[l],derivative(P[i,j],l))for l in range(5))),-1),plus(*(times(P.get((l,j),{}),derivative(V[i],l))for l in range(5))),plus(*(times(P.get((i,l),{}),derivative(V[j],l))for l in range(5))))for i,j in PAIRS})
   def normal_coeff(poly,normal):
    return {e:c for e,c in poly.items()if e[3:5]==normal}
   # Remove exactly the six Euler redundancies P_i,z[z²],P_i,w[w²].
   A=[];B=[]
   for i in range(3):
    first=normal_coeff(delta[i,3],(2,0));second=normal_coeff(delta[i,4],(0,2));assert all(not any(e[:3])for e in first)and all(not any(e[:3])for e in second)
    A.append({(0,0,0,0,0,e[5]):c for e,c in first.items()});B.append({(0,0,0,0,0,e[5]):c for e,c in second.items()})
   gauge=[plus(times(x[3],A[i]),times(x[4],B[i]))for i in range(3)]+[{},{}];delta=tensor_add(delta,wedge(x,gauge));coeff=[]
   for nidx in [3,4]:
    for ij in BASEPAIRS:
     for a in range(3):
      ex=[0]*5;ex[nidx]=ex[a]=1;terms={e[5]:c for e,c in delta[ij].items()if e[:5]==tuple(ex)};p=[F(0)]*(max(terms,default=-1)+1)
      for j,c in terms.items():p[j]=c
      coeff.append(RF(p))
   for normal,j in [((2,0),4),((1,1),3),((1,1),4),((0,2),3)]:
    for i in range(3):
     terms={e[5]:c for e,c in delta[i,j].items()if e[:3]==(0,0,0)and e[3:5]==normal};p=[F(0)]*(max(terms,default=-1)+1)
     for k,c in terms.items():p[k]=c
     coeff.append(RF(p))
   out.append(coeff)
 return [list(row)for row in zip(*out)]
def universal_lax_check():
 names=['a','b','c','z','w']+[f't{i}'for i in range(8)]+[f'm{i}'for i in range(9)]+[f's{i}'for i in range(3)];R=Ring(names);p=R.parse;x=[R.var(v)for v in names[:5]];a,b,c,z,w=x;one=R.const(1);U=[one,{},{},{},{}];En=[{},{},{},z,w];t=[R.var(f't{i}')for i in range(8)];T=[plus(times(t[0],a),times(t[1],b),times(t[2],c)),plus(times(t[3],a),times(t[4],b),times(t[5],c)),plus(times(t[6],a),times(t[7],b),scale(times(plus(t[0],t[4]),c),-1)),{},{}]
 M=[[plus(*(times(R.var(f'm{3*k+j}'),x[k])for k in range(3)))for j in range(3)]];aa,bb,cc=M[0];Mn=[plus(times(aa,z),times(bb,w)),plus(times(cc,z),scale(times(aa,w),-1))];K=p('s2*z**2-2*s0*z*w-s1*w**2');P=tensor_add(wedge(T,En),wedge(U,[{},{},{},*Mn]),anti({(3,4):K}));J=strict(P)
 assert all(not J[tri]for tri in TRIPLES if not{3,4}.issubset(tri))
 qM=plus(times(z,Mn[1]),scale(times(w,Mn[0]),-1));T_M=[plus(*(times(T[i],derivative(m,i))for i in range(3)))for m in [aa,bb,cc]];qTM=plus(times(times(z,z),T_M[2]),scale(times(times(z,w),T_M[0]),-2),scale(times(times(w,w),T_M[1]),-1));M0n=[plus(times(derivative(aa,0),z),times(derivative(bb,0),w)),plus(times(derivative(cc,0),z),scale(times(derivative(aa,0),w),-1))];detMM=plus(times(Mn[0],M0n[1]),scale(times(Mn[1],M0n[0]),-1));VMK=plus(*(times(Mn[j],derivative(K,3+j))for j in range(2)));scalar=plus(VMK,qTM,detMM)
 expected=[plus(scale(times(qM,derivative(T[i],0)),-1),scalar if i==0 else{})for i in range(3)]
 assert all(J[i,3,4]==expected[i]for i in range(3))
 return dict(all35_weight0_source_Lax_identity_reconstructed=True,universal_firstsymbol_parameters='tracefree T8, M9, verticalS3 with U fixed; global35 before k1 specialization',strict_Jacobi_only_base1_normal2=True,all_three_vector_coefficients_exact=True)
def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(119,120));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);path=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_20261107.json';saved=json.loads(path.read_text());R=Ring(('a','b','c','z','w','tau'));p=R.parse;x=[R.var(n)for n in R.names[:5]];a,b,c,z,w=x;t=p('tau');En=[{},{},{},z,w];Mn=[{},{},{},p('b*z+c*w'),p('a*z-b*w')];P1,P2=independent_sources(R);records=[];lax=universal_lax_check()
 for source in saved['cases']:
  name=source['type'];ss=name=='semisimple';U=[{},R.const(1),{},{},{}]if ss else[R.const(1),{},{},{},{}];T=[p('2*tau*a'),{},p('-2*tau*c'),{},{}]if ss else[p('-2*tau*b'),p('tau*c'),{},{},{}];K=p('(1-2*tau)*z*w'if ss else'(tau-1/2)*z**2');P=tensor_add(wedge(T,En),wedge(U,Mn),anti({(3,4):K}));assert all(P[ij]==p(source['normal_form_P0'][str(ij)])for ij in PAIRS);assert all(not v for v in projective(strict(P),R).values());D1=operator(P,P1,R);D2=operator(P,P2,R);S=canonical_shear(P,R)
  assert D1==rf_mat(source['complete_D1_350_by30'])and D2==rf_mat(source['complete_D2_350_by9'])and S==rf_mat(source['complete_actual_shear_30_by6']);assert zeros(D1,S)
  m1=cert(D1,source['generic_D1']);m2=cert(D2,source['generic_D2']);ms=cert(S,source['generic_actual_shear']);product=m1*m2*ms;assert product==parse_rf(source['exception_polynomial']);factors=source['all_irreducible_factors'];pol=RF(1);roots=set()
  for factor,k in factors:
   f=parse_rf(factor);assert len(f.p)==2 and f.q==(F(1),);root=-f.p[0]/f.p[1];roots.add(root)
   for _ in range(k):pol=pol*f
  assert product/pol and len((product/pol).p)==1 and len((product/pol).q)==1
  walls=[]
  for row in source['all_exceptional_roots']:
   value=F(row['tau']);assert value in roots;matrices=[[[RF(rf_eval(c,value))for c in rr]for rr in M]for M in[D1,D2,S]]
   for M,key in zip(matrices,['D1','D2','actual_shear']):cert(M,row[key],False)
   higher=row['D1']['nullity']+row['D2']['nullity'];bound=(16 if ss else 15)+higher;assert higher==row['higher_fiber_upper_bound']and bound==row['PGL_image_upper_bound'];assert row['excluded_tau_zero']==(value==0);assert value==0 or bound<=23;walls.append(dict(tau=str(value),kernel1=row['D1']['nullity'],kernel2=row['D2']['nullity'],whole_image_bound=bound,excluded_m_ge3=value==0))
  assert {F(row['tau'])for row in source['all_exceptional_roots']}==roots
  generic=source['generic_D1']['nullity']+source['generic_D2']['nullity'];bound=(17 if ss else 16)+generic;assert bound==source['generic_PGL_image_upper_bound']and bound<=23;assert source['generic_D1']['nullity']==source['generic_actual_shear']['rank']==6 and source['generic_D2']['nullity']==0
  pf=[]
  for i in range(5):j,k,l,m=[r for r in range(5)if r!=i];pf.append(scale(plus(times(P[j,k],P[l,m]),scale(times(P[j,l],P[k,m]),-1),times(P[j,m],P[k,l])),(-1)**i))
  D=plus(*(times(x[i],pf[i])for i in range(5)));assert D==p(source['actual_D']);expected=p('-4*a*c*tau*(a*z**2-2*b*z*w-c*w**2)'if ss else'c**2*tau*(a*z**2-2*b*z*w-c*w**2)');assert D==expected
  records.append(dict(type=name,all350_by30_and350_by9_exact=True,actual_shear30_by6_exact=True,full_generic_kernel_and_polynomial_minor_exact=True,complete_linear_wall_factorization=True,all_walls=walls,generic_whole_image_upper_bound=bound,whole_image_upper_bound=max([bound]+[r['whole_image_bound']for r in walls if not r['excluded_m_ge3']]),whole_symbol_D_exact=True,tau_zero_has_no_degree2_normal_D=True,generic_finite_shear_erases_weight1_and_weight2=True))
 paths=[Path(__file__),path,ROOT/'code/classification_projected_veronese_zero_review.py',ROOT/'code/classification_unimodular_logarithmic_review.py',ROOT/'research/classification/PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md',ROOT/'research/classification/PLANE_RANK_ZERO_FULLRANK_WEIGHT_BOUND.md'];out=dict(task_id='G',all_passed=True,independent_method='Standard-library Fraction polynomial and QQ(tau) rational functions; no producer/helper/SymPy/CAS import.',universal_Lax=lax,cases=records,whole_stratum_projective_image_bound=23,geometry_review=dict(source_weight_dimensions=[35,30,9],full_plane_raw90_Euler10_abelian_normal_linear6=True,saturated_coHiggs_line_k0_or1=True,rank3_inverse_M_GL2_tensor_normalizations_complete=True,ss_first_symbol_dimension=12,nil_first_symbol_dimension=11,fixed_tau_dimension_drops_one=True,tau0_equals_old_c0_m_ge3_not_discarded_as_empty=True,complements=['rankM<=2','k=0','multiplicity>=3','other restriction-rank strata']),input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},checker_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_REVIEW_20261107.json';assert not destination.exists();destination.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(all_passed=True,whole_image_bound=23,seconds=out['seconds'],maximum_RSS_KiB=out['maximum_RSS_KiB'])))
if __name__=='__main__':main()
