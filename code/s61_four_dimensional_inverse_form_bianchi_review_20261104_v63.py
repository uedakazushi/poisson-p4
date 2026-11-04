#!/usr/bin/env python3
"""Independent sparse-integer jet/exterior/Schouten review; no producer import."""
import argparse,ast,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_FOUR_DIMENSIONAL_INVERSE_FORM_BIANCHI_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
 c=dict(a)
 for m,v in b.items():
  x=c.get(m,0)+v
  if x:c[m]=x
  else:c.pop(m,None)
 return c
def scale(a,c):return {m:v*c for m,v in a.items()if v*c}
def mul(a,b):
 c={}
 for m,u in a.items():
  for n,v in b.items():
   key=tuple(sorted(m+n));c[key]=c.get(key,0)+u*v
 return {m:v for m,v in c.items()if v}
def total(xs):
 c={}
 for x in xs:c=add(c,x)
 return c
def const(v):return {():Q(v)}if v else{}
def var(i):return {(i,):1}
def sign(xs):return(-1)**sum(xs[i]>xs[j]for i in range(len(xs))for j in range(i+1,len(xs)))
def power(a,n):
 b=const(1)
 for i in range(n):b=mul(b,a)
 return b
def parse(text,variables):
 def visit(n):
  if isinstance(n,ast.Name):return var(variables[n.id])
  if isinstance(n,ast.Constant):assert isinstance(n.value,int);return const(n.value)
  if isinstance(n,ast.UnaryOp):
   v=visit(n.operand);return scale(v,-1)if isinstance(n.op,ast.USub)else v
  assert isinstance(n,ast.BinOp)
  a=visit(n.left)
  if isinstance(n.op,ast.Pow):assert isinstance(n.right,ast.Constant)and n.right.value>=0;return power(a,n.right.value)
  b=visit(n.right)
  if isinstance(n.op,ast.Add):return add(a,b)
  if isinstance(n.op,ast.Sub):return add(a,scale(b,-1))
  if isinstance(n.op,ast.Mult):return mul(a,b)
  if isinstance(n.op,ast.Div):assert set(b)=={()};return scale(a,1/b[()])
  raise AssertionError(ast.dump(n))
 return visit(ast.parse(text,mode='eval').body)
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/S61_FOUR_DIMENSIONAL_INVERSE_FORM_BIANCHI_20261104_v63';prp=folder/'report.json';pr=json.loads(prp.read_text());assert pr['all_passed']
 for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 source=ROOT/'code/s61_four_dimensional_inverse_form_bianchi_20261104_v63.py';assert sha(source)==pr['source_sha256']
 cp=folder/'certificate.json';assert sha(cp)==pr['evidence_sha256']['certificate.json'];cert=json.loads(cp.read_text())
 pairs=list(itertools.combinations(range(4),2));second=list(itertools.combinations_with_replacement(range(4),2));triples=list(itertools.combinations(range(4),3))
 names=[f'p{i}{j}'for i,j in pairs]+[f'p{i}{j}_{k}'for i,j in pairs for k in range(4)]+[f'p{i}{j}_{k}{l}'for i,j in pairs for k,l in second]+[f'h{i}'for i in range(4)]+[f'h{j}_{i}'for i in range(4)for j in range(4)]
 assert len(names)==110 and names==cert['jet_variable_order'];lookup={n:i for i,n in enumerate(names)}
 P=[[{}for j in range(4)]for i in range(4)];maps=[{}for i in range(4)]
 for t,(i,j)in enumerate(pairs):
  P[i][j]=var(t);P[j][i]=scale(P[i][j],-1)
  for k in range(4):
   maps[k][t]=6+4*t+k
   for l in range(4):maps[k][6+4*t+l]=30+10*t+second.index(tuple(sorted((k,l))))
 for k in range(4):
  for j in range(4):maps[k][90+j]=94+4*k+j
 def D(a,k):
  b={}
  for m,c in a.items():
   for t,v in enumerate(m):
    assert v in maps[k],('unsupported higher jet differentiated',v)
    n=list(m);n[t]=maps[k][v];key=tuple(sorted(n));b[key]=b.get(key,0)+c
  return {m:c for m,c in b.items()if c}
 def wedge(a,b):
  if not a or not b:return{}
  da=len(next(iter(a)));db=len(next(iter(b)));result={}
  for I in itertools.combinations(range(4),da+db):
   value={}
   for J in itertools.combinations(I,da):
    K=tuple(v for v in I if v not in J)
    if J in a and K in b:value=add(value,scale(mul(a[J],b[K]),sign(J+K)))
   if value:result[I]=value
  return result
 def dform(a):
  if not a:return{}
  degree=len(next(iter(a)));result={}
  for I in itertools.combinations(range(4),degree+1):
   value=total(scale(D(a.get(I[:t]+I[t+1:],{}),v),(-1)**t)for t,v in enumerate(I))
   if value:result[I]=value
  return result
 def fs(a,c):return {I:scale(v,c)for I,v in a.items()if scale(v,c)}
 def fadd(a,b):
  result=dict(a)
  for I,v in b.items():
   value=add(result.get(I,{}),v)
   if value:result[I]=value
   else:result.pop(I,None)
  return result
 def times(a,p):return {I:mul(v,p)for I,v in a.items()if mul(v,p)}
 Pform={I:P[I[0]][I[1]]for I in pairs};pf=scale(wedge(Pform,Pform)[(0,1,2,3)],Q(1,2))
 A=[[{}for j in range(4)]for i in range(4)]
 for i,j in pairs:
  k,l=[q for q in range(4)if q not in(i,j)];A[i][j]=scale(P[k][l],-sign((i,j,k,l)));A[j][i]=scale(A[i][j],-1)
 checks={}
 checks['own_Hodge_adjugate_both_products']=all(total(mul(P[i][k],A[k][j])for k in range(4))==(pf if i==j else{})and total(mul(A[i][k],P[k][j])for k in range(4))==(pf if i==j else{})for i in range(4)for j in range(4))
 dp=[D(pf,k)for k in range(4)];X=[total(D(P[i][j],i)for i in range(4))for j in range(4)]
 J={I:total(mul(P[i][l],D(P[j][k],l))for i,j,k in[I,(I[1],I[2],I[0]),(I[2],I[0],I[1])]for l in range(4))for I in triples}
 dual=[scale(J[tuple(j for j in range(4)if j!=i)],(-1)**i)for i in range(4)]
 checks['own_J_dual_sign']=all(dual[i]==add(scale(dp[i],-1),scale(total(mul(A[i][j],X[j])for j in range(4)),-1))for i in range(4))
 formA={I:A[I[0]][I[1]]for I in pairs};formdp={(i,):dp[i]for i in range(4)};dA=dform(formA)
 E=fadd(times(dA,pf),fs(wedge(formdp,formA),-1));edual=[scale(E[tuple(j for j in range(4)if j!=i)],(-1)**i)for i in range(4)]
 checks['own_E_dual_equals_minus_P_J']=all(edual[i]==scale(total(mul(P[i][j],dual[j])for j in range(4)),-1)for i in range(4))
 de=dform(E);checks['own_exterior_dE_sign']=de==fs(wedge(formdp,dA),2)
 checks['own_cleared_deRham_zero']=not fadd(times(de,pf),fs(wedge(formdp,E),-2))
 h=[var(90+i)for i in range(4)];B=add(total(mul(X[j],h[j])for j in range(4)),scale(total(mul(P[i][j],D(h[j],i))for i in range(4)for j in range(4)),-1))
 def contraction_theta(form,k):
  return {I[:I.index(k)]+I[I.index(k)+1:]:scale(v,(-1)**I.index(k))for I,v in form.items()if k in I}
 def spatial(form,k):return {I:D(v,k)for I,v in form.items()if D(v,k)}
 Hform={tuple(j for j in range(4)if j!=i):scale(h[i],(-1)**i)for i in range(4)}
 def bracket2(form):
  ans={}
  for k in range(4):ans=fadd(ans,fs(fadd(wedge(contraction_theta(Pform,k),spatial(form,k)),wedge(spatial(Pform,k),contraction_theta(form,k))),-1))
  return ans
 checks['own_Grassmann_Schouten_P_P_equals2J']=bracket2(Pform)==fs(J,2)
 checks['own_Grassmann_Schouten_P_H_equals_B']=bracket2(Hform).get((0,1,2,3),{})==B
 actualB=add(total(mul(X[j],dual[j])for j in range(4)),scale(total(mul(P[i][j],D(dual[j],i))for i in range(4)for j in range(4)),-1));checks['own_actual_Bianchi_all_second_jets']=not actualB
 minusPh=[scale(total(mul(P[i][j],h[j])for j in range(4)),-1)for i in range(4)]
 deRham=add(mul(pf,total(D(minusPh[i],i)for i in range(4))),scale(total(mul(mul(dp[i],P[i][j]),h[j])for i in range(4)for j in range(4)),2))
 koszul=scale(total(mul(P[i][j],add(mul(dual[i],h[j]),scale(mul(dual[j],h[i]),-1)))for i,j in pairs),-2)
 checks['own_formal_all110_operator_identity']=deRham==add(mul(pf,B),koszul)
 checks['own_actual_skew_Koszul_zero']=not total(mul(mul(dual[i],P[i][j]),dual[j])for i in range(4)for j in range(4))
 # Compare all saved nontrivial expressions via independent AST parsing.
 for key,values in [('dp',dp),('divergence_X',X),('Jacobi_dual',dual),('E_dual',edual)]:assert [parse(x,lookup)for x in cert[key]]==values,key
 assert parse(cert['pfaffian'],lookup)==pf
 assert [[parse(x,lookup)for x in row]for row in cert['P']]==P
 assert [[parse(x,lookup)for x in row]for row in cert['inverse_numerator_A']]==A
 for key,value in [('formal_Bianchi_operator',B),('formal_deRham_operator',deRham),('difference_Koszul_operator',koszul)]:assert parse(cert[key],lookup)==value,key
 checks['all_saved_primary_expressions_independently_bound']=True
 # Generic4x3 coefficient transport, plus an index proof in the review note.
 Qmat=[[var(110+3*i+j)for j in range(3)]for i in range(4)];F=[var(122+j)for j in range(3)]
 S=[[total(mul(mul(Qmat[i][a],P[i][j]),Qmat[j][b])for i in range(4)for j in range(4))for b in range(3)]for a in range(3)]
 checks['own_generic_QtPQ_skew']=all(S[a][b]==scale(S[b][a],-1)for a in range(3)for b in range(3))
 cof=[scale(total(mul(F[a],S[a][b])for a in range(3)),-2)for b in range(3)];std=[{}for b in range(3)]
 for a,b in itertools.combinations(range(3),2):std[b]=add(std[b],scale(mul(S[a][b],F[a]),-2));std[a]=add(std[a],scale(mul(S[a][b],F[b]),2))
 checks['own_standard_generator_Koszul_transport_all_components']=cof==std and not total(mul(cof[j],F[j])for j in range(3))
 allowed={0,5,*range(6,10),*range(26,30),90};block={}
 for m,c in koszul.items():
  if 90 in m and all(v in allowed for v in m):n=list(m);n.remove(90);block=add(block,{tuple(n):c})
 expected={(0,0,27):-2};assert block==expected
 checks['blockdiag_h0_cofactor_minus2_a_squared_f1_NOT_af_divisible']=all(5 not in m for m in block)
 assert all(checks.values()),checks
 evidence={'own_polynomial_representation':'sparse dictionary sorted variable-index tuples -> exact integer/Fraction; own jet derivation, Hodge, exterior target-subset wedge and Grassmann Schouten','universal_jet_variables':110,'own_checks':checks,'primary_certificate_all_expression_comparisons':True,'generic_transport_auxiliary_variables':15,'generic_transport_is_universal_fixture_not_beta_point':True,'standard_N_source_transport_proof':'S=Q^T P Q skew; cofactor(-2 F^T S)_b = sum_{a<b} -2 S_ab F_a +sum_{b<c}2 S_bc F_c, the corresponding standard Koszul sum; applies to arbitrary N and polynomialQ',
 'nondivisibility_witness':{'specialization':'P01=a,P23=f, allotherPij=0, allotherjetfamilies=0','h0_cofactor':'-2*a^2*partial_1(f)','pfaffian':'a*f','reason':'independent indeterminate f has exponent0 in nonzero cofactor, so af does not divide it'},
 'accepted_scope':'this precise cleared inverse-form operator equals Pfaffian times standardBianchi plus explicit Koszul atcofactor level; no general contraction/saturation result',
 'source135_transport':'index theorem applies to any polynomialJacobi correspondence Q; existing acceptedfullraw operator correspondence remains dependency, not reexpanded',
 'not_claimed':['saturation M:g^infinity=M','alllocalizedkernel classes descend','remaining43 source directions excluded','full80 kernel','full239','generalbetaexhaustion']}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),source,ROOT/'research/notes/S61_FOUR_DIMENSIONAL_INVERSE_FORM_BIANCHI_20261104_v63.md',prp,cp]
 report={'id':ID,'independent_review_passed':True,'all_own_universal110_jet_identities_passed':True,'independent_Grassmann_Schouten_sign_verified':True,'all_primary_saved_formulas_AST_bound':True,'arbitrary_standard_generator_Koszul_transport_accepted':True,'nonmonomial_saturation_gate_scope_accepted':True,'no_full135_newpoint_GB_or_primarysolver':True,'no_fullkernel_negative_or_quartic_rank_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'checks':len(checks),'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
