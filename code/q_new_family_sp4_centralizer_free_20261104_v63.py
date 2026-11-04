#!/usr/bin/env python3
"""ONE frozen oldbeta/prime common-commutant certificate for B_NEW."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63';P=32003
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mod(q):
 q=Q(q);assert q.denominator%P;return q.numerator*pow(q.denominator,-1,P)%P
def mm(A,B):return [[sum(A[i][k]*B[k][j]for k in range(4))for j in range(4)]for i in range(4)]
def comm(A,B):
 AB,BA=mm(A,B),mm(B,A);return [[AB[i][j]-BA[i][j]for j in range(4)]for i in range(4)]
def det(A):
 A=[[x%P for x in row]for row in A];d=1
 for j in range(len(A)):
  k=next((i for i in range(j,len(A))if A[i][j]),None)
  if k is None:return 0
  if k!=j:A[k],A[j]=A[j],A[k];d=-d
  pivot=A[j][j];d=d*pivot%P;iv=pow(pivot,-1,P)
  for i in range(j+1,len(A)):
   t=A[i][j]*iv%P
   for k in range(j,len(A)):A[i][k]=(A[i][k]-t*A[j][k])%P
 return d%P
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 pointp=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/certificate.json.gz';pointrp=pointp.parent/'report.json';point=read(pointp);pr=read(pointrp)
 assert pr['all_passed']and point['same_existing_point_index']==0
 krylovrp=ROOT/'reports/Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63/report.json';kr=read(krylovrp);assert kr['all_passed']
 assert sha(pointp)==kr['input_sha256'][str(pointp.relative_to(ROOT))]
 spinp=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz';spinrp=spinp.parent/'report.json';spin=read(spinp);sr=read(spinrp);assert sr['all_passed']and sha(spinp)==sr['certificate_sha256']
 spinreview=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json';rv=read(spinreview);assert rv['all_passed']and rv['sym2_adjoint_multiplicity2_and_full_C_D_covariance_accepted']
 krylovp=krylovrp.parent/'certificate.json';krylov=read(krylovp);assert sha(krylovp)==kr['evidence_sha256']['certificate.json']
 assert krylov['oldpoint_index']==0 and krylov['prime']==P and krylov['beta16_QQ']==point['beta16']
 beta=list(map(Q,point['beta16']));assert len(beta)==16 and all(q.denominator%P for q in beta)
 beta_mod=list(map(mod,beta));assert beta_mod==krylov['beta16_modp']and any(beta_mod)
 mons=spin['quadratic_monomials'];assert mons==[list(m)for m in itertools.combinations_with_replacement(range(16),2)]
 products=[beta[i]*beta[j]for i,j in mons];products_mod=[beta_mod[i]*beta_mod[j]%P for i,j in mons]
 matrices=[];matrices_mod=[];coefficient_denoms=set();coefficient_count=0
 for name in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  rows=spin[name];assert len(rows)==16 and all(len(row)==136 for row in rows)
  coeffs=[[Q(c)for c in row]for row in rows]
  for row in coeffs:
   for c in row:assert c.denominator%P;coefficient_denoms.add(c.denominator);coefficient_count+=1
  values=[sum((c*v for c,v in zip(row,products)),Q())for row in coeffs]
  values_mod=[sum(mod(c)*v for c,v in zip(row,products_mod))%P for row in coeffs]
  assert list(map(mod,values))==values_mod
  matrices.append([values[4*i:4*i+4]for i in range(4)]);matrices_mod.append([values_mod[4*i:4*i+4]for i in range(4)])
 C,D=matrices;Cm,Dm=matrices_mod;assert Cm==krylov['C4_modp']and Dm==krylov['D4_modp']
 basis=[];labels=[]
 for i,j in itertools.product(range(4),repeat=2):
  if i==j:continue
  X=[[0]*4 for _ in range(4)];X[i][j]=1;basis.append(X);labels.append('E'+str(i)+str(j))
 for i in range(3):
  X=[[0]*4 for _ in range(4)];X[i][i]=1;X[3][3]=-1;basis.append(X);labels.append('E'+str(i)+str(i)+'-E33')
 assert len(basis)==15 and all(sum(X[i][i]for i in range(4))==0 for X in basis)
 columns=[]
 for X in basis:columns.append([x for A in[C,D]for row in comm(X,A)for x in row])
 M=[[columns[j][i]for j in range(15)]for i in range(32)];Mp=[[mod(x)for x in row]for row in M]
 direct=[[x%P for A in[Cm,Dm]for row in comm(X,A)for x in row]for X in basis]
 assert Mp==[[direct[j][i]for j in range(15)]for i in range(32)]
 echelon={};selected=[]
 for i,row in enumerate(Mp):
  v=row[:]
  for j in sorted(echelon):
   t=v[j]
   if t:v=[(a-t*b)%P for a,b in zip(v,echelon[j])]
  pivot=next((j for j,c in enumerate(v)if c),None)
  if pivot is not None:
   iv=pow(v[pivot],-1,P);echelon[pivot]=[c*iv%P for c in v];selected.append(i)
  if len(echelon)==15:break
 rank=len(echelon);positive=rank==15;minor=[Mp[i]for i in selected]if positive else[];delta=det(minor)if positive else 0
 if positive:assert delta
 rowlabels=[{'covariant':name,'matrix_row':i,'matrix_column':j}for name in['C','D']for i in range(4)for j in range(4)]
 evidence={'prime':P,'same_existing_point_index':0,'same_saved_source_point':point['same_source_point'],'beta16_QQ':list(map(str,beta)),'beta16_modp':beta_mod,
 'original_quadratic_C_D_source':str(spinp.relative_to(ROOT)),'quadratic_monomial_order':mons,'coefficient_denominators':sorted(coefficient_denoms),'all4352_quadratic_coefficients_p_integral':coefficient_count==4352,
 'C4_QQ':[[str(x)for x in row]for row in C],'D4_QQ':[[str(x)for x in row]for row in D],'C4_modp':Cm,'D4_modp':Dm,'exact_QQ_evaluation_and_independent_modp_evaluation_agree':True,'old_Krylov_C_D_exactly_matched':True,
 'sl4_basis_order':labels,'sl4_basis_integer_matrices':basis,'commutator_map':'X maps([X,C],[X,D]), columnmap sl4->Mat4^2, shape32x15','row_order':rowlabels,'commutator_matrix_modp':Mp,
 'rank_modp':rank,'selected15_minor_rows':selected if positive else[],'minor_all15_columns':list(range(15))if positive else[],'selected15_minor_modp':minor,'selected_minor_determinant_modp':delta,'selected_minor_QQ_polynomial_degree':30 if positive else None,
 'p_integral_minor_model':'Z_(32003)[beta0,...,beta15], original homogeneousquadratic C,D, integral sl4 basis; selectedminor polynomialDelta evaluates to a32003unit at the same savedR-point',
 'nonempty_B_NEW_open_proved':positive,'open_definition':'rank ofcommutator(sl4,C,D)=15 (invariantopen); chosenDelta!=0 is one witnesschart',
 'stabilizer_conclusion_if_positive':'Any GL4matrix commutingwithbothC,D isscalar, sincesubtracttrace/4 gives a tracelesskernelvector. Any Sp4 stabilizer ofbeta fixesC,D, hence is±I. Odd centralcharacter gives(-I)beta=-beta, andDelta!=0 forcesbeta!=0. Thus stabilizer isidentity on thisnonemptyB_NEW open.',
 'scope':'generic Sp4 action onfirstjetB_NEW only; no claim about allB_NEW points, fullPoissoncomponent25 rationality, orbitsectionor targetfieldindex',
 'rank_deficient_action':'STOP; nootherpointprimeorsearch authorized; rankdrop would beinconclusive'}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),pointp,pointrp,spinp,spinrp,spinreview,krylovp,krylovrp,ROOT/'research/notes/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW.md']
 report={'id':ID,'all_same_point_source_and_commutator_checks_passed':True,'same_old_beta_point_index':0,'prime':P,'commutator_map_shape':[32,15],'commutator_rank':rank,'nonzero15_minor_modp':delta,'common_GL4_commutant_scalars_on_B_NEW_open':positive,'generic_firstjet_B_NEW_Sp4_stabilizer_identity_proved':positive,'fullcomponent25_or_allpoint_freeness_claimed':False,'no_new_beta_point_prime_search_GB_or_orbitsection_computation':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-started,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
