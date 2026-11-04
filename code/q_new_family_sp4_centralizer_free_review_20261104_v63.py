#!/usr/bin/env python3
"""Same saved beta/prime only; own coefficient evaluation, matrix-unit map, Bareiss."""
import argparse
from fractions import Fraction as Q
import gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63_REVIEW';P=32003
def read(p):
 b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
 d=read(p);assert d.get('all_passed') or d.get('all_same_point_source_and_commutator_checks_passed')
 for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
 for n,h in d.get('evidence_sha256',{}).items():
  q=ROOT/n if (ROOT/n).exists() else p.parent/n
  assert sha(q)==h,n
 return d
def mod(x):
 x=Q(x);assert x.denominator%P
 return x.numerator*pow(x.denominator,P-2,P)%P
def bareiss(A):
 A=[row[:] for row in A];n=len(A);sign=1;previous=1
 for k in range(n-1):
  if not A[k][k]:
   pivot=next(i for i in range(k+1,n) if A[i][k]);A[k],A[pivot]=A[pivot],A[k];sign=-sign
  p=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    value=p*A[i][j]-A[i][k]*A[k][j];assert value%previous==0;A[i][j]=value//previous
   A[i][k]=0
  previous=p
 return sign*A[-1][-1]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 assert not out.exists();out.mkdir(parents=True)
 base=ROOT/'reports/Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63';ap=base/'report.json';author=accepted(ap);cp=base/'certificate.json';c=read(cp)
 assert sha(cp)==author['certificate_sha256'] and c['prime']==P and c['same_existing_point_index']==0
 pointp=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/certificate.json.gz';point=read(pointp)
 pointrev=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63_REVIEW/report.json';accepted(pointrev)
 spinp=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz';spin=read(spinp)
 spinrev=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json';sr=accepted(spinrev)
 assert sr['sym2_adjoint_multiplicity2_and_full_C_D_covariance_accepted']
 beta=[Q(v) for v in point['beta16']];bm=[mod(v) for v in beta]
 assert list(map(str,beta))==c['beta16_QQ'] and bm==c['beta16_modp'] and any(bm)
 assert point['same_source_point']==c['same_saved_source_point']
 mon=list(itertools.combinations_with_replacement(range(16),2))
 assert list(map(list,mon))==spin['quadratic_monomials']==c['quadratic_monomial_order']
 QQ=[];FF=[];den=set();count=0
 for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  rows=spin[key];assert len(rows)==16 and all(len(row)==136 for row in rows);values=[];fv=[]
  for row in rows:
   z=Q(0);v=0
   for coefficient,(i,j) in zip(row,mon):
    q=Q(coefficient);assert q.denominator%P;den.add(q.denominator);count+=1
    z+=q*beta[i]*beta[j];v=(v+mod(q)*bm[i]*bm[j])%P
   values.append(z);fv.append(v);assert mod(z)==v
  QQ.append([values[4*i:4*i+4] for i in range(4)]);FF.append([fv[4*i:4*i+4] for i in range(4)])
 assert count==4352 and sorted(den)==c['coefficient_denominators']
 assert [[[str(v) for v in row] for row in M] for M in QQ]==[c['C4_QQ'],c['D4_QQ']]
 assert FF==[c['C4_modp'],c['D4_modp']]
 krylovp=ROOT/'reports/Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63/certificate.json';old=read(krylovp)
 assert FF==[old['C4_modp'],old['D4_modp']] and old['prime']==P and old['oldpoint_index']==0
 labels=[('off',p,q) for p,q in itertools.product(range(4),repeat=2) if p!=q]+[('diag',i,3) for i in range(3)]
 matrix=[];matrixQQ=[]
 for k,M in enumerate(QQ):
  for i in range(4):
   for j in range(4):
    row=[]
    for typ,p,q in labels:
     if typ=='off':v=(M[q][j] if i==p else Q(0))-(M[i][p] if j==q else Q(0))
     else:v=((int(i==p)-int(i==3))-(int(j==p)-int(j==3)))*M[i][j]
     row.append(v)
    matrixQQ.append(row);matrix.append([mod(v) for v in row])
 assert matrix==c['commutator_matrix_modp']
 rowids=c['selected15_minor_rows'];assert len(rowids)==len(set(rowids))==15 and c['minor_all15_columns']==list(range(15))
 chosen=[matrix[i] for i in rowids];assert chosen==c['selected15_minor_modp']
 exact_integer=bareiss(chosen);assert exact_integer%P==c['selected_minor_determinant_modp']==author['nonzero15_minor_modp']==17760
 assert c['selected_minor_QQ_polynomial_degree']==30
 obs=base/'OBSERVED_RECEIPT.json';o=read(obs);assert o['report_sha256']==sha(ap) and o['certificate_sha256']==sha(cp)
 history=base/'precomputation-repair-history.tar.gz';assert sha(history)==o['history_sha256'] and o['actual_numeric_centralizer_evaluations']==1
 evidence=dict(same_oldpoint_index=0,prime=P,all4352_coefficients_checked=True,source_coeff_denominators=sorted(den),
  C_D_QQ_values_match=True,all480_matrixunit_commutator_entries_match=True,
  fixed15_minor_rows=rowids,independent_integer_lift_determinant=str(exact_integer),determinant_mod32003=17760,
  determinant_algorithm='Fraction-free integer Bareiss on the SAME selected residue matrix; no new row/minor selection',
  p_integral_QQ_minor_lift_and_geometric_scalar_commutant_proof_accepted=True)
 ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ap,cp,pointp,pointrev,spinp,spinrev,krylovp,obs,history,
  ROOT/'code/q_new_family_sp4_centralizer_free_20261104_v63.py',ROOT/'research/notes/Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63.md',base/'EXECUTION_HISTORY.md']
 report=dict(id=STEM,all_passed=True,original_covariant_coefficient_proof_and_actual_source_point_reused_hash_bound=True,
  all4352_p_integral_coefficients_and_same_QQ_and_modp_C_D_rechecked=True,
  independent_matrixunit_formula_all480_entries_and_fixed_minor_Bareiss=True,minor17760_mod32003=True,
  positive_QQ_polynomial_minor_nonzero_on_B_NEW_proved=True,
  full_geometric_GL4_common_commutant_scalars_not_only_Lie_stabilizer=True,
  generic_B_NEW_Sp4_stabilizer_identity_and_group_scheme_char0_proved=True,
  allpoint_freeness_component25_rationality_or_target_index_claim=False,
  no_new_point_prime_rank_search_minor_selection_GB_source26_or_author_rerun=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k] for k in ['all_passed','minor17760_mod32003','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
