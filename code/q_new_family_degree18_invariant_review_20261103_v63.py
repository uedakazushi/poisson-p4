#!/usr/bin/env python3
"""Independent integer matrix-word evaluation of the exact degree18 minor."""
import hashlib,json,math,resource,signal,sys,time
from fractions import Fraction
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def mm(A,B):
 cols=list(zip(*B))
 return [[sum(a*b for a,b in zip(row,col))for col in cols]for row in A]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));sys.set_int_max_str_digits(0)
 stem='Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63_REVIEW'
 out=ROOT/'reports'/stem;out.mkdir(exist_ok=False)
 af=ROOT/'reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63'
 ar=read(af/'report.json');assert ar['all_passed']
 for name,want in ar['input_sha256'].items():assert sha(ROOT/name)==want,name
 samples=read(af/'samples.json');exact=read(af/'exact-evaluation.json')
 oldpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63/samples.json';old=read(oldpath)
 infofile=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json';info=read(infofile)
 invreview=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json'
 assert read(invreview)['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
 dimreview=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63_REVIEW.json';assert read(dimreview)['all_passed']
 charsfile=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
 assert read(charsfile)['character_records'][18]['invariant_dimension']==9
 O=s.Matrix(info['symplectic_form16']);Oi=O.inv();reps=[s.Matrix(M)for M in info['all_ten_representation_matrices']]
 # The trace-pairing inverse is taken from the exact native input, then checked
 # against a fresh standard sp4 basis; the native evaluator is never executed.
 from weight1_structure import sp_basis
 basis=sp_basis();G=s.Matrix([[s.trace(a*b)for b in basis]for a in basis]);Gi=G.inv()
 records=[];minor=[];words=samples['selected_words'];cols=samples['selected_columns']
 assert len(set(map(tuple,words)))==len(words)==9 and all(len(w)==9 for w in words)
 for k,(row,pt)in enumerate(zip(samples['probe_rows'],samples['actual_QQ_source_points'])):
  prior=old['actual_QQ_points_and_betas'][row]
  assert pt['point']==prior['point']and pt['beta16']==prior['beta16']
  assert Fraction(prior['value_pivot'])and Fraction(prior['value_pfaffian'])
  b=s.Matrix(pt['beta16']);V=s.Matrix.hstack(*[O*R*b for R in reps]);mu=V.T*b/2
  c=Gi*mu;theta=sum((v*R for v,R in zip(c,reps)),s.zeros(16));H=theta+Oi*V*Gi*V.T
  den=math.lcm(*(int(v.q)for M in [theta,H]for v in M))
  matrices=[[[int(den*M[i,j])for j in range(16)]for i in range(16)]for M in [theta,H]]
  cache={():[[int(i==j)for j in range(16)]for i in range(16)]}
  def word(w):
   w=tuple(w)
   if w not in cache:cache[w]=mm(word(w[:-1]),matrices[w[-1]])
   return cache[w]
  values=[Fraction(sum(word(w)[i][i]for i in range(16)),den**9)for w in words]
  assert values==[Fraction(exact['evaluation'][k][j])for j in cols]
  minor.append([s.Rational(v.numerator,v.denominator)for v in values])
  records.append(dict(probe_row=row,integer_matrix_common_denominator=str(den),selected_exact9=list(map(str,values))))
 det=s.Matrix(minor).det(method='domain-ge');assert det!=0 and Fraction(str(det))==Fraction(exact['exact_selected9_minor'])
 # Full exact independence on points in the actual new image proves restriction
 # injective, since the independently audited ambient invariant dimension is 9.
 evidence=out/'integer-review.json';evidence.write_text(json.dumps(dict(records=records,exact_minor=str(det)),indent=2)+'\n')
 note=ROOT/'research/notes'/f'{stem}.md'
 files=[Path(__file__),note,af/'report.json',af/'samples.json',af/'exact-evaluation.json',oldpath,infofile,invreview,dimreview,charsfile,ROOT/'code/weight1_structure.py']
 result=dict(id=stem,all_passed=True,producer_or_native_not_reexecuted=True,independent_scaled_integer_sequential_matrix_words=True,
  all81_selected_QQ_trace9_values_equal_author=True,exact9_minor_nonzero_and_equal=True,
  all9_points_linked_to_independently_accepted_actual_QQ_degree16_sources=True,
  universal_trace_invariance_and_complete_Inv18_dimension_reused_from_accepted_reviews=True,
  complete_degree18_restriction_injective=True,defining_hypersurface_degree_lower_bound=20,
  degree20_equation_or_probe_or_generic_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in files},evidence_sha256={'integer-review.json':sha(evidence)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print('INDEPENDENT DEGREE18 INTEGER WORDS AND EXACT9MINOR PASS',result['seconds'])
if __name__=='__main__':main()
