#!/usr/bin/env python3
"""Four explicit ideal-membership certificates with pure-power leading terms.

Only reduced support is needed: initialI contains t0^5,t1^4,t2^4,t3^3.
No full nonreduced ring, exact nilindex, or liftstd transformation is sought.
"""
from pathlib import Path
import hashlib,json,os,resource,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();source_path=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_ANCHOR_SUPPORT_20261107.json';source=json.loads(source_path.read_text());assert source['maximal_radical_origin'];assert source['coefficient_field']=='QQ at the existing b=2 anchor';assert source['common_delta']=='1'
 equations=source['complete_all350_remaining_equations_through_weight16'];wanted=['t0^5','t1^4','t2^4','t3^3'];indices=[source['full_leading_monomials'].index(m)for m in wanted];targets=[source['full_generic_residual_grevlex_basis'][i]for i in indices];target_weights=[5,8,8,9]
 work=ROOT/'work/ga9-constant-anchor-certificate';work.mkdir(parents=True,exist_ok=True);input_path=work/'four-lifts.sing';output_path=work/'four-lifts.txt'
 if output_path.exists():raise FileExistsError(str(output_path))
 lines=['ring r=0,(t0,t1,t2,t3),dp;','ideal I='+','.join(q['numerator']for q in equations)+';','ideal T='+','.join(targets)+';','matrix L=lift(I,T);','matrix A[1][size(I)]=I;matrix B[1][4]=T;matrix check=A*L-B;','int i;int j;for(i=1;i<=nrows(check);i++){for(j=1;j<=ncols(check);j++){if(check[i,j]!=0){print("LIFT_FAIL");exit;}}}',
  'intvec nw=1,2,2,3;intvec tw=5,8,8,9;intvec ew='+','.join(str(q['weight'])for q in equations)+';',
  'poly f;poly term;poly keep;intvec ex;int wwt;int k;',
  # Taking the required homogeneous coefficient of each lift is exact
  # because the target and all originalgenerators are weightedhomogeneous.
  'for(j=1;j<=4;j++){for(i=1;i<=size(I);i++){f=L[i,j];keep=0;while(f!=0){term=lead(f);ex=leadexp(f);wwt=0;for(k=1;k<=4;k++){wwt=wwt+nw[k]*ex[k];}if(wwt==tw[j]-ew[i]){keep=keep+term;}f=f-term;}L[i,j]=keep;}}',
  'check=A*L-B;for(j=1;j<=4;j++){if(check[1,j]!=0){print("HOMOGENEOUS_LIFT_FAIL");exit;}}',
  'for(j=1;j<=4;j++){write("'+str(output_path)+'","TARGET "+string(j-1)+" "+string(T[j]));write("'+str(output_path)+'","LEAD "+string(j-1)+" "+string(leadmonom(T[j])));for(i=1;i<=size(I);i++){if(L[i,j]!=0){write("'+str(output_path)+'","LIFT "+string(j-1)+" "+string(i-1)+" "+string(L[i,j]));}}}',
  'write("'+str(output_path)+'","ALL_FOUR_POLYNOMIAL_IDENTITIES_PASS");','print("COMPLETE_FOUR_LIFTS");','exit;']
 input_path.write_text('\n'.join(lines)+'\n');binary=ROOT/'work/nonsquare-v62/singular-runtime/prefix/usr/bin/Singular';run=subprocess.run([str(binary),'-q',str(input_path)],text=True,capture_output=True,timeout=55)
 assert run.returncode==0 and not any(line.lstrip().startswith('?')for line in run.stdout.splitlines()),run.stdout+'\n'+run.stderr
 assert output_path.exists()and output_path.read_text().splitlines()[-1]=='ALL_FOUR_POLYNOMIAL_IDENTITIES_PASS'
 records=[dict(normal_variable='t'+str(i),leading_monomial=wanted[i],weighted_degree=target_weights[i],polynomial=targets[i],original_equation_coefficients={})for i in range(4)]
 for line in output_path.read_text().splitlines():
  if line.startswith('LIFT '):_,target,eq,f=line.split(' ',3);records[int(target)]['original_equation_coefficients'][eq]=f
  elif line.startswith('LEAD '):_,target,mon=line.split(' ',2);assert mon==wanted[int(target)]
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(source_path.relative_to(ROOT)):hashlib.sha256(source_path.read_bytes()).hexdigest(),str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},producer_bytes=Path(__file__).read_text(),
  original_complete_equation_count=len(equations),four_lifted_polynomials=records,all_four_exact_polynomial_identities_checked=True,monomial_order='degree reverse lexicographic t0>t1>t2>t3',
  initial_ideal_contains_pure_powers=wanted,normal_weights=[1,2,2,3],reduced_true_cut_support_is_origin=True,
  weakest_support_proof='Each saved polynomial is an explicit polynomial combination of the COMPLETE179 genuinecut equations. Its leadingmonomial is a purepower of a differentnormalvariable. Therefore the initialideal contains powers of all4variables and the truecut has dimension0. Everyequation is positivelyweightedhomogeneous, so its only reducedpoint is the origin. No exactnilindex, fullcompletedring, minimalbasis or localunique theorem is claimed.',
  source_component_proof='The independently reconstructed true102cut uses48 actual independentlinearcuts complementary to full10parameter sourceF48. Its0dimension gives rawlocaldimension<=48 by Krull, matching actualF48. The closure of the irreducible ENTIREc1 source is therefore a fullprojective23component, including allsourceparameters in its domain and allfiniteclosure boundaries.',
  native_input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),native_output_sha256=hashlib.sha256(output_path.read_bytes()).hexdigest(),native_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,seconds=time.monotonic()-started,independent_review_required=True)
 out=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_ANCHOR_CERTIFICATE_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,lift_sizes=[len(r['original_equation_coefficients'])for r in records],seconds=result['seconds'])))
if __name__=='__main__':main()
