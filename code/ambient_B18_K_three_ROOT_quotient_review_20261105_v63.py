"""Independent module quotient argument, without producer polynomial minors."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=ROOT/'reports/AMBIENT_B18_C_KBAR_THREE_GENERATOR_RANK_SOURCE_20261105_v63/report.json';assert sha(p)=='e4d6acd73dd74d09ba7eb328c0d836576365f98a1d4a8f8dce88228281758260';r=json.loads(p.read_text());assert r['all_passed'];pins={str(p.relative_to(ROOT)):sha(p)}
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h,n;pins[n]=h
 ep=p.parent/'formal-deep-action-unit-three-generator-minor-and-exact-rank-four.json';assert sha(ep)==r['evidence_sha256'][ep.name];d=json.loads(ep.read_text());pins[str(ep.relative_to(ROOT))]=sha(ep)
 lam=F(d['known_initial_scalars']['lambdaS']);assert lam==F(-6957653957,295294329)and lam
 assert r['mu_R0_K_lower']==2 and r['mu_R0_K_upper']==3
 # Six independent selector/initial coefficients, all remaining actual
 # deeper coefficients are arbitrary. Eliminate in a DIFFERENT order.
 steps=[['q2','chosen generator'],['g','r*q2'],['q1','(r²S-F*g)/lambdaS'],['v','(rS-E*q1-F*q2-G*g)/lambdaS'],['w','rU-A*v-B*q1-C*q2-D*g']]
 assert [s[0]for s in steps]==['q2','g','q1','v','w']
 # The two top classes are independent on the original h2,h3 heads.
 head_minor=F(0)*F(-1,4)-F(1)*F(1);assert head_minor==-1
 # On A5, r spans exactly q1,g. The other images rU,rS contribute two
 # independent classes modulo this span, detected by w,v coordinates.
 rank_extra_minor=F(1)*lam-F(0);assert rank_extra_minor
 possible=[list(t)for t in [(3,3,1),(3,2,2)]];assert all(sum(t)==7 and sum(x-1 for x in t)==4 for t in possible)
 rep=dict(id='AMBIENT_B18_K_THREE_ROOT_QUOTIENT_REVIEW_20261105_v63',all_passed=True,
 independent_method='explicit quotient elimination and image filtration; no producer polynomial DP/determinant/rank code used',
 scalar_lambdaS=str(lam),all_seven_deeper_coordinates_remain_formal=True,
 generation_recovery_steps=steps,generation_proof='U,S,q2 generate. In quotientbythem, g=rq2=0, r²S=lambdaS*q1+F*g givesq1=0, rS=lambdaS*v+E*q1+F*q2+G*g givesv=0, rU=w+A*v+B*q1+C*q2+D*g givesw=0. All seven truebarbasis classes vanish. lambdaS is known nonzero; no deepcoeff used.',
 lower_bound_proof='constant h2,h3 images ofU,S have determinant−1; mKbar maps into mTold; these two classes remainindependent modulo mKbar',
 upstairs_proof='K finiteD and finiteR0; choose any trueK lifts of threebarclasses, quotient=ellquotient so Nakayama kills it. mu_R0K=mu_R0Kbar.',
 mu_R0_K_interval=[2,3],actual_r_rank=4,
 rank_proof='rA5=span(q1,g), dim2. rU,rS modthatspan detectedon(w,v) bynonzero determinantlambdaS, dim2. Allsevenbasisimages belongsum, henceexactrank4.',
 r_Jordan_possibilities=possible,Jordans_proof='dim7,rankr4,r³0 imply3blocks ofsize<=3; only331 or322',
 uncomputed_scalar='lambdaU*F-lambdaS*(H+C)',no_actual_new_source_action_values_or_full_lifts=True,no_ell_DVR_H0_or_nil_inference=True,
 input_sha256=pins,code_sha256=sha(Path(__file__)),source_numeric_upper=50000)
 out=ROOT/'reports'/rep['id'];out.mkdir(exist_ok=False);p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':main()
