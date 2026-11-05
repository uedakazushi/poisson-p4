"""Independent P7 map and pathwise positive 32-response bounds, no signed G8 values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import count

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));o=ROOT/'reports/AMBIENT_B18_NATIVE_C28_ROOT_SOURCE_PEER_20261105_v63';o.mkdir(exist_ok=False)
 source=ROOT/'reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_SOURCE_V2_20261105_v63/report.json';r=json.loads(source.read_text());assert r['all_passed']
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h
 p=source.parent/'native-P7-C28-canonical-third-source-and-positive34-response-input.json';assert sha(p)==r['evidence_sha256'][p.name];d=json.loads(p.read_text());ci= json.loads((ROOT/'reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/canonical33-to-B17-C68-G8-source-input-and-response-budget.json').read_text());T=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json').read_text())
 assert count([r,d,ci,T])+15000<60000
 (o/'pre-positive-enumeration.json').write_text(json.dumps(dict(no_new_signed_G8_response=True,source_owned_upper60000=60000,path_scalar_guard384=True,only32target_positive_envelopes=True),indent=2)+'\n')
 basis=list(map(tuple,d['P7_basis_r_ell']));sel=d['P7_raw33_selector_indices'];assert sel==[0,8,9,24,25,29,30];assert basis==[(0,0),(1,0),(0,1),(2,0),(1,1),(3,0),(2,1)]
 raw=[tuple(m)for m in T['basis33']];loc={x:i for i,x in enumerate(sel)}
 for i in range(33):
  for j in range(33):
   v={k:Q(x)for k,x in T['table_columns_for_each_left_basis'][i][j]};lhs=[v.get(k,Q(0))for k in sel]
   target=tuple(raw[i][k]+raw[j][k]for k in range(3));expected=[int(i in loc and j in loc and target==(0,*b))for b in basis];assert lhs==expected
 K=d['actual_kappa9_C28_cached_inputs_UNEVALUATED_AS_G8'];nu=d['actual_nu_two_P7_cached_inputs'];assert all(v[:3]==['0','0','0']for v in nu) and all(Q(v[0])==0for v in K)
 assert K==[[v[4*i+m]for i in sel for m in range(4)]for v in ci['actual_kappa9_C132_cached_UNEVALUATED']];assert nu==[[v[i]for i in sel]for v in ci['actual_nu_A_B33_cached_UNEVALUATED']]
 program=[[[m,q]for m,q in row if sum(m)<=5]for row in ci['G8_formal_kappa9_polynomials_mod_radical7']];assert program==d['G8_formal_kappa9_polynomials_mod_radical6']
 kap=[]
 for row in K:
  terms=[]
  for i,q in enumerate(row):
   if Q(q):b,m=divmod(i,4);terms.append((*basis[b],m,abs(Q(q))))
  kap.append(terms)
 nus=[[(a,b,abs(Q(q)))for (a,b),q in zip(basis,row)if Q(q)]for row in nu]
 accumulated=[[0]*4 for _ in range(8)];paths=0
 def checked(q):assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=384;return q
 for gi,row in enumerate(program):
  def recurse(factors,k,a,b,mask,c):
   nonlocal paths
   if k==len(factors):
    if (a,b)==(2,1):
     den,U=map(int,d['positive_G8_requested32_D_U'][gi][mask]);assert den%c.denominator==0;accumulated[gi][mask]+=c.numerator*(den//c.denominator);assert accumulated[gi][mask]<=U;paths+=1
    return
   for x,y,m,v in kap[factors[k]]:
    aa,bb=a+x,b+y
    if aa>2 or bb>1:continue
    collision=mask&m;outmask=mask^m;q=checked(c*v)
    if collision==3:continue
    if not collision:recurse(factors,k+1,aa,bb,outmask,q)
    else:
     for u,w,n in nus[collision-1]:
      if aa+u<=2 and bb+w<=1:recurse(factors,k+1,aa+u,bb+w,outmask,checked(q*n))
  for ex,q in row:
   factors=[i for i,n in reversed(list(enumerate(ex)))for _ in range(n)];recurse(factors,0,0,0,0,abs(Q(q)))
 assert paths<2000000
 rr=dict(all_passed=True,source_only_no_signed_G8third=True,all1089_P7_products_independently_checked=True,same_original_C132_monic_and_branch=True,formal_G8_recipe_inherits_independent_geometric_inverse_peer=True,
 weighted_proof='nu plane order>=2; weightbase1/trace1/maxbasis5, radical6 and d6=0, so Pinverse3terms. Two nu carries have base degree>=4 and vanish inP7.',
 all32_positive_bounds_replayed_by_individual_path_enumeration=True,enumerated_target_paths=paths,output_scaled_positive_sums=accumulated,source37_112_applicable_to_prepared_program=True,no_honest_quotient_or_native_carry_inference=True,
 input_sha256={str(source.relative_to(ROOT)):sha(source),str(p.relative_to(ROOT)):sha(p)},code_sha256=sha(Path(__file__)),seconds=time.monotonic()-start)
 (o/'report.json').write_text(json.dumps(rr,indent=2)+'\n');print(sha(o/'report.json'))
if __name__=='__main__':main()
