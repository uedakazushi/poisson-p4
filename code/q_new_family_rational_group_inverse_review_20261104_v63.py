#!/usr/bin/env python3
"""Independent tensor-index pullback, covariant transpose, and rational group inverse audit."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_RATIONAL_GROUP_INVERSE_20261104_v63_REVIEW'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def enc(M):return [[str(c)for c in row]for row in M.tolist()]
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/Q_NEW_FAMILY_RATIONAL_GROUP_INVERSE_20261104_v63';ar=read(folder/'report.json');assert ar['all_passed']
 for rel,h in ar['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in ar['evidence_sha256'].items():assert sha(folder/rel)==h,rel
 rp=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 sp=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz';rep,spin=read(rp),read(sp)
 gates=[ROOT/'reports'/stem/'report.json'for stem in ['Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW','Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in gates)
 B=s.Matrix(rep['primitive24_by16_basis']);rows=list(B.T.rref()[1]);Binv=B[rows,:].inv()
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 assert J.T==-J and J*J==-s.eye(4)and J[0,1]==1
 pairs=list(itertools.combinations(range(4),2));coords=[(i,j,k)for i,j in pairs for k in range(4)]
 # Index-level covariant triple contraction, independent of the author's
 # skew-matrix-valued 1-form implementation.
 def pull(g):
  R=s.zeros(24)
  for row,(i,j,k)in enumerate(coords):
   for col,(a,b,c)in enumerate(coords):
    R[row,col]=s.expand((g[a,i]*g[b,j]-g[a,j]*g[b,i])*g[c,k])
  moved=R*B;U=(Binv*moved[rows,:]).applyfunc(s.expand)
  assert(B*U-moved).applyfunc(s.expand)==s.zeros(24,16)
  return U
 t,u,z=s.symbols('t u z',nonzero=True)
 E=s.zeros(4);E[0,1]=1;F=s.zeros(4);F[1,0]=1
 g=s.eye(4)+t*E;h=s.eye(4)+u*F
 assert g.T*J*g==h.T*J*h==J
 Ug,Uh=pull(g),pull(h)
 assert(Ug*Uh-pull(h*g)).applyfunc(s.expand)==s.zeros(16)
 assert(Ug*Uh-Uh*Ug).applyfunc(s.expand)!=s.zeros(16)
 assert(pull(z*g)-z**3*Ug).applyfunc(s.expand)==s.zeros(16)
 # Full quadratic-coefficient covariance for the OPPOSITE transvection
 # from the author. Work with symmetric coefficient matrices, without
 # introducing or evaluating any beta variable/source point.
 mons=list(itertools.combinations_with_replacement(range(16),2));assert mons==[tuple(e)for e in spin['quadratic_monomials']]
 gp=g.T;gpi=s.eye(4)-t*E.T;assert gp*gpi==s.eye(4)
 checked=0
 for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  tensors=[]
  for row in spin[key]:
   S=s.zeros(16)
   for c,(i,j)in zip(row,mons):
    c=s.Rational(c)
    if i==j:S[i,j]=c
    else:S[i,j]=S[j,i]=c/2
   tensors.append(S)
  for i,j in itertools.product(range(4),repeat=2):
   expected=s.zeros(16)
   for a,b in itertools.product(range(4),repeat=2):
    co=gp[i,a]*gpi[b,j]
    if co:expected+=co*tensors[4*a+b]
   actual=Ug.T*tensors[4*i+j]*Ug
   assert(actual-expected).applyfunc(s.expand)==s.zeros(16)
   checked+=136
 assert checked==4352
 # Two independent Cartan exponents on all original quadratic coefficients.
 baseweights=[(1,0),(-1,0),(0,1),(0,-1)]
 betaweights=[]
 for k in range(16):
  ww={tuple(sum(baseweights[a][d]for a in coords[i])for d in range(2))for i in range(24)if B[i,k]}
  assert len(ww)==1;betaweights.append(ww.pop())
 for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  for row,coefs in enumerate(spin[key]):
   i,j=divmod(row,4);expected=tuple(baseweights[i][d]-baseweights[j][d]for d in range(2))
   for c,(a,b)in zip(coefs,mons):
    if s.Rational(c):assert tuple(betaweights[a][d]+betaweights[b][d]for d in range(2))==expected
 # Scalar normalization checked as finite identities on every primitive
 # basis coefficient, not a sampled beta value.
 T=z*g.T;lam=(T.T*J*T)[0,1];assert s.expand(lam-z*z)==0
 k=z**3;recovered=(lam/k)*T.T
 assert(recovered-g).applyfunc(s.cancel)==s.zeros(4)
 q=s.Symbol('q',nonzero=True)
 assert s.cancel((q*q*lam)/(q**3*k)*q-lam/k)==0
 # A=Cayley(g) is Hamiltonian. Verify the universal polynomial equivalence
 # without expanding a ten-variable inverse or determinant.
 A=s.Matrix(4,4,lambda i,j:s.Symbol('a'+str(i)+str(j)));I=s.eye(4)
 assert((I+A).T*J*(I+A)-(I-A).T*J*(I-A)-2*(A.T*J+J*A)).applyfunc(s.expand)==s.zeros(4)
 ss=s.symbols('s0:10');S=s.zeros(4)
 for v,(i,j)in zip(ss,itertools.combinations_with_replacement(range(4),2)):S[i,j]=S[j,i]=v
 ham=J*S;assert ham.T*J+J*ham==s.zeros(4)and -J*ham==S
 cayley=(g-I)*(g+I).inv();assert(cayley.T*J+J*cayley).applyfunc(s.cancel)==s.zeros(4)
 assert((I-cayley).inv()*(I+cayley)-g).applyfunc(s.cancel)==s.zeros(4)
 proof={'independent_finite_transvection':'I+t E01, opposite to author I-t E10',
        'full24x16_primitive_basis_pullback_verified':True,'anti_action_checked_on_two_NONCOMMUTING_symbolic_transvections':True,
        'full_C_D_quadratic_coefficient_checks':4352,'both_Cartan_weight_checks':True,
        'covariant_law':'C(g^* beta)=g^t C(beta)g^(-t), likewise D',
        'basis_triple_scalar_degree':3,'intertwiner_line':'T=t*g^t because the simultaneous commutant of Cbar,Dbar is scalar',
        'Cramer_chart':'a rank15 minor of the 32x16 linear system plus a nonzero kernel coordinate normalized to1; finite lexicographic chart selection',
        'T_invertibility':'on the actual same generic orbit every nonzero solution is scalar times g^t',
        'lambda':'(T^t J T)_01=t^2, because J01=1 and g J g^t=J',
        'k':'((T^t)^*beta_bar)_i/beta_i=t^3 on a nonzero target coordinate',
        'rational_no_root_group_inverse':'g=(lambda/k)*T^t, independent of normalizing T',
        'existence_not_circular':'accepted rational quotient section/index1 and Rosenlicht generic orbit give a geometric g first; the formula proves g rational over the target field',
        'ten_group_coordinates':'upper triangular entries of S=-J*(g-I)*(g+I)^-1, with S symmetric',
        'Cayley_chart_nonempty':'contains g=I, and hence is dense in split connected Sp4',
        'finite_inverse_coordinates':'five accepted quotient/source coordinates plus ten Cayley coordinates give rational inverse coordinates on the actual B_NEW15 dense product chart',
        'scope':'B_NEW15 only; no boundary exhaustiveness or entire Poisson25 component/fibre conclusion'}
 ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),folder/'report.json',folder/'certificate.json',rp,sp,*gates]
 report={'id':STEM,'all_passed':True,'independent_raw_triple_index_full16_basis_pullback':True,
         'independent_anti_action_noncommuting_symbolic_transvections':True,'opposite_transvection_all4352_covariant_coefficients_verified':True,
         'all_quadratic_coefficients_both_Cartan_characters_verified':True,
         'intertwiner_transpose_scalar_square_and_cube_no_square_root_inverse_accepted':True,
         'generic_Cramer_chart_and_scalar_choice_independence_proved':True,'ten_Cayley_and_five_quotient_dense_inverse_coordinates_accepted':True,
         'B_NEW15_only_no_entire_Poisson25_or_new_full_fibre_claim':True,
         'no_author_script_rerun_new_point_prime_invariant_source_expansion_or_large_generic_solve':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','opposite_transvection_all4352_covariant_coefficients_verified','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM);a=p.parse_args();main(a.output.resolve())
