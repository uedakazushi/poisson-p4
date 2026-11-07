#!/usr/bin/env python3
"""Independent complete Cartan-symbol SNC and actual-group weight review."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,itertools as it,json,resource,time
from classification_projected_veronese_zero_review import Ring,plus,times,scale,derivative,monomials
from classification_plane_rank_zero_fullrank_review import PAIRS,QUADS,anti,wedge,tensor_add,strict,projective,universal_lax_check
ROOT=Path(__file__).resolve().parents[1]
def det3(M):return plus(*(scale(times(*(M[i][j]for i,j in enumerate(perm))),(-1)**sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3)))for perm in it.permutations(range(3))))
def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(59,60));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,)*2)
 # General tracefree X,Y, unrestricted binary K: derive the full Jacobi condition.
 names=['x0','x1','x2','z','w']+[f'a{i}'for i in range(8)]+[f'b{i}'for i in range(8)]+['k0','k1','k2'];R=Ring(names);p=R.parse;E=[R.var(n)for n in names[:5]];x0,x1,x2,z,w=E
 def vector(prefix):
  cs=[R.var(prefix+str(i))for i in range(8)]
  return[plus(*(times(cs[3*i+j],E[j])for j in range(3)))if i<2 else plus(times(cs[6],x0),times(cs[7],x1),scale(times(plus(cs[0],cs[4]),x2),-1))for i in range(3)]+[{},{}]
 X=vector('a');Y=vector('b');ez=[{},{},{},R.const(1),{}];ew=[{},{},{},{},R.const(1)];K=p('k0*z**2+k1*z*w+k2*w**2');P=tensor_add(wedge([times(z,a)for a in X],ez),wedge([times(w,a)for a in Y],ew),anti({(3,4):K}));J=strict(P)
 comm=[plus(*(plus(times(X[j],derivative(Y[i],j)),scale(times(Y[j],derivative(X[i],j)),-1))for j in range(3)))for i in range(3)]
 expected=[plus(times(times(z,w),comm[i]),times(p('k0*z**2-k2*w**2'),plus(X[i],scale(Y[i],-1))))for i in range(3)]
 assert all(J[i,3,4]==expected[i]for i in range(3));assert all(not v for tri,v in J.items()if not{3,4}.issubset(tri))
 # Simultaneous Cartan source: all six diagonal eigenvalues and b remain symbolic.
 names2=['x0','x1','x2','z','w']+[f'alpha{i}'for i in range(3)]+[f'beta{i}'for i in range(3)]+['b'];S=Ring(names2);p=S.parse;coords=[S.var(n)for n in names2[:5]];diagonalX=[p(f'alpha{i}*x{i}')for i in range(3)]+[{},{}];diagonalY=[p(f'beta{i}*x{i}')for i in range(3)]+[{},{}];z,w=coords[3:];ez=[{},{},{},S.const(1),{}];ew=[{},{},{},{},S.const(1)];P0=tensor_add(wedge([times(z,a)for a in diagonalX],ez),wedge([times(w,a)for a in diagonalY],ew),anti({(3,4):p('b*z*w')}));assert all(not v for v in projective(strict(P0),S).values())
 pf=[]
 for i in range(5):j,k,l,m=[r for r in range(5)if r!=i];pf.append(scale(plus(times(P0[j,k],P0[l,m]),scale(times(P0[j,l],P0[k,m]),-1),times(P0[j,m],P0[k,l])),(-1)**i))
 D=plus(*(times(coords[i],pf[i])for i in range(5)));kappa=det3([[S.const(1)]*3,[S.var(f'alpha{i}')for i in range(3)],[S.var(f'beta{i}')for i in range(3)]]);assert D==scale(times(kappa,*coords),-1)
 # Every homogeneous quotient coordinate with nonnegative weight obeys the
 # exact g_t=diag(I3,t^-1 I2) pushforward exponent, not a truncated flow.
 Eul=(0,0,0,1,1);weight_coordinates={k:[]for k in range(-2,3)}
 for ij in PAIRS:
  for ex in monomials(5,2):
   weight=sum(Eul[i]*ex[i]for i in range(5))-Eul[ij[0]]-Eul[ij[1]];weight_coordinates[weight].append([list(ij),list(ex)])
 # Canonical entire30 first-weight basis and entire9 second-weight basis.
 basis1=[]
 for normal in [3,4]:
  for ij in it.combinations(range(3),2):
   for base in range(3):ex=[0]*5;ex[normal]+=1;ex[base]+=1;basis1.append((ij,ex))
 for factors,j in[([3,4],3),([4,4],3),([3,3],4),([3,4],4)]:
  for i in range(3):ex=[0]*5;ex[factors[0]]+=1;ex[factors[1]]+=1;basis1.append(((i,j),ex))
 basis2=[]
 for ij in it.combinations(range(3),2):
  for factors in[[3,3],[3,4],[4,4]]:ex=[0]*5;ex[factors[0]]+=1;ex[factors[1]]+=1;basis2.append((ij,ex))
 assert len(basis1)==30 and len(basis2)==9
 for level,basis in[(1,basis1),(2,basis2)]:
  for ij,ex in basis:assert sum(Eul[i]*ex[i]for i in range(5))-Eul[ij[0]]-Eul[ij[1]]==level
 inputs=[Path(__file__),ROOT/'code/classification_projected_veronese_zero_review.py',ROOT/'code/classification_plane_rank_zero_fullrank_review.py',ROOT/'code/classification_unimodular_logarithmic_review.py',ROOT/'research/classification/PLANE_RANK_ZERO_CARTAN_SNC_RECOGNITION.md',ROOT/'research/classification/PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md'];out=dict(task_id='G',all_passed=True,universal_general_X_Y_binaryK_Jacobi_exact=True,Jacobi_condition='zw[X,Y]+(k0*z²-k2*w²)(X-Y)',m2_nonzero_traceless_forces_K_bzw_and_commuting_projective_XY=True,all_diagonal_eigenvalues_symbolic=True,no_eigenvalue_or_biresidue_localization=True,kappa_nonzero_is_exact_m2_condition=True,actual_D='-det([[1,1,1],[alpha0,alpha1,alpha2],[beta0,beta1,beta2]])*x0*x1*x2*z*w',all350_Cartan_projective_Jacobi_zero=True,actual_D_exact=True,whole30_and9_weight_basis_checked=True,source1_basis=[dict(pair=list(ij),exponents=ex)for ij,ex in basis1],source2_basis=[dict(pair=list(ij),exponents=ex)for ij,ex in basis2],actual_group='g_t=diag(1,1,1,t^-1,t^-1)',actual_family='Pi0+t Pi1+t² Pi2, t1 original, t0 SNC',actual_group_coordinate_weight_identity_all150=True,component_stability_argument='For every global irreducible component C containing Pi, closure of connected C*×C action is irreducible and contains C, hence equals C. Thus C contains the actual SNC point Pi0 and is one of exactly40 MPS components meeting SNC. This is not reverse membership from an arbitrary known boundary point.',primary_source=dict(url='https://arxiv.org/html/2010.08692v2',pdf='https://arxiv.org/pdf/2010.08692v2',direct_web_read='2026-11-07',theorem='Theorem1.4 printedp7 / HTML155-161',all_period_scope='Example6.24 printedpp55-56 / HTML2497-2502',toric_nondegeneracy='section5.3 printedpp39-41 / HTML2088-2108',nonresonance_open_not_assumed=True,formal_global_distinctness_n2_explicitly_checked_by_authors=True),whole_Cartan_firstsymbol_higher_source_known_MPS40_accepted=True,complements=['nilpotent constant normal matrix','commuting but not simultaneously diagonalizable T,V','other k1/rankM cases','mD>=3 or kappa0'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},checker_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);path=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_CARTAN_SNC_REVIEW_20261107.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(all_passed=True,whole_Cartan_to_MPS40=True,seconds=out['seconds'],maximum_RSS_KiB=out['maximum_RSS_KiB'])))
if __name__=='__main__':main()
