#!/usr/bin/env python3
"""Exact normalization of existing norm output to triangular finitefree112 cover."""
import hashlib,json,math,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_FINITE_CURVE_COVER_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def poly(p):return{tuple(m):Q(c)for m,c in p if Q(c)}
def ser(p):return[[list(m),str(c)]for m,c in sorted(p.items(),reverse=True)]
def mul(p,q):
 r={}
 for m,c in p.items():
  for n,d in q.items():
   z=tuple(x+y for x,y in zip(m,n));v=r.get(z,Q())+c*d
   if v:r[z]=v
   else:r.pop(z,None)
 return r
def scale(p,c):return{m:v*c for m,v in p.items()if v*c}
def shift(p,z):return{tuple(x+y for x,y in zip(m,z)):c for m,c in p.items()}
def main():
 t=time.monotonic();nb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_20261104_v63';nr=load(nb/'report.json');np=nb/'native-certificate.json';assert sha(np)==nr['evidence_sha256'][np.name];n=load(np);h=poly(n['remaining_h_coprime_U']);R=poly(n['det_cleared_matrix_R']);one={(0,0,0):Q(1)};a={(1,0,0):Q(1)};C={(1,0,0):Q(100),(0,2,0):Q(18),(0,1,0):Q(-9),(0,0,0):Q(-2)}
 assert len(n['localizer_removal_steps'])==31
 for i,s in enumerate(n['localizer_removal_steps']):assert poly(s['removed_gcd_with_U'])==(scale(shift(C,(1,0,0)),Q(1,100))if i<4 else a)
 den=math.lcm(*(v.denominator for v in h.values()));cc=Q(math.gcd(*(int(v*den)for v in h.values())),den);assert cc==100000000;prim=scale(h,1/cc);C4=mul(mul(C,C),mul(C,C));rhs=scale(mul(shift(C4,(31,0,0)),prim),Q(n['R_QQ_content']));assert rhs==R
 assert all(m[2]==0 for m in h);emax=max(m[1]for m in h);assert emax==28;leading={m:c for m,c in h.items()if m[1]==28};assert len(leading)==1 and (0,28,0)in leading;lc=leading[(0,28,0)];E=scale(h,1/lc)
 hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';hp=hb/'certificate.json';assert sha(hp)==load(hb/'report.json')['evidence_sha256'][hp.name];H=poly(load(hp)['monic_H_laurent']);assert {(m[0],m[1],0):c for m,c in H.items()if m[2]==4}==one
 basis=[[i,j]for j in range(4)for i in range(28)];assert len(basis)==112 and len(set(map(tuple,basis)))==112
 cert={'base_ring':'QQ[a,a^-1]','coefficient_variable_order':['a','e','d'],'monic_E_e28':ser(E),'monic_H_d4':ser(H),'E_from_h_QQ_unit':str(1/lc),'E_lce_equals1':True,'H_lcd_equals1':True,'basis_e_i_d_j':basis,'basis_convention':'[i,j] means e^i*d^j','exact_R_factorization_without_irreducible_factorization':'R=R_QQ_content*a^31*C^4*(h/100000000)','R_QQ_content':n['R_QQ_content'],'h_QQ_content':'100000000','h_primitive_integer':ser(prim),'all31_removed_factors_exactly_checked':True,'whole_R_factorization_reconstructed':True,'finite_free_rank112_proof':'For R0=QQ[a,a^-1], E monic e28 gives R0[e]/E free28. Over this ring H remains monic d4, so R0[e,d]/(E,H) is free112, irrespective of reducedness/factorization.','actual_B_bad_containment':'At each FULLU B-rankbad point all maximalminors and hence f,g vanish. H is their a-unit linearcombination, so H=0; R=uf*f+vg*g=0 and the knownunit a^31*C^4 implies h=0. Therefore every B-bad point lies in V(E,H).','FULLU_inverse_NOT_adjoined_to_cover_yet':True,'not_badexistence_or_globalB3':True,'no_new_determinant_resultant_GB_point_prime_jobs':True}
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),nb/'report.json',np,hb/'report.json',hp]
 report={'id':ID,'E_terms':len(E),'E_degree_e':28,'E_degree_a':max(m[0]for m in E),'H_terms':len(H),'H_degree_d':4,'finite_free_cover_rank':112,'actual_B_bad_contained':True,'FULLU_localization_finite_NOT_asserted':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-t,'frozen':True,'independent_review_pending':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
