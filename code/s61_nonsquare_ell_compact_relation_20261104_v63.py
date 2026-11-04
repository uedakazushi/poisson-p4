#!/usr/bin/env python3
"""Extract the compact saved relation and check its literal original319 identity."""
import argparse,hashlib,json,resource,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(p,m,c):
 z=p.get(m,Q())+c
 if z:p[m]=z
 else:p.pop(m,None)
def wp(f,p):f.write(str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in sorted(p.items(),reverse=True))+'\n')
def main(out):
 started=time.monotonic();deadline=started+60;resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 pre=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63';border=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63'
 pr=read(pre/'report.json');br=read(border/'report.json');assert pr['all_source_bound_coordinate_and_polynomial_input_checks_passed']and br['generic_rank_exactly10_proved']
 for name in ['genuine16x40.txt','L-coefficients.txt','native.json']:assert sha(pre/name)==pr['artifact_sha256'][name]
 n=read(pre/'native.json');assert n['genuine_C_factors_cancelled'][1]==n['genuine_C_factors_cancelled'][6]==1
 qpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';q=read(qpath)
 spath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';src=read(spath);m3=src['normal_monomials']
 cpath=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';C319=read(cpath);defs=[C319['source_rows'][i]for i in C319['retained_row_indices']];m4=C319['weight4_columns'];m2=[row[2]for row in defs[:29]];mi2={tuple(m):j for j,m in enumerate(m2)}
 assert defs[:174]==[['quadric',i,m]for i in range(6)for m in m2]
 basispath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';basis=read(basispath);ids=basis['remaining_original_row_indices'];assert ids[42]==43 and ids[47]==49
 cubicids=list(dict.fromkeys(row[1]for row in defs if row[0]=='cubic_normalized'));assert cubicids==ids[41:]
 assert src['labels'][43]==['cubic',20]and src['labels'][49]==['cubic',61]
 meta=out/'source-metadata.txt'
 with meta.open('w')as f:
  f.write('29 91 239\n')
  for ms in [m2,m3,m4]:
   for m in ms:f.write(' '.join(map(str,m))+'\n')
  f.write('319\n')
  for row in defs:
   if row[0]=='quadric':v=[0,row[1],mi2[tuple(row[2])]]
   elif row[0]=='cubic_normalized':v=[1,row[1],row[3]]
   else:assert row[0]=='quartic';v=[2,row[1],0]
   f.write(' '.join(map(str,v))+'\n')
 qfile=out/'six-original-quadrics.txt'
 with qfile.open('w')as f:
  f.write('6 29\n')
  for row in q['quadratic_rows_sparse']:
   tab={tuple(m):{tuple(ex):Q(c)for ex,c in ts}for m,ts in row}
   for m in m2:wp(f,tab.get(tuple(m),{}))
 anchorpath=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json';anchor=read(anchorpath);ell=[{}for _ in range(8)];old=[0,1,2,4,7]
 for m,ts in anchor['ell']:
  assert sum(m)==1;j=m.index(1)
  images=[(old[j],1)]if j<5 else([(0,3),(1,-1),(6,1)]if j==5 else[(3,1)])
  for v,c in images:
   for ex,k in ts:add(ell[v],tuple(ex),Q(k)*c)
 ellfile=out/'ell-original8.txt'
 with ellfile.open('w')as f:
  f.write('8 1\n')
  for p in ell:wp(f,p)
 # Only literal small accepted residual4 coefficients are read; no operator expansion.
 rbase=ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63';R4=read(rbase/'certificate.json');assert read(rbase/'report.json')['new_known_unit_directions']==4
 assert read(ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW/report.json')['all_passed']
 cubiccols=R4['cubic_source319_indices'];look={(defs[g][1],defs[g][3]):g for g in cubiccols};piv=[row['source319_column']for row in R4['unit_pivot_steps']]
 assert [defs[i]for i in piv]==[['cubic_normalized',42,18,2],['cubic_normalized',42,18,3],['cubic_normalized',52,86,0],['cubic_normalized',43,20,0]]
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();frac=lambda c:K.from_expr(s.sympify(c));ix={g:j for j,g in enumerate(cubiccols)};mat=R4['post_pivot_cubic_matrix']
 assert all(frac(mat[i][ix[piv[j]]])==(K.one if i==j else K.zero)for i in range(4)for j in range(4))
 gc20=[look[(43,j)]for j in [1,6]];assert sum((frac(mat[3][ix[k]])for k in gc20),K.zero)==K.zero
 bparent=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json';Bp=read(bparent);assert Bp['actual319_source_definitions']==defs
 breview=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW/report.json';assert read(breview)['all11_final_vectors_have_literal_identity_quartic_block']
 source=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt';assert sha(source)==basis['remaining_matrix_sha256']
 cpp=ROOT/'code/s61_nonsquare_ell_compact_relation_20261104_v63.cpp';exe=out/'compact-relation';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(.1,deadline-time.monotonic()))
 with(out/'native.log').open('w')as f:subprocess.run([str(exe),str(source),str(meta),str(qfile),str(ellfile),str(pre/'genuine16x40.txt'),str(pre/'L-coefficients.txt'),str(out)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(.1,deadline-time.monotonic()))
 exe.unlink();native=read(out/'native.json');assert native['literal_actual319_identity_checked_all239']
 proof={'source319_definitions':defs,'source_genuine_basis_positions':[42,47],'source_original_normalized_ids':[43,49],'source_actual135_cubic_rows':[20,61],
 'coordinate_names':'Z=3*v0-v1+v6; W=v3; L is independent of Z,W in (v0,v1,v2,v4,v7)',
 'delta':'-64000*C^2','E_definition':'G6^L-G1^L=L^2*(d0*L+d1*x0+d2*y1+d3*y2+d4*y3); E is this linear form pulled back by L=sum L_j*x_j',
 'literal_full_cubic_identity':'delta*(G6-G1)-8100*q1*E=Z*A+W*B',
 'literal_source_change':['a*ell*Z=(10*a+3*e^2-3*e)*q1-3*(e-1)*q3-10*a*q4','ell*W=-(d*e+8)*q1+8*q2+d*q3'],
 'full319_kernel_recipe':{'q1':'-8100*a*ell*E-(10*a+3*e^2-3*e)*A+a*(d*e+8)*B','q2':'-8*a*B','q3':'3*(e-1)*A-a*d*B','q4':'10*a*A','other_quadric_cofactors':'0','G20_normalized43':'-a*delta*ell','G61_normalized49':'a*delta*ell','all33_quartic_cofactors':'0'},
 'new_coefficients_are_polynomials_over_QQ':True,'unknown_core_factors_or_heavy_minors_inverted':False,
 'known_unit_separator':{'source319_indices':gc20,'sum':'-1024000*a*C^3*J/243','domain':'FULLU=a*d*F*B*C*J!=0'},
 'outside_known37_gate':{'quartic_projection_Bianchi15_rank11':'accepted global11 unit pivots; every quartic-zero combination is in residual4',
 'residual4_first3_unit_pivot_indices':piv[:3],'new_vector_at_first3_pivots':[0,0,0],
 'residual4_last_G20_v1_plus_v6':0,'new_G20_v1_plus_v6':'known unit above',
 'argument':'If c*new is a combination of old22 and direct15, project to quartics then reduce the11 unit pivots. At first3 cubic pivot coordinates new=0 forces residual coefficients0,0,0; last residual has zero G20(v1+v6), hence c=0. Old37 are independent, so adjoining new gives38 pointwise independent directions on FULLU.',
 'no_full239_or_all_kernel_claim':True},
 'uniform38_primary_proof_complete_independent_review_pending':True,'generic_canonical11_rank_exact10_reused':True,
 'no_new_point_prime_GB_solver_pivot_or_quartic_source_expansion':True,
 'native':native}
 (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/'code/s61_nonsquare_ell_Lframe_preflight_20261104_v63.cpp',ROOT/'research/notes'/(ID+'.md'),pre/'report.json',pre/'native.json',pre/'genuine16x40.txt',pre/'L-coefficients.txt',border/'report.json',source,qpath,spath,cpath,basispath,anchorpath,rbase/'report.json',rbase/'certificate.json',rbase.parent/(rbase.name+'_REVIEW')/'report.json',bparent,breview]
 report={'id':ID,'all_literal_source_identity_and_outside37_gate_checks_passed':True,'uniform_pointwise_independent_kernel_directions_primary':38,'remaining_for80':42,
 'independent_review_pending':True,'all_new_denominators_only_QQ_constants':True,'no_new_core_factor_inverse':True,'no_full239_or_allpoint_cubic57_claim':True,
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'artifact_sha256':{str(p.relative_to(out)):sha(p)for p in out.iterdir()if p.is_file()},
 'native':native,'seconds':time.monotonic()-started,'combined_artifact_bytes':sum(p.stat().st_size for p in out.iterdir()if p.is_file()),'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','artifact_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
