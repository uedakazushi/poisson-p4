#!/usr/bin/env python3
"""Independent source/frame/monic-square audit; no Hilbert or rank computation."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def weighted(degree):
 weights=[1]*7+[2];ans=[]
 def rec(j,left,m):
  if j==8:
   if left==0:ans.append(m)
   return
  for k in range(left//weights[j]+1):rec(j+1,left-weights[j]*k,m+[k])
 rec(0,degree,[]);return sorted(ans,reverse=True)
def wp(f,p):
 p={tuple(m):Q(v)for m,v in p.items()if v};f.write(str(len(p))+''.join(' '+str(v)+' '+' '.join(map(str,m))for m,v in sorted(p.items()))+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';r=read(b/'report.json');p=read(b/'preflight.json');cir=read(b/'circuit.json')
 assert r['frozen']and r['all_exact_joint_coordinate_and_quadric_checks_passed']and r['invalid_exact_four_monomial_ideal_rejected']
 for f,h in r['input_sha256'].items():assert sha(ROOT/f)==h,f
 for f,h in r['evidence_sha256'].items():assert sha(b/f)==h,f
 cp=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';c319=read(cp);defs=[c319['source_rows'][j]for j in c319['retained_row_indices']];oldm2=[v[2]for v in defs[:29]]
 assert defs[:174]==[['quadric',i,m]for i in range(6)for m in oldm2]
 m2,m3,m4=[weighted(D)for D in[2,3,4]];assert list(map(len,[m2,m3,m4]))==[29,91,239]
 with(out/'independent-monomials.txt').open('w')as f:
  f.write('29 91 239\n')
  for mons in[oldm2,m2,m3]:
   for m in mons:f.write(' '.join(map(str,m))+'\n')
 assert(out/'independent-monomials.txt').read_text().split()==(b/'monomial-input.txt').read_text().split()
 qp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';qd=read(qp)
 with(out/'original-six-quadrics.txt').open('w')as f:
  f.write('6 29\n')
  for row in qd['quadratic_rows_sparse']:
   lookup={tuple(m):{tuple(e):Q(v)for e,v in terms}for m,terms in row}
   for m in oldm2:wp(f,lookup.get(tuple(m),{}))
 assert p['new_light_order']==['L','ell','x0','x3','r','Z','W']and p['old_light_order']==['x0','x1','x2','x3','x4','Z','W']
 lead=[[2,0,0,0,0,0,0,0],[1,0,1,0,0,0,0,0],[0,1,0,0,0,1,0,0],[0,1,0,0,0,0,1,0]]
 # Independent support reader. Scalar coefficients are fully checked by
 # inverse-substitution native; here their supports verify each triangular row.
 toks=iter((b/'six-prime-quadrics-numerators.txt').read_text().split());assert[int(next(toks)),int(next(toks))]==[6,29];supports=[]
 for i in range(6):
  sup=[]
  for j in range(29):
   n=int(next(toks))
   if n:sup.append(m2[j])
   for _ in range(4*n):next(toks)
  supports.append(sup)
 assert next(toks,None)is None
 independent=[]
 for D,mons,st in zip([3,4],[m3,m4],p['stages']):
  cofs=weighted(D-2);rows=[[i,m]for i in range(4)for m in cofs];first={}
  for row,(i,cof)in enumerate(rows):
   lm=tuple(x+y for x,y in zip(lead[i],cof));first.setdefault(lm,row)
  monindex={tuple(m):j for j,m in enumerate(mons)};pivotcols=[monindex[m]for m in sorted(first,reverse=True)];selected=[first[tuple(mons[j])]for j in pivotcols]
  othercols=[j for j in range(len(mons))if j not in pivotcols];overlap=[j for j in range(len(rows))if j not in selected]
  assert st['weighted_degree']==D and st['target_monomials']==mons and st['four_quadric_product_rows']==rows
  assert st['monic_pivot_source_rows']==selected and st['monic_pivot_target_columns']==pivotcols and st['retained_target_columns']==othercols
  assert st['all_unselected_overlap_rows_retained']==overlap and st['remaining_q0_q5_products']==[[i,m]for i in[4,5]for m in cofs]
  for row,col in zip(selected,pivotcols):
   i,cof=rows[row];lm=mons[col]
   assert [x+y for x,y in zip(lead[i],cof)]==lm
   for term in supports[i]:
    target=tuple(x+y for x,y in zip(term,cof));assert target in monindex and monindex[target]>=col
  independent.append(dict(degree=D,pivots=len(selected),target_remaining=len(othercols),all_overlap_rows=overlap,quad_pair_rows=2*len(cofs),all_unit_triangular_supports_checked=True))
 assert[independent[i]['pivots']for i in range(2)]==[26,98]and independent[0]['all_overlap_rows']==[7,26]
 assert len(independent[1]['all_overlap_rows'])==18
 # Row26=p4*Z is the exact-zero cubic overlap; row7=p2*L is the S overlap.
 assert p['stages'][0]['four_quadric_product_rows'][26]==[3,[0,0,0,0,0,1,0,0]]
 assert p['stages'][0]['four_quadric_product_rows'][7]==[1,[1,0,0,0,0,0,0,0]]
 cubic_after=2+14+16-1;quartic_after=18+58+112+33
 assert cubic_after==31 and quartic_after==221 and 319-98==221 and 239-98==141
 cpp=ROOT/'code/s61_nonsquare_joint_Lell_preflight_review_20261104_v63.cpp';utility=ROOT/'code/s61_nonsquare_ell_Lframe_review_20261104_v63.cpp';exe=out/'review-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();rem=lambda:max(.1,60-(time.monotonic()-began))
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=rem())
 pre=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63';run=subprocess.run([str(exe),str(out),str(b),str(pre)],capture_output=True,text=True,timeout=rem());exe.unlink();(out/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0,run.stderr;native=json.loads(run.stdout);assert native['all_passed']
 deps=[ROOT/'reports'/j/'report.json'for j in['S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW','S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW','S61_NONSQUARE_CUBIC_QUARTIC_STRUCTURAL_REDUCTIONS_20261103_v63_REVIEW']]
 assert all(read(x)['all_passed']for x in deps)
 (out/'independent-certificate.json').write_text(json.dumps({'stages':independent,'source319_quad_only_reframe':True,'retained_old_genuine_light_multipliers_and_all33_quartics':True,'cubic31_target65':True,'fullquartic221_target141':True,'not_a_Hilbert_function_or_four_monomial_ideal_quotient':True},indent=2)+'\n')
 note=ROOT/'research/notes'/(ID+'.md')
 result=dict(id=ID,all_passed=True,full_original_q6_and_joint_L_ell_frame_rebuilt=True,
  all174_prime_quadric_inverse_polynomial_identities_verified=True,
  all98_coordinate_inverse_identities_and_literal_Jacobian_unit_verified=True,
  actual_ZW_corrections82_and_nonzero_S72_retained=True,
  exact26_cubic98_quartic_unit_triangular_pivots_verified=True,
  all2_cubic18_quartic_overlaps_accounted_for_and_only_exact_cubic_zero_discarded=True,
  actual_fullquartic319_to221_source_and239_to141_target_Schur_recipe_accepted=True,
  no_four_monomial_ideal_Hilbert_rank141_full239_or_free183_claim=True,
  no_G16_quartic_coefficient_expansion_newpoint_GB_solver_or_author_rerun=True,native=native,
  input_sha256={str(x.relative_to(ROOT)):sha(x)for x in[Path(__file__),cpp,utility,note,b/'report.json',b/'preflight.json',b/'circuit.json',qp,cp,*deps]},
  evidence_sha256={x.name:sha(x)for x in out.iterdir()if x.is_file()},seconds=time.monotonic()-began,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
