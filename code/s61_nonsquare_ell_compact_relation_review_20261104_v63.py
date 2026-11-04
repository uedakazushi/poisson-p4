#!/usr/bin/env python3
"""Independent literal full319 kernel and outside-known37 polynomial audit."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def writepoly(f,p):
 p={tuple(e):Q(a)for e,a in p.items()if a};f.write(str(len(p))+''.join(' '+str(a)+' '+' '.join(map(str,e))for e,a in sorted(p.items()))+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';r=read(folder/'report.json');c=read(folder/'certificate.json')
 assert r['all_literal_source_identity_and_outside37_gate_checks_passed']and r['frozen']and r['uniform_pointwise_independent_kernel_directions_primary']==38
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in r['artifact_sha256'].items():assert sha(folder/p)==h,p
 cp=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';d=read(cp);defs=[d['source_rows'][i]for i in d['retained_row_indices']]
 assert defs==c['source319_definitions']and len(defs)==319
 m2=[x[2]for x in defs[:29]];assert defs[:174]==[['quadric',i,m]for i in range(6)for m in m2]
 sp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';src=read(sp);m3=src['normal_monomials'];m4=d['weight4_columns']
 bp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';basis=read(bp);ids=basis['remaining_original_row_indices']
 assert ids.index(43)==42 and ids.index(49)==47 and src['labels'][43]==['cubic',20]and src['labels'][49]==['cubic',61]
 raw=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt';assert sha(raw)==basis['remaining_matrix_sha256']
 for mons,degree in[(m2,2),(m3,3),(m4,4)]:
  assert len(mons)==len(set(map(tuple,mons)))and all(len(m)==8 and sum(m)+m[5]==degree for m in mons)
 im2={tuple(m):j for j,m in enumerate(m2)}
 with(out/'independent-source-metadata.txt').open('w')as f:
  f.write('29 91 239\n')
  for mons in[m2,m3,m4]:
   for m in mons:f.write(' '.join(map(str,m))+'\n')
  f.write('319\n')
  for row in defs:
   if row[0]=='quadric':v=[0,row[1],im2[tuple(row[2])]]
   elif row[0]=='cubic_normalized':v=[1,row[1],row[3]]
   else:assert row[0]=='quartic';v=[2,row[1],0]
   f.write(' '.join(map(str,v))+'\n')
 assert(out/'independent-source-metadata.txt').read_text().split()==(folder/'source-metadata.txt').read_text().split()
 qp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';qd=read(qp)
 with(out/'original-six-quadrics.txt').open('w')as f:
  f.write('6 29\n')
  for row in qd['quadratic_rows_sparse']:
   tab={tuple(m):{tuple(e):Q(a)for e,a in terms}for m,terms in row}
   for m in m2:writepoly(f,tab.get(tuple(m),{}))
 # Reuse only accepted actualD15/R4 source proof; bind the exact files
 # to their independent receipts before inspecting24 selected scalar entries.
 r4p=ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/certificate.json';r4=read(r4p)
 p15=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json';assert read(p15)['actual319_source_definitions']==defs
 deps=[ROOT/'reports'/z/'report.json'for z in['S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW','S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW','S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW','S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW','S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in deps)
 for review,p in[(deps[2],r4p),(deps[3],p15)]:assert read(review)['input_sha256'][str(p.relative_to(ROOT))]==sha(p)
 assert r4['parent_certificate_sha256']==sha(p15)
 cols=r4['cubic_source319_indices'];index={v:j for j,v in enumerate(cols)};piv=[v['source319_column']for v in r4['unit_pivot_steps']]
 assert piv==[176,177,237,181]and c['outside_known37_gate']['residual4_first3_unit_pivot_indices']==piv[:3]
 assert [defs[j]for j in piv]==[['cubic_normalized',42,18,2],['cubic_normalized',42,18,3],['cubic_normalized',52,86,0],['cubic_normalized',43,20,0]]
 assert c['known_unit_separator']['source319_indices']==[182,186]and[defs[j]for j in[182,186]]==[['cubic_normalized',43,20,1],['cubic_normalized',43,20,6]]
 a,e,dd=s.symbols('a e d');loc={'a':a,'e':e,'d':dd};selected=[]
 with(out/'residual4-selected-ratios.txt').open('w')as f:
  f.write('4 6\n')
  for i in range(4):
   row=[]
   for j in piv+[182,186]:
    text=r4['post_pivot_cubic_matrix'][i][index[j]];num,den=s.fraction(s.sympify(text,locals=loc));pn=s.Poly(num,a,e,dd,domain=s.QQ);pd=s.Poly(den,a,e,dd,domain=s.QQ)
    writepoly(f,{m:Q(str(v))for m,v in pn.terms()});writepoly(f,{m:Q(str(v))for m,v in pd.terms()});row.append(text)
   selected.append(row)
 cpp=ROOT/'code/s61_nonsquare_ell_compact_relation_review_20261104_v63.cpp';utility=ROOT/'code/s61_nonsquare_ell_Lframe_review_20261104_v63.cpp';exe=out/'review-native'
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 def remaining():return max(.1,60-(time.monotonic()-start))
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=remaining())
 pre=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63';run=subprocess.run([str(exe),str(out),str(folder),str(raw),str(pre)],check=True,capture_output=True,text=True,timeout=remaining());exe.unlink()
 (out/'native.log').write_text(run.stdout+run.stderr);own=json.loads(run.stdout);assert own['all_passed']
 proof={'original_genuine_rows':{'normalized':[43,49],'actual135':[20,61],'raw57_positions':[42,47]},
 'kernel_polynomial_terms':own['kernel_terms'],'separator_source319':[182,186],'separator':'-1024000*a*C^3*J/243',
 'residual4_selected_columns':piv+[182,186],'residual4_selected_entries':selected,
 'functional_recipe':'First subtract acceptedI11 quartic lifts. Then S=sum182+186; subtract S(R4_i)*coordinate[pivot_i] fori0..2. R4_last hasS=0. Thisannihilatesold22+D15; onnewkernel equals knownunitseparator.',
 'pointwise_all_FULLU_rank38':True,'scope':'FULLU=adFBCJ nonzero; no newfactor inversions, no allpointcubic57/fullquartic239, no wholekernel claim',
 'free281_followon_not_computed_by_review':True}
 (out/'independent-certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 note=ROOT/'research/notes'/(ID+'.md')
 result=dict(id=ID,all_passed=True,actual_original_q6_cubic43_49_and_source319_bindings_checked=True,
  direct_literal_full319_kernel_identity_all239_QQ_polynomial_coefficients=True,
  compact_original_cubic91_and_six_quad_source_changes_checked=True,
  all_cofactors_QQ_core_polynomials_no_unknown_factor_inversions=True,
  all33_quartic_cofactors_zero_and_source123_nonzeros19694_terms=True,
  separate_outside_known37_gate_and_FULLU_unit_separator_accepted=True,
  uniform_pointwise_independent_kernel_directions=38,remaining_for80=42,
  generic_rank10_already_accepted_not_recomputed=True,
  no_allpoint_cubic57_fullquartic239_wholekernel_or_free281_construction_claim=True,
  no_new_point_prime_GB_solver_pivot_choice_or_author_rerun=True,native=own,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),cpp,utility,note,folder/'report.json',folder/'certificate.json',folder/'full319-kernel.txt',folder/'compact-polynomials.txt',qp,sp,cp,bp,raw,r4p,p15,*deps]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
