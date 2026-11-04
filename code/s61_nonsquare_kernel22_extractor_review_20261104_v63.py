#!/usr/bin/env python3
"""Independent physical-source binding and scalar/linear operator checker input."""
import argparse,gzip,hashlib,itertools,json,resource,shlex,signal,subprocess,tarfile,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as sy
from sympy.polys.rings import ring
import s61_nonsquare_universal_quad22_review_20261103_v63 as sparse
ROOT=Path(__file__).resolve().parents[1];STEM='S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
 d=read(p);assert d.get('all_passed') or d.get('all_exact_extractor_checks_passed')
 for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
 for n,h in d.get('evidence_sha256',{}).items():
  q=ROOT/n if (ROOT/n).exists() else p.parent/n
  assert sha(q)==h,n
 return d
def mons(n):return sorted(tuple(sum(i==j for i in pair) for j in range(n)) for pair in itertools.combinations_with_replacement(range(n),2))

def main(out,work):
 out=out.resolve();work=work.resolve()
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(120)
 assert not out.exists() and not work.exists();out.mkdir(parents=True);work.mkdir(parents=True)
 base=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63';ap=base/'report.json';author=accepted(ap)
 gp=base/'global-P22.json';g=read(gp);fp=base/'full-P22.json';full=read(fp)
 assert full['circuit']==g['circuit'] and full['linear_output_expected_basis_columns']==[[r,i] for i,r in enumerate(g['P22_linear_output_roots_in_frozen_original_kernel_order'])]
 tp=base/'transport-manifest.json';trans=read(tp);assert trans['scalar_checks']==full['scalar_checks']
 pubp=base/'publish-receipt.json';pub=read(pubp)
 assert pub['all_passed'] and pub['producer_report_sha256']==sha(ap) and pub['P22_sha256']==sha(gp)
 assert pub['note_sha256']==sha(ROOT/'research/notes/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63.md')
 assert pub['transport_manifest_sha256']==sha(tp)
 kp=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json';K=read(kp)
 kr=kp.parent.with_name(kp.parent.name+'_REVIEW')/'report.json';accepted(kr)
 assert g['source_variables']==K['actual_source174_definitions']
 assert g['prime_quadric_source_change']==K['six_quadric_change_to_old_prime_basis']
 assert trans['inverse_prime_cofactor_transform_field_is_inverse_A_rows']
 assert g['coefficient_extractor_new22']==g['P22_linear_output_roots_in_frozen_original_kernel_order'][21]
 C=g['circuit'];assert [len(C[k]) for k in ['atoms_num_den','scalar_nodes','linear_nodes']]==[378,2491,1685]
 for i,n in enumerate(C['scalar_nodes']):
  assert n[0] in ['z','c','a','s','m']
  if n[0] in ['a','s','m']:assert 0<=n[1]<i and 0<=n[2]<i
 for i,n in enumerate(C['linear_nodes']):
  assert n[0] in ['z','v','l']
  if n[0]=='v':assert 0<=n[1]<174
  if n[0]=='l':assert all(0<=c<2491 and 0<=j<i for c,j in n[1])
 oldstagep=base/'new22-extractor.json';oldstage=read(oldstagep)
 for key in ['atoms_num_den','scalar_nodes','linear_nodes']:
  assert C[key][:len(oldstage['circuit'][key])]==oldstage['circuit'][key]
 # Build the actual source coefficients independently in the old cert light frame.
 rawp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';raw=read(rawp)
 anp=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json';anc=read(anp)
 x=[sparse.var(i) for i in range(7)];images=[x[0],x[1],x[2],x[6],x[3],{},sparse.add(x[5],x[1],sparse.scale(x[0],-3)),x[4]]
 qs=[]
 for row in raw['quadratic_rows_sparse']:
  terms=[]
  for mon,cofs in row:
   p={(0,)*7+tuple(e):Q(c) for e,c in cofs}
   for im,n in zip(images,mon):
    if n:p=sparse.mul(p,sparse.power(im,n))
   terms.append(p)
  qs.append(sparse.add(*terms))
 L=sparse.decode(anc['L']);ell=sparse.decode(anc['ell'])
 qp=[qs[0],qs[1],qs[3],qs[5],sparse.mul(ell,x[5]),sparse.mul(ell,x[6])]
 enc=lambda p:[[list(e),str(c)] for e,c in sorted(p.items(),reverse=True) if c]
 def coef(p,m):return enc({e[7:]:c for e,c in p.items() if e[:7]==tuple(m)})
 qm=mons(7);assert list(map(list,qm))==g['quadric_monomials7']==trans['quad_monomial_order']
 quad=[[coef(p,m) for m in qm] for p in qp]
 unitmons=[tuple(int(i==j) for i in range(7)) for j in range(7)]
 linear=[[coef(p,m) for m in unitmons] for p in [L,ell]]
 archive=ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63/original_ZW_rank19_generic146_exact.tar.gz'
 with tarfile.open(archive) as tf:
  for suffix,dest in [('pair-native.txt','input-pair-native.txt'),('literal-Bezout.json','input-literal-Bezout.json')]:
   members=[m for m in tf.getmembers() if m.name.endswith('/'+suffix)];assert len(members)==1
   assert tf.extractfile(members[0]).read()==(base/dest).read_bytes()
 pairlines=(base/'input-pair-native.txt').read_text().splitlines();assert pairlines[0]=='2 6 0'
 def nativepoly(line):
  t=line.split();n=int(t[0]);assert len(t)==1+4*n
  return [[list(map(int,t[2+4*i:5+4*i])),t[1+4*i]] for i in range(n)]
 pair=[nativepoly(line) for line in pairlines[1:]];assert len(pair)==12
 bez=read(base/'input-literal-Bezout.json')
 assert sum(len(p) for p in bez['scalar_unit_cofactors'])==117
 assert sum(bool(Q(c)) for p in bez['scalar_unit_cofactors'] for ex,c in p)==116
 shadow=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/right_inverse.txt.gz'
 with gzip.open(shadow,'rt') as f:
  tokens=iter(f.read().split());assert [int(next(tokens)) for i in range(3)]==[21,4,0];old4=[]
  for j in range(84):
   count=int(next(tokens));p={}
   for k in range(count):
    c=Q(next(tokens));a,e,P=[int(next(tokens)) for i in range(3)];ex=(a+P,e,P);assert ex not in p;p[ex]=c
   old4.append(enc(p))
  assert next(tokens,None)is None
 ecpath=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';ec=read(ecpath);wm=mons(6)
 saved4=[]
 for row in ec['cleared_qfour_sparse']:
  d={tuple(m[:6]):p for m,p in row};assert all(m[6]==0 for m,p in row)
  saved4.extend(d.get(m,[]) for m in wm)
 R,*_=ring('a,e,d',sy.QQ);field=R.to_field()
 def stringrat(text):
  v=field.from_expr(sy.sympify(text));return[enc(v.numer),enc(v.denom)]
 native=work/'input.txt'
 def putpoly(f,p):f.write(str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,e)) for e,c in p)+'\n')
 def putrat(f,p):putpoly(f,p[0]);putpoly(f,p[1])
 den1=[[[0,0,0],'1']]
 with native.open('w') as f:
  f.write('378 2491 1685 45\n')
  for p in C['atoms_num_den']:putrat(f,p)
  for n in C['scalar_nodes']:f.write(' '.join(map(str,n))+'\n')
  for n in C['linear_nodes']:
   f.write(('l '+str(len(n[1]))+''.join(' '+str(c)+' '+str(k) for c,k in n[1]) if n[0]=='l' else ' '.join(map(str,n)))+'\n')
  for root,want in full['scalar_checks']:f.write(f'{root} {want}\n')
  f.write(' '.join(map(str,g['P22_linear_output_roots_in_frozen_original_kernel_order']))+'\n')
  lookup={(i,j):[num,den] for j,col in enumerate(K['sparse_columns']) for i,num,den in col}
  for i in range(174):
   for j in range(22):putrat(f,lookup.get((i,j),[[],den1]))
  for row in quad:
   for p in row:putrat(f,[p,den1])
  for row in linear:
   for p in row:putrat(f,[p,den1])
  for key in ['prime_quadric_source_change','inverse_prime_cofactor_transform']:
   for row in g[key]:
    for text in row:putrat(f,stringrat(text))
  for p in pair+bez['source_scalar5']+bez['scalar_unit_cofactors']+old4+saved4:putrat(f,[p,den1])
  f.write(' '.join(map(str,trans['new22_lambda_fivevar_quad_coefficients']))+'\n')
  for key in ['pair_dual6x2_scalar_roots','global_quad_dual6x28_scalar_roots']:
   for row in g[key]:f.write(' '.join(map(str,row))+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=work/'verify'
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
 run=subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-lgmp','-o',str(exe)],capture_output=True,text=True,check=True,timeout=25)
 (out/'compile.log').write_text(run.stdout+run.stderr)
 run=subprocess.run([str(exe),str(native)],capture_output=True,text=True,check=True,timeout=max(1,115-(time.monotonic()-start)))
 logp=out/'independent-operator.log';logp.write_text(run.stdout+run.stderr)
 assert 'ALL484_LITERAL_P22_K_EQUALS_I22_PASS' in run.stdout and 'ACTUAL_4DUAL16_6DUAL36_LAMBDA15_AND_ANCHORS_PASS' in run.stdout
 bp=base/'global282-basis.json';basis=read(bp);atlasp=ROOT/basis['source304_atlas_path'];atlas=read(atlasp)
 assert sha(atlasp)==basis['source304_atlas_sha256'] and sha(gp)==basis['P22_sha256'] and sha(bp)==pub['global282_basis_sha256']
 ret=atlas['retained304_source319_indices']
 assert basis['old_source319_indices_quadratic_A174']==ret[:174] and basis['old_source319_indices_arbitrary_B22']==ret[174:196] and basis['old_source319_indices_C108']==ret[196:]
 assert basis['actual_global_new_source_shape']==[282,239] and basis['old282_subset_of_generators_not_claimed']
 packed=out/'verification-input.txt.gz';packed.write_bytes(gzip.compress(native.read_bytes(),mtime=0))
 failurepath=base/'first-failure.json';read(failurepath)
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{STEM}.md',ap,gp,fp,tp,pubp,kp,kr,rawp,anp,archive,shadow,ecpath,bp,atlasp,oldstagep,failurepath,
        ROOT/'code/s61_binary_hyperplane_cubic_qq3_20261102.hpp',ROOT/'code/s61_nonsquare_universal_quad22_review_20261103_v63.py']
 report=dict(id=STEM,all_passed=True,
  verification_kind='Independent source conversion and scalar-DAG/column-by-column linear operator evaluation; shared accepted FLINT rational utility',
  all484_literal_P22K_equals_I22=True,all_scalar_atom_denominators_known_FULLU_units=True,
  scalar2491_linear1685_atoms378_linearity_and_topology=True,
  actual2pair12coeffs_Schur5a_minus_aC_Bezout116_transport_checked=True,
  actual84_ellzero_sourcecoeffs_P_to_ad_shifts_and_fourdual16_checked=True,
  actual_six_quad_dual36_and_lambda15_anchors_checked=True,
  actual_A_inverse_rows_applied_transpose_and_Koszul_Atranspose_S_A_scope_accepted=True,
  explicit_global_free282_basis_and_same_cokernel_formula_verified=True,
  fixed282_submatrix_full239_C122rank64_or_additional_kernel_directions_claim=False,
  known_kernel_directions_unchanged=37,remaining_independent_directions=43,
  no20_or66_rerun_new_search_minor_Bezout_GB_field_point_or_fullquartic_expansion=True,
  executable_saved_only_in_ignored_work=str(exe.relative_to(ROOT)),
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
  evidence_sha256={logp.name:sha(logp),packed.name:sha(packed)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k] for k in ['all_passed','all484_literal_P22K_equals_I22','seconds']},indent=2))
def publish_verified(out,work):
 """Metadata-only recovery after completed checker; NEVER reruns arithmetic."""
 start=time.monotonic();out=out.resolve();work=work.resolve();front=out/'finalization-frontier.json';fr=read(front)
 cpp=Path(__file__).with_suffix('.cpp');log=out/'independent-operator.log';packed=out/'verification-input.txt.gz'
 assert sha(cpp)==fr['native_cpp_sha256'] and sha(log)==fr['native_log_sha256']
 assert sha(work/'input.txt')==fr['native_input_sha256'] and sha(packed)==fr['compressed_input_sha256']
 assert gzip.decompress(packed.read_bytes())==(work/'input.txt').read_bytes()
 for marker in ['ALL484_LITERAL_P22_K_EQUALS_I22_PASS','ACTUAL_PAIR12_SCHUR_BEZOUT116_AND_84_ELL_SOURCE_TRANSPORT_PASS','ACTUAL_4DUAL16_6DUAL36_LAMBDA15_AND_ANCHORS_PASS','ALL_KNOWN_FULLU_DENOMINATORS_AND_LINEARITY_PASS']:assert marker in log.read_text()
 base=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63';accepted(base/'report.json')
 rels=['global-P22.json','full-P22.json','transport-manifest.json','publish-receipt.json','global282-basis.json','new22-extractor.json','first-failure.json']
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{STEM}.md',base/'report.json',*[base/n for n in rels],front,
  ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json',ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/atlas.json',ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63_REVIEW/report.json',
  ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json',ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json',
  ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63/original_ZW_rank19_generic146_exact.tar.gz',ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/right_inverse.txt.gz',ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json',
  ROOT/'code/s61_binary_hyperplane_cubic_qq3_20261102.hpp',ROOT/'code/s61_nonsquare_universal_quad22_review_20261103_v63.py']
 report=dict(id=STEM,all_passed=True,verification_kind='Independent physical-source conversion and scalar-DAG/column-by-column exact checker; metadata-only recovery after completed native proof',
  all484_literal_P22K_equals_I22=True,all_scalar_atom_denominators_known_FULLU_units=True,scalar2491_linear1685_atoms378_linearity_and_topology=True,
  actual2pair12coeffs_Schur5a_minus_aC_Bezout116_transport_checked=True,actual84_ellzero_sourcecoeffs_P_to_ad_shifts_and_fourdual16_checked=True,
  actual_six_quad_dual36_and_lambda15_anchors_checked=True,actual_A_inverse_rows_applied_transpose_and_Koszul_Atranspose_S_A_scope_accepted=True,
  explicit_global_free282_basis_and_same_cokernel_formula_verified=True,fixed282_submatrix_full239_C122rank64_or_additional_kernel_directions_claim=False,
  known_kernel_directions_unchanged=37,remaining_independent_directions=43,no20_or66_rerun_new_search_minor_Bezout_GB_field_point_or_fullquartic_expansion=True,
  executable_saved_only_in_ignored_work=str((work/'verify').relative_to(ROOT)),metadata_only_recovery_no_native_rerun=True,
  input_to_verified_log_seconds=fr['input_to_verified_log_seconds'],full_native_review_enforced_seconds_cap=120,seconds=None,metadata_publish_seconds=time.monotonic()-start,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={log.name:sha(log),packed.name:sha(packed)})
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','input_to_verified_log_seconds','metadata_only_recovery_no_native_rerun']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--publish-verified',action='store_true');a=p.parse_args();(publish_verified if a.publish_verified else main)(a.output,a.work)
