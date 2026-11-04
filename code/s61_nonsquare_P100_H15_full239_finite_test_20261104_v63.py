#!/usr/bin/env python3
"""ONE point-free full-H15 source/NF/rank job with streamed coefficients."""
import argparse,gzip,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63';P=32003
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists();work.mkdir(parents=True,exist_ok=False);phase=read(out/'phase.json');assert phase['phase']=='implementation_persisted_before_any_coefficient_or_rank_job';paths=[]
 def load(stem,name='certificate.json'):
  p=ROOT/'reports'/stem/name;paths.append(p);return read(p)
 pr=load('S61_NONSQUARE_P100_H15_FULL239_SOURCE_PREFLIGHT_20261104_v63','report.json');assert sha(paths[-1])=='958a02292e9a88ea4151eb3c077f4198b8e5aaf84dcc041fdcabb1cfea691109'
 factor=load('S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63','native.json');assert factor['all_passed']and factor['Hbar_squarefree']
 oracle=load('S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63','oracle.json');O=oracle['arithmetic'];ON=O['nodes']
 fd=load('S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63');nf=fd['full_target_projection239_to122'];LN=fd['linear_nodes'];SN=fd['scalar_circuit']['nodes'];SA=fd['scalar_circuit']['atoms_num_den']
 h=load('S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63');ri=load('S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63','native.json');L=ri['right_inverse112x23']
 frame=load('S61_NONSQUARE_P100_G18_PURE7_20261104_v63','joint-cubic-frame.json');ad=load('S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63','adapter.json')
 cat=load('S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63');catalogue=[cat['source_rows'][i]for i in cat['retained_row_indices']];assert len(catalogue)==319
 qp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';paths.append(qp);quads=read(qp)
 sp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';paths.append(sp);oldmons=read(sp)['normal_monomials'];assert len(oldmons)==91
 bp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';paths.append(bp);br=read(bp);mp=bp.with_name('remaining.txt.gz');paths.append(mp)
 def code(m):return sum(e<<(3*i)for i,e in enumerate(m))
 # Integer metadata is the union of both frames, using ordinary degree<=4.
 # All polynomial VALUES remain sparse in the correct weighted frame (<=239).
 weights=[1]*8;mons=[m for m in itertools.product(range(5),repeat=8)if sum(m)<=4];assert len(mons)==495
 mons.sort(key=lambda m:code(m));moncodes=[code(m)for m in mons]
 joint=h['parent_joint_target239_monomials'];assert len(joint)==239
 # No extra coefficient query: exact dependency closures of frozen DAGs.
 need=set();atoms=set();todo=oracle['genuine_quartic_roots'][:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);n=ON[i]
  if n[0]=='m':todo+=n[1:]
  elif n[0]=='l':
   for c,j in n[1]:atoms.add(c);todo.append(j)
  else:assert n[0]in['z','v']
 oi={i:j for j,i in enumerate(sorted(need))};oa={i:j for j,i in enumerate(sorted(atoms))}
 lneed=set();sroot=set();todo=nf['output122_linear_roots'][:]
 while todo:
  i=todo.pop()
  if i in lneed:continue
  lneed.add(i);n=LN[i]
  if n[0]=='l':
   for c,j in n[1]:sroot.add(c);todo.append(j)
  else:assert n[0]in['z','v']
 sneed=set();sat=set();todo=list(sroot)
 while todo:
  i=todo.pop()
  if i in sneed:continue
  sneed.add(i);n=SN[i]
  if n[0]=='c':sat.add(n[1])
  elif n[0]in['a','s','m']:todo+=n[1:]
  else:assert n[0]=='z'
 si={i:j for j,i in enumerate(sorted(sneed))};li={i:j for j,i in enumerate(sorted(lneed))};sa={i:j for j,i in enumerate(sorted(sat))}
 classes=h['target122_Hdegree_classes'];assert[len(classes[str(k)])for k in [0,1,2]]==[98,23,1]
 actual16=[x['original_equation']for x in ad['basis57_source_definitions'][41:]];assert len(actual16)==16
 c2=code
 with(work/'config.txt').open('w')as f:
  f.write(str(len(mons))+'\n');f.writelines(f'{code(m)} {sum(a*b for a,b in zip(m,weights))}\n'for m in mons)
  for seq in [[code(m)for m in joint],[code(m)for m in oldmons],classes['0'],classes['1'],classes['2']]:f.write(str(len(seq))+'\n'+' '.join(map(str,seq))+'\n')
  f.write('319\n')
  for i,r in enumerate(catalogue):
   if r[0]=='quadric':assert i<174;f.write(f'0 {r[1]} {code(r[2])}\n')
   elif r[0]=='cubic_normalized':assert 174<=i<286;g=actual16.index(r[2]);assert i==174+7*g+[0,1,2,3,4,6,7].index(r[3]);f.write(f'1 {g} {r[3]}\n')
   else:assert r[0]=='quartic'and r[1]==oracle['genuine_quartic_original_rows'][i-286];f.write(f'2 {i-286} 0\n')
  f.write(str(len(factor['factors']))+'\n')
  for x in factor['factors']:assert x['multiplicity']==1;f.write(str(x['degree'])+' '+' '.join(map(str,x['coefficients']))+'\n')
 checks=0
 def qmod(c):
  nonlocal checks
  c=Q(c);assert c.denominator%P;checks+=1;return c.numerator*pow(c.denominator,-1,P)%P
 def poly(f,p):
  f.write(str(len(p))+'\n')
  for m,c in p:f.write(f'{qmod(c)} '+' '.join(map(str,m))+'\n')
 def crf(f,nd):
  for p in nd:poly(f,p)
 def rf(f,r):
  for p in [r['numerator'],r['denominator']]:
   f.write(str(len(p))+'\n')
   for m,c in p:f.write(f'{m[0]} {qmod(c)}\n')
 def wrnode(f,n,mp,ma=None,varmap=None):
  op=n[0]
  if op in ['m','a','s']:f.write(f'{op} {mp[n[1]]} {mp[n[2]]}\n')
  elif op=='c':f.write(f'c {ma[n[1]]}\n')
  elif op=='v':f.write(f'v {varmap[n[1]]if varmap is not None else n[1]}\n')
  elif op=='l':f.write('l '+str(len(n[1]))+' '+' '.join(f'{ma[c]} {mp[j]}'for c,j in n[1])+'\n')
  else:assert op=='z';f.write('z\n')
 with(work/'shared.txt').open('w')as f:
  f.write(str(len(mons))+'\n');f.writelines(f'{code(m)} {sum(a*b for a,b in zip(m,weights))}\n'for m in mons)
  f.write(str(len(oa))+'\n')
  for i in sorted(atoms):crf(f,O['atoms'][i])
  f.write(str(len(oi))+'\n')
  for i in sorted(need):wrnode(f,ON[i],oi,oa)
  f.write('33\n');f.writelines(f'{l} {oi[r]}\n'for l,r in zip(oracle['genuine_quartic_original_rows'],oracle['genuine_quartic_roots']))
  f.write(str(len(sa))+'\n')
  for i in sorted(sat):crf(f,SA[i])
  f.write(str(len(si))+'\n')
  for i in sorted(sneed):wrnode(f,SN[i],si,sa)
  f.write(str(len(li))+'\n');varmap={v:j for j,v in enumerate(nf['input_variables'])}
  for i in sorted(lneed):wrnode(f,LN[i],li,si,varmap)
  f.write('122\n'+' '.join(str(li[r])for r in nf['output122_linear_roots'])+'\n')
  for row in frame['physical_light_inverse7x7']:
   for r in row:rf(f,r)
  hs=h['scalar_circuit'];pair=h['H2_source_pair_change_scalar_roots'][0]
  for root in pair+[h['H2_known_unit_root']]:
   nd=hs['nodes'][root];assert nd[0]=='c';crf(f,hs['atoms_num_den'][nd[1]])
  for row in quads['quadratic_rows_sparse']:
   f.write(str(len(row))+'\n')
   for m,p in row:f.write(str(code(m))+'\n');poly(f,p)
  active=[i for i,row in enumerate(L)if any(r['numerator']for r in row)];assert len(active)==42
  f.write(str(len(active))+'\n')
  for i in active:
   f.write(str(i)+'\n')
   for r in L[i]:rf(f,r)
 hh=hashlib.sha256();source_core_terms=0
 with(work/'cubics.txt').open('w')as f,gzip.open(mp,'rb')as zf:
  f.write('16\n');head=zf.readline();hh.update(head);assert head.strip()==b'57 91 57'
  for row in range(57):
   for col in range(91):
    line=zf.readline();hh.update(line)
    if row>=41:
     s=line.decode().split();n=int(s[0]);assert len(s)==1+4*n;source_core_terms+=n;f.write(str(n)+'\n')
     for k in range(n):f.write(str(qmod(s[1+4*k]))+' '+' '.join(s[2+4*k:5+4*k])+'\n')
  assert not zf.read()
 assert hh.hexdigest()==br['remaining_matrix_sha256']
 phase.update(phase='all_source_inputs_bound_before_single_native_job',QQ_coefficients_punit_checked=checks,streamed_all16_cubic91_core_terms=source_core_terms,source_native_input_hashes={p.name:sha(p)for p in work.iterdir()if p.is_file()})
 (out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=work/'finite-test';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,54-(time.monotonic()-start)))
 remaining=max(0.1,56-(time.monotonic()-start))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(work),str(work),str(out),str(remaining)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,58-(time.monotonic()-start)))
 summary=read(out/'native-summary.json')if(out/'native-summary.json').exists()else{};front=read(out/'frontier.json')if(out/'frontier.json').exists()else{}
 passed=r.returncode==0 and summary.get('all_passed')is True
 result=dict(id=ID,all_passed=passed,frozen=True,independent_review_pending=True,source_preflight_SHA256=phase['source_preflight_SHA256'],prime=P,
  native_summary=summary,frontier=front,no_new_source_point_prime_or_chart_factor_omitted=True,
  positive_full239_on_finite_H15=passed,global_off_H15_full239_NOT_proved=True,
  coefficient_punit_checks=checks,source_original16_cubic91_core_terms=source_core_terms,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths+[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md']},
  generated_source_inputs={p.name:sha(p)for p in work.iterdir()if p.is_file()and p!=exe and not p.name.startswith('Gamma-')},
  temporary_Gamma_source_hashes={p.name:sha(p)for p in work.glob('Gamma-*.bin')},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 assert sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in work.iterdir()if p.is_file())<32*1024*1024
 signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','generated_source_inputs','temporary_Gamma_source_hashes','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output.resolve(),a.work.resolve())
 except Exception as e:
  if a.output.exists():(a.output/'driver-failure.json').write_text(json.dumps(dict(all_passed=False,error=repr(e),no_auto_retry=True),indent=2)+'\n')
  raise
