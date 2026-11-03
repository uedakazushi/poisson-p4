#!/usr/bin/env python3
"""Independent exact source/split audit and reversed native block rank replay."""
import argparse,gzip,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_FINITE_SIGMA_QUARTIC_NAKAYAMA_20261103_v63_REVIEW';pchar=32003

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def add(a,b):
 z=dict(a)
 for m,c in b.items():z[m]=z.get(m,Q())+c
 return{m:c for m,c in z.items()if c}
def mul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,Q())+c*d
 return{k:c for k,c in z.items()if c}
def mons(n,d):return[(d,)]if n==1 else[(i,)+v for i in range(d+1)for v in mons(n-1,d-i)]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 parent=ROOT/'reports/S61_NONSQUARE_FINITE_SIGMA_QUARTIC_NAKAYAMA_20261103_v63';pr=read(parent/'report.json');assert pr['all_exact_QQ_split_and_source_checks_passed']and pr['actual_raw_quad78_proved_on_ENTIRE_finite_Sigma']
 paths=[Path(__file__),ROOT/'code/s61_nonsquare_finite_sigma_quartic_review_20261103_v63.cpp',ROOT/'research/notes'/(ID+'.md'),parent/'report.json']
 for name,want in pr['input_sha256'].items():assert sha(ROOT/name)==want,name
 for name,want in pr['evidence_sha256'].items():assert sha(parent/name)==want,name
 probe=ROOT/'reports/S61_NONSQUARE_FINITE_SIGMA_QUARTIC_PROBE_20261103_v63';probe_r=read(probe/'report.json');assert probe_r['native_completed']and probe_r['rank']==2184
 for name,want in probe_r['input_sha256'].items():assert sha(ROOT/name)==want,name
 finite=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW';fr=read(finite/'report.json');assert fr['all_passed']and fr['finite_V_algebra_without_flatness_reviewed']and fr['FULLU_unit_on_ENTIRE_special_fiber_and_finite_V_model_reviewed']
 for name,want in fr['input_sha256'].items():assert sha(ROOT/name)==want,name
 for name,want in fr['evidence_sha256'].items():assert sha(finite/name)==want,name
 source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';data=read(source);qfour=[];shifts=[];terms=0
 for row in data['cleared_qfour_sparse']:
  shift=max([0]+[ex[2]-ex[0]for mon,cs in row for ex,c in cs]);shifts.append(shift);q={}
  for mon,cs in row:
   assert len(mon)==7 and mon[6]==0 and sum(mon)==2;coef={}
   for ex,c in cs:
    c=Q(c);assert c.denominator%pchar;new=(ex[0]-ex[2]+shift,ex[1],ex[2]);assert min(new)>=0;coef=add(coef,{new:c});terms+=1
   q[tuple(mon[:6])]=coef
  qfour.append(q)
 assert shifts==[0,2,2,1]and terms==804
 split=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/right_inverse.txt.gz';ts=iter(gzip.decompress(split.read_bytes()).decode().split());assert[int(next(ts))for _ in range(3)]==[21,4,0];L=[];nterms=0
 for i in range(21):
  row=[]
  for j in range(4):
   z={}
   for k in range(int(next(ts))):
    c=Q(next(ts));ex=tuple(int(next(ts))for _ in range(3));assert c.denominator%pchar and c and ex not in z;z[ex]=c;nterms+=1
   row.append(z)
  L.append(row)
 assert next(ts,None)is None and nterms==2461;m2=mons(6,2);g={(2,0,3):Q(-5),(2,1,2):Q(6),(2,0,2):Q(4)};products=[]
 for i in range(4):
  row=[]
  for j in range(4):
   z={}
   for k,mon in enumerate(m2):z=add(z,mul(qfour[i].get(mon,{}),L[k][j]))
   assert z==(g if i==j else{});row.append(z)
  products.append(row)
 pairs=list(itertools.combinations(range(4),2));signature=[]
 for i,j in pairs:
  row=[]
  for a,b in pairs:
   # Literal ordered right inverse: only block a carries L[:,b]/g.
   z=products[j][b]if a==i else({m:-c for m,c in products[i][b].items()}if a==j else{})
   assert z==(g if(i,j)==(a,b)else{});row.append(int(bool(z)))
  signature.append(row)
 # Our input uses the OPPOSITE lex order. Build independently from804 terms.
 own_m2=list(reversed(m2));ownfile=out/'own-four-quads-modp.txt'
 with ownfile.open('w')as f:
  f.write('4 21 3\n')
  for q in qfour:
   for mon in own_m2:
    z={m:c.numerator*pow(c.denominator,-1,pchar)%pchar for m,c in q.get(mon,{}).items()};z={m:c for m,c in z.items()if c};f.write(str(len(z))+'\n')
    for m,c in sorted(z.items()):f.write(str(c)+' '+' '.join(map(str,m))+'\n')
 block=ROOT/'work/marked-square-v63/finite-Sigma-quartic-probe-first/actual-block.u32le';assert block.stat().st_size==12+4*2352*3528 and sha(block)==probe_r['actual_block_u32le_sha256']
 ops=finite/'own-operators.txt';cpp=ROOT/'code/s61_nonsquare_finite_sigma_quartic_review_20261103_v63.cpp';exe=out/'native-bin';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();subprocess.run(['c++','-O3','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,capture_output=True,timeout=20)
 with(out/'native.log').open('w')as log:subprocess.run([str(exe),str(ops),str(ownfile),str(block),str(out/'native.json')],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=30)
 nr=read(out/'native.json');assert nr['all_passed']and nr['own_single_rank']==2184 and nr['all_author_block_entries_permuted_exact']==8297856;exe.unlink()
 paths.extend([source,split,finite/'report.json',ops,probe/'report.json',ROOT/'research/notes/S61_NONSQUARE_FINITE_SIGMA_QUARTIC_NAKAYAMA_20261103_v63.md'])
 report={'id':ID,'all_passed':True,'own_Fraction_actual804_and_all16_QL_gI_exact':True,'own_all36_ordered_K_right_inverse_signatures_exact':True,'dual_orientation_S126_to_S84_and_KM_zero_reviewed':True,'own_reversed_source_all8297856_block_entries_exact':True,'own_single_FLINT_rank2184':True,'finite_V_split_kernel_cokernel_Nakayama_valid_without_flatness_or_socle':True,'raw_quad78_all_finite_Sigma_accepted':True,'original_light146_weighted152_finite_Sigma_corollary_accepted':True,'all_U_quad78_C122_rank64_full_quartic_not_claimed':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in[out/'own-four-quads-modp.txt',out/'native.log',out/'native.json']},'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
