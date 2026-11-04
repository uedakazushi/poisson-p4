#!/usr/bin/env python3
"""Independent additional-minor norm/pencil/curve48+finite-residual review."""
import argparse,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_ADDITIONAL_NORM_PENCIL_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def poly(p):return{tuple(m):Q(c)for m,c in p if Q(c)}
def plus(p,q,c=Q(1)):
 z=p.copy()
 for m,v in q.items():
  z[m]=z.get(m,Q())+c*v
  if not z[m]:del z[m]
 return z
def mul(p,q):
 z={}
 for m,c in p.items():
  for n,d in q.items():
   k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,Q())+c*d
 return{k:v for k,v in z.items()if v}
def shifted(p,m):return{tuple(x+y for x,y in zip(k,m)):v for k,v in p.items()}
def scale(p,c):return{k:c*v for k,v in p.items()if c*v}
def divide(f,H):
 r=f.copy();q={};assert{m:v for m,v in H.items()if m[2]==4}=={(0,0,4):Q(1)}
 while r and max(m[2]for m in r)>=4:
  deg=max(m[2]for m in r);part={(m[0],m[1],deg-4):v for m,v in r.items()if m[2]==deg}
  q=plus(q,part);r=plus(r,mul(part,H),Q(-1))
 assert plus(mul(q,H),r)==f
 return q,r
class Tasks:
 def __init__(self):self.nodes=[];self.cache={};self.checks=[];self.zero=self.atom({});self.one=self.atom({(0,0,0):Q(1)})
 def atom(self,p):
  key=tuple(sorted((m,str(v))for m,v in p.items()))
  if key not in self.cache:self.cache[key]=len(self.nodes);self.nodes.append(('P',p))
  return self.cache[key]
 def op(self,k,*args):self.nodes.append((k,*args));return len(self.nodes)-1
 def add(self,x,y):return self.op('A',x,y)
 def mul(self,x,y):return self.op('M',x,y)
 def sub(self,x,y):return self.op('S',x,y)
 def det(self,A):return self.op('D',len(A),*[x for row in A for x in row])
 def total(self,vals):
  z=self.zero
  for v in vals:z=self.add(z,v)
  return z
 def eq(self,x,y):self.checks.append(('E',x,y))
 def gcd(self,x,y,z):self.checks.append(('G',x,y,z))
 def write(self,p):
  with p.open('w')as f:
   f.write(str(len(self.nodes))+'\n')
   for n in self.nodes:
    if n[0]=='P':
     p=n[1];amin=min([m[0]for m in p]+[0]);f.write('P '+str(amin)+' '+str(len(p)))
     for m,c in sorted(p.items()):f.write(' '+str(c)+' '+str(m[0]-amin)+' '+str(m[1])+' '+str(m[2]))
     f.write('\n')
    else:f.write(' '.join(map(str,n))+'\n')
   f.write(str(len(self.checks))+'\n')
   for n in self.checks:f.write(' '.join(map(str,n))+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(55);out.mkdir(parents=True,exist_ok=True)
 assert not(out/'report.json').exists()
 stems=['S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM','S61_NONSQUARE_RAW_H_B_NORM_PENCIL','S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT']
 bases=[ROOT/'reports'/(x+'_20261104_v63')for x in stems];reports=[read(b/'report.json')for b in bases]
 paths=[Path(__file__),ROOT/'code/s61_nonsquare_norm_pencil_review_20261104_v63.cpp',ROOT/'research/notes'/(ID+'.md')]
 for b,r in zip(bases,reports):
  assert r['frozen']
  for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in r['evidence_sha256'].items():assert sha(b/p)==h
  paths.extend([b/'report.json',*(b/p for p in r['evidence_sha256'])])
 peers=[ROOT/'reports'/n/'report.json'for n in ['S61_NONSQUARE_RAW_H_B_MINORS_GCD_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_PAIR_NORM_FINITE112_20261104_v63_REVIEW']]
 for p in peers:assert read(p)['all_passed'];paths.append(p)
 nb,pb,rb=bases;nn=read(nb/'native-certificate.json');pn=read(pb/'native-certificate.json');res=read(rb/'certificate.json')
 sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';sc=read(sb/'certificate.json');sr=read(sb/'native-certificate.json');paths.extend([sb/'certificate.json',sb/'native-certificate.json'])
 kb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_INPUT_20261104_v63';kc=read(kb/'certificate.json');gb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_INPUT_20261104_v63';gc=read(gb/'certificate.json');paths.extend([kb/'certificate.json',gb/'certificate.json'])
 assert sc['minor35_to_distinct_generator_index'][33]==20 and kc['selected_generator_index']==20 and kc['source_minor35_indices']==[33]
 k=poly(kc['source_k']);assert k==poly(sc['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'][20])
 prim=poly(sr['records'][33]['primitive_integer_polynomial']);assert scale(prim,1/prim[max(prim)])==k
 H=poly(kc['source_H']);assert H==poly(gc['H']);g=poly(gc['g'])
 matrices=[]
 for source,cert in[(g,gc),(k,kc)]:
  M=[[{}for _ in range(4)]for _ in range(4)]
  for j in range(4):
   quotient,remainder=divide(shifted(source,(0,0,j)),H)
   assert quotient==poly(cert['column_quotients_by_H'][j])
   for i in range(4):M[i][j]={(m[0],m[1],0):v for m,v in remainder.items()if m[2]==i}
  assert M==[[poly(p)for p in row]for row in cert['matrix_rows']];matrices.append(M)
 D=Q(kc['cost_bounds_only']['matrix_uniform_QQ_denominator']);assert kc['cost_bounds_only']['matrix_uniform_a_clear_power']==11
 Dg=Q(gc['determinant_cost_bounds_NOT_value']['matrix_uniform_QQ_denominator']);assert D/Dg==5000
 N=[[[scale(shifted(p,(11,0,0)),D)for p in row]for row in M]for M in matrices]
 assert all(v.denominator==1 and min(m)>=0 for mat in N for row in mat for p in row for m,v in p.items())
 t=Tasks();A,B=[[[t.atom(p)for p in row]for row in mat]for mat in N];detk=t.det(B);t.eq(detk,t.atom(poly(nn['det_cleared_matrix_R'])))
 # Compute all16 adjugate cofactors with last-row Laplace, not permutation.
 adj=[]
 for i in range(4):
  row=[]
  for j in range(4):
   minor=[[B[r][c]for c in range(4)if c!=i]for r in range(4)if r!=j];v=t.det(minor)
   if(i+j)%2:v=t.sub(t.zero,v)
   row.append(v)
  adj.append(row)
 for i in range(4):
  t.eq(adj[i][0],t.atom(poly(nn['adjugate_first_column'][i])))
  for j in range(4):t.eq(t.total([t.mul(B[i][z],adj[z][j])for z in range(4)]),detk if i==j else t.zero)
 t.eq(t.total([t.mul(t.atom(poly(nn['literal_membership_H_cofactor'])),t.atom(H)),t.mul(t.atom(poly(nn['literal_membership_k_cofactor'])),t.atom(k))]),detk)
 U=poly(sc['FULLU_product']);one=t.one;Uid=t.atom(U)
 def strip(record,raw):
  current=t.atom(poly(record['primitive']if 'primitive'in record else record['R_primitive']))
  content=Q(record['QQ_content']if 'QQ_content'in record else record['R_QQ_content']);t.eq(t.mul(current,t.atom({(0,0,0):content})),raw)
  for st in record['localizer_removal_steps']:
   f=t.atom(poly(st['removed_gcd_with_U']));after=t.atom(poly(st['remaining_after']));t.eq(t.mul(f,after),current);t.eq(t.mul(f,t.atom(poly(st['U_quotient']))),Uid);current=after
  result=poly(record['remaining_h_coprime_U']);t.eq(current,t.atom(result));t.gcd(current,Uid,one);return result
 K=strip(nn,detk)
 oldnorm=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_20261104_v63';on=read(oldnorm/'native-certificate.json');paths.extend([oldnorm/'report.json',oldnorm/'native-certificate.json'])
 oldh=poly(on['remaining_h_coprime_U']);G=poly(nn['gcd_newK_oldh']);h1=poly(nn['oldh_quotient_by_gcd']);k1=poly(nn['newK_quotient_by_gcd'])
 t.eq(t.mul(t.atom(G),t.atom(h1)),t.atom(oldh));t.eq(t.mul(t.atom(G),t.atom(k1)),t.atom(K));t.gcd(t.atom(oldh),t.atom(K),t.atom(G));t.gcd(t.atom(h1),t.atom(k1),one)
 # All formal lambda coefficients by row multilinearity: sixteen fixed4det.
 coeff=[t.zero]*5
 for mask in range(16):
  mat=[B[i]if(mask>>i)&1 else A[i]for i in range(4)];j=mask.bit_count();coeff[j]=t.add(coeff[j],t.det(mat))
 for c,p in zip(coeff,pn['pencil_coefficients_lambda0_to4']):t.eq(c,t.atom(poly(p)))
 t.eq(coeff[0],t.mul(t.atom(poly(on['det_cleared_matrix_R'])),t.atom({(4,0,0):Q(5000**4)})));t.eq(coeff[4],detk)
 prev=G;mixed=[]
 for i,record in enumerate(pn['mixed3_records']):
  z=strip(record,coeff[i+1]);mixed.append(z);now=poly(record['curve_gcd_after_this'])
  t.eq(t.mul(t.atom(now),t.atom(poly(record['previous_curve_quotient_by_gcd']))),t.atom(prev));t.eq(t.mul(t.atom(now),t.atom(poly(record['mixed_quotient_by_gcd']))),t.atom(z));t.gcd(t.atom(prev),t.atom(z),t.atom(now));prev=now
 final=poly(pn['final_common_curve_gcd']);assert prev==final==G and len(G)==49
 t.eq(t.mul(t.atom(final),t.atom(poly(pn['original_G_quotient_by_final']))),t.atom(G))
 for x,q in zip(mixed,pn['mixed3_quotients_by_final']):t.eq(t.mul(t.atom(final),t.atom(poly(q))),t.atom(x))
 # Exact cover and Bezout bounds, no resultant or projective-degree claim.
 assert G==poly(res['common_G'])and h1==poly(res['oldh_quotient_h1'])and k1==poly(res['newK_quotient_k1'])and H==poly(res['curve_cover_monic_H_d4'])
 assert max(m[1]for m in G)==12;lc={m:v for m,v in G.items()if m[1]==12};assert len(lc)==1 and(0,12,0)in lc
 assert poly(res['common_G_monic_e12'])==scale(G,1/lc[0,12,0])
 assert max(sum(m)for m in h1)==16 and max(sum(m)for m in k1)==25
 assert res['curve_cover_free_basis_e_i_d_j']==[[i,j]for j in range(4)for i in range(12)]
 assert res['isolated_projection_points_bound_Bezout']==400 and res['isolated_B_bad_points_bound_after_H4']==1600
 # Native tasks are fixed assertions and determinants; no author worker run.
 native=out/'independent-native-input.txt';t.write(native);assert native.stat().st_size<28*1024*1024
 exe=out/'checker';cpp=paths[1];flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-began)))
 with(out/'native.log').open('w')as f:subprocess.run([str(exe),str(native),str(out/'native-result.json')],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-began)))
 nr=read(out/'native-result.json');assert nr['all_passed'];exe.unlink()
 proof={'source_k_is_exact_gen20_minor33':True,'all8_monic_divisions_independently_rebuilt':True,
  'additionalnorm_full16_adjugate_and_H_k_membership':True,'pencil_all5_coefficients_via16_row_multilinearity_determinants':True,
  'all130_known_unit_strips_reconstructed':len(nn['localizer_removal_steps'])+sum(len(r['localizer_removal_steps'])for r in pn['mixed3_records']),
  'G49_unchanged_by_all3_pencil_coefficients':True,
  'necessary_set_cover':'Bbad subset V(H,G) union V(H,h1,k1) on FULLU; both pieces retained.',
  'curve48_meaning':'QQ[a,a^-1,e,d]/(monicG_e12,monicH_d4) is finite FREE RANK48 over QQ[a,a^-1]. Not an ordinary projective-degree48 statement.',
  'finite_residual_bound':'h1,k1 coprime oftotaldegrees16,25 -> atmost400 geometricplane intersections; monicH givesatmost4d roots perpoint ->atmost1600 geometricresidual points.',
  'normpencil_root_logic':'Over each algebraicallyclosed characteristic0 point a!=0, the norm is product_r(g(r)+lambda*k(r)) withmultiplicity. It is identicallyzero iffoneHroot annihilatesbothg,k. No etalenessassumption.',
  'no_bad_existence_or_all_minor_sufficiency':'Only g andk are used by thepencil; other20generators/FULLUlocalizer stillmatter.',
  'no_FULLU_localization_finiteness_flatness_Nakayama_or_globalH7_cubic57_quartic239_claim':True}
 (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 result={'id':ID,'all_passed':True,'source_monic_norm_pencil_and_residual_cover_accepted':True,
  'common_necessary_curve_terms':49,'finite_free_cover_rank_over_a_line':48,'finite_residual_geometric_point_bound':1600,
  'no_curve_existence_badpoint_or_globalrank_claim':True,'native':nr,
  'no_author_worker_GB_resultant_point_prime_or_207expansion':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
  'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
  'preserved_compile_failure':{str(p.relative_to(out)):sha(p)for p in(out/'failed-compile').glob('*')},'seconds':time.monotonic()-began,'frozen':True}
 assert sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256','preserved_compile_failure']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
