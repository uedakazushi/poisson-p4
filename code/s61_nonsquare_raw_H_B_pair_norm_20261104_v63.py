#!/usr/bin/env python3
"""ONE exact native4 norm/adj/Bezout, independent integer products, U-factor stripping."""
import hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_PAIR_NORM_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def poly(p):return{tuple(m):Q(c)for m,c in p if Q(c)}
def intp(p):
 den=math.lcm(*(c.denominator for c in p.values()));return {m:int(c*den)for m,c in p.items()},den
# Independent integer-denominator-cleared multiplication; no FLINT/Python symbolic product reuse.
def check(products,expect):
 rows=[]
 for A,B in products:
  a,da=intp(A);b,db=intp(B);rows.append((a,b,da*db))
 ep,de=intp(expect);den=math.lcm(de,*(x[2]for x in rows));out={m:-c*(den//de)for m,c in ep.items()}
 for a,b,d in rows:
  scale=den//d
  for m,c in a.items():
   for n,e in b.items():
    z=tuple(x+y for x,y in zip(m,n));v=out.get(z,0)+c*e*scale
    if v:out[z]=v
    else:out.pop(z,None)
 assert not out
 return True
def scalar(p,c):return {m:x*c for m,x in p.items()if x*c}
def shift(p,x):return {tuple(m[i]+x[i]for i in range(3)):c for m,c in p.items()}
def add(p,q):
 out=dict(p)
 for m,c in q.items():
  v=out.get(m,Q())+c
  if v:out[m]=v
  else:out.pop(m,None)
 return out
def line(p):return str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in sorted(p.items(),reverse=True))+'\n'
def lline(p):
 sh=min([m[0]for m in p],default=0);return str(sh)+' '+line(shift(p,(-sh,0,0)))
def stats(p):return {'terms':len(p),'totaldegree':max([sum(m)for m in p],default=-1),'degrees':[max([m[j]for m in p],default=-1)for j in range(3)],'max_coefficient_numerator_bits':max([abs(c.numerator).bit_length()for c in p.values()],default=0)}
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True)
 nb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_INPUT_20261104_v63';np=nb/'certificate.json';nr=load(nb/'report.json');assert sha(np)==nr['evidence_sha256'][np.name];c=load(np);H=poly(c['H']);g=poly(c['g']);M=[[poly(p)for p in row]for row in c['matrix_rows']];Qcol=list(map(poly,c['column_quotients_by_H']));bc=c['determinant_cost_bounds_NOT_value'];aclear=bc['matrix_uniform_a_clear_power'];D=Q(bc['matrix_uniform_QQ_denominator']);assert aclear==10
 pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_COPRIME_PAIR_PREFLIGHT_20261104_v63';pp=pb/'certificate.json';pc=load(pp);assert sha(pp)==load(pb/'report.json')['evidence_sha256'][pp.name];f=poly(pc['selected_pair_source']['polynomials'][0]);mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';mp=mb/'certificate.json';mc=load(mp);assert sha(mp)==load(mb/'report.json')['evidence_sha256'][mp.name];q=poly(mc['q'])
 sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';sp=sb/'certificate.json';sc=load(sp);assert sha(sp)==load(sb/'report.json')['evidence_sha256'][sp.name];U=poly(sc['FULLU_product']);N=[[scalar(shift(p,(aclear,0,0)),D)for p in row]for row in M];assert all(all(m[0]>=0 and x.denominator==1 for m,x in p.items())for row in N for p in row)
 ip=out/'native-input.txt'
 with ip.open('w')as stream:
  for row in N:
   for p in row:stream.write(line(p))
  for p in [H,g,f,q]+Qcol:stream.write(lline(p))
  stream.write(line(U));stream.write(str(aclear)+' '+str(D)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-pair-norm'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','norm4-native',str(exe),str(ip),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');R=poly(n['det_cleared_matrix_R']);v=list(map(poly,n['adjugate_first_column']));assert R and len(v)==4 and max(m[2]for m in R)==0
 for i in range(4):check([(N[i][j],v[j])for j in range(4)],R if i==0 else{})
 one={(0,0,0):Q(1)};sunit={(aclear,0,0):D};V=poly(n['monic_V']);QQ=poly(n['monic_Qsum']);uf=poly(n['literal_membership_f_cofactor']);vg=poly(n['literal_membership_g_cofactor'])
 check([(g,shift(V,(aclear,0,0))),(scalar(H,-1),shift(QQ,(aclear,0,0)))],scalar(R,1/D));check([(uf,f),(vg,g)],R)
 primitive=poly(n['R_primitive']);check([(primitive,{(0,0,0):Q(n['R_QQ_content'])})],R);current=primitive
 for step in n['localizer_removal_steps']:
  t=poly(step['removed_gcd_with_U']);after=poly(step['remaining_after']);check([(t,after)],current);check([(t,poly(step['U_quotient']))],U);current=after
 h=poly(n['remaining_h_coprime_U']);assert current==h and poly(n['final_gcd_h_U'])==one
 RR,a,e,d=ring('a,e,d',s.QQ);rp=lambda p:RR.from_dict({m:s.QQ(str(c))for m,c in p.items()});cg=rp(U).gcd(rp(h));assert cg/cg.LC==RR.one
 assert stats(R)['terms']<=bc['det_cleared_support_box_upper'] and stats(R)['max_coefficient_numerator_bits']<=bc['det_cleared_integer_abs_coefficient_bits_upper']
 nonconstant=any(any(m)for m in h);classification='necessary_nonzero_bad_projection_curve_only'if nonconstant else'pair_ideal_unit_on_FULLU_and_B3_global'
 proof={'classification':classification,'independent_integer_all4_matrix_adj_identities':True,'independent_integer_full_monic_and_original_pair_memberships':True,'independent_integer_QQcontent_and_ALL_localizer_removal_reconstructions':True,'different_reverse_SymPy_exact_gcd_h_U_is1':True,'known_inverses_before_strip':['a','QQconstants'],'removed_gcd_factors_unit_on_FULLU_by_literal_U_quotients':True,'actual_B_bad_projection_contained_in_V_h':True,'nonconstant_h_NOT_badexistence':True,'no_factorization_GB_primarydec_newpoints_or_primes':True,'monic_clear_scalar_D_a10':{'D':str(D),'a_power':aclear},'literal_identities':['N*adj(N)*e0=R*e0','s*g*V-s*H*Qsum=R','uf*f+vg*g=R'],'R_stats':stats(R),'primitive_R_stats':stats(primitive),'remaining_h_stats':stats(h)}
 (out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),nb/'report.json',np,pb/'report.json',pp,mb/'report.json',mp,sb/'report.json',sp]
 report={'id':ID,'classification':classification,'R_stats':stats(R),'remaining_h_stats':stats(h),'adj4_terms':list(map(len,v)),'membership_cofactor_terms':[len(uf),len(vg)],'localizer_removal_steps':len(n['localizer_removal_steps']),'native_seconds':n['native_seconds'],'peak_single_native_support':n['peak_single_polynomial_support'],'GB_newpoint_prime_factorization_jobs':0,'all_exact_native_and_independent_integer_memberships_checked':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
