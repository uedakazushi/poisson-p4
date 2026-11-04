#!/usr/bin/env python3
"""ONE prepared17x7FFLU: firstunitcofactor vector/globalP80two-minor criterion."""
import hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_P80_REPLACEMENT_SOLVE_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def poly(p):return{int(m):Q(v)for m,v in p if Q(v)}
def integerize(p):
 den=math.lcm(*(Q(v).denominator for v in p.values()));return{m:int(Q(v)*den)for m,v in p.items()},den
def check(products,expected):
 rows=[]
 for p,q in products:
  a,da=integerize(p);b,db=integerize(q);rows.append((a,b,da*db))
 ep,de=integerize(expected);den=math.lcm(de,*(x[2]for x in rows));out={m:-v*(den//de)for m,v in ep.items()}
 for p,q,d in rows:
  for m,c in p.items():
   for n,v in q.items():
    z=m+n;t=out.get(z,0)+c*v*(den//d)
    if t:out[z]=t
    else:out.pop(z,None)
 assert not out
 return True
def mul(p,q):
 a,da=integerize(p);b,db=integerize(q);r={}
 for m,c in a.items():
  for n,v in b.items():r[m+n]=r.get(m+n,0)+c*v
 return{m:Q(v,da*db)for m,v in r.items()if v}
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True);work.mkdir(exist_ok=False,parents=True)
 pb=ROOT/'reports/S61_NONSQUARE_THETA_P80_REPLACEMENT_PREFLIGHT_20261104_v63';pp=pb/'certificate.json';assert sha(pp)==load(pb/'report.json')['evidence_sha256'][pp.name];c=load(pp);A=[[poly(p)for p in row]for row in c['prepared_A17x17']];B=[[poly(p)for p in row]for row in c['prepared_B17x7']];E=poly(c['expected_prepared_det_from_saved_source_NO_newdet']);M=[[poly(p)for p in row]for row in c['prepared_integer_matrix17x24']]
 sb=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_GENERIC_20261104_v63';sp=sb/'substituted17x24-before-rank.json';assert sha(sp)==load(sb/'report.json')['evidence_sha256'][sp.name];original=[[poly(p)for p in row]for row in load(sp)['integer_matrix17x24']];rows=c['source_row_order_for_A_B'];cols=c['A_source_columns'];extras=c['B_source_columns'];rf=list(map(poly,c['full_expanded_row_units']));cf=list(map(poly,c['full_expanded_column_units']));one={0:Q(1)}
 for i in range(17):
  for j in range(24):check([(mul(rf[i],cf[j]),M[i][j])],original[i][j])
 assert A==[[M[i][j]for j in cols]for i in rows]and B==[[M[i][j]for j in extras]for i in rows]
 ub=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_MINOR_UNITS_20261104_v63';up=ub/'exact-known-unit-strip.json';assert sha(up)==load(ub/'report.json')['evidence_sha256'][up.name];Roriginal=poly(load(up)['remaining_polynomial']);assert all(v.denominator==1 for v in Roriginal.values());g=math.gcd(*(int(v)for v in Roriginal.values()));R={m:v/g for m,v in Roriginal.items()};check([(R,{0:Q(g)})],Roriginal)
 detold=poly(load(sb/'native-rank.json')['nonzero_selected_minor']);unitdet=poly(c['whole_det_unit_multiplier']);check([(unitdet,E)],detold)
 def line(p):
  assert all(v.denominator==1 for v in p.values());return str(len(p))+'\n'+''.join(str(m)+' '+str(v.numerator)+'\n'for m,v in sorted(p.items(),reverse=True))
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for p in[p for row in A for p in row]+[p for row in B for p in row]+[E,R]:f.write(line(p))
 cpp=ROOT/'code/s61_nonsquare_Theta_P80_replacement_solve_20261104_v63.cpp';exe=work/'P80-replacement-solve';allcp=work/'all119-completed.json';gcp=work/'first-unit-before-Bezout.json'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','P80replace-native',str(exe),str(inp),str(out/'native-certificate.json'),str(allcp),str(gcp)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');den=poly(n['solver_denominator']);ratio=Q(n['solver_denominator_over_expected_QQunit']);assert ratio;check([(E,{0:ratio})],den);full=load(allcp);assert poly(full['solver_denominator'])==den;proved=False
 if n['selected_Bcolumn']>=0:
  j=n['selected_Bcolumn'];i=n['selected_row'];v=list(map(poly,n['selected17_vector']));assert len(v)==17 and v==[poly(full['all119_solution_cofactors'][k][j])for k in range(17)]
  for row in range(17):check([(A[row][k],v[k])for k in range(17)],mul(den,B[row][j]))
  assert n['gcd_is_QQunit']
  if n['Bezout_saved']:
   S,T=poly(n['Bezout_S_integer']),poly(n['Bezout_T_integer']);constant=int(n['Bezout_nonzero_integer_constant']);assert constant;check([(S,v[i]),(T,R)],{0:Q(constant)});proved=True
 proof={'all408_original_source_preconditioner_identities_independent_integer_checked':True,'solver_denominator_matches_expected_up_nonzero_QQunit':str(ratio),'selected17_vector_AX_deltaB_independent_integer_checked':n['selected_Bcolumn']>=0,'selected_entry_coprime_R46_Bezout_whole_coefficients_verified':proved,'two_actual_source_minors_cover_ALL_P80':proved,'argument':'OffR46oldselecteddet=knownchartunit*R46 isunit. OnR46thechosenreplacementminor equalsXij/QQunit andiscoprimeR46byintegerBezout; originalrow/columnscalingsarechartunits. Thussource17x24fullrank17everywhere; acceptedsourcecriteriongivesfullcubic57allP80.','all119_zero_case_only_sufficient_criterion_inconclusive':n['all119_divisible_R46'],'no_P100_or_fullquartic239_inference':True,'no_K46_Gaussian_new_point_prime_GB':True,'all119_checkpoint_retained_locally':{'path':str(allcp.relative_to(ROOT)),'sha256':sha(allcp),'bytes':allcp.stat().st_size},'first_preBezout_checkpoint':{'path':str(gcp.relative_to(ROOT)),'sha256':sha(gcp),'bytes':gcp.stat().st_size}if gcp.exists()else None}
 (out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),pb/'report.json',pp,sb/'report.json',sp,sb/'native-rank.json',ub/'report.json',up,ROOT/'reports/S61_NONSQUARE_THETA_P80_FACTOR_INPUT_20261104_v63/report.json']
 report={'id':ID,'ALL_P80_full_cubic57_proved_primary':proved,'selected_row':n['selected_row'],'selected_Bcolumn':n['selected_Bcolumn'],'all119_coefficient_slots':n['all119_coefficient_slots'],'Bezout_saved':n['Bezout_saved'],'all119_divisible_R46_inconclusive':n['all119_divisible_R46'],'native_seconds':n['native_seconds'],'all_independent_source_AX_Bezout_checks':True,'checkpoint_all119_local':proof['all119_checkpoint_retained_locally'],'native_executable_SHA_not_committed':exesha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())+allcp.stat().st_size+(gcp.stat().st_size if gcp.exists()else 0)<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);p.add_argument('--work',type=Path,default=ROOT/'work/P80-replacement-solve-first');a=p.parse_args();main(a.output.resolve(),a.work.resolve())
