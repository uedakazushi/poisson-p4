#!/usr/bin/env python3
"""One source-bound L-coordinate polynomial-input preflight; no rank or determinant."""
import argparse,hashlib,itertools,json,os,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def treebytes(p):return sum(f.stat().st_size for f in p.rglob('*')if f.is_file())
def monomials(n):
 return [tuple(c.count(i)for i in range(5))for c in itertools.combinations_with_replacement(range(5),n)]
def main(out):
 start=time.monotonic();deadline=start+60;resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 old=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63';oldreport=old/'report.json';r=json.loads(oldreport.read_text());assert r['frozen']and r['status']=='blocked_no_identity_conclusion'
 names=['frozen57-source.txt','quad20x35.txt','selection.txt']
 for name in names:assert sha(old/name)==r['artifact_sha256'][name]
 circp=ROOT/'reports/S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63/circuit.json';circ=json.loads(circp.read_text());assert sha(old/names[0])==circ['source_uncompressed_sha256']
 pair=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63/input-pair-native.txt'
 cpp=ROOT/'code/s61_nonsquare_ell_Lframe_preflight_20261104_v63.cpp';exe=out/'Lframe-native'
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(.1,deadline-time.monotonic()))
 stop=None;maxrss=0
 with(out/'process.log').open('w')as f:
  proc=subprocess.Popen([str(exe),*(str(old/name)for name in names),str(pair),str(out)],stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  while proc.poll()is None:
   if time.monotonic()>=deadline:stop='total60s_deadline';break
   if treebytes(out)>32*1024*1024:stop='artifact32MiB_guard';break
   try:
    for line in Path('/proc',str(proc.pid),'status').read_text().splitlines():
     if line.startswith('VmRSS:'):maxrss=max(maxrss,int(line.split()[1])*1024)
   except FileNotFoundError:pass
   time.sleep(.05)
  if stop:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
 exe.unlink(missing_ok=True)
 nativep=out/'native.json';ok=proc.returncode==0 and nativep.exists();n=json.loads(nativep.read_text())if ok else{}
 if ok:
  assert n['actual_source_q1_q3_identities_passed']and n['all12_pair_coefficients_match_frozen_global_dual_source']
  assert n['new_quad10x26']['rows']==10 and n['new_genuine16x31']['columns']==31
 pure=monomials(3);kept=list(range(15,35))+list(range(9,15));light=['L','x0','y1','y2','y3']
 c={'old_five_lights':['x0=v0','x1=v1','x2=v2','x3=v4','x4=v7'],'old_restriction_v8':circ['restriction_v8'],
 'new_five_lights':light,'new_heavy':'H=v5','new_y1_y2_y3':['old_x2','old_x3','old_x4'],
 'L_formula_from_actual_q1':'L=-(405/4)*d(q1)/d(old_x1)/C','C':'100*a+18*e^2-9*e-2','L_old_x1_coefficient':'-40*C',
 'inverse':'old_x1=(L-L_0*x0-L_2*y1-L_3*y2-L_4*y3)/(-40*C)',
 'inverse_denominator_only_known_unit':'-40*C','actual_source_normalizations':['q1=L^2/8100','q3=L*(e*L-400*a*x0)/8100'],
 'normalized_monomial_quadric_generators':['L^2','L*x0'],'normalization_operation':'L^2=8100*q1; L*x0=(8100*e*q1-8100*q3)/(400*a)',
 'qfour4x15_row_scaling':['1','(-40*C)^2','(-40*C)^2','(-40*C)^2/C^q5_C_factors_cancelled'],
 'genuine16_row_scaling':'(-40*C)^3/C^genuine_C_factors_cancelled[row]','only_exact_polynomial_C_divisions_used':True,
 'full35_cubic_monomials_new_light_order':pure,'nine_monomial_pivot_columns':list(range(9)),
 'remaining26_full35_column_indices':kept,'remaining26_order':['20 cubics in (x0,y1,y2,y3)','6 L*quadrics in (y1,y2,y3)'],
 'genuine31_columns':['remaining26 pure cubic monomials','H*L','H*x0','H*y1','H*y2','H*y3'],
 'quad10_rows':['q0*x0','q0*y1','q0*y2','q0*y3','scaled_q5*x0','scaled_q5*y1','scaled_q5*y2','scaled_q5*y3','q0*L','scaled_q5*L'],
 'quad10_block_shape':'[[A(8x20),B(8x6)],[0(2x20),Cpair(2x6)]]',
 'pair_source':str(pair.relative_to(ROOT)),'pair_column_reversal':[5,4,3,2,1,0],
 'pair_source_identity':'newPair[0,j]=oldPair[0,5-j]; C*newPair[1,j]=((-40*C)^2/C^k)*oldPair[1,5-j]',
 'pair_source_old_row1':'C*q5 on L=x0=0','pair_global_rank2_previously_proved_not_recomputed':True,
 'actual_genuine_rows':list(range(41,57)),'polynomial_native_record_schema':'rows columns; each entry: n followed by n rational,aexp,eexp,dexp quadruples',
 'no_arithmetic_elimination_of_the_remaining10_quad_rows':True,'canonical_R11_not_constructed':True,
 'no_determinant_rank_GB_field_newpoint_or_quartic_computation':True,'allpoint_cubic57_full239_and_kernel_count_unchanged':True,
 'bounded_next_gate':'Use the already proved pair dual to split Cpair, then a structured 8x20 block; formulate exact canonical11 before any new border identity attempt. No cost guarantee or automatic job.'}
 (out/'certificate.json').write_text(json.dumps(c,indent=2)+'\n')
 inputs=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),oldreport,*(old/name for name in names),circp,pair,
 ROOT/circ['quad_source'],ROOT/circ['parametric_source'],ROOT/circ['normal_monomials_source'],ROOT/'reports/S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63/report.json']
 report={'id':ID,'all_source_bound_coordinate_and_polynomial_input_checks_passed':ok,'native_returncode':proc.returncode,'budget_stop_reason':stop,
 'no_new_rank_or_determinant_computation':True,'no_allpoint57_full239_or_new_kernel_claim':True,'canonical_R11_not_yet_computed':True,
 'source_and_recipe_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'artifact_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},
 'native_profiles':n,'seconds':time.monotonic()-start,'observed_peak_native_RSS_bytes':maxrss,'combined_artifact_bytes':treebytes(out),'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['source_and_recipe_sha256','artifact_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
