#!/usr/bin/env python3
"""Source-bound joint L,ell leading-term preconditioner, with all overlap rows."""
import argparse,hashlib,itertools,json,resource,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pure(n,D):
 return [tuple(c.count(i)for i in range(n))for c in itertools.combinations_with_replacement(range(n),D)]
def weighted(D):
 return sorted([m+(h,)for h in range(D//2+1)for m in pure(7,D-2*h)],reverse=True)
def wp(f,p):f.write(str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in sorted(p.items(),reverse=True))+'\n')
def main(out):
 began=time.monotonic();deadline=began+60;resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 compact=ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';pre=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63'
 cr=read(compact/'report.json');pr=read(pre/'report.json');assert cr['all_literal_source_identity_and_outside37_gate_checks_passed']and pr['all_source_bound_coordinate_and_polynomial_input_checks_passed']
 qfile=compact/'six-original-quadrics.txt';assert sha(qfile)==cr['artifact_sha256'][qfile.name]
 Lfile=pre/'L-coefficients.txt';assert sha(Lfile)==pr['artifact_sha256'][Lfile.name]
 anchorpath=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/certificate.json';anchor=read(anchorpath)
 ell=[{}for _ in range(7)]
 for m,ts in anchor['ell']:assert sum(m)==1;ell[m.index(1)]={tuple(ex):Q(c)for ex,c in ts}
 ellfile=out/'ell-x7.txt'
 with ellfile.open('w')as f:
  f.write('7 1\n')
  for row in ell:wp(f,row)
 C319path=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';C319=read(C319path);defs=[C319['source_rows'][i]for i in C319['retained_row_indices']];oldm2=[r[2]for r in defs[:29]]
 m2,m3,m4=map(weighted,[2,3,4]);assert list(map(len,[m2,m3,m4]))==[29,91,239]
 meta=out/'monomial-input.txt'
 with meta.open('w')as f:
  f.write('29 91 239\n')
  for mons in[oldm2,m2,m3]:
   for m in mons:f.write(' '.join(map(str,m))+'\n')
 # Exact fixed monic-pivot choices, not a Groebner/Hilbert-function assertion.
 leads=[(2,0,0,0,0,0,0,0),(1,0,1,0,0,0,0,0),(0,1,0,0,0,1,0,0),(0,1,0,0,0,0,1,0)]
 stages=[]
 for D,mons in[(3,m3),(4,m4)]:
  cofs=weighted(D-2);rows=[[i,m]for i in range(4)for m in cofs];lookup={(i,tuple(m)):j for j,(i,m)in enumerate(rows)};selected=[];piv=[];remaining=[]
  for c,m in enumerate(mons):
   valid=[i for i,lm in enumerate(leads)if all(x>=y for x,y in zip(m,lm))]
   if not valid:remaining.append(c);continue
   i=valid[0];cof=tuple(x-y for x,y in zip(m,leads[i]));selected.append(lookup[i,cof]);piv.append(c)
  assert len(selected)==len(set(selected))
  overlaps=[j for j in range(len(rows))if j not in selected]
  stage={'weighted_degree':D,'target_monomials':mons,'four_quadric_product_rows':rows,'monic_pivot_source_rows':selected,'monic_pivot_target_columns':piv,'retained_target_columns':remaining,
   'all_unselected_overlap_rows_retained':overlaps,'unit_triangular_argument':'Each selected product has precisely its distinct leading monomial with coefficient1 and only lexicographically lower terms. Descending weighted-degree/lex order gives a unit triangular square, with inverse a finite sum (I-U)^k. No S row is set0 by monomial-divisibility alone.',
   'remaining_q0_q5_products':[[i,m]for i in[4,5]for m in cofs]}
  stages.append(stage)
 assert[len(x['monic_pivot_source_rows'])for x in stages]==[26,98]
 assert[len(x['retained_target_columns'])for x in stages]==[65,141]
 assert[len(x['all_unselected_overlap_rows_retained'])for x in stages]==[2,18]
 (out/'preflight.json').write_text(json.dumps({'new_light_order':['L','ell','x0','x3','r','Z','W'],'r':'7*x2/2+3*x4','heavy':'H','old_light_order':['x0','x1','x2','x3','x4','Z','W'],
  'coordinate_Jacobian':'3200*a*C*d^2/81','stages':stages,'original_four_exact_monomial_quadric_proposal':'REJECT unlessfull corrections vanish; actual source checks below preserve corrections',
  'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[qfile,Lfile,anchorpath,C319path]}},indent=2)+'\n')
 cpp=ROOT/'code/s61_nonsquare_joint_Lell_preflight_20261104_v63.cpp';exe=out/'joint-Lell';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(.1,deadline-time.monotonic()))
 with(out/'native.log').open('w')as f:subprocess.run([str(exe),str(qfile),str(meta),str(Lfile),str(ellfile),str(out)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(.1,deadline-time.monotonic()))
 exe.unlink();n=read(out/'native.json');assert n['four_monic_leading_terms_checked']and n['false_four_monomial_ideal_rejected_by_literal_nonzero_correction']
 circ={'coefficient_domain':'QQ[a,e,d,1/(a*d*F*B*C*J)]','coordinate_inverse':'old_x=N*new_y/delta; all49+49 identities exact','delta':'3200*a*C*d^2/81',
  'six_new_generators':['p1=8100*q1','p2=8100*(e*q1-q3)/(400*a)','p3=ell*Z','p4=ell*W','p5=q0','p6=q5'],
  'full_source_inverse':['q0=p5','q1=p1/8100','q3=e*p1/8100-400*a*p2/8100','q5=p6','q2=(p4+(d*e+8)*q1-d*q3)/8','q4=((10*a+3*e^2-3*e)*q1-3*(e-1)*q3-a*p3)/(10*a)'],
  'source_change_known_units_only':['a','QQ_nonzero_constants'],'coordinate_denominator_factors_only':['a','d','C'],
  'actual_four_monomial_ideal_claim_rejected':True,'actual_generator_form':'p1=L^2+ZWcorrection; p2=L*x0+ZWcorrection; p3=ellZ,p4=ellW',
  'nonzero_overlap':'x0*p1-L*p2; full numerator/knownunitdenominator saved, must be retained after the26 pivot projection',
  'exact_zero_overlap':'W*p3-Z*p4=0; this is the ONE universal cubic overlap that may be discarded',
  'cubic_source_recipe':'Four products28 +q0/q5 products14+genuine16: pivot26 and discardonly exactzero overlap leaves31->65, same fullcubic source span. No allpoint rank31 claim.',
  'quartic_source_recipe':'Four products116+q0/q5 products58+genuine112+quartic33: pivot98 leaves18+58+112+33=221->141. All18 overlaps retained; no Hilbert function or rank141 claim.',
  'coefficient_projection_recipe':'For each degree use selected rows R and columns C frompreflight. U=P[R,C] is unittriangular from saved6 primequadrics. NF(row)=row[D]-row[C]*U^(-1)*P[R,D]. U^(-1) is the finite triangular arithmetic circuit, all source factors known units.',
  'G16_lazy_transport_recipe':'Read the samefrozen57 polynomial source rows41..56; replace oldv=(x0,x1,x2,W,x3,H,Z+x1-3*x0,x4) using inverseN/delta and retain all91 weighteddegree3 coefficients. No expanded G16 core coefficients computed inthis initial quad preflight.',
  'all33_quartic_lazy_transport_recipe':'Use acceptedTn/Qhat319 coefficient oracle; same invertible light map, Hunchanged. No new33quartic expansion.',
  'free183_followon_not_yet_constructed':'After98 target/source pivots, known15 survivein nonquad coordinates; descendK22/P22 toremainingquad76, useBnew22/Cnew excluding182/186,thenknownunitseparator. Authorized separately asfinitecircuit stage.',
  'no_rank_job_newpoint_prime_GB_or_281_expansion':True,'native':n}
 (out/'circuit.json').write_text(json.dumps(circ,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/'code/s61_nonsquare_ell_Lframe_preflight_20261104_v63.cpp',ROOT/'research/notes'/(ID+'.md'),qfile,Lfile,anchorpath,C319path,ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json',ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW/report.json']
 report={'id':ID,'all_exact_joint_coordinate_and_quadric_checks_passed':True,'invalid_exact_four_monomial_ideal_rejected':True,
  'validated_leading_monic_pivots_degree3_degree4':[26,98],'exact_target_complement_counts_not_Hilbert_function':[65,141],
  'all_overlap_source_rows_retained_degree3_degree4':[2,18],'exact_full_quartic_Schur_source_target_shape':[221,141],
  'G16_core_expansion_not_performed_initial_input_stage':True,'finite_circuit_transport_preferred_over_new_parameter_expansion':True,
  'no_new_rank_GB_point_field_or_full239_claim':True,'native':n,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
  'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'artifact_bytes':sum(p.stat().st_size for p in out.iterdir()if p.is_file()),'seconds':time.monotonic()-began,'frozen':True}
 assert report['artifact_bytes']<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
