#!/usr/bin/env python3
"""ONE15x6 raw H-linear preflight and individual known-unit Gaussian; no H1 retry."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 hb=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';hp=hb/'certificate.json';assert sha(hp)==read(hb/'report.json')['certificate_sha256'];hdata=read(hp);hr=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW/report.json';assert read(hr)['all_passed']
 pb=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63';pp=pb/'certificate.json';pair=read(pp);assert sha(pp)==read(pb/'report.json')['certificate_sha256']
 R,a,e,d=ring('a,e,d',s.QQ);dec=lambda p:R.from_dict({tuple(m):s.QQ(c)for m,c in p});enc=lambda p:[[list(m),str(c)]for m,c in p.items()]
 F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;J=5*a*d-6*e-4;known=[a,d,F,B,C,J]
 # Literal112 frozen raw H-linear entries; no scalar DAG evaluation.
 original=[]
 for row in hdata['genuine_cubic16_Hlinear_old_scalar_roots']:
  rr=[]
  for root in row:
   if not root:rr.append(R.zero);continue
   n=hdata['scalar_circuit']['nodes'][root];assert n[0]=='c';num,den=hdata['scalar_circuit']['atoms_num_den'][n[1]];assert dec(den)==R.one;rr.append(dec(num))
  original.append(rr)
 assert len(original)==16 and all(len(row)==7 for row in original)
 H=pair['constant_heavy_coefficient_change7x7'];Hi=pair['constant_heavy_inverse7x7'];assert all(sum(H[i][k]*Hi[k][j]for k in range(7))==int(i==j)for i in range(7)for j in range(7));assert all(sum(Hi[i][k]*H[k][j]for k in range(7))==int(i==j)for i in range(7)for j in range(7))
 X=[[sum((H[j][k]*row[k]for k in range(7)),R.zero)for j in range(7)]for row in original]
 delta=-64000*C**2;zz=delta*(X[6][5]-X[1][5]);ww=delta*(X[6][6]-X[1][6]);assert all(X[6][j]==X[1][j]for j in range(5))
 g1=(3*e-2)*d+20;g2=(6*e+1)*d+40;hh=C*g1-B;z=-a*F*B*g1*g2/10000000;w=-F*B*hh/10000000;assert zz==delta*z and ww==delta*w
 cZ=dec(pair['Bezout_cZ']);cW=dec(pair['Bezout_cW']);eta=dec(pair['known_unit_eta']);assert cZ*zz+cW*ww==eta and eta==s.QQ(192,125)*a*F*B**2*C**3
 gamma=delta*F*B/10000000;assert gamma==-s.QQ(4,625)*F*B*C**2
 rows=[i for i in range(16)if i!=6];M=[];clear=[]
 for i in rows:
  last=hh*X[i][5]-a*g1*g2*X[i][6];assert gamma*last==-ww*X[i][5]+zz*X[i][6];M.append(X[i][:5]+[last]);clear.append(cZ*X[i][5]+cW*X[i][6])
 terms=sum(len(p)for row in M for p in row);maxentry=max(len(p)for row in M for p in row);degree=max(sum(mon)for row in M for p in row for mon in p)if terms else 0
 pre={'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','physical_light_order':[0,1,2,3,4,6,7],'heavy_variable':5,'raw16x7':[[enc(p)for p in row]for row in original],
  'constant_coefficient_column_frame7x7':H,'constant_inverse7x7':Hi,'frameX_order':['x0','x1','x2','x3','x4','Z','W'],
  'source_delta_change':{'replace_genuine_index':6,'formula':'Cdelta=delta*(C6-C1)','delta':enc(delta),'determinant_knownunit':'-64000*C^2'},
  'pair_values':{'z_unscaled':enc(z),'w_unscaled':enc(w),'delta_z':enc(zz),'delta_w':enc(ww),'cZ':enc(cZ),'cW':enc(cW),'eta':enc(eta),'gamma':enc(gamma)},
  'coefficient_column_GL2':[['cZ/eta','cW/eta'],['-delta*w/gamma','delta*z/gamma']],'coefficient_column_GL2_inverse':[['delta*z','-gamma*cW/eta'],['delta*w','gamma*cZ/eta']],
  'target_GL2_determinant':'1/gamma (knownunit)','pair_identity_and_GL2_inverse_cleared_QQ_verified':True,'unit_Cdelta_row':'(0,0,0,0,0,1,0)',
  'other15_full_source_clearings_numerators_over_eta':[enc(p)for p in clear],'full_source_clearings':'Ci -> Ci - (numerator_i/eta)*Cdelta; actual pure coefficients retained lazy.',
  'remaining15_genuine_indices':rows,'remaining6_coefficient_columns':[0,1,2,3,4,6],'matrix15x6':[[enc(p)for p in row]for row in M],
  'matrix_terms':terms,'max_entry_terms':maxentry,'max_totaldegree':degree,'support_at_most10000_required':terms<=10000,
  'no_NF117_Hquadratic_quartic_newpoint_or_GB':True}
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');native=out/'native-input.txt';cpp=ROOT/'code/s61_nonsquare_raw_H_linear_units_20261104_v63.cpp';exe=ROOT/'work/raw-H-linear-units'
 def line(p):return str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in p.items())+'\n'
 stop='preflight_support_above10000';k=0;checkpoint=None;audit=None;exesha=None
 if terms<=10000:
  with native.open('w')as f:
   f.write('15 6 6\n')
   for p in known:f.write(line(p))
   for row in M:
    for p in row:f.write(line(p))
  command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
  with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-began)))
  exesha=sha(exe)
  with(out/'native.log').open('w')as f:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','rawHnative',str(exe),str(native),str(out/'checkpoint.json')],stdout=f,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-began)))
  assert run.returncode==0;audit=read(out/'checkpoint.json.audit.json');k=audit['completed_pivots'];stop=audit['stop_reason'];checkpoint=read(out/'checkpoint.json');assert checkpoint['completed_pivots']==k;exe.unlink(missing_ok=True)
  # Verify the selected-original-minor unit formula using only scalar constants/exponent vectors.
  scalars=[[s.QQ(1),[0]*6]for i in rows];det=[s.QQ(1),[0]*6]
  def u(q):return[s.QQ(q['constant']),q['a_d_F_B_C_J_powers']]
  def combine(v,w,divide=False):return[v[0]/w[0]if divide else v[0]*w[0],[a-b if divide else a+b for a,b in zip(v[1],w[1])]]
  for i,unit in enumerate(checkpoint['initial_row_divisors']):scalars[i]=combine(scalars[i],u(unit),True)
  for step in checkpoint['steps']:
   r=step['pivot_row'];p=u(step['pivot_unit']);assert scalars[r]==u(step['row_scale_before']);det=combine(det,combine(p,scalars[r],True))
   for elim in step['eliminations']:i=elim['row'];scalars[i]=combine(combine(scalars[i],p),u(elim['row_divisor_after']),True)
  assert det==u(checkpoint['original_selected_minor_unit'])
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),hb/'report.json',hp,hr,pb/'report.json',pp]
 result={'id':ID,'actual16_raw_H_inputs_source_bound':True,'delta_eta_GL7_and_knownunit_column_rescale_verified':True,'preflight_matrix_shape':[15,6],'preflight_terms':terms,'preflight_max_entry_terms':maxentry,
  'additional_global_H_linear_unit_pivots':k,'total_proved_H_linear_unit_lower':1+k,'full_H7_rowspan_ALL_FULLU_proved_primary':k==6,'if_fullH7_then_all_H_quadratics28_and_H1_projection23_surjective':True,
  'stop_reason':stop,'native_audit':audit,'no_full57_full239_newkernel_count_or_automatic_free156_claim':True,'pure_cubic_coefficients_and_allquartics_retained_lazy':True,
  'native_executable_sha256_not_committed':exesha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
