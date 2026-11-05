#!/usr/bin/env python3
"""One projected natural-map value, from existing source2 caches only."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_honest_nil_T_gauge120_cached_20261105_v63 as G
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,nqq,invq,mv
from ambient_B18_Told_Tnil_head4_z5_source_preflight_20261105_v63 import PATHS
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63'
GATE='AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_SOURCE_PREFLIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 gr=ROOT/'reports'/GATE/'report.json';gate=json.loads(gr.read_text());assert gate['all_passed'] and gate['whole_worker_QQ_upper']==35700
 for p,h in gate['input_sha256'].items():assert sha(ROOT/p)==h,p
 gp=ROOT/'reports'/GATE/'head4-z5-projected-source-and-cache-gate.json';assert sha(gp)==gate['evidence_sha256'][gp.name];H=json.loads(gp.read_text())
 pins={str(gr.relative_to(ROOT)):sha(gr),str(gp.relative_to(ROOT)):sha(gp),**gate['input_sha256']}
 (out/'phase.json').write_text(json.dumps({'phase':'whole SOURCE gate and fixed source cache hashes verified before value','input_sha256':pins},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 P=D['physics'];MON=[tuple(m)for m in P['base25_monomials']];MI={m:i for i,m in enumerate(MON)};G.MON=MON;G.IDX=MI;a=G.Finite(2,[[Q(c)for c in v]for v in P['nu_A_B_base25']])
 def c4parts(parts):return tuple(Q(parts[k][j])for j in range(25)for k in range(4))
 def scalar(v):return c4parts([v,[Q(0)]*25,[Q(0)]*25,[Q(0)]*25])
 def unit(mask):return tuple(Q(i==mask)for i in range(100))
 def add(x,y,c=Q(1)):return tuple(v+c*w for v,w in zip(x,y))
 pairs=P['physical150_tracepairs_base25'];rows=[18,19,28,29,60,67,72,63,64,68,73,74];columns=[]
 for row in rows:
  ca,cb=pairs['CA'][row],pairs['CB'][row];delta=[Q(x)-Q(y)for x,y in zip(ca[0],cb[0])];f=c4parts([delta,list(map(Q,ca[1])),[-Q(x)for x in cb[1]],[Q(0)]*25])
  for mask in range(4):columns.append(a.mul(f,unit(mask)))
 # All five requested old prefix firsts use only four already-saved pairs.
 # gA is the SAME physical f63 source function, not a guessed zero germ.
 jet=D['source2'];gmon=jet['source2_monomials'];gjet=jet['both57_source2']['CA']['pure_function2jets']['7790'];g=[Q(0)]*25
 for mon,c in zip(gmon,gjet):
  m=tuple(sorted(7 if i==8 else i for i in mon))
  if m in MI:g[MI[m]]+=Q(c)
 for i in range(8):
  f=[Q(0)]*25
  if i<7 and(i,)in MI:f[MI[(i,)]]=Q(-1)
  elif i==7:f=[-x for x in g]
  columns.append(a.mul(scalar(f),unit(1)))
 assert len(columns)==56
 cold=dict((i,Q(c))for i,c in H['actual_old_c0_sparse']);body=tuple(Q(0)for _ in range(100))
 for j,c in cold.items():body=add(body,columns[j],c)
 # m0*cold0=0 fixes all requested old first coefficients. Sector-B's
 # z+c*r² substitution preserves the mixed zi*r coefficient exactly.
 oldadj=list(map(Q,H['actual_old_m1_r_adjoint24']));prefixes=[];oldterm=Q(0)
 for row,w in enumerate(oldadj):
  if not w:continue
  sec,i=divmod(row,8);assert i<7;mon=tuple(sorted((i,7)));value=body[4*MI[mon]+[1,2,3][sec]];prefixes.append([row,value,w]);oldterm+=w*value
 N=D['natural']['actual_origin_data'];c0new=list(map(Q,N['all_five_new_cycle_origins44'][4]));oldn=[[Q(c)for c in r]for r in N['old_n0_56_by24']]
 X=D['X1'];xp=X['X0_pivot_columns'];raw=[{i:rec(v)for i,v in col}for col in X['all780_source_columns44_sparse_dual9']];a0=list(map(Q,H['all18_first_origin_cleanup_a0']))
 # Reconstruct the two requested exact Z1 rows from the accepted A1 and
 # trace multiplication, retaining nu-first carries. The lowerI8 is literal.
 A1=X['new_A12_by9_C36'];Z20z6=Q(0);Z42r=Q(0)
 for oldcol,c in cold.items():
  if oldcol>=48:continue
  row,mask=divmod(oldcol,4);v=[Q(0)]*100
  for j in range(4):v[j]=Q(A1[row][5][j])
  for var in[6,7]:
   for j in range(4):v[4*MI[(var,)]+j]=Q(A1[row][5][4*(var+1)+j])
  p=a.mul(tuple(v),unit(mask));Z20z6+=c*p[4*MI[(6,)]]
 assert all(Q(c)==0 for c in H['actual_top0_functional44'][:42])
 xp42r=sum((a0[k]*raw[j].get(42,{}).get(8,Q(0))for k,j in enumerate(xp)),Q(0));xp20z6=sum((a0[k]*raw[j].get(20,{}).get(7,Q(0))for k,j in enumerate(xp)),Q(0))
 Jr=Z42r+oldterm-xp42r+Q(20,3)*(Z20z6-xp20z6)
 # Mandatory actual first-relation lift and its source-bound quadratic part.
 rv=D['rightvalues'];labels=[tuple(x)for x in rv['convolution_labels']];lower=[[Q(c)for c in v]for v in rv['all35_projected_M2_lower_representatives']];lift=list(map(Q,H['exact_R35_constant_lift35']));correction=[sum((lift[p]*lower[p][labels.index((f,7))]for p in range(35)),Q(0))for f in range(3)]
 assert H['lifted_R35_U_sensitivity3']==H['lifted_R35_V_sensitivity3']==['0']*3
 factors=list(map(Q,H['quadratic_projection_factors_times_Jr']));wimage=[c*Jr-v for c,v in zip(factors,correction)]
 pair=D['filtration']['extra_actual_gr2_pair_columns'];u=[Q(pair[0].get(str(i),'0'))for i in range(3)];q2=[Q(pair[1].get(str(i),'0'))for i in range(3)];iv,det=invq([[u[i],q2[i]]for i in[1,2]]);coeff=mv(iv,[wimage[i]for i in[1,2]])
 assert all(coeff[0]*u[i]+coeff[1]*q2[i]==wimage[i]for i in range(3))
 acoef,bcoef=coeff;test=Q(1)+acoef*Q(11767,3652);forced=bool(test)
 data={'actual_old_head4_c0_support':[[j,c]for j,c in cold.items()],'all5_actual_old_prefix_r_values_with_adjoint':prefixes,'actual_old_projected_first_contribution':oldterm,
  'actual_Z20_z6_cold0':Z20z6,'literal_Z42_r_cold0_zero':Z42r,'actual_Xp42_r_a0':xp42r,'actual_Xp20_z6_a0':xp20z6,
  'actual_projected_natural_head4_first_target0_r':Jr,'actual_R35_quadratic_lift_correction3':correction,'actual_head4_z5_quadratic_image3':wimage,
  'known_actual_u_equals_e3r_image3':u,'known_actual_q2_equals_e4r2_image3':q2,'fixed_pair_inverse2':iv,'fixed_pair_determinant':det,'actual_w_equals_a_u_plus_b_q2_coefficients':[acoef,bcoef],
  'all3_quadratic_F2_image_coefficients_verified':True,'actual_qC_coefficient_of_r_kernel_vector':test,'kernel_vector':'w-a*v-b*q2; rimage=(1+a*11767/3652)*qC-(a*alpha+b)*g, qC=e4*z5*r,g=e4*r^3',
  'actual_FULL_HF55321_forced':forced,'actual_FULL_HF_if_forced':[5,5,3,2,1]if forced else None,
  'scope':'If test0 the deep oldalpha is still uncomputed; no higher callback or nilcube inference. This job evaluates only existing physical source2 and fixed matrices; every omitted coefficient is removed by the frozen exact adjoint.',
  'actual_R35_UV_sensitivities_zero_from_pinned_source_gate':True,'new_source_gauge_jet_or_higher_callback':False}
 assert nqq(data)+20000<35700
 p=out/'actual-head4-z5-map-and-F2-kernel-coefficient.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_cache_only_actual_head4_z5_quadratic_image','actual_qC_coefficient_nonzero':forced,'actual_full_HF_if_forced':[5,5,3,2,1]if forced else None,'all_source_and_three_F2_image_checks':True,'new_source_gauge_jet_or_higher_callback':False,'original_nilcube_or_generators_decided':False,'QQ_slots':nqq(data),'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
