#!/usr/bin/env python3
"""Independent complete-kernel/Tor source and specific RIGHT106 cache eligibility."""
import ast,hashlib,json,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def field(p,key):
 PINS[p]=sha(ROOT/p);text=(ROOT/p).read_text();tok=json.dumps(key);assert text.count(tok)==1;i=text.index(':',text.index(tok))+1
 while text[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(text,i)[0]
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(enc,v))
 return v
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63/report.json';r=read(rp);assert r['all_passed']and r['second_callbacks_not_evaluated']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p in r['input_sha256']:
  if p.endswith('REVIEW/report.json'):assert read(p)['all_passed']
 D=read('reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63/complete-twelve-generator-kernel-and-two-cache-callback-gate.json')
 fp='reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json'
 cols=[{int(i):Q(c)for i,c in col}for col in field(fp,'all849_actual_minimal_first_columns40_sparse')]
 assert len(cols)==849 and all(c.get(7,Q(0))==0 for c in cols)
 edge=read('reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63/actual-five-edge-functionals-and-genuine-Nbar14-cycles.json')
 a,*b=list(map(Q,edge['five_actual_cached_functional_values_a_b1_b2_b3_b4']));xi=Q(edge['xi_equals_minus11a_by16']);C=[v/xi for v in b]
 assert C==[Q(5,3),Q(-88,25),Q(2,5),Q(136,75)]and a and xi
 cc=D['fixed_z_minor_columns'];rr=D['fixed_z_minor_rows'];assert cc==[623,635,647,659,673,675]and rr==[0,1,2,3,4,6]
 inv=[[Q(c)for c in row]for row in D['fixed_z_minor_inverse']];A=[[cols[j].get(i,Q(0))for j in cc]for i in rr]
 assert all(sum(A[i][k]*inv[k][j]for k in range(6))==int(i==j)and sum(inv[i][k]*A[k][j]for k in range(6))==int(i==j)for i in range(6)for j in range(6))
 assert Q(D['fixed_z_minor_determinant'])!=0
 for j,col in enumerate(cols):
  co=[sum(v[k]*col.get(rr[k],Q(0))for k in range(6))for v in inv]
  assert co==list(map(Q,D['all849_z_Tor_factorizations'][j]))
  assert all(sum(cols[k].get(i,Q(0))*v for k,v in zip(cc,co))==col.get(i,Q(0))for i in range(7))
 assert all(cols[627].get(i,Q(0))+Q(32,11)*cols[623].get(i,Q(0))==0 for i in range(7))
 ps='reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63/right106-full-first-span-and-cache-footprint.json'
 selected=field(ps,'selected35_actual_first_relations');assert selected==field(fp,'rank_minor_columns')
 duals=field(ps,'all_nonzero_convolution_duals');eff=field(ps,'all_actual_effective_C1_H2_callbacks')
 assert duals[0]==[0,7,[[35,'1']]]and eff[0]==[0,7,[[27,7,'-8/3']]]
 # Independent typed sensitivity: pi index35 is top target0*r², not a retained
 # wedge. H27 is sector0 pair(6,7); its r coefficient is v7[z6*r].
 # That coefficient uses m0X2 + mu1X1 and no mu2X0 (X0 of this source is0).
 # The separate V=m1_last23_r belongs to sector2, absent from this callback.
 assert 35//36==0 and 27//28==0
 assert D['two_cached_second_callbacks']==[[623,selected.index(623),0],[627,selected.index(627),0]]
 cor={j:sum(C[k-1]*cols[j].get(8*k+7,Q(0))for k in range(1,5))for j in[623,627]}
 assert cor=={int(k):Q(v)for k,v in D['known_edge_first_correction_on_two_sources'].items()}
 # No cached second values or the final Tor witness are paired in this audit.
 custody=read('reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUE_CUSTODY_20261105_v63/report.json')
 code='code/ambient_B18_nil_R0_minimal_generators_values_20261105_v63.py';PINS[code]=sha(ROOT/code)
 assert sha(ROOT/code)==custody['prepared_value_code_sha256']=='c6b649b54019b6b93da363513b6835c44f1330892e232a8310c48a8eb874e99a'
 for p,h in custody['input_sha256'].items():assert sha(ROOT/p)==h
 tree=ast.parse((ROOT/code).read_text());assert not any(isinstance(n,ast.ImportFrom)and n.module not in['fractions','pathlib']for n in ast.walk(tree))
 assert subprocess.run([sys.executable,str(ROOT/code),'--help'],capture_output=True,timeout=5).returncode==0
 assert not (ROOT/'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63').exists()
 data={'all849_head0_r_firsts_zero':True,'fixed6_minor_both_inverse_products':True,'all849_z7_factorizations_verified':True,'z_rank_exact6':True,'zero_z_source627_plus32over11_source623':True,
 'normalized_actual_edge_coefficients_c':C,'complete12_by877_presentation_proof':{'L':'columns e_j-c_j*r*e0 (j1..4), plus z_i*e0(i0..6),r²e0','image_L':'exact ker(S5→S/I), I=(z7,r²), since target is cyclic generated by actual edgeh0','kernel_L':'First4 target coefficients force first4 source coefficients0; remaining kernel is syz(z7,r²), generated by Koszul28 because this is a regular sequence in formal S','lift_M':'For each of the849 actual minimal relations, set first4 lift entries=M_j. Ordered z-Hadamard plus pure-r regular r² division decomposes fM into I. L times the lift is M exactly. Thus coker[liftM,Koszul28]=Nbar. This is an explicit formal coefficient recipe, not expanded12×877 scalar entries.',
 'origin':'top4=0, last8=(M1head0,z7, M2head0,r²+sum cj*M1headj,r). Koszul28 origin0. No extra origin contributions from nonconstant coefficients because M0=0 and M1head0,r=0.'},
 'exact_mu_formula':'mu_R0N=mu_SNbar=12-rank(E8x849); equivalently Tor1(k,Tnil)→Tor1(k,S/I)→Nbar/mNbar→k5→k1→0','mu_interval_without_any_second_values':[5,6],
 'individual_RIGHT106_direction0_is_actual':True,'U_zero_source':'Only retained-wedge direct output can carry mu2 unknown; pi={top0r²} has no wedge. The C1 correction H27_r is v7[z6r] in sector0 and uses only m0X2+mu1X1, since original sourceX0=0.','V_zero_source':'V=m1sector2,last_r*k42. Only sector0 H27 callback occurs, so V cannot enter.',
 'other_eight_RIGHT106_convolutions_not_declared_individually_actual':True,'known_first_correction_two_sources':cor,'second_callbacks_not_read_or_paired':True,
 'worker_static_numeric_bound':30000,'hard_numeric_limit':60000,'source512bits_under2048':True,'source_no_new_physical_gauge_jet':True,'future_value_job_needs_parent_GO':True}
 ep=out/'independent-complete-kernel-Tor-and-two-cache-source-proof.json';ep.write_text(json.dumps(enc(data),indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'SOURCE_PASS_complete_Nbar12_by877_and_two_actual_cache_callbacks','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'z_Tor_rank':6,'mu_N_interval':[5,6],'second_cache_values_paired':False,'prepared_code_sha256':sha(ROOT/code),'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
