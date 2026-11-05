#!/usr/bin/env python3
"""Independent two-cache/carry witness: actual R0-module generator count, not D0 count."""
import hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def field(p,key):
 PINS[p]=sha(ROOT/p);s=(ROOT/p).read_text();tok=json.dumps(key);assert s.count(tok)==1;i=s.index(':',s.index(tok))+1
 while s[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(s,i)[0]
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(enc,v))
 return v
def det_bareiss(a):
 den=math.lcm(*(v.denominator for row in a for v in row));m=[[int(v*den)for v in row]for row in a];n=len(m);prev=1;sign=1;transcript=[]
 for k in range(n-1):
  p=next(j for j in range(k,n)if m[j][k])
  if p!=k:m[k],m[p]=m[p],m[k];sign=-sign
  pivot=m[k][k];checks=0
  for i in range(k+1,n):
   for j in range(k+1,n):
    numerator=pivot*m[i][j]-m[i][k]*m[k][j];assert numerator%prev==0;m[i][j]=numerator//prev;checks+=1
   m[i][k]=0
  transcript.append([k,p,pivot,prev,checks]);prev=pivot
 return F(sign*m[-1][-1],den**n),transcript
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 source=read('reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63_REVIEW/report.json');assert source['all_passed']and not source['second_cache_values_paired']
 S=read('reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63/complete-twelve-generator-kernel-and-two-cache-callback-gate.json')
 V=read('reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63/actual-two-cache-witness-and-five-minimal-bar-generators.json')
 scope=read('reports/AMBIENT_B18_NIL_R0_GENERATION_SCOPE_20261105_v63/report.json');assert scope['all_passed']
 first='reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json'
 cols=[{int(i):F(c)for i,c in col}for col in field(first,'all849_actual_minimal_first_columns40_sparse')]
 cache=read('reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63/phase1-all35-nine-convolution-projections.json')
 assert cache['convolution_labels'][0]==[0,7]
 q={}
 for rel,index,direction in S['two_cached_second_callbacks']:
  assert cache['all35_source_input'][index]['source_relation']==rel and direction==0
  q[rel]=F(cache['all35_projected_M2_lower_representatives'][index][direction])
 assert q=={623:F(6875,2072),627:F(-2500,259)}
 C=[F(5,3),F(-88,25),F(2,5),F(136,75)]
 carry={rel:sum(C[j-1]*cols[rel].get(j*8+7,F(0))for j in range(1,5))for rel in[623,627]}
 qsum=q[627]+F(32,11)*q[623];csum=carry[627]+F(32,11)*carry[623];w=qsum+csum
 assert qsum==0 and csum==w==F(-75,11)
 assert V['two_cached_second_contributions_sum']==str(qsum)and V['edge_first_contributions_sum']==str(csum)and V['nonzero_last_Tor_row_witness']==str(w)
 assert [[j,str(x)]for j,x in sorted(carry.items())]==V['two_known_edge_first_corrections']
 cc=S['fixed_z_minor_columns'];rr=S['fixed_z_minor_rows'];A=[[cols[j].get(i,F(0))for j in cc]for i in rr];det6,trace=det_bareiss(A)
 assert det6==F(S['fixed_z_minor_determinant'])and det6
 assert all(cols[627].get(i,F(0))+F(32,11)*cols[623].get(i,F(0))==0 for i in range(7))
 det7=det6*w;assert det7==F(V['fixed_rank7_certificate']['block_triangular_determinant'])and det7
 assert V['actual_complete_origin_relation_rank']==7 and V['actual_mu_S_Nbar']==5 and V['actual_mu_R0_Nil_D0']==5
 data={'actual_two_oldcache_q623_q627':q,'independent_known_first_carry623_627':carry,'q_sum':qsum,'carry_sum':csum,'actual_last_Tor_witness':w,
 'fixed6_minor_det':det6,'fixed6_integer_Bareiss_exact_division_transcript':trace,'fixed7_actual_source_columns':cc+[627],'fixed7_actual_rows':rr+[7],'fixed7_determinant':det7,
 'fixed7_source_proof':'Replace final source627 column by source627+(32/11)*source623, a determinant-one column operation because623 is already selected. All7 z entries of final column vanish; bottom is w=-75/11. Therefore det7=det6*w without needing the six other bottom entries. Full7z-row rank6 bounds entire8row Tor rank≤7; thus exact rank7.',
 'actual_mu_R0_N':5,'actual_mu_S_Nbar':5,'minimal_five_bar_generators':V['actual_minimal_bar_generators'],
 'completeness_and_minimality':'Complete formal12x877 presentation has zero top4 origin rows. Seven independent origin constraints among I*h0 leave one slot. Fixed7 rows exclude z5; full z-factorization shows e_z5 is nonzero in the quotient. Thus correctedheads4 and z5*h0 form a minimal S-generator set. They lift to N by exact connecting; any such lifts generate finiteR0-module N by Nakayama along ell.',
 'same_regular_R0_module_count_not_D0_ideal_count':True,'mu_D_N_unchanged':14,'D0_ideal_generator_upper_bound':5,'D0_ideal_minimal_generator_count':'NOT_DETERMINED','full_formal_nil_lifts':'not constructed; existence only','new_source_gauge_jet_or_second_coefficients_evaluated':False,'nilcube_or_index_improvement':False}
 ep=out/'independent-two-cache-Tor7-and-actual-R0-five-generators.json';ep.write_text(json.dumps(enc(data),indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_independent_actual_R0_module_nil_generator_count5','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'actual_mu_R0_N':5,'actual_mu_S_Nbar':5,'D0_ideal_minimal_count_determined':False,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
