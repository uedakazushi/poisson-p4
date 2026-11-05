#!/usr/bin/env python3
"""Independent existing-cache dots, canonical-beta-safe kernel and source operators."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63_REVIEW'
PARENT='AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63'
PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(path):
 p=ROOT/path;PINS[path]=sha(p);return json.loads(p.read_text())
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(enc,v))
 return v
def sp(v):return {int(i):Q(c)for i,c in v if Q(c)}
def slots(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except ValueError:return 0
 if isinstance(v,dict):return sum(slots(x)for x in v.values())
 if isinstance(v,(tuple,list)):return sum(map(slots,v))
 return 0
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp=f'reports/{PARENT}/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 p=f'reports/{PARENT}/actual-five-edge-functionals-and-genuine-Nbar14-cycles.json';v=read(p);assert PINS[p]==r['evidence_sha256'][Path(p).name]
 sr=read('reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_20261105_v63_REVIEW/report.json');assert sr['all_passed']and not sr['actual_five_edge_values_paired']
 S=read('reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_20261105_v63_REVIEW/independent-five-functional-source-and-symbolic-kernel-proof.json')
 f=read('reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json')
 m=[[Q(x)for x in row]for row in f['m0_24_by44']];iv=[[Q(x)for x in row]for row in f['selected_minor_inverse']];n=[[Q(0)]*24 for _ in range(44)]
 for j,i in enumerate(f['selected_source_columns24']):n[i]=iv[j]
 assert all(sum(m[i][k]*n[k][j]for k in range(44))==int(i==j)for i in range(24)for j in range(24))
 assert all(row[43]==0 for row in m)
 c0=[Q(i==42)-sum(n[i][k]*m[k][42]for k in range(24))for i in range(44)]
 first=read('reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63/phase1-actual-fullnine-first-V-columns.json')
 L=[sp(x)for x in first['all44_actual_linear_columns27_sparse']]
 V=[[[col.get(t*9+var,Q(0))for col in L]for t in range(3)]for var in range(9)]
 q=read('reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase4-all-twenty-actual-quadratic-values.json')
 B2=[[Q(x)for x in row]for row in q['all44_actual_V2_pureplane_nine_values']]
 g=read('reports/AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63/actual-Cprime9-actions-and-genuine-G3-target.json')
 A=[[[Q(x)for x in col]for col in mat]for mat in g['actual_R0_actions9_by9_columns']]
 W=[[-A[var][target][7]for target in range(3)]for var in range(9)]
 a0=[Q(-1,75),Q(0),Q(-14,55)];b0=[Q(92,20625),Q(0),Q(3756,75625)];c0dual=[Q(0),Q(0),Q(1)]
 cz=[Q(2,25)*W[z][0]+W[z][2]/75 for z in range(7)]
 rho=[sum(b0[t]*B2[j][3*t]+a0[t]*B2[j][3*t+1]+c0dual[t]*B2[j][3*t+2]-W[7][t]*V[8][t][j]for t in range(3))+sum(cz[z]*V[z][0][j]for z in range(7))for j in range(44)]
 assert dict((i,x)for i,x in enumerate(rho)if x)==sp(S['independent_full_unnormalized_rho44'])
 s1=[[-Q(1375,288)*n[i][23]for i in range(44)]]+[[-n[i][j]for i in range(44)]for j in[0,4,5]]
 rows=[V[8][1]]+[rho]*4;cols=[c0]+s1
 ownterms=[[[i,x,col[i],x*col[i]]for i,x in enumerate(row)if x]for row,col in zip(rows,cols)]
 vals=[sum(term[-1]for term in terms)for terms in ownterms]
 assert vals==[Q(-1,10),Q(11,96),Q(-121,500),Q(11,400),Q(187,1500)]
 assert vals==list(map(Q,v['five_actual_cached_functional_values_a_b1_b2_b3_b4']))
 for saved,own in zip(v['all_five_exact_termwise_cache_products'],ownterms):
  assert [[int(a),Q(b),Q(c),Q(d)]for a,b,c,d in saved]==own
 a=vals[0];bs=vals[1:];xi=-Q(11,16)*a;assert xi==Q(11,160)==Q(v['xi_equals_minus11a_by16'])
 free=[1,2,3,4,5,7,8,9,10,11,12,13,14,15]
 K=[[Q(i==row)for row in free]for i in range(16)]
 for j in range(4):K[6][j]=-bs[j]/xi
 assert K==[[Q(x)for x in row]for row in v['actual_edge_kernel_basis16_by14_rows']]
 assert K[0]==[0]*14 and all(K[i][j]==int(k==j)for k,i in enumerate(free)for j in range(14))
 knownlong=[Q(0)]+bs+[Q(0),xi]+[Q(0)]*9;short=[a]+[Q(0)]*15
 assert all(sum(short[i]*K[i][j]for i in range(16))==0 and sum(knownlong[i]*K[i][j]for i in range(16))==0 for j in range(14))
 # Canonical beta is not supplied. Because K[0]=0, all beta coefficients vanish.
 assert v['canonical_long_coordinate_beta_of_head0']=='FORMAL_UNEVALUATED'
 assert v['actual_edge_pivot2symbolic_rows']==[[str(a),'0'],['beta',str(xi)]]
 assert Q(v['pivot_determinant_independent_beta'])==a*xi==Q(-11,1600)
 norm=[[Q(i==0)for i in range(16)],[Q(0)]+[b/xi for b in bs]+[Q(0),Q(1)]+[Q(0)]*9]
 assert all(sum(row[i]*K[i][j]for i in range(16))==0 for row in norm for j in range(14))
 op=v['genuine_Nbar14_barcycle_operator'];assert op['shapes']==[[44,16],[16,14],[44,14]]and op['root']==2
 assert op['nodes'][0]==['actual_external_matrix','CONNECTING16','matrix_roots.actual_ell0_cycles44_by16']
 assert op['nodes'][1]==['QQ_matrix_rows',enc(K)]and op['nodes'][2]==['matrix_product',[0,1]]
 binding=op['namespace_bindings']['CONNECTING16'];assert sha(ROOT/binding['path'])==binding['sha256']
 con=read(binding['path']);assert con['matrix_roots']['actual_ell0_cycles44_by16']==58 and con['genuine_kernel_dimension']==14
 assert v['retained_unit_free_rows14']==free and v['genuine_Nbar_dimension']==14
 data={'independent_n0_576_products':True,'independent_five_cache_products':ownterms,'independent_actual_a_b1_b2_b3_b4':vals,'xi':xi,
       'canonical_beta':'FORMAL; not evaluated or zero-filled','canonical_symbolic_pivot_determinant':a*xi,'normalized_target_basis':['edge(head0)','edge(rhead0)'],
       'actual_normalized_edge2_by16':norm,'actual_kernel16_by14':K,'all28_edge_kernel_products_zero_independent_of_beta':True,'retained_I14_rows':free,
       'all14_actual_barcycle_source_operator':op,'genuine14dim_kernel_equals_Nbar':'accepted exact connecting0→N/ellN→Tnil16→genuineG3→0 and rank2',
       'no_full_R0_nil_representative_or_nilideal_product_evaluated':True}
 owned=slots([r,sr,S,f,m,iv,n,c0,first,L,V,q,B2,g,A,W,rho,s1,ownterms,vals,K,norm,data,v,con]);assert owned<60000,owned
 ep=out/'independent-five-values-full-symbolic-edge-and-fourteen-barcycles.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_independent_actual_genuine_Nbar14_barcycle_basis','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'actual_five_functionals':list(map(str,vals)),'actual_xi':str(xi),'canonical_beta_evaluated':False,'genuine_dimension':14,'owned_numeric_slots':owned,'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps(result))
if __name__=='__main__':main()
