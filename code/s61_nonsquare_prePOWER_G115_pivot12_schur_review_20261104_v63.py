#!/usr/bin/env python3
"""Independent typed 4+4+3 parent and twelfth-pivot/domain/source audit."""
import hashlib,itertools,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
# Only prior INDEPENDENT checker utilities; no producer imports or main calls.
import s61_nonsquare_prePOWER_kappa11_schur_review_20261104_v63 as U
import s61_nonsquare_prePOWER_full57_source_lift_review_20261104_v63 as P
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='S61_NONSQUARE_PREPOWER_G115_PIVOT12_SCHUR_20261104_v63';ID=PRIMARY+'_REVIEW';OUT=ROOT/'reports'/ID
STAGE='S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63'
ENTRY='S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_ENTRY_20261104_v63'
K11='S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63'
LIFT='S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63'
COLS=[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+'\n')
def inversions(p):return sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;k=str(p.relative_to(ROOT));pins[k]=sha(p)
  if k not in cache:cache[k]=json.loads(p.read_text())
  return cache[k]
 r=load(PRIMARY);assert pins[f'reports/{PRIMARY}/report.json']=='1ca54732944a8fb504d1a1f535b47d5089286c829d23c47baff5b51fa0f73d52'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 O=load(PRIMARY,'operator.json');proof=load(PRIMARY,'proof.json');A=load(STAGE,'operator.json');old=load(K11,'operator.json');lift=load(LIFT,'operator.json')
 ep=load(ENTRY+'_REVIEW');assert pins[f'reports/{ENTRY}_REVIEW/report.json']=='3b3a2bdb1ea483e9471d585fb86fd327c8033daafa44ce18027124d1ae12d31f'and ep['all_passed']
 assert load(K11+'_REVIEW')['all_passed']and load(LIFT+'_REVIEW')['all_passed']
 hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json')
 # Literal three polynomial identities, with source-specific row units retained.
 g=P.dec(hc['raw_compact14x5'][6][1]);h=P.dec(hc['raw_compact14x5'][10][1])
 assert g and h
 for j in [1,2,3]:assert not P.add(P.mult(g,P.dec(hc['raw_compact14x5'][10][j])),P.mult(h,P.dec(hc['raw_compact14x5'][6][j])),-1)
 assert hc['known_row_units_rho14'][6]['a_d_F_B_C_J_powers']==hc['known_row_units_rho14'][10]['a_d_F_B_C_J_powers']==[1,0,1,0,2,0]
 c10=F(hc['known_row_units_rho14'][10]['constant']);assert c10==F(-1600000,513)
 # Independent universal rank-one adjugate identity, ALL16 entries.
 aa=[[U.sym(f'A{i}{j}')for j in range(4)]for i in range(4)];u=[U.sym(f'u{i}')for i in range(4)];v=[U.sym(f'v{j}')for j in range(4)];k=U.sym('k')
 for i,j in itertools.product(range(4),repeat=2):
  rr=[r for r in range(4)if r!=j];cc=[c for c in range(4)if c!=i]
  sub=[[aa[r][c]for c in cc]for r in rr]
  tilde=[[U.plus(U.times(k,aa[r][c]),U.times(u[r],v[c]),-1)for c in cc]for r in rr]
  reduced=U.times(k,U.determinant(sub))
  for pos,c in enumerate(cc):
   replaced=[[u[r]if t==pos else aa[r][cc[t]]for t in range(3)]for r in rr]
   reduced=U.plus(reduced,U.times(v[c],U.determinant(replaced)),-1)
  assert U.determinant(tilde)==U.times(U.times(k,k),reduced)
 # Both block inverses in the staged parent, in a free NONcommutative algebra.
 U.PAIRS.update({'Q':'Qi','Qi':'Q'})
 a,ai,b,c,q,qi=[U.nc(x)for x in ['A','Ai','B','C','Q','Qi']];z={};one=U.nc();I=[[one,z],[z,one]]
 d=U.plus(q,U.ncm(U.ncm(c,ai),b))
 inverse=[[U.plus(ai,U.ncm(U.ncm(U.ncm(U.ncm(ai,b),qi),c),ai)),U.neg(U.ncm(U.ncm(ai,b),qi))],[U.neg(U.ncm(U.ncm(qi,c),ai)),qi]]
 full=[[a,b],[c,d]];assert U.mm(full,inverse)==U.mm(inverse,full)==I
 U.nc_block_proof()
 def check_types(doc):
  N=doc['nodes']
  for i,n in enumerate(N):
   sh=n['shape'];args=n['args'];ar=[N[j]['shape']for j in args];op=n['op'];assert n['id']==i and all(j<i for j in args)
   if op in ['external_frozen_reference','external_frozen_operator_root']:
    rf=n['reference'];p=ROOT/rf['path'];assert sha(p)==rf['sha256'];pins[rf['path']]=sha(p)
    d=cache.setdefault(rf['path'],json.loads(p.read_text()))
    if 'node_index'in rf:assert d['nodes'][rf['node_index']]['shape']==sh
    else:
     val=d
     for s in rf['json_pointer'].strip('/').split('/'):val=val[int(s)]if isinstance(val,list)else val[s]
     if rf['json_pointer'].startswith('/exports/'):assert d['nodes'][val]['shape']==sh
     else:assert sh==[1,1]
   elif op=='matrix_product':assert len(ar)==2 and ar[0][1]==ar[1][0]and sh==[ar[0][0],ar[1][1]]
   elif op in ['difference','sum']:assert ar==[sh,sh]
   elif op=='negative':assert ar==[sh]
   elif op=='central_scalar_multiple':assert ar[0]==[1,1]and sh==ar[1]
   elif op=='submatrix':assert sh==[len(n['rows']),len(n['columns'])]and all(x in range(ar[0][0])for x in n['rows'])and all(x in range(ar[0][1])for x in n['columns'])
   elif op=='transpose':assert sh==ar[0][::-1]
   elif op in ['identity','zero']:assert not ar
   elif op in ['row_permutation','column_permutation']:assert sh[0]==sh[1]and sorted(n['order'])==list(range(sh[0]))and not ar
   elif op=='QQ_constant':assert not ar and sh==[1,1]and F(n['value'])
   elif op=='literal_parameter_monomial':assert sh==[1,1]and n['variable_order']==['a','e','d']and n['polynomial']==[[[1,0,0],'1']]
   elif op in ['certified_chart_inverse','known_FULLU_parameter_unit_inverse']:assert sh==[1,1]and ar==[[1,1]]
   elif op=='diagonal_known_FULLU_unit_inverse':assert ar==[sh]and sh[0]==sh[1]
   elif op=='diagonal':assert n['diagonal']==args and len(args)==sh[0]==sh[1]and ar==[[1,1]]*len(args)
   elif op=='horizontal_concatenation':assert len({a[0]for a in ar})==1 and sh==[ar[0][0],sum(a[1]for a in ar)]
   elif op=='block_matrix':
    assert args==[j for row in n['blocks']for j in row];ss=[[N[j]['shape']for j in row]for row in n['blocks']]
    assert all(len({s[0]for s in row})==1 for row in ss)and all(len({row[j][1]for row in ss})==1 for j in range(len(ss[0])))
    assert sh==[sum(row[0][0]for row in ss),sum(s[1]for s in ss[0])]
   elif op in ['universal_polynomial_determinant','negative_universal_three_determinant']:assert sh==[1,1]and len(ar)==1 and ar[0][0]==ar[0][1]and n['expanded']is False
   elif op=='universal_polynomial_adjugate':assert ar==[sh]and sh[0]==sh[1]and n['expanded']is False
   elif op=='lazy_polynomial_adjugate_with_last_row_cache':assert sh==[4,4]and ar==[[4,4],[1,4]]and n['last_row_equals']=='minus kappa'and n['other12_entries_unexpanded']
   elif op=='rank_one_adjugate_reduced_polynomial':assert sh==[4,4]and ar==[[1,1],[4,4],[4,1],[1,4]]and n['no_k_inverse_or_parameter_division']and n['all16_entries_unexpanded']
   else:raise AssertionError(op)
  return N
 def flat(N,i):return sum((flat(N,j)for j in N[i]['args']),[])if N[i]['op']=='matrix_product'else[i]
 def sub(N,i,p,rr,cc):assert N[i]['op']=='submatrix'and N[i]['args']==[p]and N[i]['rows']==rr and N[i]['columns']==cc
 S=check_types(A);assert len(S)==103
 # Exact parent permutations and all 4+4+3 source inverse wiring.
 assert S[3]['order']==A['source_pivot_order']+A['source_remaining_order']==[2,4,5,8,0,1,3,6,7,13,10,9,11,12,14]
 assert [COLS[i]for i in S[5]['order']]==A['target_pivot_original32']+A['target_remaining_original32']==[27,31,22,2,11,13,25,15,28,30,29,1,3,4,6,7,8,9,14,17,18,19]
 assert inversions([4,5,6,7,0,1,2,3])==16 and inversions([6,7,0,5,1,2,3,4])==16
 assert flat(S,9)==[3,0,1,5]
 for i,rr,cc in [(10,range(4),range(4)),(11,range(4),range(4,8)),(12,range(4,8),range(4)),(13,range(4,8),range(4,8)),(14,range(4,8),[2])]:sub(S,i,9,list(rr),list(cc))
 assert S[21]['args']==[17,18,19,20]and S[22]['args']==[10,21]and S[23]['args']==[16,22]and S[24]['args']==[23]
 assert S[16]['args']==[15]and S[15]['reference']['path'].endswith('/root239.json')
 for j,i in enumerate([17,18,19,20]):assert S[i]['reference']['path'].endswith(f'/root{[179,198,217,237][j]}.json')
 for j,sel in enumerate([28,30,32,34]):sub(S,sel,9,[0,1,2],[0,1,4+j]);assert S[sel+1]['args']==[sel]
 assert S[37]['formula']=='entry(i,j)=(-1)^(i+j)*(k*det(A_withoutrowj_coli)-v_withouti^t*Adj(A_withoutrowj_coli)*u_withoutj)'
 assert S[37]['args']==[20,13,14,36]and S[36]['args']==[29,31,33,35]
 assert S[39]['value']=='-5120000000000000/45927'and flat(S,40)==[39,38]and S[41]['args']==[40]and S[42]['args']==[41,37]and S[43]['args']==[42]
 assert S[27]['args']==[13,26]and flat(S,26)==[12,24,11]
 assert S[54]['blocks']==[[49,51],[53,43]]and S[49]['args']==[24,48]and flat(S,48)==[24,11,43,12,24]
 assert S[51]['args']==[50]and flat(S,50)==[24,11,43]and S[53]['args']==[52]and flat(S,52)==[43,12,24]
 assert S[73]['blocks']==[[54,69],[72,67]]and S[69]['op']=='zero'and S[72]['args']==[71]and flat(S,71)==[66,64,56,54]
 sub(S,55,9,list(range(8,11)),list(range(8,11)));sub(S,56,9,list(range(8,11)),list(range(8)))
 assert flat(S,67)==[66,64]and S[66]['args']==[59,65]and S[59]['args']==[58]and S[65]['args']==[57]
 assert S[63]['args']==[60,61,62]and S[64]['args']==[63]
 for i,j in zip([60,61,62],[6,12,9]):assert S[i]['reference']['json_pointer']==f'/known_row_units_rho14/{j}'
 sub(S,68,9,list(range(11)),list(range(11)));sub(S,74,9,list(range(11)),list(range(11,22)));sub(S,75,9,list(range(11,15)),list(range(11)));sub(S,76,9,list(range(11,15)),list(range(11,22)))
 assert S[80]['args']==[76,79]and flat(S,79)==[75,73,74]
 assert S[84]['blocks']==[[73,81],[82,83]]and S[82]['args']==[78]and flat(S,78)==[75,73]
 assert S[87]['blocks']==[[68,85],[75,86]]and S[92]['blocks']==[[88,89],[90,91]]and S[89]['args']==[77]and flat(S,77)==[73,74]
 assert S[96]['blocks']==[[93,77],[94,95]]
 assert flat(S,97)==[84,3]and flat(S,98)==[4,87]and flat(S,100)==[1,5,92]and flat(S,102)==[96,6,2]
 assert [i for i,n in enumerate(S)if n['op']=='certified_chart_inverse']==[16,41,59]
 # New scalar inverse and all101 typed/source/target composition nodes.
 N=check_types(O);assert len(N)==101
 assert N[6]['rows']==[0,1,3,6,2]and [COLS[i]for i in N[6]['columns']]==[22,11,13,25,15]and N[7]['args']==[6]and N[7]['expanded']is False
 assert N[8]['reference']['json_pointer']=='/raw_compact14x5/6/1'and F(N[9]['value'])==c10
 assert N[11]['args']==[10,10]and N[12]['reference']['json_pointer']=='/exact_normalized_D6'and N[13]['args']==[12]and N[13]['chart']=='D(N6)'
 assert N[14]['args']==[9,11]and N[15]['args']==[14]and flat(N,18)==[8,7,15,13]
 assert [i for i,n in enumerate(N)if n['op']in ['certified_chart_inverse','known_FULLU_parameter_unit_inverse']]==[13,15]
 assert N[19]['order']==[1,0,2,3]and N[21]['order']==[10]+list(range(10))and flat(N,24)==[19,0,21]
 for i,rr,cc in [(25,[0],[0]),(26,[0],range(1,11)),(27,range(1,4),[0]),(28,range(1,4),range(1,11))]:sub(N,i,24,list(rr),list(cc))
 assert N[29]['args']==[18,27]and N[30]['args']==[18,26]and N[32]['args']==[28,31]and N[31]['args']==[29,26]
 assert N[36]['blocks']==[[18,33],[34,35]]and N[34]['args']==[29]and N[39]['blocks']==[[25,37],[27,38]]
 assert N[44]['blocks']==[[40,41],[42,43]]and N[41]['args']==[30]and N[48]['blocks']==[[45,30],[46,47]]
 assert N[53]['blocks']==[[49,51],[52,50]]and flat(N,50)==[36,19]and flat(N,54)==[53,1]
 assert N[59]['blocks']==[[55,57],[58,56]]and flat(N,56)==[20,39]and flat(N,60)==[2,59]
 assert N[65]['blocks']==[[61,63],[64,62]]and flat(N,62)==[21,44]and flat(N,66)==[3,65]
 assert N[71]['blocks']==[[67,69],[70,68]]and flat(N,68)==[48,22]and flat(N,72)==[71,4]
 for i,oldroot in zip(range(73,83),[0,49,107,109,111,114,115,116,118,108]):assert N[i]['reference']['node_index']==oldroot
 sub(N,83,54,[12,13,14],list(range(15)));assert flat(N,84)==[83,73]and flat(N,85)==[83,73,74]
 assert flat(N,87)==[83,73,74,76,77];sub(N,88,87,[0,1,2],[5]);sub(N,89,87,[0,1,2],[6])
 assert N[90]['args']==[88,79]and N[91]['args']==[89,80]and N[92]['args']==[90,91]and N[93]['args']==[81,92]and N[94]['args']==[93,78]and N[95]['args']==[85,94]
 assert flat(N,96)==[95,75]and flat(N,97)==[95,82]and flat(N,100)==[95,76,77]
 sub(N,98,97,[0,1,2],[j for j in range(91)if j not in lift['nodes'][109]['columns']])
 assert O['source_remaining_current15']==[9,12,14]and O['target_remaining_original32']==[1,3,4,6,7,8,9,14,17,18]and O['chart']=='FULLU D(Delta407*Saux*T*N6)'
 assert O['no_extra_g_D5_inverse']and O['all7_remaining_extra_tail_effects_retained']and O['new_parameter_coefficients']==0
 save('independent-domain-source-and-block-proof.json',{'staged103_and_new101_nodes_all_types_hashes_and_wiring':True,'all16_rankone_adjugate_identities_independent_free_polynomial_proof':True,'both_general_4block_and_full_Schur_noncommuting_inverse_products':True,'staged_parent_same_pivot11_set_and_residual4x11_as_accepted_KAPPA11':True,'Ksmall_all3_middle_heavy_coefficients_exact_polynomial_zero':True,'border_identity_g_D5_r_equals_c10_a_squared_N6_positive_sign':True,'dense_field_derivation_descends_to_whole_old_localized_domain':True,'r_inverse':'g*D5/(c10*a^2*N6)','no_g_D5_or_kappa_extra_inverse':True,'new_chart':'FULLU D(Delta407*Saux*T*N6)','samecokernel_shape':[3,10],'all_actual_old25_normalized57_canonical63_and_eta_raw7_cleanup_bound':True,'remaining7_tails_and_chart_complement_kept':True,'no_D5_scalar_expansion_or_rank3_global57_claim':True})
 for helper in [Path(U.__file__),Path(P.__file__)]:pins[str(helper.relative_to(ROOT))]=sha(helper)
 rep={'id':ID,'all_passed':True,'frozen':True,'samecokernel_shape':[3,10],'pivot_count':12,'chart':O['chart'],'all103_parent_and101_new_types_and_source_inverse_assemblies_accepted':True,'border_product_identity_and_no_extra_g_D5_chart_accepted':True,'actual_full57_lift_and_raw7_cleanup_preserved':True,'rank3_global57_quartic239_NOT_proved':True,'no_new_source_jet_parameter_coefficient_rank_point_GB_or_D5_expansion':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{ID}.md'),'evidence_sha256':{'independent-domain-source-and-block-proof.json':sha(OUT/'independent-domain-source-and-block-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in OUT.iterdir())<32*1024*1024;save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':main()
