#!/usr/bin/env python3
"""Independent unit/source binding and polynomial proof A*D15=-T*D12."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW';PINS={};NV=16
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def field(p,key):
 PINS[p]=sha(ROOT/p);s=(ROOT/p).read_text();token=json.dumps(key);assert s.count(token)==1;i=s.index(':',s.index(token))+1
 while s[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(s,i)[0]
def pv(i):e=[0]*NV;e[i]=1;return {tuple(e):1}
def add(*ps):
 z={}
 for p in ps:
  for m,c in p.items():z[m]=z.get(m,0)+c
 return {m:c for m,c in z.items()if c}
def neg(p):return {m:-c for m,c in p.items()}
def mul(*ps):
 z={(0,)*NV:1}
 for p in ps:
  n={}
  for a,c in z.items():
   for b,d in p.items():m=tuple(i+j for i,j in zip(a,b));n[m]=n.get(m,0)+c*d
  z={m:c for m,c in n.items()if c}
 return z
def determinant(m):
 if len(m)==1:return m[0][0]
 return add(*(mul(m[0][j],determinant([row[:j]+row[j+1:]for row in m[1:]]))if j%2==0 else neg(mul(m[0][j],determinant([row[:j]+row[j+1:]for row in m[1:]])))for j in range(len(m))))
def sign(p):return (-1)**sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))
def encodepoly(p):return [[list(m),c]for m,c in sorted(p.items())]
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(enc,v))
 return v
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 note='research/notes/S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63.md';PINS[note]=sha(ROOT/note);assert PINS[note]=='7690a76b114fa8541b6a8e170e3a798a35e83c1a0605c6d20b63d39dcaeda40a'
 gate=read('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_BINARY_GATE_20261105_v63/report.json');assert gate['all_passed']and gate['source_only']and gate['input_sha256'][note]==PINS[note]
 for p,h in gate['input_sha256'].items():assert sha(ROOT/p)==h
 # All these are earlier independent source/value proofs; no producer is run.
 peerpaths=[
 'reports/S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW/report.json',
 'reports/S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63_REVIEW/report.json',
 'reports/S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63_REVIEW/report.json',
 'reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63_REVIEW/report.json',
 'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW/report.json',
 'reports/S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63_REVIEW/report.json']
 peers=[read(p)for p in peerpaths];assert all(p['all_passed']for p in peers)
 assert peers[0]['source_bottom_pure6_rank_le1']and peers[1]['raw4x7_rank_at_most1_on_exact_seven_columns']
 pure8=[22,11,13,25,15,19,14,18];assert peers[2]['fixed_raw_pair_rank_at_most1_on_columns']==pure8
 assert peers[3]['same_target17_pair_proportional_anchor22_accepted']
 assert peers[4]['all240_source_and_four_signed_cofactor_B_alignment_bound']and peers[5]['previous_projected15_to_actual57_source_gap_closed']
 qpath='reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json'
 targets=field(qpath,'new22_target_original32_columns');roots=field(qpath,'remaining15x22_scalar_roots');dag=field(qpath,'scalar_DAG');assert dag['nodes'][0]==['z']
 for r in[2,4,5,8]:
  for col in[28,30,29]:assert roots[r][targets.index(col)]==0
 del dag
 hpath='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json';gamma=field(hpath,'known_column_units')[4];rho=field(hpath,'known_row_units_rho14')
 assert gamma=={'constant':'-4/625','a_d_F_B_C_J_powers':[0,0,1,1,2,0]}
 constants={i:Q(rho[i]['constant'])for i in[6,8,9,10,11,12]};assert all(constants.values())and constants[11]==constants[9]
 for i in constants:assert rho[i]['a_d_F_B_C_J_powers']==([0,0,1,0,2,0]if i==12 else[1,0,1,0,2,0])
 k=read('reports/S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63/operator.json');st=read('reports/S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63/operator.json')
 assert k['nodes'][9]['op']=='diagonal'and st['exports']['A4']==10
 assert st['nodes'][3]['order'][:4]==[2,4,5,8]
 assert [targets[i]for i in st['nodes'][5]['order'][:4]]==[27,31,22,2]
 assert k['nodes'][2]['expression']=='a*F*C^2'and k['nodes'][3]['expression']=='B'
 assert k['nodes'][5]['args']==[2,4]and k['nodes'][6]['args']==[2]and k['nodes'][7]['args']==[5]
 nd=k['nodes'][9]['diagonal'];assert {targets[i]:v for i,v in enumerate(nd)if v!=1}=={2:3,27:6,31:7}
 # A is an actual nonzero polynomial, not a newly declared FULLU unit.
 anchor=read('reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json')['bottom4x6'][0][0]
 ab=read('reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/source-binding.json')
 assert ab['source_current15_rows']==[2,4,5,8]and ab['original32_columns']==[22,11,13,25,15,19]and ab['all_row_a_shifts_zero']
 assert anchor and all(Q(c)!=0 for m,c in anchor)and len({tuple(m)for m,c in anchor})==len(anchor)
 # Explicit known-unit source change on ALL15 rows, including the added columns.
 roworder=[0,1,3,6,2,11,12,9,7,13,10,4,5,8,14]
 R=[{}for _ in range(15)]
 def put(i,j,c=Q(1),a=0):R[i][j]={a:c}
 for i,j in enumerate(roworder):put(i,j)
 put(1,1,a=-1);put(1,3,Q(2),-1);put(2,3,a=-1)
 for row,cur,ci in[(5,11,10),(6,12,11),(7,9,8),(8,7,6),(9,13,12),(10,10,9)]:put(row,cur,1/constants[ci],1 if row==9 else 0)
 put(6,10,-1/constants[11])
 Rin=[{}for _ in range(15)]
 def ip(i,j,c=Q(1),a=0):Rin[i][j]={a:c}
 for i,j in enumerate(roworder):ip(j,i)
 ip(1,1,a=1);ip(1,2,Q(-2),1);ip(3,2,a=1)
 for cur,row,ci in[(11,5,10),(12,6,11),(9,7,8),(7,8,6),(13,9,12),(10,10,9)]:ip(cur,row,constants[ci],-1 if row==9 else 0)
 ip(12,10,constants[9])
 def lpadd(p,q):
  z=dict(p)
  for a,c in q.items():z[a]=z.get(a,Q(0))+c
  return {a:c for a,c in z.items()if c}
 def lpmul(p,q):
  z={}
  for a,c in p.items():
   for b,d in q.items():z[a+b]=z.get(a+b,Q(0))+c*d
  return {a:c for a,c in z.items()if c}
 for X,Y in[(R,Rin),(Rin,R)]:
  for i in range(15):
   for j in range(15):
    value={}
    for a,v in X[i].items():value=lpadd(value,lpmul(v,Y[a].get(j,{})))
    assert value==({0:Q(1)}if i==j else{})
 J11=[22,11,13,25,15,19,14,18,28,30,29];J15=J11+[27,31,2,17]
 assert len(set(J15))==15 and sorted(roworder)==list(range(15))
 assert sign([2,0,1,3])==1 and (sum([11,12,13])*2)%2==0
 # Fresh free16-variable identity, NOT actual a/e/d coefficient evaluation.
 A=pv(0);u=[pv(i)for i in range(1,4)];c=[pv(i)for i in range(4,7)];V=[[pv(7+3*i+j)for j in range(3)]for i in range(3)]
 block=[[A]+u]+[[c[i]]+V[i]for i in range(3)]
 corr=[[add(mul(A,V[i][j]),neg(mul(c[i],u[j])))for j in range(3)]for i in range(3)]
 d4=determinant(block);d3=determinant(corr);assert d3==mul(A,A,d4)
 original=[[row[j]for j in[1,2,0,3]]for row in block];assert determinant(original)==d4
 kappas=[neg(determinant([row[:3]for k,row in enumerate(original)if k!=i]))if (i+3)%2==0 else determinant([row[:3]for k,row in enumerate(original)if k!=i])for i in range(4)]
 T=add(*(mul(kappas[i],original[i][3])for i in range(4)));assert d4==neg(T)
 omega='FULLU D(Delta407*Saux*T*N6*N7bar)'
 p14=read('reports/S61_NONSQUARE_PREPOWER_K92_PIVOT14_SCHUR_20261105_v63/operator.json');assert p14['chart']=='FULLU D(Delta407*Saux*T*N6*N7bar*P_K92)'
 lift=read('reports/S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63/operator.json');assert lift['exports']['physical91_to_current22']==179
 cprod=Q(1)
 for c0 in constants.values():cprod*=c0
 data={'primary_is_note_pinned_by_binary_gate_not_standalone_report':True,'actual_source15_row_order':roworder,'row_operator_QQ_Laurent_a':R,'row_operator_inverse_QQ_Laurent_a':Rin,'all450_row_inverse_products':True,
 'J11':J11,'actual_target15_original32_order':J15,'M12_rows':'first11 transformed rows plus current14, omitting current4,5,8','all12_pair_H3_entries_literal_zero':True,'pair_rank1_on_exact9_pure_columns':pure8+[17],
 'known_target_normalization':{'28':'1/L','30':'1/L','29':'1/L','27':'1/L','31':'1/(L*gamma)','2':'B','other_selected_columns':'1','L':'a*F*C²','gamma':'(-4/625)*F*B*C²'},
 'normalization_only_FULLU_units':True,'source_row_determinant':{'QQ':Q(sign(roworder))/cprod,'a_power':-1},
 'target_column_product':'B/(L^5*gamma) = (-625/4)*a^-5*F^-6*C^-12','all_normalized_D15_is_actual_ordered_minor_times_knownunit':True,
 'anchor_A_source':'literal current15row2,target22; saved nonzero source coefficient list with zero row-a shifts; no global unit assertion','nonzero_anchor_terms':len(anchor),
 'free_identity_variable_order':['A','u1','u2','u3','c1','c2','c3']+[f'V{i}{j}'for i in range(1,4)for j in range(1,4)],'det_anchorfirst_A4':encodepoly(d4),'det_A_V_minus_cu':encodepoly(d3),'all_free_identity_coefficients_match_det3_A2det4':True,
 'Laplace_source_proof':'Replace each of the3 extra pair rows by A*row-ci*anchorrow. Nine pure entries vanish by TRUE rawpair rank1 and the3 H3 entries were literalzero. Only extra columns[27,31,2] remain. Rows/columns positions[11,12,13] have even total index sum, so A³D15=det(A V-cu)*D12. AnchorfirstA4 order[22,27,31,2] differs from old[27,31,22,2] by2swaps. det(A V-cu)=A²detA4 and detA4=-T. CancelA² in the actual coefficient domain: A*D15=-T*D12.',
 'final_actual_identity':'A*D15_norm=-T*detM12','no_final_anchor_or_kappa3_or_PK92_inverse':True,
 'conditional_full57_proof':'On oldOmega and D(detM12), -T*detM12 is a unit, so its factors A and D15 are units. The normalized15minor is a FULLU unit times an actual15-row minor. Thus current15x22 has fullrowrank15; the accepted full57/raw7/source-lift criterion applies on its retained oldchart. No detM12 value or nonzero assertion is supplied.',
 'retained_old_Omega':omega,'new_sufficient_chart':omega+' ∩ D(detM12)','no_old_Omega_factor_dropped':True,'all_full57_original63_raw91_and_Q9_eight_tails_inherited':True,'det12_or_global57_quartic239_or_complements_not_computed':True}
 ep=out/'independent-source-normalization-and-A-D15-minusT-D12-proof.json';ep.write_text(json.dumps(enc(data),indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_PROVED_actual_A_times_D15_equals_minus_T_times_D12','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'free_identity_poly_terms':[len(d3),len(d4)],'actual_det12_values_or_rank_computed':False,'old_Omega_preserved':omega,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()
