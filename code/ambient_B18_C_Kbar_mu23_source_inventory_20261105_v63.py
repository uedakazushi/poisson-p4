#!/usr/bin/env python3
"""SOURCE only: fourteen z-action obstructions as common37 unknown M2 dots."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
import ambient_B18_C_Kbar_r2_initial_source_20261105_v63 as I
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_MU23_SOURCE_INVENTORY_20261105_v63'
PATHS={
 'old40':I.PATHS['old40'],'old40_peer':I.PATHS['old40_peer'],
 'linear':I.PATHS['linear'],'exact':I.PATHS['exact'],'exact_peer':I.PATHS['exact_peer'],
 'seven':I.PATHS['seven'],'seven_peer':'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63_REVIEW/report.json',
 'three':'reports/AMBIENT_B18_C_KBAR_THREE_GENERATOR_RANK_SOURCE_20261105_v63/formal-deep-action-unit-three-generator-minor-and-exact-rank-four.json',
 'three_peer':'reports/AMBIENT_B18_K_THREE_ROOT_QUOTIENT_REVIEW_20261105_v63/report.json',
 'linear_rank':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
 'linear_actual':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'linear_source':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'operator':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'operator_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pscale(p,c):return {k:I.mulq(v,c)for k,v in p.items()if v and c}
def pplus(a,b,c=Q(1)):return I.polyadd(a,b,c)
def colvec(v):return{(i//8,I.exponent(i%8)):I.parse(c)for i,c in v}
def dot(a,b):
 out=Q(0)
 for x,y in zip(a,b):out=I.addq(out,I.mulq(x,y))
 return out
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'old40':['homogeneous_module_generators40','standard_degree1','standard_degree2'],
 'linear':['linear_basis37'],'exact':['exact_length','exact_Loewy_HF'],
 'seven':['actual_correction_coefficients'],'three':['known_initial_scalars'],
 'linear_rank':['complete_linear_rank_certificate'],
 'linear_actual':['actual40_by657_linear_relation_columns'],
 'linear_source':['actual_selected37_original768_constant_combinations'],
 'operator':['actual_minimal_generator_presentation_shape','actual_relation_matrix_node','actual_relation_order']}
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 prior=S.count(d)+35000;assert prior<150000
 assert d['exact']['exact_length']==11 and d['exact']['exact_Loewy_HF']==[5,3,2,1]
 assert d['operator']['actual_minimal_generator_presentation_shape']==[5,825]
 cert=d['linear_rank']['complete_linear_rank_certificate'];chosen=cert['selected_original_columns'];assert len(chosen)==37
 norm=[colvec(v)for v in cert['basis_normalized_columns']]
 gs=[{(t,tuple(e)):I.parse(c)for t,e,c in row}for row in d['old40']['homogeneous_module_generators40']]
 assert norm==gs[:37]
 # Derive true constant lifts of all normalized first37 from the existing
 # elimination receipts, never substituting homogeneous rows for M2.
 T=[];original=[colvec(d['linear_actual']['actual40_by657_linear_relation_columns'][j])for j in chosen]
 record={r['column']:r for r in cert['fixed_elimination_receipts']if 'divide'in r}
 for pos,j in enumerate(chosen):
  rec=record[j];assert colvec(rec['normalized'])==norm[pos]
  row=[Q(int(k==pos))for k in range(37)]
  for b,c in rec['subtract']:
   assert b<pos;c=I.parse(c)
   row=[I.addq(x,-I.mulq(c,y))for x,y in zip(row,T[b])]
  pivot=I.parse(rec['divide']);row=[I.divq(x,pivot)for x in row];T.append(row)
  out={}
  for c,o in zip(row,original):out=pplus(out,o,c)
  assert out==norm[pos]
 # Complete gr2 covector; only SOURCE initial reductions, no M2 queried.
 q1=(4,I.exponent(5,7));q2=(4,I.exponent(7,7));v=(4,I.exponent(5));w=(3,I.exponent(7));r=I.exponent(7)
 pi={};nfhash=hashlib.sha256()
 m2=[tuple(e)for e in __import__('itertools').combinations_with_replacement(range(8),2)]
 for t in range(5):
  for k,inds in enumerate(m2):
   inp={(t,I.exponent(*inds)):Q(1)};nf,tr=I.reduce(inp,gs);assert set(nf)<={q1,q2}
   if nf.get(q2,Q(0)):pi[36*t+k]=nf[q2]
   nfhash.update((json.dumps([I.enc(nf),tr],separators=(',',':'))+'\n').encode())
 coeff=d['seven']['actual_correction_coefficients'];a3,c3,as_,cs=[I.parse(coeff[k])for k in['a3=b3/b1','c3=(p3-a3*p1)/d','aS=(b2-b3/4)/b1','cS=(p2-p3/4-aS*p1)/d']]
 U0={(3,I.exponent()):Q(1),(1,I.exponent()):-a3};S0={(2,I.exponent()):Q(1),(3,I.exponent()):Q(-1,4),(1,I.exponent()):-as_}
 rU={};rS={}
 for out,g in[(rU,U0),(rS,S0)]:
  for(t,e),c in g.items():out[(t,r)]=c
 nU,_=I.reduce(rU,norm);nS,_=I.reduce(rS,norm);assert nU.get(w)==1 and set(nU)<={w,v}and set(nS)<={v}
 A=nU.get(v,Q(0));lam=nS[v];assert lam
 forms=[];proofs=[]
 for X,X0,cx in[('U',U0,c3),('S',S0,cs)]:
  for zi in range(7):
   z=I.exponent(zi);linear={(t,z):c for(t,e),c in X0.items()};nf,tr=I.reduce(linear,norm);assert set(nf)<={v,w}
   av=nf.get(v,Q(0));aw=nf.get(w,Q(0));k=I.divq(I.addq(av,-I.mulq(aw,A)),lam)
   F1=pplus(pplus(linear,rU,-aw),rS,-k);zero,tr=I.reduce(F1,norm);assert zero=={}
   basisco=[Q(0)]*37
   for j,shift,s in tr:assert shift==[0]*8;basisco[j]=I.addq(basisco[j],I.parse(s))
   actualco=[dot(basisco,[T[b][j]for b in range(37)])for j in range(37)]
   rebuilt={}
   for q,g in zip(actualco,original):rebuilt=pplus(rebuilt,g,q)
   assert rebuilt==F1
   F2={(4,I.exponent(zi,7)):-cx};F2=pplus(F2,{q2:I.addq(I.mulq(aw,c3),I.mulq(k,cs))})
   nf2,tr2=I.reduce(F2,gs);constant=nf2.get(q2,Q(0))
   forms.append([constant]+[-x for x in actualco]);proofs.append({'source':[zi,X],'first_NF_v_w':[str(av),str(aw)],'subtract_rU_rS':[str(aw),str(k)],'degree1_original_source':I.enc(F1),'degree2_explicit_source':I.enc(F2),'actual_selected37_M2_weights':[str(x)for x in actualco],'degree2_known_NF':I.enc(nf2),'linear37_normalized_membership':tr,'quadratic_initial_membership':tr2})
 # Rank of the FUNCTIONAL FAMILY only, not the values of the actions.
 bb={};selected=[]
 for j,row in enumerate(forms):
  row=row.copy()
  for p,b in bb.items():
   a=row[p]
   if a:row=[I.addq(x,-I.mulq(a,y))for x,y in zip(row,b)]
  non=[i for i,x in enumerate(row)if x]
  if non:
   p=non[0];pivot=row[p];bb[p]=[I.divq(x,pivot)for x in row];selected.append(j)
 evidence={
 'scope':'SOURCE only: fourteen unknown z-actions are affinefunctionals of SAME37 actualOLDM2scalars. NoM2/actionvalueor generatorchoicecomputed.',
 'input_sha256':pins,'old_operator':[5,825,d['operator']['actual_relation_matrix_node']],
 'same_Rbar_frame':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
 'actual_selected37_sources':[[j,168+j,combo]for j,combo in zip(chosen,d['linear_source']['actual_selected37_original768_constant_combinations'])],
 'normalized37_actual37_constant_lift':[[str(x)for x in row]for row in T],
 'all1480_lift_multiplybacks':True,'actual_pi_q2_degree2_sparse':[[i,str(x)]for i,x in sorted(pi.items())],
 'all180_initial_reduction_membership_stream_sha256':nfhash.hexdigest(),
 'known_rU_rS_firsts':{'A':str(A),'lambdaS':str(lam)},
 'all14_z_U_S_full_source_and_coherent_M2_functionals':proofs,
 'all14_affine_rows_constant_then37_unknown_t':[[str(x)for x in row]for row in forms],
 'unpaired_affine_family_rank':len(bb),'fixed_independent_functional_rows':selected,
 'actual_unknowns':'t_j=pi_q2(M2ofTRUEoldselected37source_j), exactOLD5x825root29/168+original657index, fullSchur/Eulertransport. NOTnativehomogeneousquadratic.',
 'mu2_iff_test':'Kbar/rKbar hasQQbasisU,S,q2. AllziU/ziS maponlyq2, allziq2areoldm3=QQg⊂rKbar. Hence mu_R0K=2 iff at leastone14functional_nonzero; mu3 iffallzero. FiniteellNakayama, notDgeneratorcount.',
 'no_M3_required':'Onlygr2q2coeffofdegree1relationsaftercoherentM2clearing; oldm3tailsdrop inKbar/rKbar sinceg=rq2. NoH/G38cubicremainder/DeltaM3query.',
 'next_cache_frontier':['Readactualpi_q2supportbeforeassumingliteralhead4r²','BindactualOLDPi0/Pi1raw56+wedge84domain','Forresidualoriginal768constantssourceX0zero, keeprawP_g2 plusoldfirstSchurcarry; noPi2X0 shortcutwithoutsourcebinding','Inventorytrueg-cofactor_r²leaves againstuniv5eightvalues/copiedCBsource; keepallopaqueunknownsecondsformal'],
 'original_H0_nil_or_new_response_claimed':False,'maximum_SOURCE_pre_bits':I.MAXPRE,'selective_receipts':receipts}
 owned=S.count(d)+S.count(T)+S.count(forms)+S.count(evidence)+25000;assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'fourteen-z-actions-common37-M2-source-functionals-and-minimal-generator-criterion.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_MU2_OR3_FOURTEEN_AFFINE_ACTIONS','unpaired_affine_family_rank':len(bb),'actual_pi_q2_support':[[i,str(x)]for i,x in sorted(pi.items())],'whole_SOURCE_owned_upper':owned,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'input_sha256':pins,'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
