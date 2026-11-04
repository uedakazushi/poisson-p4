#!/usr/bin/env python3
"""Typed actual11pivot from kappa3*D8=cprod*Saux*T, no D8 coefficient expansion."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 dr='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63';old=load(dr,'operator.json');assert load(dr+'_REVIEW','report.json')['all_passed']
 es='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63';er=load(es,'report.json');en=load(es,'native.json');eb=load(es,'source-binding.json');assert er['all_passed']and not er['output_is_zero']and er['output_terms']==7208
 assert eb['common_B_clearing_power']==1 and eb['source_current15_rows']==[2,4,5,8]and eb['cofactor_roots']==[179,198,217,237]
 assert load('S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW','report.json')['all_passed']
 assert load('S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW','report.json')['all_passed']
 hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');rank=load('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63')
 for i in[1,3,4,7]:assert hc['known_row_units_rho14'][i]['a_d_F_B_C_J_powers']==[1,0,1,0,2,0]
 assert all(not hc['raw_compact14x5'][i][j]for i in[1,3,4,7]for j in[1,2,3])
 cprod=Q(hc['known_row_units_rho14'][1]['constant'])*Q(hc['known_row_units_rho14'][3]['constant'])*Q(hc['known_row_units_rho14'][4]['constant'])
 nodes=[]
 def n(op,shape,args=(),**kw):
  assert all(0<=i<len(nodes)for i in args);i=len(nodes);nodes.append(dict(id=i,op=op,shape=shape,args=list(args),**kw));return i
 def ref(st,name,pointer,shape):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return n('external_frozen_input',shape,reference=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),json_pointer=pointer))
 def oldref(exp,shape):return ref(dr,'operator.json',f"/exports/{exp}",shape)
 def prod(a,b):assert nodes[a]['shape'][1]==nodes[b]['shape'][0];return n('matrix_product',[nodes[a]['shape'][0],nodes[b]['shape'][1]],[a,b])
 def scale(a,b):return n('central_scalar_multiple',nodes[b]['shape'],[a,b])
 def neg(a):return n('negative',nodes[a]['shape'],[a])
 def select(a,rr,cc):return n('submatrix',[len(rr),len(cc)],[a],rows=rr,columns=cc)
 def eye(k):return n('identity',[k,k])
 def zero(m,k):return n('zero',[m,k])
 def block(rows):
  h=[nodes[r[0]]['shape'][0]for r in rows];w=[nodes[x]['shape'][1]for x in rows[0]]
  for i,r in enumerate(rows):
   for j,x in enumerate(r):assert nodes[x]['shape']==[h[i],w[j]]
  return n('block_matrix',[sum(h),sum(w)],[x for r in rows for x in r],blocks=rows)
 def inv(a,proof):return n('certified_parameter_unit_inverse',[1,1],[a],proof=proof)
 M=oldref('input15x22',[15,22]);one=n('one',[1,1]);L=n('named_FULLU_unit',[1,1],expression='a*F*C^2');B=n('named_FULLU_unit',[1,1],expression='B');gamma=ref('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json','/known_column_units/4',[1,1]);Lg=prod(L,gamma);iL=inv(L,'knownFULLUunit');iLg=inv(Lg,'knownFULLUunit');iB=inv(B,'knownFULLUunit')
 cols=[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31]
 diag=[B if j==1 else iL if j==17 else iLg if j==21 else one for j in range(22)];idiag=[iB if j==1 else L if j==17 else Lg if j==21 else one for j in range(22)]
 Dcol=n('diagonal',[22,22],diag,diagonal=diag);Dcoli=n('diagonal',[22,22],idiag,diagonal=idiag)
 rows8=[0,1,3,6,2,4,5,8];Hrows=[7,13,10];prow=rows8+Hrows;rest=[i for i in range(15)if i not in prow];targets=[22,11,13,25,15,2,27,31,28,30,29];pc=[cols.index(j)for j in targets];rc=[j for j in range(22)if j not in pc]
 assert rest==[9,11,12,14]and[cols[j]for j in rc]==[1,3,4,6,7,8,9,14,17,18,19]
 Pr=n('row_permutation',[15,15],order=prow+rest);Pri=n('transpose',[15,15],[Pr]);Pc=n('column_permutation',[22,22],order=pc+rc);Pci=n('transpose',[22,22],[Pc]);MM=prod(prod(Pr,prod(M,Dcol)),Pc)
 A8=select(MM,list(range(8)),list(range(8)));B3=select(MM,list(range(8,11)),list(range(8,11)));C38=select(MM,list(range(8,11)),list(range(8)))
 Sax=ref('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63','native.json','/exact_Saux_Laurent_polynomial',[1,1]);T=ref(es,'root239.json','/polynomial',[1,1]);kap=ref(es,'root237.json','/polynomial',[1,1]);c=n('QQ_constant',[1,1],value=str(cprod));rhs=prod(prod(c,Sax),T);irhs=inv(rhs,'NewchartD(Saux*T); cprod nonzeroQQ. No additionalkappa3chart.');D8=n('universal_determinant',[1,1],[A8],expanded=False);D8inv=prod(kap,irhs)
 A8i=scale(D8inv,n('universal_polynomial_adjugate',[8,8],[A8],expanded=False))
 N3=oldref('H5_B3',[3,3]);Delta=ref('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63','certificate.json','/nonzero_selected_minor_Delta',[1,1]);iDelta=inv(Delta,'NewchartD(Delta407). No p orDelta5inverse.')
 rhos=[ref('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json',f'/known_row_units_rho14/{i}',[1,1])for i in[6,12,9]];irhos=[inv(x,'knownFULLUunit')for x in rhos];RD=n('diagonal',[3,3],rhos,diagonal=rhos);RDi=n('diagonal',[3,3],irhos,diagonal=irhos)
 B3i=prod(scale(iDelta,n('universal_polynomial_adjugate',[3,3],[N3],expanded=False)),RDi)
 K=select(MM,list(range(11)),list(range(11)));Ki=block([[A8i,zero(8,3)],[neg(prod(prod(B3i,C38),A8i)),B3i]])
 Bb=select(MM,list(range(11)),list(range(11,22)));Cc=select(MM,list(range(11,15)),list(range(11)));Dd=select(MM,list(range(11,15)),list(range(11,22)));KB=prod(Ki,Bb);CK=prod(Cc,Ki);res=n('difference',[4,11],[Dd,prod(Cc,KB)])
 LS=block([[Ki,zero(11,4)],[neg(CK),eye(4)]]);LSi=block([[K,zero(11,4)],[Cc,eye(4)]]);UT=block([[eye(11),neg(KB)],[zero(11,11),eye(11)]]);UTi=block([[eye(11),KB],[zero(11,11),eye(11)]])
 fs=prod(LS,Pr);fsi=prod(Pri,LSi);ft=prod(prod(Dcol,Pc),UT);fti=prod(prod(UTi,Pci),Dcoli);sv=select(fs,list(range(11,15)),list(range(15)));old25=prod(sv,oldref('source_current15_to_original25',[15,25]))
 proof=dict(cprod=str(cprod),kappa_last_row_coefficient_root237=True,column2_Bclear_is_known_FULLUunit=True,
  polynomial_identity='kappa3*D8tilde=cprod*Saux*T; T=B*kappa*Mcol2',
  row_replacement_determinant='Replace lastbottomrow by kappa0*row4+kappa1*row5+kappa2*row6+kappa3*row7; determinantmultiplier kappa3.',
  last_row_expansion_sign=1,last_row_index=7,last_nonzero_column_index=5,
  D7_source_identity='D7=c1*c3*c4*Saux, since rho_i=c_i*L and Bhat_u4=delta2*G_u4/c4-mu0*G_pair1/c1-mu1*G_pair3/c3',
  D7_row_swap_parity=8,D7_column_swap_parity=10,before_p_inversion='Derive atdelta2nonzero then clearknownFULLUdenominators and descend integralpolynomial identity; no p inverse in final operator.',
  inverse_D8_formula='kappa3/(cprod*Saux*T)',kappa3_additional_chart_NOT_needed=True,
  upper8x3_heavy_literal_zeros=True,B3actual_equals='diag(rho6,rho12,rho9)*acceptedN3',
  whole_chart='FULLU D(Delta407*Saux*T)',all3_chart_factors_nonzero=True,
  no_D8_or_11minor_coefficients_expanded=True,projected15_only_fullrawH7_eta_Cdelta_HZ_HW_lift_SEPARATE=True)
 op=dict(nodes=nodes,exports=dict(input15x22=M,normalized_column_units=Dcol,actual_A8=A8,unexpanded_D8=D8,certified_D8_inverse=D8inv,actual_B3=B3,pivot11=K,pivot11_inverse=Ki,residual4x11=res,full_source15=fs,full_source15_inverse=fsi,full_target22=ft,full_target22_inverse=fti,surviving4_source_original25=old25),
  source_pivot11=prow,source_remaining4=rest,target_pivot11_original32=targets,target_remaining11_original32=[cols[j]for j in rc],
  all_original8_Q9_tail_effects_retained=True,tail2_used_as_unit_pivot=True,remaining7_extra_tails=[1,3,4,6,7,8,9],
  identities=['kappa3*detA8=cprod*Saux*T','detA8*certified_D8_inverse=1 withoutkappa3inverse','K11=[[A8,0],[C38,B3]], detK11=detA8*prod(rho)*Delta407','K11*Ki=Ki*K11=I11','fullsource*M15*fulltarget=diag(I11,residual4x11); bothinversepairs exact'],
  chart='FULLU D(Delta407*Saux*T)',raw_H7_lift_NOT_computed=True,new_scalar_coefficients=0,rank4_or_global57_NOT_proved=True)
 dump(out/'operator.json',op);dump(out/'proof.json',proof);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,independent_review_pending=True,typed_nodes=len(nodes),actual_samecokernel_shape=[4,11],source_projected15_only=True,
  localizers=['Delta407','Saux','T7208'],no_extra_kappa3_p_R_h22_Delta5_inverse=True,all8_tail_effects_retained=True,new_parameter_coefficients=0,
  no_D8_minor_or_residual_scalar_expansion=True,no_remaining4rank_or_global57_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
