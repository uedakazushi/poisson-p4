#!/usr/bin/env python3
"""Source-only thirteenth pivot with the actualN7 monomial/content normalization."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K120_PIVOT13_SCHUR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 prev='S61_NONSQUARE_PREPOWER_G115_PIVOT12_SCHUR_20261104_v63';po=load(prev,'operator.json');assert load(prev)['all_passed'];assert load(prev+'_REVIEW')['all_passed']
 en='S61_NONSQUARE_PREPOWER_K120_T14_NESTED_ENTRY_20261104_v63';er=load(en);nv=load(en,'native.json');assert er['all_passed']and not er['normalized_N7_is_zero']and er['normalized_N7_terms']==14986
 globalD=int(nv['normalized_N7']['global_denominator']);pol=nv['normalized_N7']['integer_coefficients'];assert len(pol)==14986 and globalD>0
 content=__import__('functools').reduce(__import__('math').gcd,(abs(int(c))for m,c in pol));amin=min(m[0]for m,c in pol);assert content>0 and amin==9
 bar=[[[m[0]-amin,m[1],m[2]],str(int(c)//content)]for m,c in pol];assert all(int(c)%content==0 for m,c in pol);assert all(int(c)*content==int(oldc) and m[0]+amin==oldm[0]and m[1:]==oldm[1:]for(m,c),(oldm,oldc)in zip(bar,pol));assert min(m[0]for m,c in bar)==0
 primitive=dict(polynomial=bar,coefficient_field='ZZ',content_positive=str(content),original_global_denominator=str(globalD),removed_known_a_power=amin,rational_scalar=str(__import__('fractions').Fraction(content,globalD)),exact_identity='N7=scalar*a^9*N7bar',known_FULLU_units_only=True)
 dump(out/'primitive-N7bar.json',primitive)
 N=[]
 def n(op,sh,args=(),**kw):
  assert all(0<=i<len(N)for i in args);i=len(N);N.append(dict(id=i,op=op,shape=sh,args=list(args),**kw));return i
 def ref(st,name,ptr,sh):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return n('external_frozen_reference',sh,reference=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),json_pointer=ptr))
 def root(st,node,sh):
  p=ROOT/'reports'/st/'operator.json';pins[str(p.relative_to(ROOT))]=sha(p);x=read(p);assert x['nodes'][node]['shape']==sh;return n('external_frozen_operator_root',sh,reference=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),node_index=node))
 def eye(m):return n('identity',[m,m])
 def zero(m,k):return n('zero',[m,k])
 def prod(a,b):assert N[a]['shape'][1]==N[b]['shape'][0];return n('matrix_product',[N[a]['shape'][0],N[b]['shape'][1]],[a,b])
 def sub(a,b):assert N[a]['shape']==N[b]['shape'];return n('difference',N[a]['shape'],[a,b])
 def add(a,b):assert N[a]['shape']==N[b]['shape'];return n('sum',N[a]['shape'],[a,b])
 def neg(a):return n('negative',N[a]['shape'],[a])
 def scale(s,a):assert N[s]['shape']==[1,1];return n('central_scalar_multiple',N[a]['shape'],[s,a])
 def sel(a,r,c):return n('submatrix',[len(r),len(c)],[a],rows=r,columns=c)
 def block(rr):
  h=[N[r[0]]['shape'][0]for r in rr];w=[N[x]['shape'][1]for x in rr[0]]
  for i,r in enumerate(rr):
   for j,x in enumerate(r):assert N[x]['shape']==[h[i],w[j]]
  return n('block_matrix',[sum(h),sum(w)],[x for r in rr for x in r],blocks=rr)
 def diag(a,b):return block([[a,zero(N[a]['shape'][0],N[b]['shape'][1])],[zero(N[b]['shape'][0],N[a]['shape'][1]),b]])
 R=ref(prev,'operator.json','/exports/residual3x10',[3,10]);Ls=ref(prev,'operator.json','/exports/fullsource15',[15,15]);Lsi=ref(prev,'operator.json','/exports/fullsource15_inverse',[15,15]);Us=ref(prev,'operator.json','/exports/fulltarget22',[22,22]);Usi=ref(prev,'operator.json','/exports/fulltarget22_inverse',[22,22])
 n6=ref('S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_ENTRY_20261104_v63','native.json','/exact_normalized_D6',[1,1]);pbar=out/'primitive-N7bar.json';barroot=n('external_integral_polynomial',[1,1],reference=dict(path=str(pbar.relative_to(ROOT)),sha256=sha(pbar),json_pointer='/polynomial'));ibar=n('certified_chart_inverse',[1,1],[barroot],chart='D(N7bar)',actual_source_nonzero_primary=en,independent_nonzero_gate_pending=True)
 hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json');c11=__import__('fractions').Fraction(hc['known_row_units_rho14'][11]['constant']);qscalar=__import__('fractions').Fraction(content,globalD);unit=n('literal_known_FULLU_parameter_unit',[1,1],QQ_scalar=str(c11*qscalar),polynomial_monomial=[9,0,0],variable_order=['a','e','d']);iunit=n('known_FULLU_parameter_unit_inverse',[1,1],[unit]);ir=prod(prod(n6,iunit),ibar)
 rp=[1,0,2];cp=[7,0,1,2,3,4,5,6,8,9];Pr=n('row_permutation',[3,3],order=rp);Pri=n('transpose',[3,3],[Pr]);Pc=n('column_permutation',[10,10],order=cp);Pci=n('transpose',[10,10],[Pc]);RR=prod(prod(Pr,R),Pc);s=sel(RR,[0],[0]);B=sel(RR,[0],list(range(1,10)));C=sel(RR,[1,2],[0]);D=sel(RR,[1,2],list(range(1,10)));ci=scale(ir,C);ib=scale(ir,B);S=sub(D,prod(ci,B))
 LL=block([[ir,zero(1,2)],[neg(ci),eye(2)]]);LLi=block([[s,zero(1,2)],[C,eye(2)]]);UU=block([[eye(1),neg(ib)],[zero(9,1),eye(9)]]);UUi=block([[eye(1),ib],[zero(9,1),eye(9)]]);W=sel(prod(LL,Pr),[1,2],list(range(3)))
 Fsrc=prod(diag(eye(12),prod(LL,Pr)),Ls);Fsrci=prod(Lsi,diag(eye(12),prod(Pri,LLi)));Ftgt=prod(Us,diag(eye(12),prod(Pc,UU)));Ftgti=prod(diag(eye(12),prod(UUi,Pci)),Usi)
 s25=ref(prev,'operator.json','/exports/surviving3_source_original25',[3,25]);s57=ref(prev,'operator.json','/exports/surviving3_corrected_full57',[3,57]);s63=ref(prev,'operator.json','/exports/surviving3_corrected_original63',[3,63]);p84=ref(prev,'operator.json','/exports/surviving3_actual_pure84',[3,84]);r7=ref(prev,'operator.json','/exports/surviving3_raw7',[3,7]);new25=prod(W,s25);new57=prod(W,s57);new63=prod(W,s63);new84=prod(W,p84);new7=prod(W,r7)
 proof=dict(actual_residual_identity='s*N6=c11*N7=c11*qscalar*a^9*N7bar',
  known_unit_logic='c11andqscalararenonzeroQQunits,aunitinFULLU,N6alreadyunitonoldchart;onD(N7bar),sinverse=N6/(c11*qscalar*a^9*N7bar).No newg/D5/k/p/R/h22unitclaim.',
  exact_content_and_monomial_identity_checked_all14986_integer_coefficients=True,no_polynomial_products_or_factorization=True,
  row_column_order=dict(old_remaining_current15=[9,12,14],new_pivot_current15=12,remaining_current15=[9,14],old_targets=[1,3,4,6,7,8,9,14,17,18],new_pivot_original32=14,new_targets=[1,3,4,6,7,8,9,17,18]),
  generic_block_identity='LL*[s B;C D]*UU=diag(1,D-C*sinverse*B),bothinverseoperatorsliteralblockidentities',
  actual_source_faithfulness='W timespreviouscorrectedactual57vectors;raw7zero preservedbylinearity,canonical63andpure84followthesameW. AllQ9taileffectsremaininoldfulloperators.',
  remaining2rank2_suffices_for_full_cubic57=True,no_remainingrank2_or_global57_claim=True,
  peer_scope='PrimaryN7arithmeticexact;forthcomingindependentmodsymbolicpeerwouldproveactualsourceNONZERO,char0,notfull14986QQcoefficientidentity.')
 dump(out/'proof.json',proof)
 operator=dict(nodes=N,exports=dict(residual2x9=S,pivot_scalar=s,pivot_scalar_inverse=ir,N7bar=barroot,fullsource15=Fsrc,fullsource15_inverse=Fsrci,fulltarget22=Ftgt,fulltarget22_inverse=Ftgti,surviving2_source_original25=new25,surviving2_corrected_full57=new57,surviving2_corrected_original63=new63,surviving2_actual_pure84=new84,surviving2_raw7=new7),source_pivot_order_current15=po['source_pivot_order_current15']+[12],source_remaining_current15=[9,14],target_pivot_original32=po['target_pivot_original32']+[14],target_remaining_original32=[1,3,4,6,7,8,9,17,18],
  chart='FULLU D(Delta407*Saux*T*N6*N7bar)',all7_remaining_extra_tail_effects_retained=True,new_expanded_parameter_coefficients=0,primitive_integer_content_normalization_only=True)
 dump(out/'operator.json',operator);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,independent_peer_pending=True,typed_nodes=len(N),samecokernel_shape=[2,9],actual_pivot_count=13,N7bar_terms=len(bar),removed_known_a_power=9,positive_integer_content_bits=content.bit_length(),nonzero_QQ_scalar=str(qscalar),no_new_scalar_coefficient_evaluation=True,no_new_unknown_localizer_beyond_N7bar=True,actual_full57_source_lift_by_linearity=True,no_global57_or_remainingrank_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
