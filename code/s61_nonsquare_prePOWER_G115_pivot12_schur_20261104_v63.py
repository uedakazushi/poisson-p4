#!/usr/bin/env python3
"""Source-only twelfth pivot from the accepted literal-pair bordered6."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_G115_PIVOT12_SCHUR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 st='S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63';so=load(st,'operator.json');assert load(st)['all_passed']
 en='S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_ENTRY_20261104_v63';er=load(en);assert er['all_passed']and er['normalized_D6_terms']==9047 and not er['normalized_D6_is_zero'];assert pins['reports/'+en+'/report.json']=='2a727a0a37b0d89ef04b18f95ae3aebf2c82146f14c073789f27008690680e42'
 lift='S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63';lo=load(lift,'operator.json');assert load(lift+'_REVIEW')['all_passed'];hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json')
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
 R=ref(st,'operator.json','/exports/residual4x11',[4,11]);Ls=ref(st,'operator.json','/exports/fullsource15',[15,15]);Lsi=ref(st,'operator.json','/exports/fullsource15_inverse',[15,15]);Us=ref(st,'operator.json','/exports/fulltarget22',[22,22]);Usi=ref(st,'operator.json','/exports/fulltarget22_inverse',[22,22])
 M=ref('S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63','operator.json','/exports/input15x22',[15,22]);labels=[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31];literal5=sel(M,[0,1,3,6,2],[labels.index(j)for j in[22,11,13,25,15]]);D5=n('universal_polynomial_determinant',[1,1],[literal5],expanded=False,order=5)
 g=ref('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json','/raw_compact14x5/6/1',[1,1]);c10=n('QQ_constant',[1,1],value=hc['known_row_units_rho14'][10]['constant']);a=n('literal_parameter_monomial',[1,1],variable_order=['a','e','d'],polynomial=[[[1,0,0],'1']]);a2=prod(a,a)
 d6=ref(en,'native.json','/exact_normalized_D6',[1,1]);iD6=n('certified_chart_inverse',[1,1],[d6],chart='D(N6)',nonzero_primary=en,peer_pending=True);ic10a2=n('known_FULLU_parameter_unit_inverse',[1,1],[prod(c10,a2)]);ir=prod(prod(prod(g,D5),ic10a2),iD6)
 rp=[1,0,2,3];cp=[10]+list(range(10));Pr=n('row_permutation',[4,4],order=rp);Pri=n('transpose',[4,4],[Pr]);Pc=n('column_permutation',[11,11],order=cp);Pci=n('transpose',[11,11],[Pc]);RR=prod(prod(Pr,R),Pc);r=sel(RR,[0],[0]);B=sel(RR,[0],list(range(1,11)));C=sel(RR,list(range(1,4)),[0]);D=sel(RR,list(range(1,4)),list(range(1,11)));ci=scale(ir,C);ib=scale(ir,B);S=sub(D,prod(ci,B))
 LL=block([[ir,zero(1,3)],[neg(ci),eye(3)]]);LLi=block([[r,zero(1,3)],[C,eye(3)]]);UU=block([[eye(1),neg(ib)],[zero(10,1),eye(10)]]);UUi=block([[eye(1),ib],[zero(10,1),eye(10)]])
 Fsrc=prod(diag(eye(11),prod(LL,Pr)),Ls);Fsrci=prod(Lsi,diag(eye(11),prod(Pri,LLi)));Ftgt=prod(Us,diag(eye(11),prod(Pc,UU)));Ftgti=prod(diag(eye(11),prod(UUi,Pci)),Usi)
 # Explicit accepted full57 lift and Cdelta cleanup for all three survivors.
 Q25=root(lift,lo['exports']['source15_to_original25'],[15,25]);Q57=root(lift,lo['exports']['lift25_to_full57'],[25,57]);Q63=root(lift,lo['exports']['full57_to_original63'],[57,63]);H7=root(lift,109,[57,7]);Xt=root(lift,111,[7,7]);cdelta=root(lift,114,[1,57]);bz=root(lift,115,[1,1]);bw=root(lift,116,[1,1]);ieta=root(lift,118,[1,1]);F57=root(lift,108,[57,91]);surv15=sel(Fsrc,[12,13,14],list(range(15)));surv25=prod(surv15,Q25);unc=prod(surv25,Q57);H=prod(prod(unc,H7),Xt);hz=sel(H,list(range(3)),[5]);hw=sel(H,list(range(3)),[6]);corr=prod(scale(ieta,add(prod(hz,bz),prod(hw,bw))),cdelta);corrected=sub(unc,corr);orig63=prod(corrected,Q63);full91=prod(corrected,F57);physicalHcols=lo['nodes'][109]['columns'];pure84=sel(full91,list(range(3)),[j for j in range(91)if j not in physicalHcols]);raw7=prod(prod(corrected,H7),Xt)
 proof=dict(nonempty_chart_bound_to_frozen_entry=True,scalar_product_identity='g*D5_G20*r=c10*a^2*N6 in FULLU[Delta407^-1,Saux^-1,T^-1]',
  why_identity_global_on_old_chart='Ratio proved on dense extraD(gD5); both sides are elements of the localized domain, so equality holds in the whole old localized ring.',
  no_extra_g_D5_chart='On D(N6), product g*D5*r is a unit; each factor is then a unit. Formula r^-1=g*D5/(c10*a^2*N6) uses only existingFULLU units plusN6 inverse.',
  no_D5_expansion=True,forward_inverse_blockproof='LL*[r B;C D]*UU=diag(1,D-C*r^-1*B), LL*LLinverse=I4, UU*UUinverse=I11; permutations and prior inverses compose.',
  all_remaining3_source_raw7_zero='Existing HLprefix/5projectedheavy zero plus exact etaCdelta cleanup; the accepted all7 cleanup identity applies to each new source covector.',
  source_original25_and_canonical63_explicit=True,canonical_source135_normalization_scope_unchanged=True,remaining_kernel_rank3_sufficient_for_full_cubic57=True,
  no_remaining3rank_or_global57_claim=True,all7_remaining_extra_tail_effects_preserved=True)
 dump(out/'proof.json',proof)
 operator=dict(nodes=N,exports=dict(residual3x10=S,pivot_scalar=r,pivot_scalar_inverse=ir,D5_G20_unexpanded=D5,fullsource15=Fsrc,fullsource15_inverse=Fsrci,fulltarget22=Ftgt,fulltarget22_inverse=Ftgti,surviving3_source_original25=surv25,surviving3_corrected_full57=corrected,surviving3_corrected_original63=orig63,surviving3_actual_pure84=pure84,surviving3_raw7=raw7),
  source_pivot_order_current15=so['source_pivot_order']+[11],source_remaining_current15=[9,12,14],target_pivot_original32=so['target_pivot_original32']+[19],target_remaining_original32=[1,3,4,6,7,8,9,14,17,18],chart='FULLU D(Delta407*Saux*T*N6)',no_extra_g_D5_inverse=True,all7_remaining_extra_tail_effects_retained=True,new_parameter_coefficients=0,all_actual_source_coefficients_external=True)
 dump(out/'operator.json',operator);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,independent_peer_pending=True,typed_nodes=len(N),samecokernel_shape=[3,10],pivot_count=12,no_new_scalar_coefficients=True,no_D5_expansion=True,no_extra_g_D5_chart=True,actual_full57_lift_and_raw7_cleanup_explicit=True,no_global57_or_residual_rank_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
