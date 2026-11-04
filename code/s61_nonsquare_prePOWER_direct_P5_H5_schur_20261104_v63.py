#!/usr/bin/env python3
"""Actual direct P5/H5 typed Schur operator: no R434 inverse or new scalar values."""
import argparse, hashlib, json, math, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_DIRECT_P5_H5_SCHUR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
    def load(st,name):
        p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    dr='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63';old=load(dr,'operator.json');assert load(dr+'_REVIEW','report.json')['all_passed']
    fi='S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63';fifth=load(fi,'report.json');assert fifth['all_passed']and fifth['actual_new_scalar_coefficients']==0
    en='S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63';bind=load(en,'source-binding.json');assert bind['Mark_entry']==[0,9]
    pf='S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63';pre=load(pf,'operator-and-cost.json');assert pre['source_current15_rows']==[5,2,4]
    assert pre['known_rho_LCM_powers']==[1,0,1,0,2,0]
    assert load('S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW','report.json')['all_passed']
    nn=old['nodes']
    def closure(root):
        seen=set();todo=[root]
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo+=nn[i]['args']
        return sorted(seen)
    forbidden=[]
    for key in ['input15x22','hybrid_row_unit','heavy_column_unit','H5_inverse','source_current15_to_original25']:
        for i in closure(old['exports'][key]):
            x=nn[i]
            if x['op']=='certified_chart_unit_inverse'and ('R434'in x.get('unit_proof','')or'h22'in x.get('unit_proof','')):forbidden.append((key,i))
    assert not forbidden
    nodes=[]
    def n(op,shape,args=(),**kw):
        assert all(0<=i<len(nodes)for i in args);i=len(nodes);nodes.append(dict(id=i,op=op,shape=shape,args=list(args),**kw));return i
    def ref(exp,shape):return n('accepted_external_export',shape,reference=dict(path=f'reports/{dr}/operator.json',sha256=pins[f'reports/{dr}/operator.json'],export=exp,node=old['exports'][exp]))
    def prod(a,b):assert nodes[a]['shape'][1]==nodes[b]['shape'][0];return n('matrix_product',[nodes[a]['shape'][0],nodes[b]['shape'][1]],[a,b])
    def sel(a,rr,cc):return n('submatrix',[len(rr),len(cc)],[a],rows=rr,columns=cc)
    def neg(a):return n('negative',nodes[a]['shape'],[a])
    def scale(s,a):return n('central_scalar_multiple',nodes[a]['shape'],[s,a])
    def eye(k):return n('identity',[k,k])
    def zero(m,k):return n('zero',[m,k])
    def block(rows):
        h=[nodes[r[0]]['shape'][0]for r in rows];w=[nodes[x]['shape'][1]for x in rows[0]]
        for i,r in enumerate(rows):
            for j,x in enumerate(r):assert nodes[x]['shape']==[h[i],w[j]]
        return n('block_matrix',[sum(h),sum(w)],[x for r in rows for x in r],blocks=rows)
    M=ref('input15x22',[15,22]);src=[]
    for i in pre['source_coefficients3']:
        src.append(n('accepted_external_scalar_DAG_node',[1,1],reference=dict(path=f'reports/{pf}/operator-and-cost.json',sha256=pins[f'reports/{pf}/operator-and-cost.json'],node=i),semantics='Literal source cofactor in the accepted768 operation DAG.'))
    cu,c1,c3=src;icu=n('certified_chart_scalar_inverse',[1,1],[cu],proof='cu=rhoLCM*delta2/rho_u; rhoLCM=a*F*C^2; delta2 divides Delta5; all rho known FULLU units.',chart='D(Delta5)')
    R=n('single_row_linear_replacement',[15,15],src,replaced_row=5,source_rows=[5,2,4],coefficients=src,all_other_rows='identity')
    nc1=neg(scale(icu,c1));nc3=neg(scale(icu,c3))
    Ri=n('single_row_linear_replacement',[15,15],[icu,nc1,nc3],replaced_row=5,source_rows=[5,2,4],coefficients=[icu,nc1,nc3],all_other_rows='identity')
    # Inverse rows are recorded explicitly and checked by the universal elementary formula.
    assert nodes[Ri]['coefficients']==nodes[Ri]['args']
    S=ref('hybrid_row_unit',[15,15]);Si=n('diagonal_known_unit_inverse',[15,15],[S],proof='Only the five selected H5 source rows divided by accepted known rho units.')
    T=ref('heavy_column_unit',[22,22]);Ti=n('diagonal_known_unit_inverse',[22,22],[T],proof='Only heavy target31 divided by accepted known gamma unit.')
    prows=[0,1,3,6,5,2,4,7,13,10];rrows=[8,9,11,12,14]
    pcols=old['context']['pivot9_target_current22_order'][:4]+[11]+old['context']['pivot9_target_current22_order'][4:]
    rcols=[j for j in range(22)if j not in pcols];assert len(set(pcols))==10
    cols=[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31]
    assert [cols[j]for j in pcols]==[22,11,13,25,15,27,31,28,30,29]
    assert [cols[j]for j in rcols]==[1,2,3,4,6,7,8,9,14,17,18,19]
    Pr=n('row_permutation',[15,15],order=prows+rrows);Pri=n('transpose',[15,15],[Pr]);Pc=n('column_permutation',[22,22],order=pcols+rcols);Pci=n('transpose',[22,22],[Pc])
    MM=prod(prod(prod(Pr,prod(S,prod(R,M))),T),Pc)
    P5=sel(MM,list(range(5)),list(range(5)));H5=sel(MM,list(range(5,10)),list(range(5,10)));C55=sel(MM,list(range(5,10)),list(range(5)))
    sax=n('accepted_external_polynomial',[1,1],reference=dict(path=f'reports/{en}/native.json',sha256=pins.setdefault(f'reports/{en}/native.json',sha(ROOT/f'reports/{en}/native.json')),json_pointer='/exact_Saux_Laurent_polynomial'),terms=7022)
    isax=n('certified_chart_scalar_inverse',[1,1],[sax],chart='D(Saux)',proof='Actual polynomial detP5=Saux by universal bordered determinant identity; existing exact source768 replay, and independent symbolicFp nonzero peer.')
    adj5=n('universal_polynomial_cofactor_adjugate',[5,5],[P5],expanded=False);P5i=scale(isax,adj5);H5i=ref('H5_inverse',[5,5])
    Ki=block([[P5i,zero(5,5)],[neg(prod(prod(H5i,C55),P5i)),H5i]])
    K=sel(MM,list(range(10)),list(range(10)));B=sel(MM,list(range(10)),list(range(10,22)));C=sel(MM,list(range(10,15)),list(range(10)));D=sel(MM,list(range(10,15)),list(range(10,22)))
    KB=prod(Ki,B);CK=prod(C,Ki);res=n('difference',[5,12],[D,prod(C,KB)])
    L=block([[Ki,zero(10,5)],[neg(CK),eye(5)]]);Li=block([[K,zero(10,5)],[C,eye(5)]])
    U=block([[eye(10),neg(KB)],[zero(12,10),eye(12)]]);Ui=block([[eye(10),KB],[zero(12,10),eye(12)]])
    fs=prod(prod(prod(L,Pr),S),R);fsi=prod(prod(prod(Ri,Si),Pri),Li)
    ft=prod(prod(T,Pc),U);fti=prod(prod(Ui,Pci),Ti)
    surviving=sel(fs,list(range(10,15)),list(range(15)));source25=prod(surviving,ref('source_current15_to_original25',[15,25]))
    op=dict(coefficient_ring='FULLU additionally localize ONLY at Delta5*Saux; R434,D4,h22 never inverted.',nodes=nodes,
      exports=dict(row_replacement15=R,row_replacement15_inverse=Ri,P5=P5,H5=H5,pivot10=K,pivot10_inverse=Ki,residual5x12=res,
       full_source15=fs,full_source15_inverse=fsi,full_target22=ft,full_target22_inverse=fti,surviving5_source_original25=source25),
      context=dict(pivot_current15=prows,pivot_original32=[cols[j]for j in pcols],remaining_current15=rrows,
       remaining_original25=[18,19,21,22,24],remaining_original32=[cols[j]for j in rcols],all8_Q9_tails_retained=True,
       hybrid='H5 source rows only scaled by rho; u4 replaced BEFORE that scaling; all other source rows unchanged.',
       determinant='detK10=detP5*Delta5=Saux*Delta5',no_R434_D4_h22_inverse=True,
       complements_retained=['Delta5=0','Saux=0','other H5 charts including p=0'],new_scalar_values=0),
      identities=['R15*R15inverse=R15inverse*R15=I15 from cu unit and unchanged pair source rows',
       'row5(R15*M) has heavy5 all zero and pure17 equal accepted Bhat_u4',
       'detP5=D4*Bhat15-Bhatpiv4*AdjP4*P4col15=Saux as a polynomial identity, no D4 inverse',
       'P5*AdjP5=AdjP5*P5=Saux I5; H5*H5inverse=I5 independently of R434',
       'fullsource*M*fulltarget=diag(I10,residual5x12); both full inverse pairs are identities',
       'on overlap D(R434), this is row/column equivalent to the accepted fifth-pivot5x12',
       'Source full rank implies actualfullcubic57; converse is not asserted'])
    dump(out/'operator.json',op)
    proof=dict(all_shapes_and_accepted_source_orders_checked=True,
      actual_source3cofactor_external_node_ids=pre['source_coefficients3'],
      heavy_row_zero_proof=bind['source_Hpair_identity'],
      determinant_polynomial_proof='Universal bordered determinant identity, valid first at generic invertibleP4 and then in the polynomial/Laurent integral domain; Saux is precisely that expression.',
      H5_inverse_import_closure=closure(old['exports']['H5_inverse']),
      no_R434_inverse_in_H5_or_prefix=True,
      row_replacement_inverse_proof='Solve newrow5=cu*oldrow5+c1*oldrow2+c3*oldrow4; the other14 rows remain unchanged, so oldrow5=(newrow5-c1*newrow2-c3*newrow4)/cu.',
      canonical_rank='On FULLU D(Delta5*Saux), rank(projected15x22)=10+rank(residual5x12). No rank5 computation.',
      chart_extension_requires_no_claim_about_nonempty_VR_intersection=True)
    dump(out/'proof.json',proof)
    # Metadata-only R factor/gcd preflight, no multiplication or factorization.
    rs='S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63';rp=load(rs,'native.json')['remaining_polynomial'];qq=[Q(c)for m,c in rp];den=math.lcm(*(x.denominator for x in qq));ints=[x.numerator*(den//x.denominator)for x in qq];content=math.gcd(*ints)
    assert len(rp)==434 and max(sum(m)for m,c in rp)==20
    expbox=[[min(m[j]for m,c in rp),max(m[j]for m,c in rp)]for j in range(3)]
    prim=[x//content for x in ints];factorpre=dict(source=str((ROOT/'reports'/rs/'native.json').relative_to(ROOT)),input434_terms=True,total_degree=20,
      exponent_box=expbox,common_QQ_denominator=str(den),integer_content=str(content),primitive_coefficient_max_bits=max(abs(x).bit_length()for x in prim),
      no_factor_or_gcd_evaluated=True,previous_whole_p_division_not_exact=True,proper_p_factor_sharing_unknown=True,
      proposed_one_native_factor_or_gcd_budget=dict(seconds=60,address_space=2000000000,artifact_bytes=33554432,total_explicit_terms=150000),
      conservative_factor_dense_terms_upper=20*math.comb(23,3),
      full_factor_internal_growth_NOT_certified=True,authorization_required_before_any_new_scalar_job=True)
    dump(out/'R434-factor-preflight.json',factorpre)
    note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
    result=dict(id=ID,all_passed=True,frozen=True,actual_samecokernel_shape=[5,12],typed_nodes=len(nodes),
      parameter_localizers=['Delta5','Saux7022'],R434_h22_D4_not_inverted=True,new_scalar_values=0,
      all8_Q9_tails_retained=True,no_new_point_prime_minor_rank_GB=True,global_FULLU57_NOT_proved=True,
      input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},raw_artifact_bytes_before_report=size,seconds=time.monotonic()-start)
    dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
