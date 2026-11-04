#!/usr/bin/env python3
"""Typed direct4+5 unit Schur reduction, no h22 inverse or scalar expansion."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63'
def read(p):
    with p.open()as f:return json.load(f)
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def pol(p):return{tuple(m):Q(c)for m,c in p if Q(c)}
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(55);out.mkdir(parents=True,exist_ok=False);R=ROOT/'reports';pins={}
    def load(stem,name='certificate.json'):
        p=R/stem/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    qstem='S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63';pstem='S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63';cstem='S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63';hstem='S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63';gstem='S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';sstem='S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63'
    q=load(qstem);pair=load(pstem);cp=load(cstem,'checkpoint.json');hc=load(hstem);rank=load(gstem);strip=load(sstem,'native.json');h9=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json')
    for stem in [qstem,pstem,cstem,hstem,gstem,'S61_NONSQUARE_PREPOWER_H22_FIXED_MINORS_20261104_v63']:
        peer=load(stem+'_REVIEW','report.json');assert peer['all_passed']
    assert pins['reports/S61_NONSQUARE_PREPOWER_H22_FIXED_MINORS_20261104_v63_REVIEW/report.json']=='f2d235534778c7bc868c6e013ee85c9bd1f81e45c02575fdaf4968178c9b4037'
    roots=q['remaining15x22_scalar_roots'];oldlabels=q['new15_source_original25_labels'];cols=q['new22_target_original32_columns'];pure=cols[:17];heavy=cols[17:];assert oldlabels==[9]+hc['final14_source_labels']and heavy==[27,28,29,30,31]
    Prows=[0,1,3,6];Pcols=[22,11,13,25];I=[1,3,6,12,9];pi=[0,4,1,3,2];Hrows=[i+1 for i in I];Hcols=[heavy[j]for j in pi];pivotrows=Prows+Hrows;pivotcols=[cols.index(j)for j in Pcols+Hcols];restrows=[i for i in range(15)if i not in pivotrows];restcols=[i for i in range(22)if i not in pivotcols]
    assert len(set(pivotrows))==9 and Hrows==[2,4,7,13,10]and restrows==[5,8,9,11,12,14]and [cols[i]for i in restcols]==[1,2,3,4,6,7,8,9,14,15,17,18,19]
    # Resolve only the structural zero/leaf identities in the accepted DAG.
    QN=q['scalar_DAG']['nodes'];cache={}
    def resolve(i):
        if i in cache:return cache[i]
        n=QN[i];op=n[0]
        if op=='z':v=0
        elif op=='one':v=1
        elif op=='checkpoint':v=('checkpoint',n[1],n[2]) if cp['matrix25x32'][n[1]][n[2]][0] else 0
        elif op=='accepted_pair_second':v=('accepted_pair_second',n[1])
        elif op=='new_unit_atom':v=('new_unit_atom',n[1])
        else:
            a,b=resolve(n[1]),resolve(n[2])
            if op=='m':v=0 if a==0 or b==0 else b if a==1 else a if b==1 else('m',a,b)
            elif op=='s':v=0 if a==b else a if b==0 else('s',a,b)
            else:assert op=='a';v=b if a==0 else a if b==0 else('a',a,b)
        cache[i]=v;return v
    for i in Prows:
        for j in range(17,22):assert resolve(roots[i][j])==0
    for lab,col,trace in hc['all70_literal_external_root_traces']:
        i=oldlabels.index(lab);want=0 if trace==[0]else tuple(trace[-1]);assert resolve(roots[i][cols.index(col)])==want
    N=hc['raw_compact14x5'];assert all(not N[i][j]for i in I[:2]for j in pi[2:])
    assert [rank['input14_remaining_rows'][i]for i in I]==rank['nonzero_minor_of14x5_rows']and [rank['input5_columns'][j]for j in pi]==rank['nonzero_minor_of14x5_cols']
    assert [rank['input14_remaining_rows'][i]for i in I[:2]]==rank['sparse_pair_rows']and [rank['input5_columns'][j]for j in pi[:2]]==rank['sparse_pair_columns']
    for i,br in zip(I[2:],rank['selected_rows7']):
        for j,bc in zip(pi[2:],rank['selected_cols3']):assert pol(N[i][j])==pol(rank['actual_block7x3'][br][bc])
    assert strip['named_unit_powers_a_d_F_B_C_J']==[4,0,0,0,3,0]and len(strip['exact_unscaled_pure4_minor'])==818 and len(strip['remaining_polynomial'])==434
    nodes=[]
    def n(op,shape,args=(),**kw):
        assert len(nodes)<150000 and len(shape)==2 and all(k>0 for k in shape)and all(0<=k<len(nodes)for k in args);i=len(nodes);nodes.append(dict(id=i,op=op,shape=list(shape),args=list(args),**kw));return i
    def ref(stem,name,pointer,shape,**kw):
        p=R/stem/name;pins[str(p.relative_to(ROOT))]=sha(p);return n('external_frozen_operator',shape,reference=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),json_pointer=pointer),**kw)
    def prod(a,b):sa,sb=nodes[a]['shape'],nodes[b]['shape'];assert sa[1]==sb[0];return n('matrix_product',[sa[0],sb[1]],[a,b])
    def neg(a):return n('negative',nodes[a]['shape'],[a])
    def sub(a,b):assert nodes[a]['shape']==nodes[b]['shape'];return n('difference',nodes[a]['shape'],[a,b])
    def scale(c,a):assert nodes[c]['shape']==[1,1];return n('central_scalar_multiple',nodes[a]['shape'],[c,a])
    def sel(a,rows,columns):return n('submatrix',[len(rows),len(columns)],[a],rows=rows,columns=columns)
    def ident(k):return n('identity',[k,k])
    def zero(m,k):return n('zero',[m,k])
    def adj(a):s=nodes[a]['shape'];assert s[0]==s[1];return n('universal_polynomial_cofactor_adjugate',s,[a],expanded=False)
    def inv(a,proof):assert nodes[a]['shape']==[1,1];return n('certified_chart_unit_inverse',[1,1],[a],unit_proof=proof,expanded=False)
    def block(rows):
        heights=[nodes[r[0]]['shape'][0]for r in rows];widths=[nodes[k]['shape'][1]for k in rows[0]]
        for i,row in enumerate(rows):
            assert len(row)==len(widths)
            for j,k in enumerate(row):assert nodes[k]['shape']==[heights[i],widths[j]]
        return n('block_matrix',[sum(heights),sum(widths)],[k for row in rows for k in row],block_rows=rows)
    M=ref(qstem,'certificate.json','/remaining15x22_scalar_roots',[15,22],scalar_interpreter='/scalar_DAG',source_scope='Accepted actualprePOWER15x22; externalroot IDs stayforeign.')
    rho=[ref(hstem,'certificate.json',f'/known_row_units_rho14/{i}',[1,1],scalar_semantics='constant times named a,d,F,B,C,J powers')for i in range(14)]
    irho=[inv(v,'Only nonzeroQQ and knownFULLU factors.')for v in rho];gamma=ref(hstem,'certificate.json','/known_column_units/4',[1,1],scalar_semantics='knownFULLU monomial');igamma=inv(gamma,'Only knownFULLU units.')
    one=n('one',[1,1]);Sdiag=n('diagonal',[15,15],[one]+[irho[i]if i in I else one for i in range(14)],diagonal=[one]+[irho[i]if i in I else one for i in range(14)],hybrid_only_H5_rows_scaled=True)
    Tdiag=n('diagonal',[22,22],[igamma if j==21 else one for j in range(22)],diagonal=[igamma if j==21 else one for j in range(22)])
    SdiagInv=n('diagonal',[15,15],[one]+[rho[i]if i in I else one for i in range(14)],diagonal=[one]+[rho[i]if i in I else one for i in range(14)]);TdiagInv=n('diagonal',[22,22],[gamma if j==21 else one for j in range(22)],diagonal=[gamma if j==21 else one for j in range(22)])
    Prow=n('row_permutation',[15,15],order=pivotrows+restrows);Pcol=n('column_permutation',[22,22],order=pivotcols+restcols)
    ProwInv=n('transpose',[15,15],[Prow]);PcolInv=n('transpose',[22,22],[Pcol])
    MM=prod(prod(prod(Prow,prod(Sdiag,M)),Tdiag),Pcol)
    P4=sel(MM,list(range(4)),list(range(4)));C54=sel(MM,list(range(4,9)),list(range(4)));H5=sel(MM,list(range(4,9)),list(range(4,9)))
    pair2=sel(H5,[0,1],[0,1]);C32=sel(H5,[2,3,4],[0,1]);B3=sel(H5,[2,3,4],[2,3,4]);A2,A3=adj(pair2),adj(B3)
    delta2=ref(gstem,'certificate.json','/literal_pair_determinant',[1,1]);delta3=ref(gstem,'certificate.json','/nonzero_selected_minor_Delta',[1,1]);Delta=ref(gstem,'certificate.json','/nonzero_minor_of14x5',[1,1]);iDelta=inv(Delta,'NEW chartD(Delta5). Delta5=delta2*delta3, with unchanged selectedrow/col ordering.')
    Hnum=block([[scale(delta3,A2),zero(2,3)],[neg(prod(prod(A3,C32),A2)),scale(delta2,A3)]]);Hinv=scale(iDelta,Hnum)
    D4=ref(sstem,'native.json','/exact_unscaled_pure4_minor',[1,1]);R434=ref(sstem,'native.json','/remaining_polynomial',[1,1]);iD4=inv(D4,'D4=a^4*C^3*R434; a,C knownFULLUunits, NEW chartD(R434). No h22 inverse.');P4inv=scale(iD4,adj(P4))
    K9=sel(MM,list(range(9)),list(range(9)));Kinv=block([[P4inv,zero(4,5)],[neg(prod(prod(Hinv,C54),P4inv)),Hinv]])
    B=sel(MM,list(range(9)),list(range(9,22)));C=sel(MM,list(range(9,15)),list(range(9)));D=sel(MM,list(range(9,15)),list(range(9,22)));KB=prod(Kinv,B);CK=prod(C,Kinv);residual=sub(D,prod(C,KB))
    L15=block([[Kinv,zero(9,6)],[neg(CK),ident(6)]]);Linv=block([[K9,zero(9,6)],[C,ident(6)]]);U22=block([[ident(9),neg(KB)],[zero(13,9),ident(13)]]);Uinv=block([[ident(9),KB],[zero(13,9),ident(13)]])
    fullsource=prod(prod(L15,Prow),Sdiag);fulltarget=prod(prod(Tdiag,Pcol),U22);fullsourceInv=prod(prod(SdiagInv,ProwInv),Linv);fulltargetInv=prod(prod(Uinv,PcolInv),TdiagInv);survivingsource=sel(fullsource,list(range(9,15)),list(range(15)))
    # Actual covector transport to the original25 source coordinates; no old H5adj.
    S16=ref(qstem,'certificate.json','/source16_roots',[16,16],scalar_interpreter='/scalar_DAG');S17=ref(pstem,'certificate.json','/source17_roots',[17,17],scalar_interpreter='/scalar_DAG')
    assert oldlabels==q['old16_source_original25_rows'][1:]and q['old16_source_original25_rows']==pair['old17_source_original25_rows'][1:]==cp['remaining_rows'][1:]
    E15=n('coordinate_row_selector',[15,16],rows=list(range(1,16)));E16=n('coordinate_row_selector',[16,17],rows=list(range(1,17)));embed17=n('coordinate_row_embedding',[17,25],columns=cp['remaining_rows']);prefix=ident(25);stepnodes=[]
    for j,st in enumerate(cp['all_completed_unit_steps']):
        p=ref(cstem,'checkpoint.json',f'/all_completed_unit_steps/{j}/pivot',[1,1]);ip=inv(p,st['known_unit']);step=n('accepted_forward_Gaussian_source_step',[25,25],[p,ip],pivot_row=st['pivot_row'],targets=[t for t,cof in st['eliminations']],cofactor_ref=dict(path=f'reports/{cstem}/checkpoint.json',json_pointer=f'/all_completed_unit_steps/{j}/eliminations'),formula='rowr:=rowr/p; rowt:=rowt-c_t*rowr. This is FORWARDsource operator, not its inverse.',inverse_formula='oldrowr=p*newrowr; oldrowt=newrowt+c_t*newrowr.');prefix=prod(step,prefix);stepnodes.append(step)
    T15=prod(prod(prod(prod(prod(E15,S16),E16),S17),embed17),prefix);source6old25=prod(survivingsource,T15)
    assert nodes[residual]['shape']==[6,13]and nodes[source6old25]['shape']==[6,25]
    extra8=[j for j in pure if j not in h9['original_target_labels']];assert extra8==[1,2,3,4,6,7,8,9]and all(j in [cols[k]for k in restcols]for j in extra8)
    operator={'coefficient_ring':'FULLU=QQ[a,e,d,1/(a*d*F*B*C*J)], additionally localize at R434*Delta5 ONLY.','nodes':nodes,'new_expanded_parameter_coefficients':0,'adjugates_unexpanded':'Only universal adj4,adj2,adj3 are used; no freshadj5 or scalarcofactor expansion.','exports':dict(input15x22=M,hybrid_row_unit=Sdiag,heavy_column_unit=Tdiag,row_permutation=Prow,column_permutation=Pcol,permuted15x22=MM,P4=P4,H5=H5,H5_pair2=pair2,H5_B3=B3,H5_cofactor_numerator=Hnum,H5_inverse=Hinv,D4=D4,R434=R434,Delta5=Delta,pivot9=K9,pivot9_inverse=Kinv,source_left15=L15,source_left15_inverse=Linv,target_right22=U22,target_right22_inverse=Uinv,full_source15=fullsource,full_source15_inverse=fullsourceInv,full_target22=fulltarget,full_target22_inverse=fulltargetInv,residual6x13=residual,surviving6_source_current15=survivingsource,source_current15_to_original25=T15,surviving6_source_original25=source6old25),'context':dict(P4_current15_rows=Prows,P4_original25_sources=[oldlabels[i]for i in Prows],P4_original32_targets=Pcols,H5_compact14_rows=I,H5_current15_rows=Hrows,H5_original25_sources=[oldlabels[i]for i in Hrows],H5_compact5_column_order=pi,H5_original32_targets=Hcols,pivot9_source_order=pivotrows,pivot9_target_current22_order=pivotcols,remaining6_current15=restrows,remaining6_original25=[oldlabels[i]for i in restrows],remaining13_original32_seed_labels=[cols[i]for i in restcols],Q9_selected_h9_labels=h9['original_target_labels'],Q9_eight_other_pure_tails=[dict(original32=j,actual_root=roots[0][cols.index(j)],remaining13_position=[cols[k]for k in restcols].index(j))for j in extra8],all13_Q9_Bblock_entries_retained=True,pure4_determinant='D4=a^4*C^3*R434 (unscaled P4rows)',H5_determinant='Delta5=-(6/5)*B*p^2*Delta407',hybrid_scaling='ONLY selectedH5rows dividedbyrho; P4rows unchanged; lastheavycol dividedbygamma.',chart='FULLU intersectD(R434*Delta5); h22 is never inverted.',source_target_order='row presentation: L15*Prow*Sdiag*M15*Tdiag*Pcol*U22=diag(I9,residual6x13)',uncomputed_operators='All matrix/cofactor/Schur values remain finite typed expressions; no scalar residual entry was evaluated.'),'identities':['H5=[[Pair2,0],[C32,B3]]','H5*Hnum=Hnum*H5=Delta5*I5 by existing2+3adj identities','P4*Adj4(P4)=Adj4(P4)*P4=D4*I4','pivot9=[[P4,0],[C54,H5]],detpivot9=D4*Delta5','pivot9*pivot9_inverse=pivot9_inverse*pivot9=I9','L15*Linv=Linv*L15=I15; U22*Uinv=Uinv*U22=I22','L15*permuted15x22*U22=diag(I9,residual6x13)','coker(input15x22)=coker(residual6x13) onthe statedchart','rank(input15x22)=9+rank(residual6x13); fullrowrank6 NOTcomputed'],'not_claimed':['R434 orDelta5globalFULLUunits','h22inverse orh22chartrestriction','A freebasis of a projectivekernel','rank6 orfullcubic57 orquartic239','Removal ofR434=0 orDelta5=0 complements','Raw16 H7 kernelbasis withoutthe separateCdelta lift']}
    dump(out/'operator.json',operator);proof={'P4_and_heavy5_zero_source_checks':True,'all70_actual_heavy_root_traces_preserved':True,'all_staged_B3_and_pair_source_indices_match':True,'all9_pivot_sources_disjoint':True,'hybrid_row_units_keep_literal_D4':True,'H5_two_plus_three_inverse_does_not_invert_p_separately':True,'full8_Q9_other_pure_tails_retained':True,'all_block_shapes_checked':True,'cokernel_theorem':'The displayed explicit invertible source/target transforms give diag(I9,S6x13), so cancelling identity9 gives samecokernel with free6source and13target presentation; this doesnot assert the originalmap kernel isfree.','new_scalar_coefficients':0,'chart_complements_preserved':['R434=0','Delta5=0']};dump(out/'proof.json',proof)
    note=ROOT/f'research/notes/{ID}.md';assert note.exists();result={'id':ID,'all_structural_source_and_type_checks_passed':True,'frozen':True,'actual_samecoker_shape':[6,13],'source9_pivot_shape':[9,9],'no_h22_inverse':True,'chart':'FULLU D(R434*Delta5)','all8_other_Q9_pure_tails_retained':True,'typed_nodes':len(nodes),'new_expanded_parameter_coefficients':0,'no_new_minor_rank_or_NF_job':True,'no_rank6_or_full57_claim':True,'input_sha256':pins,'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note]},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);dump(out/'report.json',result);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','source_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
