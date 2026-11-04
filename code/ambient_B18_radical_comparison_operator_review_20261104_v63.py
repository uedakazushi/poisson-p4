#!/usr/bin/env python3
"""Independent polynomial, full-chain and typed-expression comparison review."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW';ZERO=(0,)*9
def read(p):
    with p.open()as f:return json.load(f)
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def dec(rows):
    r={}
    for m,c in rows:assert len(m)==9 and tuple(m)not in r and Q(c);r[tuple(m)]=Q(c)
    return r
def add(a,b):
    r=dict(a)
    for m,c in b.items():
        r[m]=r.get(m,Q(0))+c
        if not r[m]:del r[m]
    return r
def mul(a,b):
    r={}
    for m,c in a.items():
        for n,d in b.items():
            k=tuple(x+y for x,y in zip(m,n));r[k]=r.get(k,Q(0))+c*d
    return {m:c for m,c in r.items()if c}
def sm(a,c):return {m:v*c for m,v in a.items()if v*c}
def var(j):return {tuple(int(i==j)for i in range(9)):Q(1)}
def const(c):return {ZERO:Q(c)}if c else{}
def embed4(rows):
    r={}
    for m,c in rows:ex=[0]*9;ex[6],ex[7],ex[8],ex[5]=m;r[tuple(ex)]=Q(c)
    return r
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);R=ROOT/'reports';b=R/'AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63';rep=read(b/'report.json')
    assert sha(b/'report.json')=='79b94d1cf15d3ad05555f78c6cc06afc2c90833cd04134930f3ce72d4a6a65bf'and rep['all_passed'];ep=b/'actual-radical-normalizer-comparison-operator.json';assert sha(ep)=='376d5fd2b8f2ab8a45651b0374d5b983ebe2df8f8a507108675942e853adbba5';e=read(ep)
    for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
    oldctx=read(R/'AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json');B=read(ROOT/e['old_actual_chain_context']['path']);H=read(ROOT/e['honest_C0_P_context']['path']);nil=read(R/'AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json')
    for p in rep['input_sha256']:
        if p.endswith('_REVIEW/report.json'):assert read(ROOT/p)['all_passed']
    finitegate=R/'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json';assert read(finitegate)['all_passed'];finitenote=ROOT/'research/notes/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63.md';assert 'tensor後もSpec上surjective' in finitenote.read_text()
    assert e['accepted_native_C_coefficient_interface']['path']=='reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json'
    interface=oldctx['source_bound_interface'];adapt=oldctx['source_adapter_T'];assert interface['operator_has_no_unspecified_coefficients_or_nonunit_quotients']and adapt['d_never_inverted']
    assert adapt['core4_indices']==[0,1,2,6]and adapt['core4_order']==['a','d','v2','v6']and adapt['full_N']=='T(kappa)*R12'and adapt['full_N_inverse']=='R12^-1*T(kappa)^-1'
    assert B['all640_UP_equals_QB_verified']and B['scalar_references_keep_actual_original_r_g_functions_opaque']and B['physical_target_rows12']==[18,19,28,29,60,67,72,63,64,68,73,74]
    assert sha(ROOT/e['old_actual_chain_context']['path'])==e['old_actual_chain_context']['sha256']and sha(ROOT/e['honest_C0_P_context']['path'])==e['honest_C0_P_context']['sha256']
    assert B['scalar_source_artifacts']['C0']['path']==e['honest_C0_P_context']['path']and B['scalar_source_artifacts']['C0']['sha256']==e['honest_C0_P_context']['sha256']
    assert B['operator_scalar_nodes'][:2]==H['DAG_scalar_nodes'][:2]==[['c',0,1],['c',1,1]]
    # Complete native tensor complex; check every defining column independently.
    variables=[var(i)for i in range(9)];regular=variables[:6];c,r,l=variables[6:];quadrics=[mul(c,r),mul(c,l),mul(r,l)];R9=regular+quadrics;d=add(add(c,r),l)
    assert [dec(p)for p in e['R9_generators_polynomials9']]==R9 and dec(e['d_polynomial'])==d
    labels=[]
    for total in range(9):
        expected=[]
        for hb in range(3):
            k=total-hb
            if 0<=k<=6:
                for wedge in itertools.combinations(range(6),k):
                    for j in range([1,3,2][hb]):expected.append((tuple(wedge),hb,j))
        supplied=[(tuple(x['regular_wedge']),x['HB_degree'],x['HB_basis'])for x in e['native_free_basis_labels'][total]];assert supplied==expected;labels.append(supplied)
    ranks=[len(x)for x in labels];assert ranks==[1,9,35,77,105,91,49,15,2]==e['native_Betti'];assert e['restricted_R0_resolution_ranks']==[4*x for x in ranks]
    hb2=[{0:l,1:sm(r,-1)},{0:l,2:sm(c,-1)}];maps=[];native_terms=0;checked_columns=0
    for degree,block in enumerate(e['native_tensor_Koszul6_HB2_resolution'],1):
        assert block['source_degree']==degree and block['shape']==[ranks[degree-1],ranks[degree]];lookup={v:i for i,v in enumerate(labels[degree-1])};columns=[]
        for label,supplied in zip(labels[degree],block['columns']):
            I,bhb,j=label;want={}
            for k,v in enumerate(I):
                target=lookup[(I[:k]+I[k+1:],bhb,j)];want[target]=sm(variables[v],(-1)**k)
            if bhb:
                terms={0:quadrics[j]}if bhb==1 else hb2[j]
                for dest,p in terms.items():
                    target=lookup[(I,bhb-1,dest)];want[target]=add(want.get(target,{}),sm(p,(-1)**len(I)))
            got={i:dec(p)for i,p in supplied};assert got=={i:p for i,p in want.items()if p};native_terms+=sum(map(len,got.values()));columns.append(got);checked_columns+=1
        assert len(columns)==ranks[degree];maps.append(columns)
    chains=0
    for j in range(1,8):
        for col in maps[j]:
            result={}
            for k,p in col.items():
                for i,q in maps[j-1][k].items():result[i]=add(result.get(i,{}),mul(q,p))
            assert all(not p for p in result.values());chains+=1
    assert chains==374
    for col in maps[1]:
        p={}
        for i,q in col.items():p=add(p,mul(mul(d,R9[i]),q))
        assert not p
    D=[[dec(p)for p in row]for row in e['G8_polynomial_factor_D8_by9']];assert len(D)==8 and all(len(row)==9 for row in D)
    expectedD=[[{}for _ in range(9)]for _ in range(8)]
    for i in range(5):expectedD[i][i]=const(1)
    for i in range(3):
        for j in range(4):expectedD[5+i][5+j]=embed4(nil['source_factorization_coefficients4_by3'][j][i])
    assert D==expectedD
    oldCI=adapt['CI_generators3']
    for i in range(3):assert {tuple(m):Q(n,den)for m,n,den in oldCI[i]}=={tuple(m):Q(v)for m,v in nil[['old_source_Q0','old_source_Q1','old_source_f'][i]]}
    G=variables[:5]+[embed4(p)for p in nil['actual_core_generators3']]
    for i in range(8):
        actual={}
        for j in range(9):actual=add(actual,mul(D[i][j],R9[j]))
        assert actual==G[i]
    # Semantic fingerprints distinguish local C/R0 nodes from foreign IDs.
    descriptors={};counter=Counter();fingerprints=[];types=[]
    def key(tag,*payload):
        value=hashlib.sha256(json.dumps([tag,*payload],separators=(',',':')).encode()).hexdigest();descriptors.setdefault(value,(tag,payload));return value
    def cq(q):return key('Q',str(Q(q)))
    zz,oo=cq(0),cq(1)
    def cs(xs):
        xs=sorted(x for x in xs if x!=zz)
        return zz if not xs else xs[0]if len(xs)==1 else key('Csum',xs)
    def cm(xs):
        if zz in xs:return zz
        xs=sorted(x for x in xs if x!=oo)
        return oo if not xs else xs[0]if len(xs)==1 else key('Cproduct',xs)
    def comp(x,i):
        tag,payload=descriptors[x]
        if tag=='Q':return x if i==0 else zz
        if tag=='embedOld':return key('oldR0',payload[0])if i==0 else zz
        if tag=='basis':return oo if i==payload[0]else zz
        return key('component',x,i)
    new=e['new_Cscalar_and_R0component_nodes'];assert len(new)==53610
    foreignB=len(B['operator_scalar_nodes']);foreignH=len(H['DAG_scalar_nodes'])
    for i,n in enumerate(new):
        tag=n[0];counter[tag]+=1
        if tag=='C_QQ':f=cq(Q(n[1],n[2]));typ='both'
        elif tag in ['actual_proper_u','actual_core_GL4','C_basis']:
            if tag=='actual_proper_u':assert 0<=n[1]<5;f=key('u',n[1])
            elif tag=='actual_core_GL4':assert n[1]in ['h','c','r','l'];f=key('core',n[1])
            else:assert 1<=n[1]<4;f=key('basis',n[1])
            typ='C'
        elif tag=='accepted_native_unit_inverse_N12_entry':assert 0<=n[1]<12 and 0<=n[2]<8;f=key('Ninv',n[1],n[2]);typ='C'
        elif tag in ['C_sum','C_product']:
            assert n[1]and all(0<=x<i and types[x]in ['C','both']for x in n[1]);f=(cs if tag=='C_sum'else cm)([fingerprints[x]for x in n[1]]);typ='C'
        elif tag=='C_embed_old_R0':assert 1<n[1]<foreignB;f=key('embedOld',n[1]);typ='C'
        elif tag=='R0_component_of_C':assert 0<=n[1]<i and types[n[1]]in ['C','both']and 0<=n[2]<4;f=comp(fingerprints[n[1]],n[2]);typ='R0'
        elif tag=='R0_old_operator':assert 0<=n[1]<foreignB;f=key('oldR0',n[1]);typ='R0'
        elif tag=='R0_original_honest_C0_scalar':assert 1<n[1]<foreignH;f=key('honestR0',n[1]);typ='R0'
        else:assert tag=='R0_neg'and 0<=n[1]<i and types[n[1]]in ['R0','both'];f=key('R0neg',fingerprints[n[1]]);typ='R0'
        fingerprints.append(f);types.append(typ)
    assert fingerprints[:2]==[zz,oo]
    cv=[key('u',i)for i in range(5)]+[key('core','h')]+[key('core',x)for x in ['c','r','l']]
    def cp(p):return cs([cm([cq(v),*[cv[j]for j,k in enumerate(m)for _ in range(k)]])for m,v in p.items()])
    A=e['native_new_A12_by9_nodes'];expectedA=[]
    for i in range(12):
        row=[]
        for j in range(9):
            want=cs([cm([key('Ninv',i,k),cp(D[k][j])])for k in range(8)]);assert fingerprints[A[i][j]]==want and types[A[i][j]]in ['C','both'];row.append(want)
        expectedA.append(row)
    basis=[oo,*[key('basis',i)for i in range(1,4)]];dn=cp(d);qn=[cm([dn,cp(p)])for p in R9];Q36=e['actual_Qrad4_by36_columns']
    assert len(Q36)==36
    for j in range(9):
        for a in range(4):assert [fingerprints[x]for x in Q36[4*j+a]]==[comp(cm([qn[j],basis[a]]),k)for k in range(4)]
    S140=e['actual_complete_SyzQrad36_by140_sparse_columns'];assert len(S140)==140
    for j,col in enumerate(maps[1]):
        for a in range(4):
            want={}
            for i,p in col.items():
                x=cm([cp(p),basis[a]])
                for k in range(4):
                    y=comp(x,k)
                    if y!=zz:want[4*i+k]=y
            assert {i:fingerprints[v]for i,v in S140[4*j+a]}==want
    def old(x):return [zz,oo][x]if x in [0,1]else key('embedOld',x)
    Res=e['actual_Res4_Atranspose36_by48_columns'];assert len(Res)==48
    for i in range(12):
        for a in range(4):assert [fingerprints[v]for v in Res[4*i+a]]==[comp(cm([expectedA[i][j],basis[a]]),k)for j in range(9)for k in range(4)]
    oldcols=B['B48_by640_sparse_columns'];BC=e['actual_Brad9_C_values_by640_columns'];BR=e['actual_Brad36_R0_values_by640_columns'];assert len(oldcols)==len(BC)==len(BR)==640
    for col,bc,br in zip(oldcols,BC,BR):
        grouped={}
        for i,v in col:assert 0<=i<48 and 0<=v<foreignB;grouped.setdefault(i//4,[]).append(cm([old(v),basis[i%4]]))
        group={i:cs(v)for i,v in grouped.items()};want=[cs([cm([expectedA[i][j],v])for i,v in group.items()])for j in range(9)]
        assert [fingerprints[v]for v in bc]==want;assert [fingerprints[v]for v in br]==[comp(v,k)for v in want for k in range(4)]
    P=H['honest_D0_relation_columns'];P8=e['actual_P8_by640_columns'];assert len(P)==len(P8)==640
    for col,p in zip(P,P8):assert [fingerprints[v]for v in p]==[[zz,oo][i]if i in [0,1]else key('honestR0',i)for i in col[1:]]
    U=[[comp(old(v),0)for v in col]for col in B['exact_existing_U9_C4']];assert U[0]==[oo,zz,zz,zz]
    V44=e['actual_V3_by44_rows_after_unit1_cancellation'];assert len(V44)==3
    for i,row in enumerate(V44,1):assert [fingerprints[v]for v in row]==[fingerprints[Q36[j][i]]for j in range(36)]+[key('R0neg',U[j][i])if U[j][i]!=zz else zz for j in range(1,9)]
    X=e['actual_X44_by780'];assert X['unreduced_shapes']=={'V':[4,45],'X':[45,780]}and X['final_shapes']=={'V':[3,44],'X':[44,780]}
    assert 'stack(actual_Brad36_R0_values_by640_columns,actual_P8_by640_columns)'==X['first640_columns']and 'stack(actual_complete_SyzQrad36_by140_sparse_columns,zero8)'==X['last140_columns']
    assert e['actual_C_variable_rules']['actual_core_GL4']==nil['forward_GL4']and e['actual_C_variable_rules']['GL4_variable_input_order']==['a','d','x','t']
    assert e['NilD0_exactly_H1_of_actual_complex']and e['no_actual_D0_nil_generator_list_cube_or_depth_claimed']
    audit={'all53610_typed_nodes_and_local_foreign_reference_bounds_checked':True,'foreign_oldB_operator_nodes':foreignB,'foreign_honestC0_scalar_nodes':foreignH,'all8_actual_G8_D_R9_polynomial_identities':True,'all108_A_entries_match_same_Ninverse_E8_D':True,'native_all_basis_and_tensor_differentials_checked':checked_columns,'all_native_chain_columns_checked':chains,'native_Betti':ranks,'actual_Res4_Atranspose_all48_columns_checked':True,'no_inner_mu_transpose':True,'all640_actual_Brad_Cvalues_and_R0components_checked':True,'actual_complete_Srad140_checked':True,'all640_P8_and_V3x44_entries_checked':True,'source_variable_rules_and_unit_cancellation_checked':True,'all_old_r_g_foreign_roots_preserved':True,'new_polynomial_terms':native_terms+sum(len(p)for row in D for p in row),'typed_tags':dict(counter),'no_scalar_series_source_jet_point_GB_or_old_function_division':True}
    assert audit['new_polynomial_terms']==1654;dump(out/'arithmetic-and-typed-audit.json',audit)
    proof={'native_exactness':'Tensor a Koszul resolution on(u5,h) with the HB resolution of(cr,cl,rl). The six variables are regular on the three-axis quotient. Multiplying only d1 by ambient regular d preserves higher kernels. Allentries liein maximalideal, so native C-resolution isminimal.','Res4':'The SAME C isfree4 over SAME physicalR0; restriction exact. F_old=A(dR9) gives Qold=Qrad Res4(A^t), block(j,i)=mu(Aij), not mu(Aij)^t.','chain':'All640 old source-bound UP=QoldB identities transport via A^t. Every native35 syzygy restricts to4actualcolumns, giving COMPLETE Srad140, not apartialkernel.','homology':'For Qrad*b=U*p map(b,p)to[p]inD0. Ifp=Pz subtract(Brad*z,Pz); the remainder(b-Brad*z,0) is killed exactly bycompleteSrad140. This proves H1=ker(D0->Ered).','unit_cancel':'U0=e0; eliminate row0 using sourcecolumn-U0. Chainidentity forces removedXrowzero; otherXrows unchanged, yielding44x780->3x44.','nil':'D0->E4 finite integral and Spec-surjective comesfromthe actualtwo finite normalizer covers and their tensor. Passing to Ered retains integral Spec-surjectivity. Its kernel is contained in everyprime ofD0, while reducedEred kills every nilpotent, so ker=Nil(D0).','scope':'Exact actualhomology model only; no D0nilgeneratorlist, principality, cubeindex, socle, depth or finiteTlength computed.'};dump(out/'proof.json',proof)
    note=ROOT/f'research/notes/{ID}.md';assert note.exists();result={'id':ID,'all_passed':True,'frozen':True,'complete_source_bound_radical_comparison':True,'NilD0_equals_H1':True,'actual_X_shape':[44,780],'actual_V_shape':[3,44],'no_actual_nil_generators_or_depth_claim':True,'source_report_sha256':sha(b/'report.json'),'source_evidence_sha256':sha(ep),'integral_Spec_surjection_gate_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [finitegate,finitenote]},'arithmetic':audit,'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note]},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);dump(out/'report.json',result);print(json.dumps({'id':ID,'all_passed':True,'seconds':result['seconds'],'report_sha256':sha(out/'report.json'),'typed_nodes':len(new),'native_chain_columns':chains},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
