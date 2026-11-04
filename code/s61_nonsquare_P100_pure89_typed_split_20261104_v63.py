#!/usr/bin/env python3
"""Source-bound finite typed global9 splitting; no scalar coefficient expansion."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE89_TYPED_SPLIT_20261104_v63'
def read(p):
    with p.open()as f:return json.load(f)
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(55);out.mkdir(parents=True,exist_ok=False);R=ROOT/'reports';pins={}
    def load(stem,name='report.json'):
        p=R/stem/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    pp='S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63_REVIEW';peer=load(pp);assert pins[f'reports/{pp}/report.json']=='1125cf14c2594385395aee6a397e66acd9616230a3c4989ce709fcfdf6df1076'and peer['all_passed']and peer['positive_full239_on_entire_H15_only']
    cp='S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63_REVIEW';coverpeer=load(cp);assert pins[f'reports/{cp}/report.json']=='2b4abe5029157c3aa08dbdd9ba26f256cb9b49b03bfa77146f88cf91d1e48bcd'and coverpeer['all_passed']
    typed=load('S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63','operator.json');ctx=typed['context'];oldsrc=load('S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63','source-binding.json');old=load('S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63','native.json');g=load('S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63','native.json');fac=load('S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63','native.json')
    assert g['gcd_degree']==15 and g['all5_literal_gcd_times_quotient_equals_input_verified']and g['oldprime_cofactor_gcd_is1']and g['cofactor_reduction_degrees_preserved']and g['literal_modular_Bezout_verified']
    C=ctx['selected_local98_columns8']+[25];T=[j for j in range(98)if j not in C];assert len(C)==9 and len(T)==89 and C==[0,1,2,3,4,5,6,9,25]
    assert all(c['nine_local98_target_columns']==C for c in oldsrc['five_cases'])and len(oldsrc['five_cases'])==len(old['determinants5'])==5
    assert ctx['H1_source_frame'].startswith('UsePiActual=delta7*PiOrig')and typed['nodes'][typed['exports']['full_Gamma319x98']]['shape']==[319,98]
    nodes=[]
    def n(op,shape,args=(),**kw):
        i=len(nodes);assert len(shape)==2 and all(isinstance(k,int)and k>0 for k in shape);assert all(0<=j<i for j in args);nodes.append(dict(id=i,op=op,shape=list(shape),args=list(args),**kw));return i
    def ext(stem,name,pointer,shape):
        p=R/stem/name;pins[str(p.relative_to(ROOT))]=sha(p);return n('external_exact_operator',shape,reference=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),json_pointer=pointer))
    def tr(a):s=nodes[a]['shape'];return n('transpose',[s[1],s[0]],[a])
    def prod(a,b):sa,sb=nodes[a]['shape'],nodes[b]['shape'];assert sa[1]==sb[0];return n('matrix_product',[sa[0],sb[1]],[a,b])
    def add(a,b):assert nodes[a]['shape']==nodes[b]['shape'];return n('sum',nodes[a]['shape'],[a,b])
    def neg(a):return n('negate',nodes[a]['shape'],[a])
    def scale(c,a):assert nodes[c]['shape']==[1,1];return n('central_scalar_multiple',nodes[a]['shape'],[c,a])
    def adj(a):s=nodes[a]['shape'];assert s[0]==s[1];return n('universal_polynomial_adjugate',s,[a],evaluated=False)
    def det(a):s=nodes[a]['shape'];assert s[0]==s[1];return n('universal_polynomial_determinant',[1,1],[a],evaluated=False)
    def identity(k):return n('identity',[k,k])
    def scalarref(stem,name,pointer,kind='QQ_polynomial'):
        z=ext(stem,name,pointer,[1,1]);nodes[z]['coefficient_kind']=kind;return z
    def invunit(a,proof):assert nodes[a]['shape']==[1,1];return n('certified_regular_unit_inverse',[1,1],[a],proof=proof,evaluated=False)
    def sumall(xs):
        assert xs;z=xs[0]
        for x in xs[1:]:z=add(z,x)
        return z
    Gamma=ext('S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63','operator.json',f"/nodes/{typed['exports']['full_Gamma319x98']}",[319,98]);A=tr(Gamma)
    E=n('coordinate_column_embedding',[98,9],columns=C);Pi=tr(E);ET=n('coordinate_column_embedding',[98,89],columns=T);PiT=tr(ET);F=prod(Pi,A)
    J=ext('S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63','operator.json',f"/nodes/{typed['exports']['source_embedding15x319']}",[15,319]);Jt=tr(J)
    H=scalarref('S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63','native.json','/exact_gcd_monic');one=n('constant_one',[1,1])
    Q=[scalarref('S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63','native.json',f'/all5_quotients/{i}')for i in range(5)]
    beta=n('canonical_QQ_polynomial_extendedgcd_family',[5,1],Q,algorithm='Fixed order Qi: startg=0, repeated monic QQ[z] xgcd, propagate earlier weights; finalg=1.',coprime_certificate='Saved literal divisions and degree-preserved Fp32003 cofactor Bezout1, independently accepted.',evaluated=False,not_a_saved_QQ_Bezout=True)
    RHparts=[];oldbindings=[]
    for i,(case,d)in enumerate(zip(oldsrc['five_cases'],old['determinants5'])):
        assert case['old_cover_basis_index']==d['cover_basis_index'];assert Fraction(d['QQ_unit_content'])!=0
        Wi=n('known_unit_weighted_row_selector',[9,15],rows=case['nine_source_rows_in_order'],LCD_powers5=case['effective_LCD_powers5'],five_factors=['z','z+20','z+60','203*z+2700','z+12'])
        S=prod(Jt,tr(Wi));X=prod(F,S);Di=det(X)
        Ui=n('saved_known_unit_content_times_five_factor_product',[1,1],QQ_content_reference=f"/determinants5/{i}/QQ_unit_content",powers=d['known_unit_multiplicities5'],source_native_path='reports/S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63/native.json',not_unknown_factor=True)
        iUi=invunit(Ui,'Nonzero QQ content and only the five whole-P100 curve units.')
        bi=n('component_polynomial',[1,1],[beta],index=i);RHparts.append(scale(prod(bi,iUi),prod(S,adj(X))))
        oldbindings.append(dict(case=d['cover_basis_index'],source_selector=Wi,actual_projected_matrix=X,determinant=Di,unit=Ui,quotient=Q[i],certified_identity='detX=Ui*H*Qi',literal_native_determinant_pointer=f'/determinants5/{i}/exact_determinant'))
    RH=sumall(RHparts)
    finiteparts=[];gparts=[];finitebindings=[]
    Hbar=n('saved_Fp_monic_polynomial',[1,1],coefficients=fac['Hbar'],prime=32003)
    for j,f in enumerate(fac['factors']):
        rc=load('S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63',f'factor{j}.json');rows=rc['selected_source319_rows'];assert rc['rank']==98 and len(set(rows))==98
        Sj=n('coordinate_column_embedding',[319,98],columns=rows);Bj=prod(A,Sj);Dj=det(Bj);Rj=prod(prod(Sj,adj(Bj)),E)
        ej=n('canonical_integer_lift_Fp_CRT_idempotent',[1,1],[Hbar],factor_index=j,factor=f['coefficients'],prime=32003,algorithm='h_j=Hbar/f_j; v_j=canonicalmonicFpExtendedGCD(h_j,f_j).inverse; reduce(h_j*v_j,Hbar); lift every coefficient to0..32002.',evaluated=False,not_a_saved_idempotent=True)
        finiteparts.append(scale(ej,Rj));gparts.append(prod(ej,Dj));finitebindings.append(dict(factor=j,selected_source_embedding=Sj,global_selected_matrix=Bj,global_Dj=Dj,projected_adjugate_source_lift=Rj,integer_CRT_lift=ej,saved_nonzero_Dj_reduction=rc['selected98_minor_determinant'],certified_identity='A*Rj=Dj*E; Pi*A*Rj=Dj*I9'))
    G=sumall(gparts);Rfinite=sumall(finiteparts)
    MG=n('canonical_multiplication_mod_monic_polynomial',[15,15],[G,H],basis='1,z,...,z^14',algorithm='Reduce each z^k*G in QQ[z]/H. Invert its known curve denominators moduloH by canonical extended Euclid; all are coprime toH.',evaluated=False)
    norm=det(MG);inorm=invunit(norm,'A QQ constant, not a function of z. Its reduction is Norm_Sbar(Gbar); CRT gives Gbar|fj=Dj|fj nonzero in every factor, so norm is a32003-unit and nonzeroQQ.')
    e0=n('constant_standard_column',[15,1],index=0);bvec=scale(inorm,prod(adj(MG),e0));b=n('coefficient_column_to_QQ_polynomial',[1,1],[bvec],basis='1,z,...,z^14',evaluated=False)
    numerator=add(one,neg(prod(b,G)));a=n('certified_exact_division_in_curve_unit_ring',[1,1],[numerator,H],algorithm='Choose a finite common denominator D consisting only of QQ constants and five curve-unit factors. For numerator=N/D, use exact QQ[z] quotient N/H and return (N/H)/D.',divisibility_proof='MG*bvec=e0 =>bG=1 moduloH; gcd(H,D)=1 implies H divides N literally. H is never inverted as a regular unit.',evaluated=False)
    R9=add(scale(a,RH),scale(b,Rfinite));FR=prod(F,R9);AR=prod(A,R9);V=prod(PiT,AR)
    I9,I89,I98,I319=identity(9),identity(89),identity(98),identity(319)
    ordered=n('vertical_block',[98,319],[F,prod(PiT,A)]);lower=n('unit_lower_block_frame',[98,98],[V],blocks='[[I9,0],[-V,I89]]',determinant=1);upper=n('unit_lower_block_frame',[98,98],[V],blocks='[[I9,0],[V,I89]]',determinant=1)
    targetP=add(I98,neg(prod(AR,Pi)));residual=prod(prod(PiT,targetP),A);sourceP=add(I319,neg(prod(R9,F)));rowresidual=tr(residual)
    # Formal dimension assertions, no scalar source/coefficient evaluation.
    assert nodes[RH]['shape']==nodes[R9]['shape']==[319,9]and nodes[FR]['shape']==[9,9]and nodes[residual]['shape']==[89,319]
    operator={'coefficient_ring':'T100=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))]','explicit_new_parameter_polynomial_or_RF_terms':0,'typed_coefficient_operators_executed':False,'canonical_QQ_Bezout_CRT_norm_and_exact_division_are_unexpanded':True,'nodes':nodes,'exports':dict(Gamma319x98=Gamma,A98x319=A,Pi9x98=Pi,E_C=E,E_T=ET,F9x319=F,R_H=RH,G=G,norm15=norm,b=b,a=a,R9=R9,F_R9=FR,target_projection=targetP,ordered_target_matrix=ordered,ordered_unipotent_target_change=lower,ordered_unipotent_target_inverse=upper,source_projection=sourceP,residual89x319=residual,residual319x89=rowresidual),'context':dict(C9_oldlocal98=C,T89_oldlocal98=T,C9_original122=[ctx['pure98_target122_labels'][i]for i in C],all319_sources_and_pure_tails_retained=True,raw_target_frame='[A*R9,E_T] = [E_C,E_T]*[[I,0],[V,I]]; determinant is the known target permutation sign. Ordered frame determinant1.',Gamma_formula=ctx['general_full_source_formula'],old_five_bindings=oldbindings,finite_seven_bindings=finitebindings),'identities':['F*R_H=H*I9','A*Rfinite=G*E_C','b*G=1 moduloH; a*H+b*G=1 inT100','F*R9=I9','target_projection=I98-A*R9*Pi isidempotent withimagekerPi','source_projection=I319-R9*F isidempotent; A*source_projection=target_projection*A','ordered_target_change*ordered_target_inverse=I98','residual=PiT*(I98-A*R9*Pi)*A; residual*R9=0','cokerA ≅coker(residual89x319), inverse induced by E_T; all319sourcegenerators retained'],'not_claimed':['Expanded R9,QQ Bezout,CRT integer coefficients,norm15 or89-matrix entries','A chosen free basis of kerF ofrank310','Full89 surjectivity/full239 onwholeP100','Any scalar coefficient,NF,rank or newpoint computation']}
    dump(out/'operator.json',operator)
    proof={'global_projection9_surjective':True,'source_bound_typed_global_rightinverse_constructor':True,'samecoker_pure89_with_all319_sources':True,'orientation':'A=transpose of accepted row319x98 Gamma; all maps act on column coefficient vectors.','global_surjectivity_proof':['Old five9-minors all use fixedC9, and their unit-normalized polynomials have exactgcdH. Thus F surjects offH.','OnH the accepted fullGamma98 proof implies its fixed9projection surjects.','For an explicit finite constructor use typed polynomial extendedgcd of the saved coprime Qi, fiveadj9 source lifts R_H, seven actualglobal selectedmatrix adj98 source lifts Rj, CRT integerlifts, norm15 inverse, exactdivision byH.','The written formula R9=aR_H+bRfinite has FR9=(aH+bG)I9=I9.','Split targets inorderedC,T frame, keep all319redundant sourcegenerators; no sourcekernel freeness assumption.'],'regularity':['Only nonzeroQQ constants and the five certified curve units are inverted globally.','The scalar norm is nonzeroQQ, certified by fullFpCRT; no Dj or H is inverted as a global regular function.','ExactH division is a polynomial divisibility operation after clearing known-unit denominators.'],'uncomputed_coefficients_honestly_recorded':True,'no_new_arithmetic_job':True}
    dump(out/'proof.json',proof)
    note=ROOT/f'research/notes/{ID}.md';assert note.exists();paths=[Path(__file__),note]
    report={'id':ID,'all_structural_and_binding_checks_passed':True,'frozen':True,'global_fixed9_projection_surjective':True,'finite_typed_R9_constructor_saved':True,'samecoker319x89':True,'new_expanded_parameter_coefficients':0,'canonical_coefficient_operators_not_executed':True,'typed_nodes':len(nodes),'input_sha256':pins,'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'no_full239_or_free310_claim':True,'seconds':time.monotonic()-start,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);dump(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','source_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
