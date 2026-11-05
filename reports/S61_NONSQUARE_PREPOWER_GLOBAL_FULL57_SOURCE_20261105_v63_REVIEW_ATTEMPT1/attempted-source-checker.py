#!/usr/bin/env python3
"""Independent global source/domain audit. No new matrix values or rank job."""
import collections,hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW'
P='S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63'
PRIMARY='f56e3d125e18ff1aa1c3d8720dc574062ff5a54759ec0bf13d17fc17267b7499'

def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def dump(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def field(t,k):
    j=t.index(':',t.index(json.dumps(k)))+1
    while t[j].isspace():j+=1
    return json.JSONDecoder().raw_decode(t,j)[0]
def spans(t,start):
    assert t[start]=='[';i=start+1
    while True:
        while t[i].isspace() or t[i]==',':i+=1
        if t[i]==']':return
        a=i;dep=0;s=False;escape=False
        while i<len(t):
            c=t[i]
            if s:
                if escape:escape=False
                elif c=='\\':escape=True
                elif c=='"':s=False
            elif c=='"':s=True
            elif c in '[{':dep+=1
            elif c in ']}':
                if dep==0:break
                dep-=1
                if dep==0:i+=1;break
            elif c==',' and dep==0:break
            i+=1
        yield a,i
def poly(x):return {tuple(m):Q(c) for m,c in x if Q(c)}
def add(a,b,k=Q(1)):
    d=dict(a)
    for m,c in b.items():
        d[m]=d.get(m,Q(0))+k*c
        if not d[m]:del d[m]
    return d
def mul(a,b):
    d={}
    for m,c in a.items():
        for n,e in b.items():
            k=tuple(x+y for x,y in zip(m,n));d[k]=d.get(k,Q(0))+c*e
    return {m:c for m,c in d.items() if c}
def mon(n,i):return {tuple(int(j==i) for j in range(n)):Q(1)}
def constant(n,c):return {(0,)*n:Q(c)}
def count(x):
    if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
    if isinstance(x,(list,tuple)):return sum(map(count,x))
    return int(isinstance(x,(int,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())

def main():
    began=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        primarypath=ROOT/'reports'/P/'report.json';assert sha(primarypath)==PRIMARY
        pr=read(primarypath);assert pr['all_passed'] and pr['source_only']
        pins={str(primarypath.relative_to(ROOT)):PRIMARY}
        for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        for p,h in pr['evidence_sha256'].items():
            pp=ROOT/'reports'/P/p;assert sha(pp)==h;pins[str(pp.relative_to(ROOT))]=h
        inv=read(ROOT/'reports'/P/'unit-inventory.json');recipe=read(ROOT/'reports'/P/'global-operator-recipe.json')
        def get(st,filename='report.json'):
            p=ROOT/'reports'/st/filename;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
        def accepted(st,flags=()):
            a=get(st);assert a['all_passed']
            for f in flags:assert a[f] is True,(st,f)
            return a
        def verify_pins(a):
            # Existing exact proofs are reused with literal content hashes, not rerun.
            for group in ['input_sha256','checked_dependency_sha256','evidence_sha256']:
                for p,h in a.get(group,{}).items():
                    if group=='evidence_sha256':p='reports/'+a['id']+'/'+p
                    assert sha(ROOT/p)==h,(group,p)
                    pins[p]=h
        product=accepted('S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63',
          ['actual_original42_by84_product_rank_exact41_on_whole_FULLU'])
        pair=accepted('S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW',
          ['raw24_module_surjective_on_FULLU_via_accepted_gcd_binary_finiteSigma_reviews'])
        bridge=accepted('S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW',
          ['pointwise_RAW24_rank24_implies_original_product_rank41'])
        for a in [product,pair,bridge]:verify_pins(a)
        assert 'all geometric points' in product['field']
        assert product['FULLU_original']=='a*d*F*B*C*J!=0, J=5*a*d-6*e-4'
        basis=accepted('S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW',
          ['exact_all6_QQ_relations_all91_columns','full_localized_row_kernel_generated_by_six_relations'])
        verify_pins(basis)
        bdata=read(ROOT/'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json')
        last=bdata['relations'][5];assert last['deleted_original_row']==14
        assert all(c['original_row']<42 for c in last['coefficients'])
        assert next(c['terms'] for c in last['coefficients'] if c['original_row']==14)==[['12/5',1,0,0]]
        assert bdata['primary_triangular_minor']=='-12*a^3/5'
        assert bdata['primary_deleted_order']==[48,55,59,57,61,14]
        # Read the literal projection identities, not the large source matrices.
        lf='S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63'
        lft=(ROOT/'reports'/lf/'certificate.json').read_text()
        assert field(lft,'inverse_denominator_only_known_unit')=='-40*C'
        assert field(lft,'actual_source_normalizations')==['q1=L^2/8100','q3=L*(e*L-400*a*x0)/8100']
        assert field(lft,'normalization_operation')=='L^2=8100*q1; L*x0=(8100*e*q1-8100*q3)/(400*a)'
        assert len(field(lft,'nine_monomial_pivot_columns'))==9 and len(field(lft,'remaining26_order'))==26
        native=(ROOT/'reports'/lf/'native.json').read_text();cancels=field(native,'genuine_C_factors_cancelled')
        assert field(native,'actual_source_q1_q3_identities_passed') is True
        assert field(native,'q5_C_factors_cancelled')==2
        del native,lft
        # Named-unit polynomials, independent small multiplication checks only.
        a=mon(3,0);d=mon(3,2);e=mon(3,1);one=constant(3,1)
        fac=[a,d,poly([[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']]),
             poly([[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]),
             poly([[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']]),
             poly([[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']])]
        def unit(c,powers):
            z=constant(3,c)
            for f,n in zip(fac,powers):
                assert n>=0
                for _ in range(n):z=mul(z,f)
            return z
        pairfile=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63/certificate.json'
        hpfile=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63/certificate.json'
        pt=pairfile.read_text();eta=poly(field(pt,'known_unit_eta'));del pt
        ht=hpfile.read_text();delta=poly(field(ht,'known_unit_delta'));del ht
        assert eta==unit(Q(192,125),[1,0,1,2,3,0]) and delta==unit(-64000,[0,0,0,0,2,0])
        fspeer=accepted('S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63_REVIEW',
          ['actual_full57_original63_raw91_source_lifts_and_eta_cleanup_accepted'])
        verify_pins(fspeer)
        full=get('S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63','operator.json')
        assert full['nodes'][169]['args']==[164,168] and full['nodes'][169]['shape']==[91,32]
        assert full['exports']['lift25_to_full57']==49
        assert full['nodes'][168]['entries']==[[5,0,167],[16,0,165]]
        assert full['nodes'][167]['op']=='negative' and full['nodes'][167]['args']==[166]
        fullraw=[5,16];del full
        augfile=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json'
        at=augfile.read_text();mu=field(at,'genuine16_source_clearing_units');labels25=field(at,'source25_labels')
        assert field(at,'duplicate_zero_in_both_pi_and_psi') is True
        assert field(at,'removed_genuine_relative_row')==6 and field(at,'kept_genuine_relative_row')==1
        assert field(at,'rawHZ_HW_indices_in_frozen91')==fullraw
        for i in range(16):assert poly(mu[i])==unit(-64000,[0,0,0,0,3-cancels[i],0])
        assert mu[1]==mu[6] and poly(mu[1])==delta;assert len(labels25)==25;del at,mu
        cpfile=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json'
        ct=cpfile.read_text();assert field(ct,'completed_pivots')==8
        remrow=field(ct,'remaining_rows');remcol=field(ct,'remaining_columns')
        stepstart=ct.index('[',ct.index('"all_completed_unit_steps"'));step_spans=list(spans(ct,stepstart))
        expected=[(8,20,'-5/4',[0]*6),(1,10,'-5/4',[0]*6),(2,16,'5',[1,0,0,0,0,0]),
          (3,12,'-5/4',[0]*6),(0,0,'25/4',[1,0,0,0,0,0]),(6,23,'-1/6',[1,0,0,0,0,0]),
          (7,24,'-1/6',[1,0,0,0,0,0]),(11,26,'2/15625',[0,0,1,2,2,0])]
        assert len(step_spans)==8
        tiny=[]
        for k,((lo,hi),(r,c,v,powr)) in enumerate(zip(step_spans,expected)):
            st=ct[lo:hi];pivot=field(st,'pivot');ku=field(st,'known_unit')
            assert field(st,'pivot_row')==r and field(st,'pivot_column')==c
            assert Q(ku['constant'])==Q(v) and ku['a_d_F_B_C_J_powers']==powr
            assert poly(pivot[0])==mul(poly(pivot[1]),unit(Q(v),powr))
            assert inv['eight_actual_pivots'][k]['unit_powers']==powr
            assert Q(inv['eight_actual_pivots'][k]['inverse_constant'])==1/Q(v)
            assert inv['eight_actual_pivots'][k]['inverse_powers']==[-i for i in powr]
            tiny.append({'step':k,'row':r,'column':c,'constant':v,'powers':powr,'exact_named_unit_identity':True})
        mstart=ct.index('[',ct.index('"matrix25x32"'));rr=list(spans(ct,mstart));assert len(rr)==25
        cols={r:list(spans(ct,rr[r][0])) for r,*_ in expected}
        # Only the 36 diagonal/lower-triangle cells are decoded, not matrix800 or its tails.
        for i,(r,c,v,powr) in enumerate(expected):
            for j in range(i+1):
                k=expected[j][1];lo,hi=cols[r][k];x=json.loads(ct[lo:hi])
                assert poly(x[0])=={m:v*int(i==j) for m,v in poly(x[1]).items() if i==j}
        assert set(remrow)==set(range(25))-{r for r,c,v,powr in expected}
        assert set(remcol)==set(range(32))-{c for r,c,v,powr in expected}
        del ct,cols,rr,st,pivot
        cppeer=accepted('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW',
          ['all8_known_unit_fixed_operations_and_all_checkpoint_entries_verified'])
        verify_pins(cppeer)
        n9name='S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63'
        n10name='S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63'
        for st,flags in [(n9name,['actual_one_new_unit_and_same_cokernel_verified']),
                         (n10name,['actual_GL3_unit_completion_and_same_cokernel_verified','all3914_DAG_nodes_and330_remaining_entries_verified'])]:
            verify_pins(accepted(st+'_REVIEW',flags))
        n9=get(n9name,'certificate.json');n10=get(n10name,'certificate.json')
        assert n9['old17_source_original25_rows']==remrow and n9['old24_target_original32_columns']==remcol
        assert n10['old16_source_original25_rows']==n9['new16_source_original25_rows']
        assert n10['old23_target_original32_columns']==n9['new23_target_original32_columns']
        nodes=n10['scalar_DAG']['nodes'];assert len(nodes)==3914 and nodes[:2138]==n9['scalar_DAG']['nodes']
        allowed={'a','s','m','one','z','checkpoint','new_unit_atom','accepted_pair_second'}
        for i,x in enumerate(nodes):
            assert x[0] in allowed
            if x[0] in ('a','s','m'):assert all(0<=j<i for j in x[1:])
        den_count=0
        for num,den in n9['actual_six_dot_and_binding']['transported_pair_second_column']+n10['scalar_DAG']['new_unit_atoms_num_den']:
            assert len(den)==1 and den[0][0][1:]==[0,0] and den[0][0][0] in [0,1,2] and Q(den[0][1]);den_count+=1
        completion_checks=[]
        for cert,skey,tkey,ikey,colskey in [(n9,'source17_roots','target24_roots','target24_inverse_roots','old24_target_original32_columns'),
                                          (n10,'source16_roots','target23_roots','target23_inverse_roots','old23_target_original32_columns')]:
            S=cert[skey];n=len(S)
            assert all(S[i][i]==1 for i in range(n))
            assert all(S[i][j]==0 for i in range(n) for j in range(n) if i!=j and j!=0)
            assert all(S[0][j]==int(j==0) for j in range(n))
            ts=cert[tkey];ti=cert[ikey];nc=len(ts)
            assert all(len(v)==nc for v in ts+ti)
            assert all(0<=q<len(nodes) for v in ts+ti for q in v)
            for h in [27,28,29,30,31]:
                j=cert[colskey].index(h)
                assert all(ts[i][j]==int(i==j) and ti[i][j]==int(i==j) for i in range(nc))
            completion_checks.append({'source_dimension':n,'target_dimension':nc,'source_unipotent_I_plus_u_e0':True,
              'source_inverse_I_minus_u_e0':True,'all_five_heavy_target_columns_literal_identity_in_both_directions':True,
              'saved_target_forward_inverse_identity_reused_from_exact_independent_peer':True})
        assert n9['tail_clear_original_columns']==[1,2,3]
        assert n10['tail_clear_original32_columns']==[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19]
        labels22=n10['new22_target_original32_columns'];roots=n10['remaining15x22_scalar_roots']
        assert len(roots)==15 and all(len(v)==22 and all(0<=r<len(nodes) for r in v) for v in roots)
        assert recipe['all330_actual_matrix_roots']==roots and recipe['new22_labels']==labels22
        assert all(roots[0][labels22.index(j)]==0 for j in [27,28,29,30,31])
        h9=[22,11,13,25,15,19,14,18,17];other=[j for j in labels22 if j not in h9+[27,28,29,30,31]]
        assert other==[1,2,3,4,6,7,8,9] and recipe['other8_Q9_tail_columns']==other
        # Raw-H cleanup identity, a fresh universal polynomial proof, independent of any evaluation.
        z,w,cz,cw,hz,hw=[mon(6,i) for i in range(6)]
        eta0=add(mul(cz,z),mul(cw,w));bn=add(mul(cz,hz),mul(cw,hw));psi=add(mul(z,hw),mul(w,hz),-1)
        assert add(add(mul(eta0,hz),mul(z,bn),-1),mul(cw,psi))=={}
        assert add(add(mul(eta0,hw),mul(w,bn),-1),mul(cz,psi),-1)=={}
        # Source row normalize + shear inverse, with generic p,c in a Laurent ring.
        pp={(1,0):Q(1)};ip={(-1,0):Q(1)};cc={(0,1):Q(1)};oo={(0,0):Q(1)}
        forward=[[ip,{}],[mul(cc,ip),oo]];forward[1][0]={m:-v for m,v in forward[1][0].items()}
        inverse=[[pp,{}],[cc,oo]]
        for X,Y in [(forward,inverse),(inverse,forward)]:
            for i,j in itertools.product(range(2),repeat=2):
                v={}
                for k in range(2):v=add(v,mul(X[i][k],Y[k][j]))
                assert v==(oo if i==j else {})
        notation=accepted('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_NOTATION_20261105_v63')
        verify_pins(notation)
        Tpeer=accepted('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW')
        verify_pins(Tpeer)
        coefpeer=accepted('S61_NONSQUARE_PREPOWER_M12_AXIS_ONE_COEFFICIENT_20261105_v63_REVIEW',
          ['single_a15_coefficient_independently_verified'])
        directpeer=accepted('S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW')
        proof={
          'ALL_FULLU_ring':'QQ[a,e,d,1/(a*d*F*B*C*J)]',
          'product41_proof_dependencies_literal_and_scope_verified':True,
          'product41_is_actual_all_point_image_rank_NOT_abstract_generic_row_rank':True,
          'RAW24_all_point_surjection':pair['id'],'pointwise_17_plus24_bridge':bridge['id'],
          'product_unit_bridge_factors':['80a','-16CJ/243','5a','aC/240','16C^2/81','160aC/81'],
          'upper41_nonzero_syzygy_coefficient':'10a',
          'selected41_span':'actual last relation has only product rows<42 and deleted14 coefficient12a/5',
          'last_product_relation_source_rows':[x['original_row'] for x in last['coefficients']],
          'six_source_relation_upper57':True,'unit_triangular_minor':'-12*a^3/5',
          'phi_definition':'Physical91 -> Z=W=0, invertible L frame, degree3 quotient by (L^2,L*x0); retain26pure+H5; append psi=-w_delta*HZ+z_delta*HW.',
          'phi_products_in_Q10_proof':'80a row identities put q2,q4 in (q1,q3) moduloZW. These lie in (L^2,Lx0); all other products map into q0/scaledq5 times5lights; products have no H, so psi=0.',
          'all16_genuine_clearing_known_units':True,'Cdelta_mu1_equals_mu6':'-64000*C^2',
          'Cdelta_actual_source':'delta*(full57row47-full57row42)',
          'Cdelta_phi_zero_both_pi_and_psi':True,'eta_unit':'192*a*F*B^2*C^3/125',
          'raw7_cleanup_fresh_universal_identities_verified':True,
          'raw7_cleanup_retains_projected_images_and_full_pure_tails':True,
          'eight_named_unit_pivots':tiny,'eight_forward_inverse_universal_check':True,
          'eight_selected_checkpoint_unit_upper_triangular_block_36_cells_verified':True,
          'upper_triangular_inverse':'sum((-strictUpper)^k,k=0..7), no new inversion',
          'completion_checks':completion_checks,'only_QQ_and_a_leaf_denominators_checked':den_count,
          'all3914_DAG_nodes_admissible_topological_and_no_inverse_opcode':True,
          'all330_root_references_and_source_target_frame_verified':True,
          'all8_other_Q9_pure_tail_labels':other,'tail2_is_one_of_eight_not_zero_assumed':True,
          'H3_targets':[28,30,29],'remaining_heavy5':[27,28,29,30,31],
          'global_rank_identity':'rank(old25)=10+rank(current15) over every residue field of FULLU',
          'full57_sufficiency_pointwise_proof':[
            'rankold25=25 makes Q10+retainedgenuine15 independent.',
            'A retainedgenuine15 relation modulo product span maps through phi to a Q10 relation, so all15 coefficients vanish.',
            'Adding Cdelta still leaves only its own coefficient after phi; eta-unit rawH pair forces it zero.',
            'Thus16 independent genuine classes modulo actualproduct41 give image rank57; source6relations give upper57.'],
          'only_named_FULLU_units_in_global_prefix':True,
          'not_inverted':['Delta407','Saux','N6','N7bar','P_K92','h22','p','R434','kappa3','A_anchor'],
          'T_notation':'literal KAPPA_EXTRA2/root239, not Q2 auxiliary900a-126e-121',
          'DIRECT_identity_reused':'A*D15_normalized=-T*D12 before A inversion',
          'corollary_after_this_peer':'On FULLU D(T*D12), RHS unit gives A and D15 unit, hence current15rank15 and actualfull57.',
          'chart_nonempty':'T and D12 are actual nonzero rational functions over the integral domain R, by existing independent T and one-coefficient peers.',
          'oldOmega_and_all_primary_frozen_files_unchanged':True,
          'no_strict_chart_inclusion_badFitting_equality_or_ALLFULLU57_claim':True,
          'no_quartic239_exhaustion_claim':True,
          'no_new_scalar_matrix_source_values_rank_determinant_or_field_computation':True,
          'max_numeric_static_upper':120000,'hard_terms':150000,'hard_RAW_bytes':32*1024*1024,'hard_seconds':60,
          'memory_bytes':2000000000,
        }
        dump(out/'independent-global-unit-source-and-rank-proof.json',proof)
        assert time.monotonic()-began<60
        report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PRIMARY,
          'all_FULLU_product41_source_scope_accepted':True,'global_rank25_equals10_plus_current15_accepted':True,
          'global_current15_fullrank_implies_actual_full57_accepted':True,
          'FULLU_D_T_D12_sufficient_nonempty_chart_accepted':True,
          'not_an_ALLFULLU_rank57_or_bad_Fitting_ideal_equality_claim':True,
          'code_sha256':sha(__file__),'input_sha256':pins,
          'evidence_sha256':{'independent-global-unit-source-and-rank-proof.json':sha(out/'independent-global-unit-source-and-rank-proof.json')},
          'no_author_producer_imports':True,'no_new_scalar_matrix_values_or_rank_job':True,'elapsed_seconds':time.monotonic()-began}
        dump(out/'report.json',report);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
    except Exception as e:
        dump(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_positive_claim':True});raise
if __name__=='__main__':main()
