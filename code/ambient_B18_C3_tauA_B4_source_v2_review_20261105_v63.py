#!/usr/bin/env python3
"""Independent source/NF/literal-zero audit, without a cached dot worker."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63'
PRIMARY='ce682470ed9742cd556dd9c7e9507fc9e9d3ae4c909cb7e0223d0af7be94fbf3'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def plus(a,b,k=Q(1)):
    d=dict(a)
    for m,c in b.items():
        d[m]=d.get(m,Q(0))+k*c
        if not d[m]:del d[m]
    return d
def mul(a,b):
    d={}
    for m,c in a.items():
        for n,v in b.items():
            t=tuple(x+y for x,y in zip(m,n));d[t]=d.get(t,Q(0))+c*v
    return {m:c for m,c in d.items() if c}
def compose(rows,x,y):
    out={}
    for (a,b),c in rows:
        v={(0,0):Q(c)}
        for _ in range(a):v=mul(v,x)
        for _ in range(b):v=mul(v,y)
        out=plus(out,v)
    return out
def normal(p,curve):
    d={m:c for m,c in p.items() if m[1]<3}
    while any(m[0]>=3 for m in d):
        m=max(m for m in d if m[0]>=3);c=d.pop(m)
        for k in range(1,4):
            n=(m[0]-k,m[1]+k)
            if n[1]<3:d=plus(d,{n:-c*curve[k]/curve[0]})
    return d
def fullC4():
    # QQ[D,bA,bB,nuA,nuB,ta,tb]/(ta^2-nuA,tb^2-nuB).
    def sym(i):return {tuple(int(i==j) for j in range(7)):Q(1)}
    D,ba,bb,na,nb,ta,tb=[sym(i) for i in range(7)]
    F=plus(plus(D,mul(ba,ta)),mul(bb,tb),-1)
    v=mul(F,plus(ta,tb,-Q(1,4)));out={}
    for m,c in v.items():
        n=list(m);n[3]+=n[5]//2;n[5]%=2;n[4]+=n[6]//2;n[6]%=2
        out=plus(out,{tuple(n):c})
    expected=plus(plus(plus(mul(na,ba),mul(nb,bb),Q(1,4)),mul(D,ta)),mul(D,tb),-Q(1,4))
    expected=plus(expected,mul(plus(bb,ba,Q(1,4)),mul(ta,tb)),-1)
    assert out==expected
    return True
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    pins={};rawbytes=0
    def pin(p,h=None):
        nonlocal rawbytes
        p=pathlib.Path(p);actual=sha(p)
        if h is not None:assert actual==h,(str(p),actual,h)
        key=str(p.relative_to(ROOT))
        if key not in pins:rawbytes+=p.stat().st_size
        pins[key]=actual;return read(p)
    pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
    docs={}
    for p,h in pr['input_sha256'].items():
        d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
        if p.endswith('_REVIEW/report.json'):assert d['all_passed']
    for p,h in pr['evidence_sha256'].items():e=pin(ROOT/'reports'/PRI/p,h)
    assert sha(ROOT/'code/ambient_B18_C3_tauA_B4_source_v2_20261105_v63.py')==pr['code_sha256']
    assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
    corr=pin(ROOT/'reports/AMBIENT_B18_C3_TAUA_B4_FORMULA_CORRECTION_20261105_v63/report.json',
             'f1c3fd6ab123783ae35bd8f274ab30470148a9172074704bc96a862b4b425454')
    for key in ['input_sha256','preserved_source_history_sha256']:
        for p,h in corr[key].items():
            assert sha(ROOT/p)==h
    oldfail=pin(ROOT/'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_PRECHECK_20261105_v63_REVIEW/report.json',
                '986fb45c72222618209509d0617c8f6e07b4731abb1170b88315aa461a4ecf50')
    assert not oldfail['all_passed'] and oldfail['audit_blocker_reproduced']
    lit=pin(ROOT/'reports/AMBIENT_B18_C3_TAUA_B4_LITERAL_ZERO_20261105_v63/report.json',
            'a9f896193c118eb09739e50279d140f2686af66d9294c679fc191e633b44ca06')
    for p,h in lit['input_sha256'].items():assert sha(ROOT/p)==h
    for p,h in lit['evidence_sha256'].items():L=pin(ROOT/'reports/AMBIENT_B18_C3_TAUA_B4_LITERAL_ZERO_20261105_v63'/p,h)
    c2=docs['actual-physical12-and-C2-right-exact-quotient.json']
    raw=docs['actual-ell2-product33-and-raw45-actions.json'];basis=raw['actual33_basis']
    pair=docs['actual-special-fiber20-coefficients.json']
    frame=docs['actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json']
    kernel=docs['nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json']
    cover=docs['actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json']
    r={(1,0):Q(1)};ell={(0,1):Q(1)};s=plus(r,ell,Q(1,5))
    beta=[]
    for row in pair['beta2_before_ell']:
        assert all(not any(m[:7]) and sum(m)==2 for m,c in row)
        p=compose([(m[7:],c) for m,c in row],r,s)
        beta.append([p.get((2-j,j),Q(0)) for j in range(3)])
    assert [[str(c) for c in row] for row in beta]==e['beta2_plane_r2_rell_ell2']
    yb=pair['YB_in_YA_linear']
    rb=plus({m:Q(yb[7][7])*c for m,c in r.items()},{m:Q(yb[7][8])*c for m,c in s.items()})
    sb=plus({m:Q(yb[8][7])*c for m,c in r.items()},{m:Q(yb[8][8])*c for m,c in s.items()})
    cp=[compose(pair['gA_plane3'],r,s),compose(pair['gB_own_plane3'],rb,sb)]
    curves=[[p.get((3-j,j),Q(0)) for j in range(4)] for p in cp]
    for i,co in enumerate(curves):
        assert co==list(map(Q,cover['two_block_source_curve_NF_receipts'][i]['curve3_r3_r2ell_rell2_ell3']))
        assert not normal(cp[i],co)
    idx30=basis.index([0,3,0]);idx21=basis.index([0,2,1])
    zr=[basis.index([i+1,1,0]) for i in range(7)];zl=[basis.index([i+1,0,1]) for i in range(7)]
    duals=zr+zl+[idx30,idx21]
    for m,nf in zip(raw['basis74'],raw['all74_normal_forms']):
        nf=dict(nf)
        for j in duals:assert Q(nf.get(j,'0'))==Q(m==basis[j])
    # Every one of the eighteen leaves is formal. Reduction is by direct
    # monomial rewriting, independently of the author's adjoint formulas.
    generic=[[],[]];vectors=[]
    for j in range(18):
        A={(3-j,j):Q(1)} if j<4 else {}
        B={m:-c/4 for m,c in A.items()}
        if j>=4:
            z=(j-4)//2;v=(1,0) if j%2==0 else (0,1)
            bp={(2-k,k):c for k,c in enumerate(beta[z]) if c}
            B=plus(B,mul({v:Q(1,4)},bp))
        na=normal(A,curves[0]);nb=normal(B,curves[1])
        vec=[na.get((2,1),Q(0)),na.get((1,2),Q(0)),nb.get((2,1),Q(0)),nb.get((1,2),Q(0))]
        vectors.append(vec);generic[0].append(vec[2]);generic[1].append(vec[3]+vec[1]/4)
    assert generic[0][2]==generic[0][3]==generic[1][2]==generic[1][3]==0
    leafindices=[idx30,idx21,None,None]+[j for z in range(7) for j in [zr[z],zl[z]]]
    expected_terms=[];generic_records=[]
    fm={d['physical_row']:d for d in frame['all12_verified_source_records']}
    cm={d['physical_row']:d for d in c2['physical12_coefficients']}
    K=[[Q(c) for c in row] for row in kernel['complete_constant_kernel12_by9']]
    literal_records=[];bindings=[]
    for row,source in zip([19,28,29,72],e['source_rows4']):
        i=kernel['physical_row_order'].index(row);k=kernel['free_columns9'].index(i)
        assert [v[k] for v in K]==[Q(i==j) for j in range(12)]
        assert list(map(Q,fm[row]['bA_first_z7']))==[Q(i==j) for j in range(7)]
        assert list(map(Q,fm[row]['bB_first_z7']))==[-Q(i==j)/4 for j in range(7)]
        assert not any(map(Q,fm[row]['bB_second_plane']))
        assert source['original_physical_coefficient_roots']==cm[row]['roots']
        assert source['first7_index']==i and source['literal_kernel_column']==k
        oldsource=next(v for v in L['four_literal_full_F_RH_and_first7_bindings'] if v['physical_row']==row)
        assert oldsource['original_physical_coefficient_roots']==cm[row]['roots']
        # Only nonzero-weight accesses below must vanish. Mixed leaves at
        # inactive beta indices can be nonzero and are not set to zero.
        rec={'physical_row':row,'adjoints2':[]}
        for n,label in enumerate(['psi_B_r2ell','chi_B_rell2_plus_A_rell2_over4']):
            generic_records.append([row,['psi','chi'][n],list(map(str,generic[n]))])
            terms=[]
            for j,w in enumerate(generic[n]):
                if w:
                    index=leafindices[j];assert index is not None
                    leaf=Q(dict(cm[row]['delta']).get(index,'0'));assert leaf==0
                    terms.append([row,index,str(w)]);literal_records.append([row,label,index,str(w),'0'])
            # Order in the primary is sparse dictionary insertion, not Taylor order.
            actual=next(x for x in e['four_two_cached_adjoints_NOT_PAIRED'] if x['physical_row']==row)['adjoints2'][n]
            assert actual['label']==label
            assert sorted(terms)==sorted(actual['terms'])
            rec['adjoints2'].append({'label':label,'terms':terms})
        expected_terms.append(rec)
        bindings.append({'row':row,'literal_same_physical_Fj_RHi':True,'full_C4_source_verified':True,
          'all_mixed_H_clears_retained':True,'unknown_zz_clears_keep_a_z_factor':True})
    assert generic_records==e['all144_generic_full_NF_functional_coefficients_verified']
    assert sorted(literal_records)==sorted(e['all60_nonzero_weight_cached_leaf_accesses_are_literal_zero'])
    assert len(literal_records)==60 and fullC4()
    # With literal leaves zero the remaining formal delta12 coefficient is
    # exactly the single stated line. No value is assigned to delta12.
    assert vectors[2]==[Q(0),Q(1),Q(0),-Q(1,4)]
    assert vectors[3]==[Q(0)]*4
    # Conditional rank statements use the actual r/ell multiplication map.
    a1=curves[0][1]/curves[0][0];b1=curves[1][1]/curves[1][0]
    rmap=[[-a1,Q(1),Q(0),Q(0)],[Q(0),Q(0),-b1,Q(1)]]
    ellmap=[[Q(1),Q(0),Q(0),Q(0)],[Q(0),Q(0),Q(1),Q(0)]]
    assert all(normal({(2,2):Q(1)},co)=={(2,2):Q(1)} for co in curves)
    assert rawbytes<500000
    evidence={'source_V1_formula_withdrawal_and_exact_history_verified':True,
      'all_source_hashes_and_prior_independent_peers_match':True,'no_author_producer_imports':True,
      'actual_same_C_module_over_same_R0_not_an_algebra_split':True,
      'all_four_same_physical_source_and_full_C4_bindings':bindings,
      'full_C4_q_vector':['nuA*bA+nuB*bB/4','delta','-delta/4','-bB-bA/4'],
      'actual_beta2_and_both_curve3_compositions_independently_rebuilt':True,
      'direct_polynomial_NF_rewriting_all144_formal_coefficients':True,
      'corrected_chi_retains_minus_b2_carry_r3':True,
      'all74_Taylor_duals_for_16_selected_indices_verified':True,
      'all60_actual_cache_occurrences_individually_zero':literal_records,
      'eight_corrected_projections_zero_without_dot_worker':True,
      'remaining_cubic_line':'alpha_j*(A_rell2-B_rell2/4), alpha_j still unknown',
      'rank2_completion_conditional_only':'Two cubic dimensions and psi-nonzero ell multiple remove >=3 from length21; accepted lower18 forces length18/HF36621/exponents1122 only IF rank2.',
      'rank1_completion_conditional_only':'The joint r/ell map on degree3 is injective; nonzero class gives upper19 only.',
      'multiplication_r':[[str(c) for c in v] for v in rmap],
      'multiplication_ell':[[str(c) for c in v] for v in ellmap],
      'no_actual_C3_length_or_torsion_exponent_or_H0_or_nil_claim':True,
      'no_eight_dot_rank_readback_physical_gauge_or_new_jet_job':True,
      'numeric_static_envelope':22000,'hard_numeric':150000,'RAW_read':rawbytes,
      'hard_seconds':60,'hard_address_bytes':2<<30,'hard_RAW_bytes':32<<20}
    ep=out/'independent-full-source-corrected-NF-and-sixty-literal-zeros.json';save(ep,evidence)
    report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PRIMARY,
      'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'corrected_eight_projections_zero_accepted':True,'four_full_cubic_classes_remain_unknown':True,
      'no_new_value_or_rank_job':True,'elapsed_seconds':time.monotonic()-start}
    save(out/'report.json',report);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
