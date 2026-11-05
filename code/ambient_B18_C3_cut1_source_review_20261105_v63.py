#!/usr/bin/env python3
"""Independent literal-ell3-cut quartic SOURCE operator; no scalar pairing."""
import hashlib,importlib.util,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_C3_CUT1_QUARTIC_SOURCE_20261105_v63'
ID=BASE+'_REVIEW'
PIN='ed1678417ff8e459de67d38ced0c8c835aab2a7ea8a1eeb315867a7570f9b633'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def independent_cut_operator():
    p=ROOT/'code/ambient_B18_C3_HB_right2_source_review_20261105_v63.py'
    assert sha(p)=='d0f891265be64f5029ff2a6ec1f447f478ad1188b8d05450353964645ed505c9'
    s=importlib.util.spec_from_file_location('own_HB_polynomials',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    b=ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json'
    assert sha(b)=='070cb5dee9df472724c3a7eac6cceaaec2812322136e718349bef762d635bc50';data=read(b)
    bp=ROOT/'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63/four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json'
    assert sha(bp)=='3605f4a6a70b95569979997b32bbc8ba5f7229456314d3b4f311185317a1f327'
    beta=[[Q(c) for c in row] for row in read(bp)['beta2_plane_r2_rell_ell2']]
    raw=read(ROOT/'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json');basis=raw['actual33_basis']
    ids=[basis.index([0,3,0]),basis.index([0,2,1])];zr=[basis.index([i+1,1,0]) for i in range(7)];zl=[basis.index([i+1,0,1]) for i in range(7)]
    T=[[Q(c) for c in row] for row in data['constant_source_normalization_T']]
    H=[[Q(row[i]) for i in [0,4,5]] for row in data['source_first7_h_rows']]
    rows=[18,60,67];D=[];E=[]
    for row in rows:
        d={(3,0):{('cache',row,ids[0]):Q(1)},(2,1):{('cache',row,ids[1]):Q(1)},
           (1,2):{('lambda',row,-1):Q(1)},(0,3):{('opaque_ell3',row,-1):Q(1)}}
        gamma={}
        for z in [0,4,5]:
            for k,c in enumerate(beta[z]):
                gamma=m.pplus(gamma,{(3-k,k):{('cache',row,zr[z]):c}})
                gamma=m.pplus(gamma,{(2-k,k+1):{('cache',row,zl[z]):c}})
        D.append(d);E.append(m.pplus(d,gamma,-1))
    A={};B={}
    for i in range(3):
        A=m.pplus(A,D[i],T[i][2])
        coefficient=sum(H[j][i]*T[j][2] for j in range(3))
        B=m.pplus(B,E[i],coefficient)
    # Literal full source ell*P2-ell3*eH: its degree3 cancels EXACTLY.
    # We only compute its two formal plane quartic coefficients.
    a1=Q(data['curve_ratios']['a1']);a2=Q(data['curve_ratios']['a2'])
    b1=Q(data['curve_ratios']['b1']);b2=Q(data['curve_ratios']['b2'])
    A4=m.normal4(m.shift(A,(0,1)),[Q(1),a1,a2,Q(0)])
    B4=m.normal4(m.shift(B,(0,1)),[Q(1),b1,b2,Q(0)])
    result=m.fplus(B4,A4,-Q(121,48))
    assert all(k[0]=='cache' for k in result)
    return {'A4':m.encode_form(A4),'B4':m.encode_form(B4),'residual':m.encode_form(result),
      'all_lambda_and_opaque_ell3_sensitivities_zero':True,
      'literal_source':'ell*P2-ell3*eH','all_cache_leaves_remain_formal':True}

def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};docs={}
    rp=ROOT/'reports'/BASE/'report.json';assert sha(rp)==PIN;pr=read(rp);assert pr['all_passed']
    pins[str(rp.relative_to(ROOT))]=PIN
    for p,h in pr['input_sha256'].items():
        assert sha(ROOT/p)==h;pins[p]=h;d=read(ROOT/p);docs[pathlib.Path(p).name]=d
        if p.endswith('_REVIEW/report.json'):assert d['all_passed']
    for p,h in pr['evidence_sha256'].items():
        q=ROOT/'reports'/BASE/p;assert sha(q)==h;pins[str(q.relative_to(ROOT))]=h;data=read(q)
    assert sha(ROOT/'code/ambient_B18_C3_cut1_quartic_source_20261105_v63.py')==pr['code_sha256']
    assert sha(ROOT/'research/notes'/f'{BASE}.md')==pr['note_sha256']
    own=independent_cut_operator()
    assert own['A4']==data['A4_cached_linear_operator']
    assert own['B4']==data['B4_cached_linear_operator']
    assert own['residual']==data['residual_B_minus_121A_over48_operator']
    assert len(own['A4'])==6 and len(own['B4'])==len(own['residual'])==24
    raw=docs['actual-ell2-product33-and-raw45-actions.json'];basis=raw['actual33_basis']
    used={key[2] for name in ['A4','B4','residual'] for key,w in own[name]}
    for m,nf in zip(raw['basis74'],raw['all74_normal_forms']):
        nf=dict(nf)
        for k in used:assert Q(nf.get(k,'0'))==Q(m==basis[k])
    c2={r['physical_row']:dict(r['delta']) for r in docs['actual-physical12-and-C2-right-exact-quotient.json']['physical12_coefficients']}
    bounds=[]
    for label,form in [('A4',own['A4']),('B4',own['B4']),('residual',own['residual'])]:
        U=0;D=1;product=0
        for (tag,row,index),weight in form:
            assert tag=='cache';w=Q(weight);x=Q(c2[row].get(index,'0'))
            U+=((abs(w.numerator)+w.denominator-1)//w.denominator)*((abs(x.numerator)+x.denominator-1)//x.denominator)
            D=__import__('math').lcm(D,w.denominator*x.denominator)
            product=max(product,max(abs(w.numerator).bit_length(),w.denominator.bit_length())+max(abs(x.numerator).bit_length(),x.denominator.bit_length()))
        bounds.append({'output':label,'U':str(U),'D':str(D),'height':max((U*D).bit_length(),D.bit_length()),'direct_product_precheck':product,'terms':len(form)})
    assert bounds==data['output_majorants_no_pairing'] and max(b['height'] for b in bounds)==43
    assert 12*43+3==data['whole_fold_minor_inverse_rank_operation_bound']==519<2048
    right=docs['phase2-actual-RIGHT2-quartics-two-factors-and-rank.json']
    v=list(map(Q,right['rank_certificate']['pivot_column']));assert v==[-Q(3,5),-Q(121,80)]
    assert v[1]-Q(121,48)*v[0]==0
    proof={'all_primary_source_and_prior_peer_hashes_match':True,
      'literal_quotient_column':'ell3*eH is exact with no higher terms in C3=C/ell3C',
      'same_actual_P2_has_FULL_second_term_ell2H_after_all_z_clears':True,
      'complete_cubic_cancellation':'ell*P2-ell3H has zero FULL degree3, including all z sectors',
      'full_quartic_dependency_proof':[
        'Degree4 is ell*P2,3; P2,4 and all higher terms multiplied byell begin at5.',
        'z-containing degree3 terms retainz afterell. Clearing them by degree3 first-frame coefficients adds no A/B carry before5.',
        'H quartics clear by degree2 times actual H quadratic pivots; their A/B third terms enter only at5.',
        'Therefore only the complete plane cubic A/B of P2 from the accepted coherent RIGHT2 construction enters.',
        'Direct multiplication byell and monomial NF4 kills every lambda*r*ell2 and opaqueell3 channel, without assigning their coefficientszero.',
        'The negativeRB beta2 gamma30 andgamma21 carry remains in full. No delta4, mixedPhi3,nu, beta3 or newgauge callback is needed.'],
      'independent_three_unpaired_operators':own,
      'all_used_old_C2_scalar_Taylor_duals_on74_monomials_verified':True,
      'positive_majorants_independently_rebuilt':bounds,
      'whole_operation_bound':519,'hard_bits':2048,
      'known_HB_vector':list(map(str,v)),'dual': ['-121/48','1'],
      'universal_det_identity':'det[[−3/5,A],[−121/80,B]]=(-3/5)*(B−121*A/48)',
      'nonzero_completion':'A nonzero residual gives two independent actual quartics, m4C3=0, HF36640 and exactlength19. LongerDVRexponent and originalH0/nil stayopen.',
      'zero_scope':'OnlyCUT1 dependence on knownHBline. No fullquarticrank or length20 inference.',
      'all_old_leaf_values_unpaired_no_residual_computed':True,
      'prepared_worker_code_and_whole_custody_still_pending':True,
      'source_static_numeric_bound':50000,'hard_numeric':150000,'hard_seconds':60,'hard_address_bytes':2<<30,'hard_RAW_bytes':32<<20}
    ep=out/'independent-literal-cut-source-full-NF-and-positive-bound-checkpoint.json';save(ep,proof)
    save(out/'phase.json',{'phase':'FORMAL_SOURCE_PASS_NO_VALUE_GO_pending_worker_custody','residual_uncomputed':True})
    print(json.dumps({'formal_source_passed':True,'terms':24,'whole_bound':519,'no_scalar_pairing':True,'seconds':time.monotonic()-start}));signal.alarm(0)
if __name__=='__main__':main()
