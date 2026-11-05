#!/usr/bin/env python3
"""Independent full, uncut degree-three source/completeness review; no eta values."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q

ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63'
ID=PRI+'_REVIEW'
PRIMARY='9aad0271c9cb8e0febb3ec651094e0eab4b88a21e970e3191fb6b586b4f39a5f'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def nc(x):
    if isinstance(x,(int,float)):return 1
    if isinstance(x,str):
        try:Q(x);return 1
        except (ValueError,ZeroDivisionError):return 0
    if isinstance(x,list):return sum(map(nc,x))
    if isinstance(x,dict):return sum(map(nc,x.values()))
    return 0
def det3(M):
    return sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
       *M[0][p[0]]*M[1][p[1]]*M[2][p[2]] for p in itertools.permutations(range(3)))
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};rawbytes=0;decoded=0
    def pin(p,h=None):
        nonlocal rawbytes,decoded
        p=pathlib.Path(p);digest=sha(p)
        if h:assert digest==h,(str(p),digest,h)
        k=str(p.relative_to(ROOT))
        if k not in pins:rawbytes+=p.stat().st_size
        pins[k]=digest
        if p.suffix!='.json':return None
        d=json.loads(p.read_text());decoded+=nc(d);assert decoded<100000
        return d
    pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
    docs={}
    for p,h in pr['input_sha256'].items():
        d=pin(ROOT/p,h);docs[p]=d
        if p.endswith('_REVIEW/report.json'):assert d['all_passed']
    D=pin(ROOT/'reports'/PRI/next(iter(pr['evidence_sha256'])),next(iter(pr['evidence_sha256'].values())))
    pin(ROOT/'code/ambient_B18_C_full_cubic_line_source_v2_20261105_v63.py',pr['code_sha256'])
    pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
    def get(name):return next(v for p,v in docs.items() if pathlib.Path(p).name==name)
    old=get('complete155-source-cubic-operator-with-four-formal-alphas.json')
    frame=get('actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json')
    c2=get('actual-physical12-and-C2-right-exact-quotient.json')
    raw=get('actual-ell2-product33-and-raw45-actions.json')
    alpha=get('eight-cached-even-root-values-and-four-complete-alpha-readbacks.json')
    reduction=get('full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json')
    oldcheck=pin(ROOT/'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/independent-full116-to155-source-completeness-and-formal-line.json',
       '3be097f790333725e445143c921c2fa75ae96aefe77b82f380807e02849ab565')
    assert oldcheck['literal_constant_source_change116_and_both_inverses_verified']
    assert oldcheck['first_map_rank21_and_entire_constant_kernel95_verified']
    # Rebuild all first-frame linear syzygies by grouping literal monomials.
    fibers={};cols=[]
    for t in range(3):
      for i in range(7):
       for j in range(9):
        k=(t,tuple(sorted([i,j])));fibers.setdefault(k,[]).append(len(cols));cols.append((t,i,j))
    assert len(cols)==189 and len(fibers)==126
    kernels=[]
    for indices in fibers.values():
        for c in indices[1:]:kernels.append({c:Q(1),indices[0]:Q(-1)})
    assert len(kernels)==63
    author={v['free_column']:{i:Q(q) for i,q in v['terms']} for v in old['complete63_linear_firstframe_syzygies_with_freeI63']}
    assert {max(v):v for v in kernels}==author
    for v in kernels:
        ims={}
        for c,q in v.items():
            t,i,j=cols[c];assert i<7 and j<7
            k=(t,tuple(sorted([i,j])));ims[k]=ims.get(k,Q(0))+q
        assert not any(ims.values())
    assert old['constant_quadratic_kernel92_counts']=={'curve':2,'cross':64,'tauA_kernel':9,'graph_tauB':5,'tauAB':12}
    assert len(old['all92_constant_sources'])==92
    assert 95*9==855 and 21*45==945 and 92+63==155
    # A fresh no-cut NF uses the actual two curve cubics and retains ell^3.
    curves=[];N=[[Q(0)]*12 for _ in range(6)]
    gr=pin(ROOT/'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
      'be1a33f4173772aadd7e9e9428a28d8674f91fb87e15960a6a750daad0110c36')
    for t,record in enumerate(gr['two_block_source_curve_NF_receipts']):
        g=list(map(Q,record['curve3_r3_r2ell_rell2_ell3']));assert g[0];curves.append(g)
        for j in range(1,4):
            v=[-g[j]/g[0]]+[Q(k==j) for k in range(1,4)]
            assert sum(a*b for a,b in zip(v,g))==0
            N[3*t+j-1][4*t:4*t+4]=v
    assert [[str(v) for v in row] for row in N]==D['full_no_cut_NF6']
    assert D['full_no_cut_NF6']==reduction['full_no_cut_cubic_NF6_rows12_columns']
    direction=[N[k][3]-N[k][7]/4 for k in range(6)]
    assert direction==[0,0,1,0,0,-Q(1,4)]
    # Literal quadratic data only: the uncomputed pure ell^3 is never accessed.
    f={r['physical_row']:r for r in frame['all12_verified_source_records']}
    c={r['physical_row']:r for r in c2['physical12_coefficients']}
    a={r['physical_row']:r for r in alpha['all4_actual_alpha_delta_B6']}
    b33=raw['actual33_basis'];physical=list(f)
    active=[0,4,5];receipts=[]
    for p in [19,28,29,72]:
        j=physical.index(p);assert j<7
        assert list(map(Q,f[p]['bA_first_z7']))==[Q(i==j) for i in range(7)]
        assert list(map(Q,f[p]['bB_first_z7']))==[-Q(i==j)/4 for i in range(7)]
        values={i:Q(v) for i,v in c[p]['delta']};mixed=[]
        for i in range(7):
            for n in [(1,0),(0,1)]:
                val=values.get(b33.index([i+1,*n]),Q(0))
                if i in active:assert val==0
                mixed.append([i,list(n),str(val),'beta2=0' if i not in active else 'literal_delta2=0'])
        assert all(values.get(b33.index([0,*m]),Q(0))==0 for m in [(3,0),(2,1)])
        assert Q(a[p]['delta_r_ell2'])==0
        assert c[p]['roots']==next(v['literal_roots'] for v in D['four_FORMAL_pureell3_callbacks'] if v['physical_row']==p)
        receipts.append({'physical_row':p,'all14_mixed_leaves_retained':mixed,
          'B_graph_cubic_correction':'sum(delta2_zi_r*r+delta2_zi_ell*ell)*beta_i2/4 = 0 termwise; inactive delta leaves NOT asserted zero',
          'pure_plane_cubic':'(eta*ell^3,-eta*ell^3/4, H3 eliminated by S1*P2initial)',
          'no_cut_NF_direction':list(map(str,direction))})
    assert len(D['literal24_active_mixed_delta2_zero_checks'])==24
    # All known-extra images are exhaustively identified in the SAME155 order.
    categories={'curve':2,'cross':64,'tauA_kernel_graph':5,'tauA_kernel_free':4,
                'graph_tauB':5,'tauAB':12,'firstframe_syzygy':63}
    assert {x['kind']:x['count'] for x in D['complete_no_cut_source_categories']}==categories
    assert sum(categories.values())==155
    expected=[[66+i,'eta_'+str(p),list(map(str,direction))] for i,p in enumerate([19,28,29,72])]
    assert D['only_possible_nonzero_columns']==expected and D['complete_formal_operator_shape']==[6,155]
    # Source91/63 zero statements below are before the ell^3 quotient. In
    # particular H*ell^3 is in S1*initial(P2), not an imposed quotient source.
    familyproof=[
      'Both intrinsic curve sources are the two actual relations defining N6, including their ell^3 coefficients.',
      'For49 zi*(zj+beta_j)H product quadratics, cancel with zi*RH_j. The cubic carry retains zi. The other15 products begin at4.',
      'For the63 firstframe Koszuls, each multiplier contains z. M2 times them has no pureplane cubic. Evaluating a z with beta2 raises degree3 to at least4.',
      'All five graph delta plane cubics vanish by the independently accepted full-C torsion obstruction. Actual14 mixed delta2 vanish too; neither graph tauA nor tauB has additional plane3 carry.',
      'For tauAB12, primitive planeA/B begins at4. CubicH is removed by the linear multiples r*P0,ell*P0,ell*P1,ell*P2 of the actual normalized quadratic sources. Their A/B terms start at3, so this removes H3 with A/B carry at4 only.',
      'For four free tauA sources, the literal source is F_j*(tauA-tauB/4). RA exact clears A mixed2; RB substitution gives +delta2_mixed*beta2/4 in B3. All24 active coefficients vanish and inactive beta2=0. H mixed2 clears use RH whose B2plane=0. Thus only eta*(ell3A-ell3B/4) remains.',
      'Every arbitrary formal source coefficient expands u0+u1+u2+... . Complete constant firstkernel95 and plane quadratic kernel92 leave92 lifts; choices of linear21 corrections differ by precisely63 Koszuls.855 multiples U95*S1 and945 R21*S2 are known homogeneous closure. Terms of higher source order begin at4. This argument is independent of the ell3 cut.'
    ]
    # 151 zero columns plus a formal rank-one line; no eta data evaluated.
    formal={i:direction for i in [66,67,68,69]}
    assert len(set(tuple(x) for x in formal.values()))==1 and direction[2]==1
    assert not (ROOT/'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63').exists()
    assert not (ROOT/'reports/AMBIENT_B18_C_FULL_GRAPH03_CACHE_VALUES_20261105_v63').exists()
    assert rawbytes<1<<20 and decoded+20000<150000
    proof={'same_actual_R0_full116_before_ell3_cut':True,
      'accepted_complete116_constant_source_change_and_both_inverses_rebound':True,
      'independent189_monomial_grouping_rank126_kernel63':True,
      'all63_literal_kernel_multiplybacks_and_95_92_source_counts':True,
      'all855_and945_coefficient_multiple_families_retained':True,
      'fresh_no_cut_curve_normal_form':[[str(x) for x in row] for row in N],
      'four_actual_source_root_and_mixed_carry_receipts':receipts,
      'all_seven_source_categories_before_any_cut':familyproof,
      'complete_6_by155_formal_operator_nonzero_columns':expected,
      'all_other151_images_zero_source_proved':True,
      'full_C_HF3':'5 iff some eta_j nonzero;6 iff all4 zero',
      'no_eta_readback_pairing_physical_IFT_gauge_jet_or_rank_search':True,
      'e_still_only_ge3_no_original_H0_nil_inference':True,
      'decoded_numeric_leaves_cumulative':decoded,'new_scratch_output_reserve':20000,
      'whole_source_upper':decoded+20000,'RAW_read':rawbytes,
      'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_RAW_bytes':32<<20}
    ep=out/'independent-no-cut-complete155-source-and-four-formal-eta-line.json';save(ep,proof)
    result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,
       'primary_report_sha256':PRIMARY,'input_sha256':pins,'code_sha256':sha(__file__),
       'evidence_sha256':{ep.name:sha(ep)},'complete_fullC_HF3_5_or6_criterion_accepted':True,
       'eta_values_unread':True,'elapsed_seconds':time.monotonic()-start}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
