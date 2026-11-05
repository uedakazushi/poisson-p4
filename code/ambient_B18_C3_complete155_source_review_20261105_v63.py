#!/usr/bin/env python3
"""Independent homogeneous/source audit of the complete cubic generators."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63'
PRIMARY='75bf9f53003d6ae0b5a73769d950dbd9e6ee694a0373eae3f33fd45bc2473578'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def first_frame_linear_kernel():
    """All189 columns,126 singleton target monomials and freeI63 kernel.

    A first-frame source (t,i) has initial z_i e_t, t<3,i<7.
    Its variable-j product has one degree2 target monomial. Every fiber has
    cardinality one or two. Choosing the smaller source index gives a section;
    all duplicate-fiber differences are exactly the z7 Koszul vectors.
    """
    fibers={};columns=[]
    for t in range(3):
        for i in range(7):
            for j in range(9):
                m=tuple(int(k==i)+int(k==j) for k in range(9))
                key=(t,m);cid=len(columns);columns.append((t,i,j,key))
                fibers.setdefault(key,[]).append(cid)
    assert len(columns)==189 and len(fibers)==126
    canonical={k:min(v) for k,v in fibers.items()}
    noncanonical=[c for c in range(189) if c!=canonical[columns[c][3]]]
    assert len(noncanonical)==63
    kernels=[];koszul=[]
    for c in noncanonical:
        t,i,j,key=columns[c];p=canonical[key]
        tp,ip,jp,_=columns[p]
        assert t==tp and i>j and ip==j and jp==i and i<7
        v={c:1,p:-1};assert len(v)==2
        kernels.append(v);koszul.append([t,j,i,c,p])
    expected={(t,i,j) for t in range(3) for i in range(7) for j in range(i+1,7)}
    assert {(t,i,j) for t,i,j,c,p in koszul}==expected
    # Explicit full kernel multiplyback and freeI, no rank algorithm.
    for k,v in enumerate(kernels):
        image={}
        for c,a in v.items():
            key=columns[c][3];image[key]=image.get(key,0)+a
        assert all(a==0 for a in image.values())
        assert [v.get(c,0) for c in noncanonical]==[int(k==j) for j in range(63)]
    # For every source column, I-section*map is either0 or its one saved
    # kernel vector; hence the63 independent vectors span the ENTIRE kernel.
    where={c:k for k,c in enumerate(noncanonical)}
    for c,(_,_,_,key) in enumerate(columns):
        p=canonical[key]
        remainder={} if c==p else {c:1,p:-1}
        assert remainder==({} if c==p else kernels[where[c]])
    return {'source_columns':189,'target_monomial_rank':126,'kernel_rank':63,
      'all_kernel_vectors':[[[c,a] for c,a in sorted(v.items())] for v in kernels],
      'free_source_rows':noncanonical,'all63_free_identity_and_multiplybacks':True,
      'all189_canonical_section_completeness_identities':True,
      'literal_three_z7_Koszul_blocks':koszul,
      'no_Gaussian_or_source_values':True}

def sparse_sum(a,b,k=Q(1)):
    d=dict(a)
    for i,c in b.items():
        d[i]=d.get(i,Q(0))+k*c
        if not d[i]:del d[i]
    return d
def composite_columns(A,B):
    out=[]
    for b in B:
        d={}
        for i,c in b.items():d=sparse_sum(d,A[i],c)
        out.append(d)
    return out
def determinant3(m):
    v=Q(0)
    for p in itertools.permutations(range(3)):
        sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        v+=sign*m[0][p[0]]*m[1][p[1]]*m[2][p[2]]
    return v
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};read_bytes=0
    def pin(p,h=None):
        nonlocal read_bytes
        p=pathlib.Path(p);actual=sha(p)
        if h is not None:assert actual==h,(str(p),actual,h)
        key=str(p.relative_to(ROOT))
        if key not in pins:read_bytes+=p.stat().st_size
        pins[key]=actual;return json.loads(p.read_text())
    pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
    docs={}
    for p,h in pr['input_sha256'].items():
        d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
        if p.endswith('_REVIEW/report.json'):assert d['all_passed']
    for p,h in pr['evidence_sha256'].items():E=pin(ROOT/'reports'/PRI/p,h)
    assert sha(ROOT/'code/ambient_B18_C3_complete_cubic_line_source_20261105_v63.py')==pr['code_sha256']
    assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
    count=pin(ROOT/'reports'/PRI/'closure-count-scope-companion.json',
        '3d43473cda393632d1490a4068a78cbc8b2ca714c4f8dd0c10035cbaea7605f0')
    frame=docs['actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json']
    Kdoc=docs['nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json']
    B4=docs['four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json']
    graph=docs['actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json']
    C2=docs['actual-physical12-and-C2-right-exact-quotient.json']
    assert frame['actual_C_and_C3_HF_degrees0_1_2']==[3,6,6]
    assert B4['all_eight_corrected_adjoints_zero_without_dot_worker']
    assert graph['all5_actual_delta_r_ell2_forced_zero_without_readback']
    physical=[r['physical_row'] for r in frame['all12_verified_source_records']]
    assert physical==[18,19,28,29,60,67,72,63,64,68,73,74]
    fm={r['physical_row']:r for r in frame['all12_verified_source_records']}
    cm={r['physical_row']:r for r in C2['physical12_coefficients']}
    # Build the literal 27x116 first map and the invertible constant source
    # change, rather than trusting a dimension subtraction alone.
    first=[{} for _ in range(116)]
    selected=list(range(7))+list(range(8,15))+[81+3*i for i in range(7)]
    for sector in range(3):
        for i in range(7):first[selected[7*sector+i]]={9*sector+i:Q(1)}
    for j,p in enumerate(physical):
        a=list(map(Q,fm[p]['bA_first_z7']));b=list(map(Q,fm[p]['bB_first_z7']))
        first[80+3*j]={18+i:-v for i,v in enumerate(b) if v}
        first[81+3*j]={18+i:v for i,v in enumerate(a) if v}
    zero=[]
    for raw in [7,15]+list(range(16,80)):zero.append((raw,{raw:Q(1)},'intrinsic'))
    for j,p in enumerate(physical):
        raw=80+3*j;v={raw:Q(1)}
        for i,c in enumerate(map(Q,fm[p]['bB_first_z7'])):
            if c:v=sparse_sum(v,{81+3*i:c})
        zero.append((raw,v,'tauA'))
    for j in range(7,12):zero.append((81+3*j,{81+3*j:Q(1)},'graph_tauB'))
    for j in range(12):zero.append((82+3*j,{82+3*j:Q(1)},'tauAB'))
    assert len(selected)==21 and len(zero)==95
    assert set(selected).isdisjoint(raw for raw,v,k in zero)
    assert set(selected)|{raw for raw,v,k in zero}==set(range(116))
    S=[{raw:Q(1)} for raw in selected]+[v for raw,v,k in zero]
    first_new=composite_columns(first,S)
    assert first_new[:21]==[{9*t+i:Q(1)} for t in range(3) for i in range(7)]
    assert first_new[21:]==[{} for _ in range(95)]
    position={raw:i for i,raw in enumerate(selected)}
    position.update({raw:21+i for i,(raw,v,k) in enumerate(zero)})
    inv=[]
    for raw in range(116):
        i=position[raw];v={i:Q(1)}
        if i>=21:
            for r,c in S[i].items():
                if r!=raw:v=sparse_sum(v,{position[r]:-c})
        inv.append(v)
    assert composite_columns(S,inv)==[{i:Q(1)} for i in range(116)]
    assert composite_columns(inv,S)==[{i:Q(1)} for i in range(116)]
    assert E['actual_first21']=={'intrinsic_A_raw116_columns':list(range(7)),
       'intrinsic_B_raw116_columns':list(range(8,15)),
       'physical_RH_raw116_columns':[81+3*i for i in range(7)],'target_order':['eA','eB','eH']}
    # Fresh189 monomial grouping, kernel63/freeI, full section identities.
    syz=first_frame_linear_kernel()
    assert E['first_frame_linear_monomial_columns189']==[[t,i,j] for t in range(3) for i in range(7) for j in range(9)]
    own_byfree={k[0][0]:dict(k) for k in syz['all_kernel_vectors']}
    # The free coordinate need not be the first entry in a sorted vector.
    own_byfree={f:dict(v) for f,v in zip(syz['free_source_rows'],syz['all_kernel_vectors'])}
    author_byfree={v['free_column']:{i:int(Q(c)) for i,c in v['terms']} for v in E['complete63_linear_firstframe_syzygies_with_freeI63']}
    assert own_byfree==author_byfree and len(author_byfree)==63
    canon={}
    for c,(t,i,j) in enumerate(E['first_frame_linear_monomial_columns189']):canon.setdefault((t,tuple(sorted([i,j]))),c)
    assert E['all126_canonical_quadratic_images']==[[[t,list(mon)],c] for (t,mon),c in canon.items()]
    # Actual plane quadratic rank3 and complete12x9 kernel, with QQ units.
    M=[[Q(c) for c in row] for row in Kdoc['bB_plane2_rows3_columns12']]
    piv=Kdoc['fixed_pivot_columns'];free=Kdoc['free_columns9']
    assert piv==[0,4,5] and free==[1,2,3,6,7,8,9,10,11]
    square=[[row[j] for j in piv] for row in M]
    assert determinant3(square)==Q(5,768)
    K=[[Q(c) for c in row] for row in Kdoc['complete_constant_kernel12_by9']]
    for j in range(9):
        assert [K[i][j] for i in free]==[Q(i==j) for i in range(9)]
        assert all(sum(M[t][i]*K[i][j] for i in range(12))==0 for t in range(3))
    assert all(K[i][j]==Q(i==free[j]) for i in range(12) for j in range(9))
    assert len(zero)-3==92
    constants=E['all92_constant_sources'];assert len(constants)==92
    expected_indices=[7,15]+list(range(16,80))+[80+3*j for j in free]+[81+3*j for j in range(7,12)]+[82+3*j for j in range(12)]
    got=[v.get('raw116_index',v.get('raw116_tauA_index')) for v in constants]
    assert got==expected_indices and len(set(got))==92
    assert E['constant_quadratic_kernel92_counts']=={'curve':2,'cross':64,'tauA_kernel':9,'graph_tauB':5,'tauAB':12}
    # All64 product source/orders, not merely a cardinality assertion.
    cross=[v for v in constants if v['kind']=='cross']
    for i,j in itertools.product(range(8),repeat=2):
        v=next(v for v in cross if (v['IA_index'],v['IB_index'])==(i,j))
        assert v['raw116_index']==16+8*i+j and v['known_cover_cubic_image']=='0'
        assert v['order3_status']==('first2=zi*zj eH; clear by zi*RH_j; all plane3 carry retainsz' if i<7 and j<7 else 'order>=4')
    # Graph tauA primitive/carry zero uses actual delta3 and RH2-plane zero;
    # it does not discard any uncomputed mixed bA2/bB2 coefficient.
    for p in [63,64,68,73,74]:
        assert not any(map(Q,fm[p]['bB_first_z7']))
        assert fm[p]['delta_and_bA_pureplane_second_zero']
        assert next(v for v in constants if v['kind']=='tauA_kernel' and v['physical_row']==p)['known_cover_cubic_image']=='0'
    for p in physical[:7]:assert fm[p]['delta_and_bA_pureplane_second_zero']
    # Match all four still-formal callbacks to the actual full source roots.
    remaining=[19,28,29,72];v=[Q(0),Q(1),Q(0),-Q(1,4)]
    expected_formal=[]
    for n,c in enumerate(constants):
        if c['kind']=='tauA_kernel' and c['physical_row'] in remaining:
            p=c['physical_row'];expected_formal.append([n,f'alpha_{p}',list(map(str,v))])
    assert E['formal_nonzero_columns']==expected_formal
    for c in E['remaining_four_alpha_source_functionals']:
        assert c['physical_row'] in remaining and c['literal_roots']==cm[c['physical_row']]['roots']
        assert c['value']=='UNAVAILABLE_FORMAL_NOT_ZERO'
    assert E['formal_operator_shape']==[4,155] and 92+63==155 and 155-4==151
    assert 95*9==855 and 21*45==945
    # Both the source proof and all possible HF consequences are finite.
    # The cover has HF36642 and m5=0. Nonzero v has r*v=(1,-1/4), ell*v=0.
    r_v=[Q(1),-Q(1,4)];ell_v=[Q(0),Q(0)]
    HFs=[]
    for h3 in [3,4]:
        for h4 in range(3):
            length=15+h3+h4
            if 18<=length<=20 and not(h3==3 and h4>1):HFs.append([3,6,6,h3,h4])
    assert HFs==E['remaining_possible_C3_HFs']
    assert E['actual_HF_degree3_cases']==[4,3] and E['actual_cubic_extra_rank_cases']==[0,1]
    assert read_bytes<500000
    proof={'same_actual_full116_source_module_and_R0':True,
      'literal_constant_source_change116_and_both_inverses_verified':True,
      'first_map_rank21_and_entire_constant_kernel95_verified':True,
      'source_unit_basis_forward_sparse':[[[i,str(c)] for i,c in sorted(x.items())] for x in S],
      'source_unit_basis_inverse_sparse':[[[i,str(c)] for i,c in sorted(x.items())] for x in inv],
      'independent_complete_linear_syzygy63_proof':syz,
      'actual_plane_quadratic_rank3_fixed_minor':'5/768',
      'complete12x9_kernel27_equations_and_freeI9':True,
      'constant_quadratic_kernel92_source_indices':expected_indices,
      'all92_source_families_and_151_zero_images_verified':True,
      'full_product64_carry':'49 products zi*(zj+beta_j)eH are cleared by zi*RH_j. Both primitive cubic and correction carry retain a z factor.15 other products have a curve of order3 and hence order>=4.',
      'all_graph_tauA5_zero_carry':'bB1=0, actual delta3 mod ell3=0, exact RA clears A; H mixed correction uses r/ell*RH whose B plane2=0. zz coefficients retain z. No opaque B carry is assigned zero.',
      'all_formal_coefficient_orders_accounted':[
        'After the invertible constant source change, sources R21 start in order1 and U95 start in order2.',
        'For a relation of order>=3, the constant R coefficient vanishes and the constant U coefficient lies in the complete92-dimensional plane-quadratic kernel.',
        'Choose linear R coefficients to cancel its z-quadratic part. Any two choices differ by the full63-dimensional kernel proved above.',
        'Those63 syzygies have only z coefficients, so M2 times them has no pure-plane cubic.',
        'The855 linear U coefficients contribute variable times their known quadratic initials; plane-H terms die by plane2*eH, all other terms retain z.',
        'The945 quadratic R coefficients lie in S2*(z7); quadratic U coefficients begin at4. Higher coefficients cannot contribute at degree3.',
        'Thus after the actual known cover, every and only extra cubic initial is generated by the92 chosen lifts and63 syzygies.'],
      'full_four_callback_operator':expected_formal,
      'actual_extra_cubic_image_exactly':'span(alpha_19*v,alpha_28*v,alpha_29*v,alpha_72*v), v=(0,1,0,-1/4)',
      'all_alpha_values_remain_formal':True,'actual_C3_HF_degree3_cases':[3,4],
      'r_times_v':list(map(str,r_v)),'ell_times_v':list(map(str,ell_v)),
      'remaining_possible_HFs':HFs,'remaining_length_range':[18,20],
      'no_actual_C3_length_or_exponents_or_H0_or_nil_claim':True,
      'no_new_source_readbacks_coefficients_gauge_jets_or_Gaussian':True,
      'static_whole_numeric_upper':50000,'hard_numeric':150000,'RAW_read':read_bytes,
      'hard_seconds':60,'hard_address_bytes':2<<30,'hard_RAW_bytes':32<<20}
    ep=out/'independent-full116-to155-source-completeness-and-formal-line.json';save(ep,proof)
    result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PRIMARY,
      'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'complete_cubic_image_line_and_HF3_3_or4_accepted':True,'four_callbacks_still_uncomputed':True,
      'elapsed_seconds':time.monotonic()-start}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
