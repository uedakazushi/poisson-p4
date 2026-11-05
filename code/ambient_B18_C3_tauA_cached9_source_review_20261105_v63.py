#!/usr/bin/env python3
"""Independent tauA source/kernel/functional/prepared-code audit; no dots."""
import ast,hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63'
PRIMARY='8234e80934cba6d761eb72666d3d4337fcfa07a1c5d12ca9b0776f85d9339bd0'
WORKER='5891d81746e4eb9ad9df0e5d93879c0d8246cadfc496e03a17e79ccf42249888'
def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def load(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def field(t,k):
    i=t.index(':',t.index(json.dumps(k)))+1
    while t[i].isspace():i+=1
    return json.JSONDecoder().raw_decode(t,i)[0]
def mm(a,b):return [[sum((x*y for x,y in zip(row,col)),Q(0)) for col in zip(*b)] for row in a]
def enc(x):return [[str(c) for c in row] for row in x]

def main():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        rp=ROOT/'reports'/PRI/'report.json';assert sha(rp)==PRIMARY;rep=load(rp);assert rep['all_passed'] and not rep['cubic_values_evaluated']
        pins={str(rp.relative_to(ROOT)):PRIMARY}
        for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        for p,h in rep['evidence_sha256'].items():
            pp=ROOT/'reports'/PRI/p;assert sha(pp)==h;pins[str(pp.relative_to(ROOT))]=h
        D=load(ROOT/'reports'/PRI/'nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json')
        for p in rep['input_sha256']:
            if p.endswith('_REVIEW/report.json'):assert load(ROOT/p)['all_passed']
        def src(suffix):return next(ROOT/p for p in rep['input_sha256'] if p.endswith(suffix))
        raw=load(src('actual-ell2-product33-and-raw45-actions.json'));basis=raw['actual33_basis']
        ix3=basis.index([0,3,0]);ix21=basis.index([0,2,1]);assert [ix3,ix21]==[29,30]
        readback=[]
        for mon,v in zip(raw['basis74'],raw['all74_normal_forms']):
            cols=dict(v);a=Q(cols.get(ix3,'0'));b=Q(cols.get(ix21,'0'))
            assert a==Q(mon==[0,3,0]) and b==Q(mon==[0,2,1])
            readback.append([mon,[str(a),str(b)]])
        assert readback==D['raw33_literal_pureplane_cubic_dual_receipts'] and len(readback)==74
        frame=load(src('actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json'))
        actual=frame['all12_verified_source_records'];order=[r['physical_row'] for r in actual]
        assert order==D['physical_row_order']
        B=[[Q(r['bB_second_plane'][i]) for r in actual] for i in range(3)]
        assert enc(B)==D['bB_plane2_rows3_columns12']
        piv=D['fixed_pivot_columns'];assert piv==[0,4,5]
        P=[[B[i][j] for j in piv] for i in range(3)];Pi=[[Q(c) for c in r] for r in D['fixed_inverse']]
        I=[[Q(i==j) for j in range(3)] for i in range(3)]
        assert mm(P,Pi)==mm(Pi,P)==I
        det=Q(0)
        for p in itertools.permutations(range(3)):
            v=Q((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)))
            for i,j in enumerate(p):v*=P[i][j]
            det+=v
        assert det==Q(5,768)==Q(D['fixed_minor'])
        K=[[Q(c) for c in row] for row in D['complete_constant_kernel12_by9']];free=D['free_columns9']
        assert free==[j for j in range(12) if j not in piv] and len(K)==12 and all(len(r)==9 for r in K)
        assert mm(B,K)==[[Q(0)]*9 for _ in range(3)]
        assert [K[i] for i in free]==[[Q(i==j) for j in range(9)] for i in range(9)]
        # With P invertible, each pivot block is uniquely -P^-1 Bfree: this proves kernel completeness.
        for i,pi in enumerate(piv):
            for j,f in enumerate(free):assert K[pi][j]==-sum((Pi[i][a]*B[a][f] for a in range(3)),Q(0))
        pt=src('actual-special-fiber20-coefficients.json').read_text();g=field(pt,'gA_plane3');del pt
        ga=[Q(0)]*4
        for (nr,ns),c in g:
            assert nr+ns==3
            for j in range(ns+1):ga[j]+=Q(c)*math.comb(ns,j)*Q(1,5**j)
        assert ga==[Q(125,24),Q(5,3),Q(1,6),Q(2,375)]
        ratio=ga[1]/ga[0];assert ratio==Q(8,25)==Q(D['gA_reduction_ratio'])
        assert ga[1]-ratio*ga[0]==0 # functional kills gA3; ell3 and r*ell2 have weight0.
        # Only form the symbolic adjoints, never pair a signed cached cubic leaf with any weight here.
        adj=[]
        for j in range(9):
            v=[]
            for i in range(12):
                if K[i][j]:v.extend([[i,'delta',ix21,str(K[i][j])],[i,'delta',ix3,str(-ratio*K[i][j])]])
            adj.append(v)
        assert adj==D['nine_cached_eA_cubic_adjoints']
        weights=[Q(x[3]) for row in adj for x in row];weightbits=max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in weights)
        C2=load(src('actual-physical12-and-C2-right-exact-quotient.json'));leaves=[]
        for i,r in enumerate(C2['physical12_coefficients']):
            assert r['physical_row']==order[i];a=dict(r['delta'])
            for idx in [ix3,ix21]:
                x=Q(a.get(idx,'0'));leaves.append([i,idx,max(abs(x.numerator).bit_length(),x.denominator.bit_length())])
        assert leaves==D['cache_leaf_size_receipts_only_no_dots']
        h=max(r[2] for r in leaves);n=max(map(len,adj));dotbound=n*(h+weightbits)+math.ceil(math.log2(n))+1
        assert (weightbits,h,n,dotbound)==(5,18,2,48)
        assert D['height_profile']['source_dot_bits_bound']==dotbound and sum(map(len,adj))==18
        custodypath=ROOT/'reports/AMBIENT_B18_C3_TAUA_CACHED9_VALUE_CUSTODY_20261105_v63/report.json'
        assert sha(custodypath)=='4a1c34a0cbd9a25028ab06ed52f0c086f80a31dad4c2f57f91d6da9c755daf52'
        custody=load(custodypath);assert custody['all_passed'];pins[str(custodypath.relative_to(ROOT))]=sha(custodypath)
        for p,h in custody['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        companion=custodypath.parent/'independent-overflow-fixture-companion.json'
        assert sha(companion)=='f8eba09d71b13e1433fcfe5518dc5e5c7aead9564f12a777d79921354f74cd3e'
        comp=load(companion);assert comp['all_passed'] and comp['blocked_before_product']
        assert comp['worker_sha256']==WORKER and sum(comp['operand_numerator_bits'])==2051>2048
        pins[str(companion.relative_to(ROOT))]=sha(companion)
        worker=ROOT/'code/ambient_B18_C3_tauA_cached9_values_20261105_v63.py';assert sha(worker)==WORKER
        code=worker.read_text();tree=ast.parse(code);functions={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
        assert {'parse','product','plus','main','numeric_count'}<=set(functions)
        seg=lambda f:ast.get_source_segment(code,functions[f])
        assert 'assert max(' in seg('product') and seg('product').index('assert max(')<seg('product').index('return a * b')
        assert 'a.denominator.bit_length() + b.denominator.bit_length()' in seg('product')
        assert seg('plus').index('assert max(')<seg('plus').index('return a + b')
        assert '+ 1,' in seg('plus') and 'a.denominator.bit_length() + b.denominator.bit_length()' in seg('plus')
        assert seg('parse').index('assert max(')<seg('parse').index('q = Q(s)')
        assert 'p = product(w, x)' in code and 'value = plus(value, p)' in code
        assert 'owned_bound = input_count + 6000' in code and 'owned_bound <= 18000 < 30000' in code
        assert 'inputs_raw + 100000 <= 400000' in code and 'out.mkdir(exist_ok=False)' in code
        assert not any(isinstance(n,ast.ImportFrom) and n.module and n.module.startswith('ambient_') for n in tree.body)
        # The value input is only source report/evidence, existing C2 cache and this compact peer report.
        worker_in= (ROOT/'reports'/PRI/'report.json').stat().st_size + (ROOT/'reports'/PRI/'nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json').stat().st_size + src('actual-physical12-and-C2-right-exact-quotient.json').stat().st_size
        assert worker_in+20000+100000<400000
        proof={'same_actual_C_and_R0_frame':True,'all_fixed_source_and_custody_hashes_verified':True,
          'complete3x12_bB2_map_and12x9_kernel_verified':True,'both3x3_inverse_products':True,
          'all27_kernel_equations':True,'free_tail_I9_and_pivot_formula_prove_completeness':True,
          'full_kernel12x9':enc(K),'free_physical_rows':[order[i] for i in free],
          'all148_raw33_normal_form_coefficients_are_literal_Taylor_functionals':True,
          'multi_z_cubic_reason':'Monic z_i(z_j+beta_j) reduces two z to beta(order2)*z, raising order and retaining z; filtered normal forms cannot introduce a pureplane cubic from those higher terms.',
          'source_first_clearing':'q_l=F_l eA+sum_i bB_l,zi RH_i: zero first, exact eA=delta_l, plane2=-bB_l2 eH.',
          'second_clearing':'For K in kerbB2 all remaining M2 has z. eA clears use exact zi eA, hence no higher eA carry. RB and RH retain their actual mixed/beta terms but have zero eA component.',
          'unknown_eB_eH_carries_retained_and_annihilated_by_eA_projection':True,
          'exact_projected_cubic_formula':'pi_eA(corrected M3_K)=sum_l K_l delta_l3',
          'gA3_coefficients_r3_r2ell_rell2_ell3':list(map(str,ga)),
          'eA_functional':'delta_r2ell-(8/25)*delta_r3; weights of unknown rell2/ell3 and all eB/H are0.',
          'gA3_and_ell3_annihilated_by_functional':True,'symbolic_nine_adjoints_match_without_pairing':True,
          'adjoint_terms18_and_leaves24':True,'prior_dot_height48_bits':True,
          'operational_parse_product_add_guards_before_each_QQ_operation_verified':True,
          'overflow_harness_correction_preserves_unrun_worker':True,
          'static_owned_worker18000_hard30000_and_RAW400000_cap_accepted':True,
          'positive_conditional_proof':'Nonzero phi gives an actual cubic outside the plain length21 initial model. Its ell multiple has nonzero eA r2ell2 in CI(3,3), so adds an independent quartic; C3 length<=19.',
          'quartic_nonzero_reason':'Modulo ell3, r*gA3 and ell*gA3 have a triangular unit2x2 block on r4,r3ell. Thus r2ell2 survives, and ell times any cubic maps to phi times this class.',
          'no_conditional_upper18_or_both_torsion_exponents2_claim':True,
          'allzero_consequence_only_nine_projections_negative':True,
          'no_actual_nine_pairings_or_new_source_jet_values_evaluated':True,
          'source_review_hard_numeric_cap':150000,'future_worker_static_numeric_cap':18000,
          'future_worker_hard_numeric_cap':30000,'future_worker_RAW_cap':400000,'seconds':60,'address_bytes':2<<30}
        save(out/'independent-nine-kernel-eA-functional-and-guard-proof.json',proof)
        result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
          'code_sha256':sha(__file__),'evidence_sha256':{'independent-nine-kernel-eA-functional-and-guard-proof.json':sha(out/'independent-nine-kernel-eA-functional-and-guard-proof.json')},
          'complete_kernel_and_coherent_eA_projection_accepted':True,'prepared_worker_sha256':WORKER,
          'no_reviewer_nine_dot_evaluation':True,'conditional_nonzero_upper19_only':True,'elapsed_seconds':time.monotonic()-started}
        save(out/'report.json',result);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}))
    except Exception as e:
        save(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_source_PASS_or_dot_result':True});raise
if __name__=='__main__':main()
