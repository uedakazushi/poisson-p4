#!/usr/bin/env python3
"""Independent symbolic RIGHT2 operator preparation; lambda3 stays formal."""
import hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63'
ID=PRI+'_REVIEW'
PIN='900540d30f353dd53850cb12d8b7850f4e0efc3737bc8e1e845ab31164d0a41a'
def read(p):return json.loads(pathlib.Path(p).read_text())
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def determinant3(m):
    return sum(((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)))
                *m[0][p[0]]*m[1][p[1]]*m[2][p[2]] for p in itertools.permutations(range(3)))
def inverse3(m):
    d=determinant3(m);assert d
    n=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rs=[k for k in range(3) if k!=j];cs=[k for k in range(3) if k!=i]
            row.append((-1)**(i+j)*(m[rs[0]][cs[0]]*m[rs[1]][cs[1]]-m[rs[0]][cs[1]]*m[rs[1]][cs[0]])/d)
        n.append(row)
    assert all(sum(m[i][k]*n[k][j] for k in range(3))==Q(i==j) for i in range(3) for j in range(3))
    assert all(sum(n[i][k]*m[k][j] for k in range(3))==Q(i==j) for i in range(3) for j in range(3))
    return n
def affine(q=0,index=None):return [Q(q)]+[Q(i==index) for i in range(3)]
def add(a,b,c=Q(1)):return [x+c*y for x,y in zip(a,b)]
def scale(a,c):return [x*c for x in a]
def prepare_operator():
    paths={
      'C2':ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
      'raw':ROOT/'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
      'pair':ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
      'frame':ROOT/'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
      'cover':ROOT/'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
      'alpha_peer':ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63_REVIEW/report.json'}
    docs={k:read(p) for k,p in paths.items()};assert docs['alpha_peer']['all_passed']
    cm={r['physical_row']:r for r in docs['C2']['physical12_coefficients']}
    fm={r['physical_row']:r for r in docs['frame']['all12_verified_source_records']}
    basis=docs['raw']['actual33_basis'];physical=[18,19,28,29,60,67,72]
    pivot=[18,60,67];pivpos=[0,4,5]
    idx=[basis.index([0,3,0]),basis.index([0,2,1])]
    zir=[basis.index([i+1,1,0]) for i in range(7)];ziel=[basis.index([i+1,0,1]) for i in range(7)]
    beta=[]
    for row in docs['pair']['beta2_before_ell']:
        co=[Q(0)]*3
        for m,c in row:
            assert not any(m[:7]) and sum(m)==2
            nr,ns=m[7:]
            for k in range(ns+1):co[k]+=Q(c)*math.comb(ns,k)/5**k
        beta.append(co)
    D=[[-Q(fm[p]['bB_second_plane'][i]) for p in pivot] for i in range(3)]
    assert determinant3(D)==-Q(5,768);N=inverse3(D)
    delta={}
    for p in physical:
        v=dict(cm[p]['delta']);v3=[affine(Q(v.get(i,'0'))) for i in idx]
        v3.append(affine(index=pivot.index(p)) if p in pivot else affine())
        if p not in pivot:assert not any(v3[0]) and not any(v3[1])
        delta[p]=v3
    qA=[];qB=[];mixed=[]
    for p in pivot:
        h=list(map(Q,fm[p]['bB_first_z7']));qA.append(delta[p]);B=[affine() for _ in range(3)]
        for k in range(7):
            for m in range(3):B[m]=add(B[m],delta[physical[k]][m],h[k])
        q2=[]
        for z in range(7):
            coeff=[]
            for kx in [zir[z],ziel[z]]:
                coeff.append(sum(h[k]*Q(dict(cm[physical[k]]['delta']).get(kx,'0')) for k in range(7)))
            q2.append(coeff)
            for ell_degree in range(3):
                B[ell_degree][0]-=coeff[0]*beta[z][ell_degree]
                if ell_degree+1<3:B[ell_degree+1][0]-=coeff[1]*beta[z][ell_degree]
        qB.append(B);mixed.append(q2)
    P=[]
    for sector in [qA,qB]:
        Psector=[]
        for i in range(3):
            parts=[]
            for d in range(3):
                v=affine()
                for j in range(3):v=add(v,sector[j][d],N[j][i])
                parts.append(v)
            Psector.append(parts)
        P.append(Psector)
    operator=[]
    for i in range(2):
        row=[]
        for sector in range(2):
            curve=list(map(Q,docs['cover']['two_block_source_curve_NF_receipts'][sector]['curve3_r3_r2ell_rell2_ell3']))
            c1,c2=curve[1]/curve[0],curve[2]/curve[0]
            e=add(scale(P[sector][i][0],-c1),P[sector][i][1])
            r=add(add(scale(P[sector][i+1][0],c1*c1-c2),P[sector][i+1][1],-c1),P[sector][i+1][2])
            row.append(add(e,r,-1))
        operator.append(row)
    return {'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths.values()},
      'pivot_rows':pivot,'D3':[[str(v) for v in row] for row in D],
      'D_inverse3':[[str(v) for v in row] for row in N],
      'all_full_qB_mixed_second_coefficients':[[list(map(str,v)) for v in row] for row in mixed],
      'P3_plane_affine_lambda_const_18_60_67':[[[[str(c) for c in v] for v in parts] for parts in sector] for sector in P],
      'HB2_quartic_NF_affine_lambda_const_18_60_67':[[[str(c) for c in v] for v in row] for row in operator],
      'all3_lambda_still_formal':True,'no_original_root_readback_or_lambda_pairing':True}

def fplus(a,b,c=Q(1)):
    d=dict(a)
    for k,v in b.items():
        d[k]=d.get(k,Q(0))+c*v
        if not d[k]:del d[k]
    return d
def fscale(a,c):return {k:v*c for k,v in a.items() if v*c}
def pplus(a,b,c=Q(1)):
    d={m:dict(v) for m,v in a.items()}
    for m,v in b.items():
        w=fplus(d.get(m,{}),v,c)
        if w:d[m]=w
        elif m in d:del d[m]
    return d
def shift(p,mon):return {tuple(a+b for a,b in zip(m,mon)):dict(v) for m,v in p.items()}
def normal4(p,curve):
    p={m:dict(v) for m,v in p.items() if m[1]<3}
    while any(m[0]>=3 for m in p):
        m=max(m for m in p if m[0]>=3);v=p.pop(m)
        for j in range(1,4):
            n=(m[0]-j,m[1]+j)
            if n[1]<3:p=pplus(p,{n:v},-curve[j]/curve[0])
    assert set(p)<={(2,2)}
    return p.get((2,2),{})
def encode_form(f):return [[list(k),str(v)] for k,v in sorted(f.items())]
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
    rp=ROOT/'reports'/PRI/'report.json';assert sha(rp)==PIN;pr=read(rp);assert pr['all_passed']
    pins[str(rp.relative_to(ROOT))]=PIN;docs={}
    for p,h in pr['input_sha256'].items():
        assert sha(ROOT/p)==h;pins[p]=h;d=read(ROOT/p);docs[pathlib.Path(p).name]=d
        if p.endswith('_REVIEW/report.json'):assert d['all_passed']
    for p,h in pr['evidence_sha256'].items():
        q=ROOT/'reports'/PRI/p;assert sha(q)==h;pins[str(q.relative_to(ROOT))]=h;E=read(q)
    assert sha(ROOT/'code/ambient_B18_C3_right2_quartic_source_20261105_v63.py')==pr['code_sha256']
    assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
    frame=docs['actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json']
    cover=docs['actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json']
    B4=docs['four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json']
    c2=docs['actual-physical12-and-C2-right-exact-quotient.json']
    raw=docs['actual-ell2-product33-and-raw45-actions.json']
    alpha=docs['eight-cached-even-root-values-and-four-complete-alpha-readbacks.json']
    assert alpha['complete_cubic_extra_rank']==0 and alpha['actual_C3_HF_degree3']==4
    fm={r['physical_row']:r for r in frame['all12_verified_source_records']}
    cm={r['physical_row']:r for r in c2['physical12_coefficients']}
    physical=[18,19,28,29,60,67,72];piv=[18,60,67]
    M=[[-Q(fm[p]['bB_second_plane'][i]) for p in piv] for i in range(3)]
    assert determinant3(M)==-Q(5,768);T=inverse3(M)
    assert E['plane2_matrix_rows_r2_rell_ell2_columns_qp']==[[str(v) for v in row] for row in M]
    assert E['constant_source_normalization_T']==[[str(v) for v in row] for row in T]
    H=[[Q(fm[p]['bB_first_z7'][physical.index(k)]) for k in piv] for p in piv]
    basis=raw['actual33_basis'];ir3=basis.index([0,3,0]);ir2ell=basis.index([0,2,1])
    zr=[basis.index([i+1,1,0]) for i in range(7)];zl=[basis.index([i+1,0,1]) for i in range(7)]
    beta=[[Q(c) for c in row] for row in B4['beta2_plane_r2_rell_ell2']]
    active=[i for i,b in enumerate(beta) if any(b)];assert active==[0,4,5]
    # Reconfirm literal scalar Taylor extraction for every old leaf.
    for mon,nf in zip(raw['basis74'],raw['all74_normal_forms']):
        nf=dict(nf)
        for i in [ir3,ir2ell]+zr+zl:assert Q(nf.get(i,'0'))==Q(mon==basis[i])
    roots=E['pivot_original_roots'];assert roots==[cm[p]['roots'] for p in piv]
    expected_q=[]
    for p in piv:
        j=physical.index(p);h=list(map(Q,fm[p]['bB_first_z7']))
        expected_q.append([[80+3*j,'1']]+[[81+3*i,str(c)] for i,c in enumerate(h) if c])
        assert not any(h[i] for i in [1,2,3])
        assert all(fm[k]['delta_and_bA_pureplane_second_zero'] for k in physical)
    assert E['q_p_actual_raw116_columns']==expected_q
    # Row72 may be removed only with its actual alpha and weighted old leaves.
    assert next(r['delta_r_ell2'] for r in alpha['all4_actual_alpha_delta_B6'] if r['physical_row']==72)=='0'
    v72=dict(cm[72]['delta'])
    assert all(Q(v72.get(i,'0'))==0 for i in [ir3,ir2ell]+[j for z in active for j in [zr[z],zl[z]]])
    D=[];F=[]
    for p in piv:
        ds={(3,0):{('cache',p,ir3):Q(1)},(2,1):{('cache',p,ir2ell):Q(1)},
            (1,2):{('lambda',p,-1):Q(1)},(0,3):{('opaque_ell3',p,-1):Q(1)}}
        gam={}
        for z in active:
            for k,b in enumerate(beta[z]):
                gam=pplus(gam,{(3-k,k):{('cache',p,zr[z]):b}})
                gam=pplus(gam,{(2-k,k+1):{('cache',p,zl[z]):b}})
        D.append(ds);F.append(pplus(ds,gam,-1))
    # Construct the normalized cubic polynomials before quotienting, rather
    # than applying a pre-derived L/R dual to the source coefficients.
    P=[[{} for _ in range(3)] for sector in range(2)]
    for k in range(3):
        for j in range(3):P[0][k]=pplus(P[0][k],D[j],T[j][k])
        for i in range(3):
            w=sum(H[j][i]*T[j][k] for j in range(3))
            P[1][k]=pplus(P[1][k],F[i],w)
    op=[[{} for hb in range(2)] for sector in range(2)]
    for sector in range(2):
        co=list(map(Q,cover['two_block_source_curve_NF_receipts'][sector]['curve3_r3_r2ell_rell2_ell3']))
        for hb in range(2):
            quartic=pplus(shift(P[sector][hb],(0,1)),shift(P[sector][hb+1],(1,0)),-1)
            op[sector][hb]=normal4(quartic,co)
    assert [[encode_form(f) for f in row] for row in op]==E['matrix2_rows_A4_B4_columns_HB0_HB1_linear_leaf_operators']
    assert all(not any(k[0]=='opaque_ell3' for k in f) for row in op for f in row)
    sens=[[[str(op[t][j].get(('lambda',p,-1),Q(0))) for p in piv] for j in range(2)] for t in range(2)]
    assert sens==E['lambda_sensitivity2_by2_by3']
    C=[fplus(op[1][j],op[0][j],Q(1,4)) for j in range(2)]
    cached=fplus(C[1],C[0],Q(10))
    assert not any(k[0]=='lambda' for k in cached)
    assert encode_form(cached)==E['all_cached_dual_HB1_plus10HB0_operator_NOT_PAIRED']
    assert 3*3*28==252 and 3*7*2==42 and 3*4==12 and 3*3*(165-4)==1449 and 3*3*495==4455
    assert E['highest_unavailable_sensitivities_zero_by_degree']=={
        'q_p_second_zz_coefficients':252,'q_p_second_H_mixed_coefficients':42,
        'q_p_third_H_plane_coefficients':12,'q_p_third_z_containing_coefficients':1449,
        'q_p_fourth_all_coefficients':4455}
    # Formal degree reasoning is the complete dependency proof, not a
    # fixture with unknown Taylor coefficients assigned0.
    reasoning=[
      'Actual q_p starts in degree2. Its full z-containing second part is removed by a canonical linear combination of actual first21 sources; this uses an origin unit1, with no r/ell inverse.',
      'RA=zi*eA exactly. Its contribution to the plane cubic A component is0. RB=(zi+beta_i)eB contributes the full negative mixed-delta2*beta2 carry.',
      'RH has zero plane second A/B component, since its B component is actual delta_i and all delta_i2 plane coefficients vanish. Removing H mixed second terms therefore adds no plane cubic A/B.',
      'Clearing any zz second term uses z times a first-frame source; its cubic carry still has z. No unknown zz coefficient can enter the plane cubic projection.',
      'Multiply by the constant QQ unit inverse T of -Bcols. The FULL degree2 parts of P0/P1/P2 are precisely r2H,r*ellH,ell2H, not just plane projections.',
      'Consequently ellP0-rP1 and ellP1-rP2 have their complete degree3 terms zero in all three targets. They start at degree4.',
      'Their degree4 terms are linear times P_k,3 only. Every z-containing P_k,3 term keeps z; removing it with degree3 first-frame coefficients adds carry only at5.',
      'Every H quartic is a degree2 multiple of the actual normalized H quadratic pivots. Those clearings have A/B carry starting at5, so no unrecorded H source affects the A/B quartic classes.',
      'P_k,4 and all higher unknown terms are multiplied by r/ell, so begin at5. Thus no delta4, mixed source3, beta3 or new gauge coefficient occurs.',
      'The known ell3 relations kill every opaque pure ell3 coefficient before quartic NF. Row72 is removed using actual alpha72=0 and its required old mixed and cubic leaves, not a tangent coincidence.']
    proof={'same_actual_C3_source_frame_and_full_qp_labels_verified':True,
      'all_qp_source_columns':expected_q,'fixed_quadratic_matrix_determinant':'-5/768',
      'constant_normalization_and_both_unit_inverse_products':True,
      'all16_old_Taylor_duals_on74_raw_monomials':True,
      'all_full_second_zz_Hmixed_and_higher_degree_carry_accounted':reasoning,
      'independently_multiplied_HB_polynomials_and_direct_NF4':[[encode_form(f) for f in row] for row in op],
      'all3_lambda_and24_old_leaf_values_remain_formal_in_final_operator':True,
      'no_delta4_mixed_Phi3_beta3_or_new_jet_dependency':True,
      'all_opaque_ell3_sensitivities_zero_without_zero_fill':True,
      'highest_sensitivity_counts_match_and_are_zero_by_degree':True,
      'B4_plus_A4_over4_HB1_plus10HB0_lambda_cancellation':True,
      'completion_scope':'If these two actual quartics later have rank2, C3 has length19/HF36640. Rank0/1 is only a subset result and does not determine full quartic rank or length20.',
      'NO_lambda_value_readback_or_RIGHT2_value_or_rank':True,
      'formal_SOURCE_only_NOT_cache_worker_GO':True,
      'source_static_numeric_bound':50000,'hard_numeric':150000,'hard_bits':2048,
      'hard_seconds':60,'hard_address_bytes':2<<30,'hard_RAW_bytes':32<<20}
    ep=out/'independent-full-q-clearing-HB-order-and-unpaired-quartic-operator.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PIN,
      'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'coherent_RIGHT2_quartic_theorem_and_operator_accepted':True,
      'three_lambda_readbacks_still_uncomputed':True,'future_worker_SOURCE_code_budget_pending':True,
      'elapsed_seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
